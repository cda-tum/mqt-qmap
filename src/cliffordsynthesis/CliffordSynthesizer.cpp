//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#include "LogicTerm/Logic.hpp"
#include "utils/logging.hpp"

#include <chrono>

namespace cs {

void CliffordSynthesizer::synthesize(const Configuration& config) {
  configuration = config;

  INFO() << "Optimization target: " << toString(configuration.target);

  const auto start = std::chrono::high_resolution_clock::now();

  // create the general configuration for the SAT encoder
  auto encoderConfig           = EncoderConfig();
  encoderConfig.initialTableau = &initialTableau;
  encoderConfig.targetTableau  = &targetTableau;
  encoderConfig.nQubits        = initialTableau.getQubitCount();
  encoderConfig.timestepLimit  = configuration.initialTimestepLimit;
  encoderConfig.targetMetric   = configuration.target;
  encoderConfig.useMaxSAT      = configuration.useMaxSAT;
  encoderConfig.nThreads       = configuration.nThreads;
  encoderConfig.useMultiGateEncoding =
      requiresMultiGateEncoding(encoderConfig.targetMetric);

  // First, determine an initial guess for the number of timesteps. This can
  // either be specified as a configuration parameter or starts at 1.
  determineInitialTimestepLimit(encoderConfig);

  // Then, determine an upper bound for the number of timesteps by solving the
  // SAT problem repeatedly with increasing timestep limits until a satisfying
  // assignment is found. This uses the general SAT encoding without any
  // objective function regardless of the configuration.
  const auto [lower, upper] = determineUpperBound(encoderConfig);

  // if the upper bound is 0, the solution does not require any gates and the
  // synthesis is done.
  if (upper == 0U) {
    INFO() << "No gates required.";
    return;
  }

  encoderConfig.timestepLimit = upper - 1U;
  if (configuration.target == TargetMetric::TWO_QUBIT_GATES) {
    encoderConfig.twoQubitGateLimit = results.getTwoQubitGates();
  }

  // Once a valid upper bound is found, the SAT problem is solved again with the
  // objective function encoded. For the MaxSAT solver, this involves a single
  // call to the solver. For the binary search approach, the SAT problem is
  // solved repeatedly until a solution with a certain timestep limit T is
  // found, but no solution with timestep limit T-1 could be determined.
  if (configuration.useMaxSAT) {
    runMaxSAT(encoderConfig);
  } else {
    runBinarySearch(encoderConfig.timestepLimit, lower, upper, encoderConfig);
  }

  if (configuration.target == TargetMetric::TWO_QUBIT_GATES) {
    // While the MaxSAT approach is guaranteed to find the optimal solution
    // within a given timestep limit, the binary search approach is not.
    // Therefore, we need to check whether the solution found by the binary
    // search approach is actually optimal. To this end, we iteratively try to
    // lower the limit on the number of two-qubit gates by one and check whether
    // a valid solution can be found.
    if (!configuration.useMaxSAT) {
      minimizeTwoQubitGatesFixedGateCount(results.getGates(), encoderConfig);
    }

    // At this point, we have found the optimal solution for the number of
    // two-qubit gates with respect to the considered gate-count limit. However,
    // it is possible that there is a solution with fewer two-qubit gates that
    // uses more gates overall. To find such a solution, we run the solver once
    // more with an increased gate count limit.
    if (configuration.tryHigherGateLimitForTwoQubitGateOptimization) {
      const auto gateLimit =
          std::max(static_cast<std::size_t>(
                       std::round(static_cast<double>(results.getGates()) *
                                  configuration.gateLimitFactor)),
                   results.getGates() + 1U);
      minimizeTwoQubitGatesFixedGateCount(gateLimit, encoderConfig);
    }

    // While the solution at this point is optimal with respect to the number of
    // two-qubit gates, it is possible that there is a solution with fewer gates
    // overall. To find such a solution, we run the solver once more with a
    // fixed limit on the number of two-qubit gates and the goal to minimize the
    // number of gates overall.
    if (configuration.minimizeGatesAfterTwoQubitGateOptimization) {
      minimizeGatesFixedTwoQubitGateCount(encoderConfig);
    }
  }

  if (configuration.target == TargetMetric::DEPTH &&
      configuration.minimizeGatesAfterDepthOptimization) {
    // While the solution at this point is guaranteed to be depth-optimal, it
    // might contain more gates than necessary. To find a solution with fewer
    // gates, we run the solver once more with a fixed depth limit and the goal
    // to minimize the number of gates.
    minimizeGatesFixedDepth(encoderConfig);
  }

  results.setSolverCalls(solverCalls);

  const auto end = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double> diff = end - start;
  INFO() << "Synthesis took " << diff.count() << " seconds";
  results.setRuntime(diff.count());
}

void CliffordSynthesizer::determineInitialTimestepLimit(EncoderConfig& config) {
  if (config.timestepLimit == 0U) {
    // Infeasibility checks are inexpensive, so we can afford to start with a
    // timestep limit of 1 and increase geometrically until we find a feasible
    // solution.
    config.timestepLimit = 1U;
  }
}

std::pair<std::size_t, std::size_t>
CliffordSynthesizer::determineUpperBound(EncoderConfig config) {
  std::size_t lowerBound = 0U;
  std::size_t upperBound = config.timestepLimit;

  INFO() << "Searching for upper bound for the number of timesteps starting "
         << "with " << upperBound;

  config.useMaxSAT = false;
  while (!results.sat()) {
    results = callSolver(config);
    if (!results.sat()) {
      lowerBound = upperBound + 1U;
      upperBound *= 2U;
      INFO() << "No solution found for " << config.timestepLimit
             << " timestep(s). Doubling timestep limit to " << upperBound;
      config.timestepLimit = upperBound;
    }
  }

  if (config.targetMetric == TargetMetric::GATES ||
      config.targetMetric == TargetMetric::TWO_QUBIT_GATES) {
    upperBound = std::min(upperBound, results.getGates());
  } else if (config.targetMetric == TargetMetric::DEPTH) {
    upperBound = std::min(upperBound, results.getDepth());
  }

  INFO() << "Found upper bound for the number of timesteps: " << upperBound;
  return {lowerBound, upperBound};
}

void CliffordSynthesizer::minimizeGatesFixedDepth(EncoderConfig config) {
  if (results.getDepth() == 0U) {
    return;
  }

  if (results.getDepth() == results.getGates()) {
    return;
  }

  INFO() << "Found a depth-optimal circuit with depth " << results.getDepth()
         << " and " << results.getGates()
         << " gate(s). Trying to minimize the number of gates.";

  config.targetMetric         = TargetMetric::GATES;
  config.timestepLimit        = results.getDepth();
  config.useMultiGateEncoding = true;
  config.useMaxSAT            = configuration.useMaxSAT;

  if (config.useMaxSAT) {
    runMaxSAT(config);
  } else {
    config.gateLimit = results.getGates();
    runBinarySearch(*config.gateLimit, results.getDepth(), results.getGates(),
                    config);
  }
  INFO() << "Found a depth " << results.getDepth() << " circuit with "
         << results.getGates() << " gate(s).";
}

void CliffordSynthesizer::minimizeTwoQubitGatesFixedGateCount(
    const std::size_t gateCount, CliffordSynthesizer::EncoderConfig config) {
  if (results.getTwoQubitGates() == 0U) {
    return;
  }

  INFO() << "Trying to find a solution with less than "
         << results.getTwoQubitGates() << " two-qubit gates and at most "
         << gateCount << " gates.";

  config.targetMetric         = TargetMetric::TWO_QUBIT_GATES;
  config.timestepLimit        = gateCount;
  config.useMultiGateEncoding = false;
  config.useMaxSAT            = true;
  config.twoQubitGateLimit    = results.getTwoQubitGates() - 1U;

  runMaxSAT(config);

  INFO() << "Found a circuit with " << results.getTwoQubitGates()
         << " two-qubit gate(s) and " << results.getGates()
         << " gate(s) overall.";
}

void CliffordSynthesizer::minimizeGatesFixedTwoQubitGateCount(
    CliffordSynthesizer::EncoderConfig config) {
  if (results.getGates() == 0U) {
    return;
  }

  if (results.getTwoQubitGates() == results.getGates()) {
    return;
  }

  INFO() << "Found a two-qubit gate-count-optimal circuit with "
         << results.getTwoQubitGates() << " two-qubit gate(s) and "
         << results.getGates()
         << " gate(s) overall. Trying to minimize the number of gates.";

  config.targetMetric         = TargetMetric::GATES;
  config.timestepLimit        = results.getGates();
  config.useMultiGateEncoding = false;
  config.useMaxSAT            = configuration.useMaxSAT;
  config.twoQubitGateLimit    = results.getTwoQubitGates();

  if (config.useMaxSAT) {
    runMaxSAT(config);
  } else {
    runBinarySearch(config.timestepLimit, results.getTwoQubitGates(),
                    results.getGates(), config);
  }
  INFO() << "Found a circuit with " << results.getTwoQubitGates()
         << " two-qubit gate(s) and " << results.getGates()
         << " gate(s) overall.";
}

void CliffordSynthesizer::runMaxSAT(const EncoderConfig& config) {
  INFO() << "Running MaxSAT scheme with timestep limit "
         << config.timestepLimit;
  const auto r = callSolver(config);
  if (r.sat()) {
    INFO() << "Found a solution.";
  } else {
    INFO() << "No solution found.";
  }
  updateResults(configuration, r, results);
}

Results CliffordSynthesizer::callSolver(const EncoderConfig& config) {
  ++solverCalls;
  auto encoder = encoding::SATEncoder(config);
  return encoder.run();
}

void CliffordSynthesizer::updateResults(const Configuration& config,
                                        const Results&       newResults,
                                        Results&             currentResults) {
  if (!newResults.sat()) {
    return;
  }

  switch (config.target) {
  case TargetMetric::GATES:
    if ((newResults.getGates() < currentResults.getGates()) ||
        ((newResults.getGates() == currentResults.getGates()) &&
         (newResults.getTwoQubitGates() < currentResults.getTwoQubitGates()))) {
      currentResults = newResults;
    }
    break;
  case TargetMetric::TWO_QUBIT_GATES:
    if ((newResults.getTwoQubitGates() < currentResults.getTwoQubitGates()) ||
        ((newResults.getTwoQubitGates() == currentResults.getTwoQubitGates()) &&
         (newResults.getGates() < currentResults.getGates()))) {
      currentResults = newResults;
    }
    break;
  case TargetMetric::DEPTH:
    if ((newResults.getDepth() < currentResults.getDepth()) ||
        ((newResults.getDepth() == currentResults.getDepth()) &&
         (newResults.getGates() < currentResults.getGates()))) {
      currentResults = newResults;
    }
    break;
  }
}
} // namespace cs
