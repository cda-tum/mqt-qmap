//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#include "LogicTerm/Logic.hpp"
#include "utils/logging.hpp"

#include <chrono>
#include <fstream>

namespace cs {

void CliffordSynthesizer::synthesize(const Configuration& config) {
  configuration = config;

  // initialize logging
  if (plog::get() == nullptr) {
    static plog::ConsoleAppender<plog::TxtFormatter> consoleAppender;
    plog::init(plog::none, &consoleAppender);
  }
  plog::get()->setMaxSeverity(configuration.verbosity);

  INFO() << "Optimization target: " << toString(configuration.target);

  const auto start = std::chrono::high_resolution_clock::now();

  // create the general configuration for the SAT encoder
  auto encoderConfig                = EncoderConfig();
  encoderConfig.initialTableau      = &initialTableau;
  encoderConfig.targetTableau       = &targetTableau;
  encoderConfig.nQubits             = initialTableau.getQubitCount();
  encoderConfig.timestepLimit       = configuration.initialTimestepLimit;
  encoderConfig.targetMetric        = configuration.target;
  encoderConfig.useMaxSAT           = configuration.useMaxSAT;
  encoderConfig.useSymmetryBreaking = configuration.useSymmetryBreaking;
  encoderConfig.nThreads            = configuration.nThreads;
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
  // Otherwise, the determined upper bound is used as an initial timestep limit.
  encoderConfig.timestepLimit = upper;

  // Once a valid upper bound is found, the SAT problem is solved again with
  // the objective function encoded.
  switch (config.target) {
  case TargetMetric::Gates:
    gateOptimalSynthesis(encoderConfig, lower, upper);
    break;
  case TargetMetric::Depth:
    if (configuration.heuristic) {
      depthHeuristicSynthesis(encoderConfig);
    } else {
      depthOptimalSynthesis(encoderConfig, lower, upper);
    }
    break;
  case TargetMetric::TwoQubitGates:
    twoQubitGateOptimalSynthesis(encoderConfig, 0U, results.getTwoQubitGates());
    break;
  }

  results.setSolverCalls(solverCalls);

  const auto end = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double> diff = end - start;
  INFO() << "Synthesis took " << diff.count() << " seconds";
  results.setRuntime(diff.count());
}

void CliffordSynthesizer::determineInitialTimestepLimit(EncoderConfig& config) {
  if (config.timestepLimit != 0U) {
    INFO() << "Using configured initial timestep limit: "
           << config.timestepLimit;
    return;
  }

  // in case no circuit was provided as input, the best guess is to start low
  // and increase the limit geometrically until a solution is found.
  if (!results.sat()) {
    config.timestepLimit = 1U;
    INFO()
        << "No initial circuit specified. Using initial timestep limit of 1.";
    return;
  }

  // If an initial circuit is specified, there already is a satisfying result.
  // We can use its gates or depth as a starting point for the number of
  // timesteps. Furthermore, no upper bound needs to be determined.
  if (requiresMultiGateEncoding(config.targetMetric)) {
    config.timestepLimit = results.getDepth();
    INFO() << "Using initial circuit's depth as initial timestep limit: "
           << config.timestepLimit;
  } else {
    config.timestepLimit = results.getGates();
    INFO() << "Using initial circuit's gate count as initial timestep limit: "
           << config.timestepLimit;
  }
}

std::pair<std::size_t, std::size_t>
CliffordSynthesizer::determineUpperBound(EncoderConfig config) {
  // In case the synthesis was started with a circuit, the upper bound is
  // inherently given by the circuit and does not need to be computed here.
  if (results.sat()) {
    return {0U, config.timestepLimit};
  }

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

  if (config.targetMetric == TargetMetric::Gates ||
      config.targetMetric == TargetMetric::TwoQubitGates) {
    upperBound = std::min(upperBound, results.getGates());
  } else if (config.targetMetric == TargetMetric::Depth) {
    upperBound = std::min(upperBound, results.getDepth());
  }

  INFO() << "Found upper bound for the number of timesteps: " << upperBound;
  return {lowerBound, upperBound};
}

void CliffordSynthesizer::gateOptimalSynthesis(EncoderConfig     config,
                                               const std::size_t lower,
                                               const std::size_t upper) {
  // Gate-optimal synthesis is achieved by determining a timestep limit T such
  // that there exists a solution with T gates, but no solution with T-1 gates.
  // This procedure uses an encoding where a single gate is allowed per timestep
  // and guarantees optimality, i.e., there is no solution with fewer gates.

  if (configuration.useMaxSAT) {
    // The MaxSAT solver can determine the optimal T with a single call by
    // minimizing over the number of applied gates.
    runMaxSAT(config);
  } else {
    // The binary search approach calls the SAT solver repeatedly with varying
    // timestep (=gate) limits T until a solution with T gates is found, but no
    // solution with T-1 gates could be determined.
    runBinarySearch(config.timestepLimit, lower, upper, config);
  }
}

void CliffordSynthesizer::depthOptimalSynthesis(
    CliffordSynthesizer::EncoderConfig config, const std::size_t lower,
    const std::size_t upper) {
  // Depth-optimal synthesis is achieved by determining a timestep limit T such
  // that there exists a solution with depth T, but no solution with depth T-1.
  // This procedure uses an encoding where multiple gates are allowed per
  // timestep (as long as they can be executed in parallel). This procedure is
  // guaranteed to produce a depth-optimal circuit. However, the number of gates
  // in the resulting circuit is not necessarily minimal, i.e., there may be a
  // solution with fewer gates and the same depth. To this end, an optimization
  // pass is provided that additionally minimizes the number of gates.

  if (configuration.useMaxSAT) {
    // The MaxSAT solver can determine the optimal T with a single call by
    // minimizing over the layers of gates (=timesteps) in the resulting
    // circuit.
    runMaxSAT(config);
  } else {
    // The binary search approach calls the SAT solver repeatedly with varying
    // timestep (=depth) limits T until a solution with depth T is found, but no
    // solution with depth T-1 could be determined.
    runBinarySearch(config.timestepLimit, lower, upper, config);
  }

  if (configuration.minimizeGatesAfterDepthOptimization) {
    // To find a solution with fewer gates, we run the solver once more with a
    // fixed depth limit and the goal to minimize the number of gates.
    minimizeGatesFixedDepth(config);
  }
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

  config.targetMetric         = TargetMetric::Gates;
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

void CliffordSynthesizer::twoQubitGateOptimalSynthesis(
    EncoderConfig config, const std::size_t lower, const std::size_t upper) {
  // Two-qubit gate-optimal synthesis is achieved by minimizing over the number
  // of two-qubit gates. This procedure uses the same encoding as gate-optimal
  // synthesis, but with a different objective function. In contrast to the
  // gate-optimal synthesis, this procedure is only guaranteed to produce a
  // two-qubit gate-optimal circuit with respect to a given timestep limit T.
  // There might be a solution with fewer two-qubit gates that requires more
  // gates overall. To this end, an optimization pass is provided that explores
  // whether increasing the timestep limit can reduce the number of two-qubit
  // gates. Furthermore, the number of gates in the resulting circuit is not
  // necessarily minimal, i.e., there may be a solution with fewer gates and the
  // same number of two-qubit gates. To this end, a further optimization pass is
  // provided that additionally minimizes the number of gates.

  if (configuration.useMaxSAT) {
    // The MaxSAT solver can determine the optimal number of two-qubit gates
    // with a single call by minimizing over the number of two-qubit gate
    // variables.
    runMaxSAT(config);
  } else {
    // The binary search approach calls the SAT solver repeatedly with varying
    // two-qubit gate count limits G until a solution with G two-qubit gates is
    // found, but no solution with G-1 two-qubit gates could be determined.
    config.twoQubitGateLimit = upper;
    runBinarySearch(*config.twoQubitGateLimit, lower, upper, config);
  }

  // To find a solution with even fewer two-qubit gates but more gates overall,
  // we run the solver once more with an increased gate count limit.
  if (configuration.tryHigherGateLimitForTwoQubitGateOptimization) {
    const auto gateLimit =
        std::max(static_cast<std::size_t>(
                     std::round(static_cast<double>(results.getGates()) *
                                configuration.gateLimitFactor)),
                 results.getGates() + 1U);
    minimizeTwoQubitGatesFixedGateCount(gateLimit, config);
  }

  // While the solution at this point is optimal with respect to the number of
  // two-qubit gates, it is possible that there is a solution with fewer gates
  // overall. To find such a solution, we run the solver once more with a
  // fixed limit on the number of two-qubit gates and the goal to minimize the
  // number of gates overall.
  if (configuration.minimizeGatesAfterTwoQubitGateOptimization) {
    minimizeGatesFixedTwoQubitGateCount(config);
  }
}

void CliffordSynthesizer::minimizeTwoQubitGatesFixedGateCount(
    const std::size_t gateCount, CliffordSynthesizer::EncoderConfig config) {
  if (results.getTwoQubitGates() == 0U) {
    return;
  }

  INFO() << "Trying to find a solution with less than "
         << results.getTwoQubitGates() << " two-qubit gates and at most "
         << gateCount << " gates.";

  config.targetMetric         = TargetMetric::TwoQubitGates;
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

  config.targetMetric         = TargetMetric::Gates;
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
  auto       encoder = encoding::SATEncoder(config);
  const auto res     = encoder.run();
  if (configuration.dumpIntermediateResults && res.sat()) {
    const auto filename = configuration.intermediateResultsPath +
                          "intermediate_" + std::to_string(solverCalls) +
                          ".qasm";
    INFO() << "Dumping circuit to " << filename;
    std::ofstream file(filename);
    file << res.getResultCircuit();
    file.close();
  }
  return res;
}

void CliffordSynthesizer::updateResults(const Configuration& config,
                                        const Results&       newResults,
                                        Results&             currentResults) {
  if (!newResults.sat()) {
    return;
  }

  switch (config.target) {
  case TargetMetric::Gates:
    if ((newResults.getGates() < currentResults.getGates()) ||
        ((newResults.getGates() == currentResults.getGates()) &&
         (newResults.getTwoQubitGates() < currentResults.getTwoQubitGates()))) {
      currentResults = newResults;
    }
    break;
  case TargetMetric::TwoQubitGates:
    if ((newResults.getTwoQubitGates() < currentResults.getTwoQubitGates()) ||
        ((newResults.getTwoQubitGates() == currentResults.getTwoQubitGates()) &&
         (newResults.getGates() < currentResults.getGates()))) {
      currentResults = newResults;
    }
    break;
  case TargetMetric::Depth:
    if ((newResults.getDepth() < currentResults.getDepth()) ||
        ((newResults.getDepth() == currentResults.getDepth()) &&
         (newResults.getGates() < currentResults.getGates()))) {
      currentResults = newResults;
    }
    break;
  }
}

void CliffordSynthesizer::depthHeuristicSynthesis(
    CliffordSynthesizer::EncoderConfig config) {
  auto optimalConfig                 = configuration;
  optimalConfig.heuristic            = false;
  optimalConfig.target               = TargetMetric::Depth;
  optimalConfig.initialTimestepLimit = configuration.split_size;

  qc::CircuitOptimizer::reorderOperations(initialCircuit.value());
  qc::QuantumComputation optCircuit{initialCircuit->getNqubits()};

  std::size_t nPartitions = initialCircuit->size() / configuration.split_size;

  for (std::size_t i = 0; i < nPartitions; ++i) {
    std::optional<Tableau> subTargetTableauOpt{};
    if (i == nPartitions - 1) {
      subTargetTableauOpt = Tableau{
          *initialCircuit, 0, std::numeric_limits<std::size_t>::max(), true};
    } else {
      subTargetTableauOpt =
          Tableau{*initialCircuit, 0, i + 1 * configuration.split_size, true};
    }

    const Tableau subInitTableau{*initialCircuit, 0,
                                 i * configuration.split_size, true};

    CliffordSynthesizer synth(subInitTableau, *subTargetTableauOpt);
    synth.synthesize(optimalConfig);
    const auto& subCircuit = synth.getResultCircuit();
    for (const auto& op : subCircuit) {
      optCircuit.emplace_back(op->clone());
    }
    results.setRuntime(results.getRuntime() + synth.results.getRuntime());
  }
}
} // namespace cs
