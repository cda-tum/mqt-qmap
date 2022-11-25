//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#include "LogicTerm/Logic.hpp"
#include "cliffordsynthesis/SATEncoder.hpp"
#include "utils/logging.hpp"

#include <chrono>

namespace cs {

void CliffordSynthesizer::synthesize(const Configuration& config) {
  configuration = config;

  INFO() << "Optimization target: " << toString(configuration.target);

  const auto start = std::chrono::high_resolution_clock::now();

  timestepLimit = configuration.initialTimestepLimit;
  if (timestepLimit == 0U) {
    // Infeasibility checks are inexpensive, so we can afford to start with a
    // timestep limit of 1 and increase geometrically until we find a feasible
    // solution.
    timestepLimit = 1U;

    // if a circuit was provided to start with, we can use its characteristics
    // as an upper bound for the number of timesteps.
    if (configuration.target == TargetMetric::GATES &&
        initialGates.has_value()) {
      timestepLimit = *initialGates;
    }
  }

  determineUpperBound();

  if (configuration.useMaxSAT) {
    runMaxSAT();
  } else {
    runBinarySearch(timestepLimit, lowerTimestepLimit, timestepLimit);
  }

  if (configuration.target == TargetMetric::DEPTH &&
      configuration.minimizeGatesAfterDepthOptimization) {
    minimizeGatesFixedDepth();
  }

  results.setSolverCalls(solverCalls);

  const auto end = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double> diff = end - start;
  INFO() << "Synthesis took " << diff.count() << " seconds";
  results.setRuntime(diff.count());
}

void CliffordSynthesizer::determineUpperBound() {
  INFO() << "Searching for upper bound for the number of timesteps starting "
         << "with " << timestepLimit;
  while (!results.sat()) {
    results = mainOptimization(false);
    if (!results.sat()) {
      lowerTimestepLimit = timestepLimit + 1U;
      INFO() << "No solution found for " << timestepLimit << " timestep(s). "
             << "Doubling timestep limit to " << 2 * timestepLimit;
      timestepLimit *= 2U;
    }
  }
  INFO() << "Found upper bound for the number of timesteps: " << timestepLimit;
}

void CliffordSynthesizer::minimizeGatesFixedDepth() {
  if (results.getDepth() == 0U) {
    return;
  }

  if (results.getDepth() == results.getGates()) {
    return;
  }

  INFO() << "Found a depth-optimal circuit with depth " << results.getDepth()
         << " and " << results.getGates()
         << " gate(s). Trying to minimize the number of gates.";
  configuration.target               = TargetMetric::GATES;
  timestepLimit                      = results.getDepth();
  configuration.useMultiGateEncoding = true;

  if (configuration.useMaxSAT) {
    runMaxSAT();
  } else {
    configuration.gateLimit = results.getGates();
    runBinarySearch(*configuration.gateLimit, timestepLimit,
                    results.getGates());
  }
  INFO() << "Found a depth " << results.getDepth() << " circuit with "
         << results.getGates() << " gate(s).";
}

void CliffordSynthesizer::runMaxSAT() {
  INFO() << "Running MaxSAT scheme with timestep limit " << timestepLimit;
  results = mainOptimization(true);
}

Results CliffordSynthesizer::mainOptimization(const bool useMaxSAT) {
  using namespace logicbase;

  ++solverCalls;
  const auto start = std::chrono::high_resolution_clock::now();

  auto encoder = encoding::SATEncoder(initialTableau.getQubitCount(),
                                      timestepLimit, useMaxSAT);
  encoder.createFormulation(initialTableau, targetTableau, configuration);
  encoder.produceInstance();
  const auto solverResult = encoder.solve();

  const auto end     = std::chrono::high_resolution_clock::now();
  const auto runtime = std::chrono::duration<double>(end - start);

  Results res{};
  res.setRuntime(runtime.count());
  res.setSolverResult(solverResult);

  if (solverResult == Result::SAT) {
    encoder.extractResultsFromModel(res);
  } else if (solverResult == Result::UNSAT) {
    encoder.cleanup();
    return res;
  } else {
    FATAL() << "Solver returned NDEF. Something went wrong.";
  }
  encoder.cleanup();

  return res;
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
