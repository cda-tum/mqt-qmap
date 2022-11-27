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
  auto encoderConfig           = encoding::SATEncoder::Configuration();
  encoderConfig.initialTableau = &initialTableau;
  encoderConfig.targetTableau  = &targetTableau;
  encoderConfig.nQubits        = initialTableau.getQubitCount();
  encoderConfig.timestepLimit  = configuration.initialTimestepLimit;
  encoderConfig.targetMetric   = configuration.target;
  encoderConfig.useMaxSAT      = configuration.useMaxSAT;
  encoderConfig.nThreads       = configuration.nThreads;
  encoderConfig.useMultiGateEncoding =
      requiresMultiGateEncoding(encoderConfig.targetMetric);

  // determine an initial guess for the number of timesteps
  determineInitialTimestepLimit(encoderConfig);

  // determine an upper bound for the number of timesteps by solving the SAT
  // problem repeatedly with increasing timestep limits until a satisfying
  // assignment is found.
  const auto [lower, upper] = determineUpperBound(encoderConfig);

  encoderConfig.timestepLimit = upper;
  // run the optimal synthesis approach
  if (configuration.useMaxSAT) {
    runMaxSAT(encoderConfig);
  } else {
    runBinarySearch(encoderConfig.timestepLimit, lower, upper, encoderConfig);
  }

  if (configuration.target == TargetMetric::DEPTH &&
      configuration.minimizeGatesAfterDepthOptimization) {
    minimizeGatesFixedDepth(encoderConfig);
  }

  results.setSolverCalls(solverCalls);

  const auto end = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double> diff = end - start;
  INFO() << "Synthesis took " << diff.count() << " seconds";
  results.setRuntime(diff.count());
}

void CliffordSynthesizer::determineInitialTimestepLimit(
    encoding::SATEncoder::Configuration& config) {
  if (config.timestepLimit == 0U) {
    // Infeasibility checks are inexpensive, so we can afford to start with a
    // timestep limit of 1 and increase geometrically until we find a feasible
    // solution.
    config.timestepLimit = 1U;
  }
}

std::pair<std::size_t, std::size_t> CliffordSynthesizer::determineUpperBound(
    encoding::SATEncoder::Configuration config) {
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
  INFO() << "Found upper bound for the number of timesteps: " << upperBound;
  return {lowerBound, upperBound};
}

void CliffordSynthesizer::minimizeGatesFixedDepth(
    encoding::SATEncoder::Configuration config) {
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

void CliffordSynthesizer::runMaxSAT(
    const encoding::SATEncoder::Configuration& config) {
  INFO() << "Running MaxSAT scheme with timestep limit "
         << config.timestepLimit;
  results = callSolver(config);
}

Results CliffordSynthesizer::callSolver(
    const encoding::SATEncoder::Configuration& config) {
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
