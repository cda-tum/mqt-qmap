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

  if (configuration.useMaxSAT) {
    runMaxSAT();
  } else {
    runBinarySearch();
  }

  const auto end = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double> diff = end - start;
  INFO() << "Synthesis took " << diff.count() << " seconds";
  results.setRuntime(diff.count());
}

void CliffordSynthesizer::runMaxSAT() {
  bool found = false;
  INFO() << "Running MaxSAT scheme with timestep limit " << timestepLimit;
  while (!found) {
    results = mainOptimization();
    if (results.sat()) {
      found = true;
    } else {
      timestepLimit *= 2;
      INFO() << "No solution found. Doubling timestep limit. New limit: "
             << timestepLimit;
    }
  }
  results.setSolverCalls(solverCalls);
}

void CliffordSynthesizer::runBinarySearch() {
  INFO() << "Running binary search scheme with timestep limit "
         << timestepLimit;
  std::size_t upper = timestepLimit;
  std::size_t lower = 0U;

  INFO() << "Searching for upper bound.";
  while (!results.sat()) {
    results = mainOptimization();
    if (!results.sat()) {
      lower = upper + 1;
      upper *= 2;
      timestepLimit = upper;
      INFO() << "No solution found. Doubling timestep limit. New limit: "
             << timestepLimit;
    } else {
      INFO() << "Found solution with timestep limit " << timestepLimit << ". "
             << "Searching for lower bound.";
    }
  }

  while (lower != upper) {
    timestepLimit = (lower + upper) / 2;
    INFO() << "Searching for solution with timestep limit " << timestepLimit
           << " (lower: " << lower << ", upper: " << upper << ")";
    const auto r = mainOptimization();
    updateResults(configuration, r, results);
    if (r.sat()) {
      upper = timestepLimit;
      INFO() << "Found solution. Adjusting upper bound to " << upper;
    } else if (r.unsat()) {
      lower = timestepLimit + 1U;
      INFO() << "No solution found. Adjusting lower bound to " << lower;
    } else {
      FATAL() << "Unexpected result from main optimization";
    }
  }
  INFO() << "Binary search finished. Found solution with timestep limit "
         << lower;
  results.setSolverCalls(solverCalls);
}

Results CliffordSynthesizer::mainOptimization() {
  using namespace logicbase;

  ++solverCalls;
  const auto start = std::chrono::high_resolution_clock::now();

  auto encoder =
      encoding::SATEncoder(initialTableau.getQubitCount(), timestepLimit);
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
