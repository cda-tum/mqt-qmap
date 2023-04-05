//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "TargetMetric.hpp"
#include "nlohmann/json.hpp"

#include <plog/Log.h>
#include <thread>

namespace cs {
struct Configuration {
  Configuration() = default;

  /// General configuration for the synthesis algorithm
  std::size_t    initialTimestepLimit    = 0U;
  bool           useMaxSAT               = false;
  TargetMetric   target                  = TargetMetric::Gates;
  bool           useSymmetryBreaking     = true;
  bool           dumpIntermediateResults = false;
  std::string    intermediateResultsPath = "./";
  plog::Severity verbosity               = plog::Severity::warning;

  /// Settings for the SAT solver
  std::size_t nThreads = 1U;

  /// Settings for depth-optimal synthesis
  bool minimizeGatesAfterDepthOptimization = false;

  /// Settings for two-qubit gate-optimal synthesis
  bool   tryHigherGateLimitForTwoQubitGateOptimization = false;
  double gateLimitFactor                               = 1.1;
  bool   minimizeGatesAfterTwoQubitGateOptimization    = false;

  // Settings for the heuristic solver
  bool        heuristic         = false;
  std::size_t splitSize         = 5U;
  std::size_t nThreadsHeuristic = std::thread::hardware_concurrency();

  [[nodiscard]] nlohmann::json json() const {
    nlohmann::json j;
    j["initial_timestep_limit"] = initialTimestepLimit;
    j["use_max_sat"]            = useMaxSAT;
    j["target_metric"]          = toString(target);
    j["use_symmetry_breaking"]  = useSymmetryBreaking;
    j["n_threads"]              = nThreads;
    j["minimize_gates_after_depth_optimization"] =
        minimizeGatesAfterDepthOptimization;
    j["try_higher_gate_limit_for_two_qubit_gate_optimization"] =
        tryHigherGateLimitForTwoQubitGateOptimization;
    j["gate_limit_factor"] = gateLimitFactor;
    j["minimize_gates_after_two_qubit_gate_optimization"] =
        minimizeGatesAfterTwoQubitGateOptimization;
    j["heuristic"]           = heuristic;
    j["split_size"]          = splitSize;
    j["n_threads_heuristic"] = nThreadsHeuristic;
    return j;
  }

  friend std::ostream& operator<<(std::ostream&        os,
                                  const Configuration& config) {
    os << config.json().dump(2);
    return os;
  }
};
} // namespace cs
