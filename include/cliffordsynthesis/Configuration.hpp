//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "TargetMetric.hpp"
#include "nlohmann/json.hpp"

namespace cs {
struct Configuration {
  Configuration() = default;

  /// General configuration for the synthesis algorithm
  std::size_t  initialTimestepLimit = 0U;
  bool         useMaxSAT            = false;
  TargetMetric target               = TargetMetric::GATES;

  /// Settings for the SAT solver
  std::size_t nThreads = 1U;

  /// Settings for depth-optimal synthesis
  bool minimizeGatesAfterDepthOptimization = false;

  /// Settings for two-qubit gate-optimal synthesis
  bool   tryHigherGateLimitForTwoQubitGateOptimization = false;
  double gateLimitFactor                               = 1.1;
  bool   minimizeGatesAfterTwoQubitGateOptimization    = false;

  [[nodiscard]] nlohmann::json json() const {
    nlohmann::json j;
    j["initial_timestep_limit"] = initialTimestepLimit;
    j["use_max_sat"]            = useMaxSAT;
    j["target_metric"]          = toString(target);
    j["n_threads"]              = nThreads;
    j["minimize_gates_after_depth_optimization"] =
        minimizeGatesAfterDepthOptimization;
    j["try_higher_gate_limit_for_two_qubit_gate_optimization"] =
        tryHigherGateLimitForTwoQubitGateOptimization;
    j["gate_limit_factor"] = gateLimitFactor;
    j["minimize_gates_after_two_qubit_gate_optimization"] =
        minimizeGatesAfterTwoQubitGateOptimization;

    return j;
  }

  friend std::ostream& operator<<(std::ostream&        os,
                                  const Configuration& config) {
    os << config.json().dump(2);
    return os;
  }
};
} // namespace cs
