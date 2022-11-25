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

  std::size_t  initialTimestepLimit = 0U;
  bool         useMaxSAT            = false;
  TargetMetric target               = TargetMetric::GATES;

  std::size_t nThreads = 1U;

  [[nodiscard]] nlohmann::json json() const {
    nlohmann::json j;
    j["initial_timestep_limit"] = initialTimestepLimit;
    j["use_max_sat"]            = useMaxSAT;
    j["target_metric"]          = toString(target);
    j["n_threads"]              = nThreads;

    return j;
  }

  friend std::ostream& operator<<(std::ostream&        os,
                                  const Configuration& config) {
    os << config.json().dump(2);
    return os;
  }
};
} // namespace cs
