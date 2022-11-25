//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include <iostream>
#include <string>

namespace cs {
enum class TargetMetric { GATES, DEPTH };

static std::string toString(const TargetMetric target) {
  switch (target) {
  case TargetMetric::GATES:
    return "gates";
  case TargetMetric::DEPTH:
    return "depth";
  }
  return "Error";
}

[[maybe_unused]] static TargetMetric
targetMetricFromString(const std::string& target) {
  if (target == "gates") {
    return TargetMetric::GATES;
  }
  if (target == "depth") {
    return TargetMetric::DEPTH;
  }
  return TargetMetric::GATES;
}
} // namespace cs
