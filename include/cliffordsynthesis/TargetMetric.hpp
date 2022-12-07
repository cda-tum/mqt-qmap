//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include <iostream>
#include <string>

namespace cs {
enum class TargetMetric { GATES, TWO_QUBIT_GATES, DEPTH };

static std::string toString(const TargetMetric target) {
  switch (target) {
  case TargetMetric::GATES:
    return "gates";
  case TargetMetric::TWO_QUBIT_GATES:
    return "two_qubit_gates";
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
  if (target == "two_qubit_gates") {
    return TargetMetric::TWO_QUBIT_GATES;
  }
  if (target == "depth") {
    return TargetMetric::DEPTH;
  }
  return TargetMetric::GATES;
}
} // namespace cs
