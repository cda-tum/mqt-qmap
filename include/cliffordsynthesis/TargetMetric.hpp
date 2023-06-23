//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include <iostream>
#include <string>

namespace cs {
enum class TargetMetric { Gates, TwoQubitGates, Depth, STDepth };

[[maybe_unused]] static inline std::string toString(const TargetMetric target) {
  switch (target) {
  case TargetMetric::Gates:
    return "gates";
  case TargetMetric::TwoQubitGates:
    return "two_qubit_gates";
  case TargetMetric::Depth:
    return "depth";
  case TargetMetric::STDepth:
    return "sTDepth";
  }
  return "Error";
}

[[maybe_unused]] static TargetMetric
targetMetricFromString(const std::string& target) {
  if (target == "gates") {
    return TargetMetric::Gates;
  }
  if (target == "two_qubit_gates") {
    return TargetMetric::TwoQubitGates;
  }
  if (target == "depth") {
    return TargetMetric::Depth;
  }
  if (target == "sTDepth") {
    return TargetMetric::STDepth;
  }
  return TargetMetric::Gates;
}
} // namespace cs
