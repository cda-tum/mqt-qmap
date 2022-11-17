/*
 * This file is part of the MQT QMAP library which is released under the MIT
 * license. See file README.md or go to
 * https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
 */

#ifndef CS_TARGETMETRIC_HPP
#define CS_TARGETMETRIC_HPP
#include <iostream>

namespace cs {
enum class TargetMetric { GATES, TWO_QUBIT_GATES, DEPTH, FIDELITY };

static std::string toString(const TargetMetric target) {
  switch (target) {
  case TargetMetric::GATES:
    return "gates";
  case TargetMetric::TWO_QUBIT_GATES:
    return "two_qubit_gates";
  case TargetMetric::DEPTH:
    return "depth";
  case TargetMetric::FIDELITY:
    return "fidelity";
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
  if (target == "fidelity") {
    return TargetMetric::FIDELITY;
  }
  return TargetMetric::GATES;
}
} // namespace cs
#endif // CS_TARGETMETRIC_HPP
