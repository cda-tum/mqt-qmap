//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/encoding/ObjectiveEncoder.hpp"

#include "LogicTerm/LogicTerm.hpp"
#include "cliffordsynthesis/encoding/MultiGateEncoder.hpp"
#include "utils/logging.hpp"

namespace cs::encoding {

using namespace logicbase;

LogicTerm
ObjectiveEncoder::collectGateCount(const bool includeSingleQubitGates,
                                   const bool includeTwoQubitGates) const {
  auto cost = LogicTerm(0);
  for (std::size_t t = 0U; t < T; ++t) {
    if (includeSingleQubitGates) {
      collectSingleQubitGateTerms(t, cost, std::plus{});
    }
    if (includeTwoQubitGates) {
      collectTwoQubitGateTerms(t, cost, std::plus{});
    }
  }
  return cost;
}

void ObjectiveEncoder::optimizeGateCount(
    const bool includeSingleQubitGates, const bool includeTwoQubitGates) const {
  const auto cost =
      collectGateCount(includeSingleQubitGates, includeTwoQubitGates);
  dynamic_cast<LogicBlockOptimizer*>(lb.get())->minimize(cost);
}

LogicTerm
ObjectiveEncoder::collectDepth(const bool includeSingleQubitGates,
                               const bool includeTwoQubitGates) const {
  auto cost = LogicTerm(0);
  for (std::size_t t = 0U; t < T; ++t) {
    auto anyGate = LogicTerm(false);
    if (includeSingleQubitGates) {
      collectSingleQubitGateTerms(t, anyGate, std::logical_or{});
    }
    if (includeTwoQubitGates) {
      collectTwoQubitGateTerms(t, anyGate, std::logical_or{});
    }
    cost = cost + LogicTerm::ite(anyGate, LogicTerm(1), LogicTerm(0));
  }
  return cost;
}

void ObjectiveEncoder::optimizeDepth(bool includeSingleQubitGates,
                                     bool includeTwoQubitGates) const {
  const auto cost = collectDepth(includeSingleQubitGates, includeTwoQubitGates);
  dynamic_cast<LogicBlockOptimizer*>(lb.get())->minimize(cost);
}

void ObjectiveEncoder::optimizeMetric(TargetMetric targetMetric) const {
  switch (targetMetric) {
  case TargetMetric::GATES:
    optimizeGateCount();
    break;
  case TargetMetric::DEPTH:
    optimizeDepth();
    break;
  default:
    FATAL() << "Unknown target metric: " << toString(targetMetric);
  }
}
} // namespace cs::encoding
