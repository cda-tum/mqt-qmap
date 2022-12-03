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
ObjectiveEncoder::collectGateCount(const bool includeSingleQubitGates) const {
  auto cost = LogicTerm(0);
  for (std::size_t t = 0U; t < T; ++t) {
    if (includeSingleQubitGates) {
      collectSingleQubitGateTerms(t, cost, std::plus{});
    }
    collectTwoQubitGateTerms(t, cost, std::plus{});
  }
  return cost;
}

void ObjectiveEncoder::optimizeGateCount(
    const bool includeSingleQubitGates) const {
  DEBUG() << "Optimizing " << (includeSingleQubitGates ? "" : "two-qubit ")
          << "gate count";
  const auto cost = collectGateCount(includeSingleQubitGates);
  dynamic_cast<LogicBlockOptimizer*>(lb.get())->minimize(cost);
}

void ObjectiveEncoder::optimizeDepth(bool includeSingleQubitGates) const {
  DEBUG() << "Optimizing " << (includeSingleQubitGates ? "" : "two-qubit ")
          << "depth";
  auto* optimizer = dynamic_cast<LogicBlockOptimizer*>(lb.get());
  for (std::size_t t = 0U; t < T; ++t) {
    auto anyGate = LogicTerm(false);
    if (includeSingleQubitGates) {
      collectSingleQubitGateTerms(t, anyGate, std::logical_or{});
    }
    collectTwoQubitGateTerms(t, anyGate, std::logical_or{});
    optimizer->weightedTerm(anyGate, 1);
  }
  optimizer->makeMinimize();
}

void ObjectiveEncoder::optimizeMetric(TargetMetric targetMetric) const {
  switch (targetMetric) {
  case TargetMetric::GATES:
    optimizeGateCount();
    break;
  case TargetMetric::TWO_QUBIT_GATES:
    optimizeGateCount(false);
    break;
  case TargetMetric::DEPTH:
    optimizeDepth();
    break;
  default:
    FATAL() << "Unknown target metric: " << toString(targetMetric);
  }
}
} // namespace cs::encoding
