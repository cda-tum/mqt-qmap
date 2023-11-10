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

void ObjectiveEncoder::optimizeDepth() const {
  DEBUG() << "Optimizing depth";
  auto* optimizer = dynamic_cast<LogicBlockOptimizer*>(lb.get());

  auto noGateIndex = singleQubitGates.gateToIndex(qc::OpType::None);
  for (std::size_t t = 0U; t < T; ++t) {
    const auto& gS     = gvars->gS[t];
    auto        noGate = LogicTerm(true);
    for (std::size_t q = 0U; q < N; ++q) {
      noGate = noGate && gS[noGateIndex][q];
    }
    optimizer->weightedTerm(!noGate, 1);
  }
  optimizer->makeMinimize();
}

void ObjectiveEncoder::optimizeMetric(TargetMetric targetMetric) const {
  switch (targetMetric) {
  case TargetMetric::Gates:
    optimizeGateCount();
    break;
  case TargetMetric::TwoQubitGates:
    optimizeGateCount(false);
    break;
  case TargetMetric::Depth:
    optimizeDepth();
    break;
  default:
    FATAL() << "Unknown target metric: " << toString(targetMetric);
  }
}
} // namespace cs::encoding
