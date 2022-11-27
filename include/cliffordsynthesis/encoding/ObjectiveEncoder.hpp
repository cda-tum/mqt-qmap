//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "cliffordsynthesis/TargetMetric.hpp"
#include "cliffordsynthesis/encoding/GateEncoder.hpp"
#include "cliffordsynthesis/encoding/TableauEncoder.hpp"

#include <cstddef>
#include <functional>
#include <memory>

namespace cs::encoding {

class ObjectiveEncoder {
public:
  ObjectiveEncoder(const std::size_t nQubits, const std::size_t timestepLimit,
                   TableauEncoder::Variables*                    tvars,
                   GateEncoder::Variables*                       gvars,
                   const std::shared_ptr<logicbase::LogicBlock>& lb)
      : N(nQubits), T(timestepLimit), tvars(tvars), gvars(gvars), lb(lb) {}

  template <class Op>
  void limitGateCount(std::size_t maxGateCount, Op op,
                      bool includeSingleQubitGates = true,
                      bool includeTwoQubitGates    = true) const {
    const auto cost =
        collectGateCount(includeSingleQubitGates, includeTwoQubitGates);
    lb->assertFormula(
        op(cost, logicbase::LogicTerm(static_cast<int>(maxGateCount))));
  }

  void optimizeMetric(TargetMetric targetMetric) const;

  void optimizeGateCount(bool includeSingleQubitGates = true,
                         bool includeTwoQubitGates    = true) const;

  void optimizeDepth(bool includeSingleQubitGates = true,
                     bool includeTwoQubitGates    = true) const;

protected:
  // number of qubits N
  std::size_t N{};
  // timestep limit T
  std::size_t T{};

  // the tableau variables
  TableauEncoder::Variables* tvars;

  // the gate variables
  GateEncoder::Variables* gvars;

  // the logic block
  std::shared_ptr<logicbase::LogicBlock> lb;

  [[nodiscard]] logicbase::LogicTerm
  collectGateCount(bool includeSingleQubitGates = true,
                   bool includeTwoQubitGates    = true) const;

  [[nodiscard]] logicbase::LogicTerm
  collectDepth(bool includeSingleQubitGates = true,
               bool includeTwoQubitGates    = true) const;

  template <class Op>
  void collectGateTerms(std::size_t pos, logicbase::LogicTerm& terms,
                        Op op) const {
    collectSingleQubitGateTerms(pos, terms, op);
    collectTwoQubitGateTerms(pos, terms, op);
  }

  template <class Op>
  void collectSingleQubitGateTerms(std::size_t pos, logicbase::LogicTerm& terms,
                                   Op op) const {
    const auto& singleQubitGates = gvars->gS[pos];
    for (std::size_t q = 0U; q < N; ++q) {
      for (const auto gate : GateEncoder::SINGLE_QUBIT_GATES) {
        if (gate == qc::OpType::None) {
          continue;
        }
        terms = op(terms, singleQubitGates[GateEncoder::gateToIndex(gate)][q]);
      }
    }
  }

  template <class Op>
  void collectTwoQubitGateTerms(std::size_t pos, logicbase::LogicTerm& terms,
                                Op op) const {
    const auto& twoQubitGates = gvars->gC[pos];
    for (std::size_t ctrl = 0U; ctrl < N; ++ctrl) {
      for (std::size_t trgt = 0U; trgt < ctrl; ++trgt) {
        terms = op(terms, twoQubitGates[ctrl][trgt]);
        terms = op(terms, twoQubitGates[trgt][ctrl]);
      }
    }
  }
};

} // namespace cs::encoding
