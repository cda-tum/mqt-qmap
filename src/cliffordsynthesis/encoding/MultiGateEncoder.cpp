//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/encoding/MultiGateEncoder.hpp"

#include "LogicTerm/LogicTerm.hpp"
#include "utils/logging.hpp"

namespace cs::encoding {

using namespace logicbase;

void encoding::MultiGateEncoder::assertConsistency() const {
  DEBUG() << "Asserting gate consistency";
  for (std::size_t t = 0U; t < T; ++t) {
    // asserting only a single gate is applied on each qubit.
    for (std::size_t q = 0U; q < N; ++q) {
      LogicVector gateVariables{};
      vars.collectSingleQubitGateVariables(t, q, gateVariables);
      vars.collectTwoQubitGateVariables(t, q, true, gateVariables);
      vars.collectTwoQubitGateVariables(t, q, false, gateVariables);

      IF_PLOG(plog::verbose) {
        TRACE() << "Gate variables at time " << t << " and qubit " << q;
        for (const auto& var : gateVariables) {
          TRACE() << var.getName();
        }
      }

      assertExactlyOne(gateVariables);
    }
  }
}

void encoding::MultiGateEncoder::assertGateConstraints() {
  DEBUG() << "Asserting gate constraints";
  for (std::size_t t = 0U; t < T; ++t) {
    TRACE() << "Asserting gate constraints at time " << t;
    rChanges = tvars->r[t];
    assertSingleQubitGateConstraints(t);
    assertTwoQubitGateConstraints(t);
    TRACE() << "Asserting r changes at time " << t;
    lb->assertFormula(tvars->r[t + 1] == rChanges);
  }
}

void encoding::MultiGateEncoder::assertSingleQubitGateConstraints(
    const std::size_t pos) {
  const auto& singleQubitGates = vars.gS[pos];
  for (std::size_t q = 0U; q < N; ++q) {
    for (const auto gate : SINGLE_QUBIT_GATES) {
      const auto changes = createSingleQubitGateConstraint(pos, q, gate);
      lb->assertFormula(
          LogicTerm::implies(singleQubitGates[gateToIndex(gate)][q], changes));

      DEBUG() << "Asserting " << toString(gate) << " on " << q;
    }
  }
}

LogicTerm encoding::MultiGateEncoder::createSingleQubitGateConstraint(
    const std::size_t pos, const std::size_t qubit, const qc::OpType gate) {
  auto changes = LogicTerm(true);

  changes = changes && (tvars->x[pos + 1][qubit] ==
                        tvars->singleQubitXChange(pos, qubit, gate));
  changes = changes && (tvars->z[pos + 1][qubit] ==
                        tvars->singleQubitZChange(pos, qubit, gate));
  rChanges =
      rChanges ^ LogicTerm::ite(vars.gS[pos][gateToIndex(gate)][qubit],
                                tvars->singleQubitRChange(pos, qubit, gate),
                                LogicTerm(0, N));

  return changes;
}

void encoding::MultiGateEncoder::assertTwoQubitGateConstraints(
    const std::size_t pos) {
  const auto& twoQubitGates = vars.gC[pos];
  for (std::size_t ctrl = 0U; ctrl < N; ++ctrl) {
    for (std::size_t trgt = 0U; trgt < N; ++trgt) {
      if (ctrl == trgt) {
        continue;
      }
      const auto changes = createTwoQubitGateConstraint(pos, ctrl, trgt);
      lb->assertFormula(LogicTerm::implies(twoQubitGates[ctrl][trgt], changes));

      DEBUG() << "Asserting CNOT on " << ctrl << " and " << trgt;
    }
  }
}

LogicTerm encoding::MultiGateEncoder::createTwoQubitGateConstraint(
    std::size_t pos, std::size_t ctrl, std::size_t trgt) {
  auto changes              = LogicTerm(true);
  const auto [xCtrl, xTrgt] = tvars->twoQubitXChange(pos, ctrl, trgt);
  const auto [zCtrl, zTrgt] = tvars->twoQubitZChange(pos, ctrl, trgt);

  changes = changes && (tvars->x[pos + 1][ctrl] == xCtrl);
  changes = changes && (tvars->x[pos + 1][trgt] == xTrgt);
  changes = changes && (tvars->z[pos + 1][ctrl] == zCtrl);
  changes = changes && (tvars->z[pos + 1][trgt] == zTrgt);

  rChanges = rChanges ^ LogicTerm::ite(vars.gC[pos][ctrl][trgt],
                                       tvars->twoQubitRChange(pos, ctrl, trgt),
                                       LogicTerm(0, N));

  return changes;
}
} // namespace cs::encoding
