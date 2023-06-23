//
// Created by Velsh Aleksei on 16.06.23.
//

#include "cliffordsynthesis/encoding/STQGatesEncoder.hpp"

#include "LogicTerm/LogicTerm.hpp"
#include "utils/logging.hpp"

namespace cs::encoding {

using namespace logicbase;

void encoding::STQGatesEncoder::assertConsistency() const {
  DEBUG() << "Asserting gate consistency";
  // T/2
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

void encoding::STQGatesEncoder::assertGateConstraints() {
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

void encoding::STQGatesEncoder::assertSingleQubitGateConstraints(
    const std::size_t pos) {
  for (std::size_t q = 0U; q < N; ++q) {
    assertZConstraints(pos, q);
    assertXConstraints(pos, q);
    assertRConstraints(pos, q);
  }
}

void STQGatesEncoder::assertRConstraints(const std::size_t pos,
                                         const std::size_t qubit) {
  for (const auto gate : SINGLE_QUBIT_GATES) {
    rChanges =
        rChanges ^ LogicTerm::ite(vars.gS[pos][gateToIndex(gate)][qubit],
                                  tvars->singleQubitRChange(pos, qubit, gate),
                                  LogicTerm(0, static_cast<std::int16_t>(S)));
  }
}

void encoding::STQGatesEncoder::assertTwoQubitGateConstraints(
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

LogicTerm encoding::STQGatesEncoder::createTwoQubitGateConstraint(
    std::size_t pos, std::size_t ctrl, std::size_t trgt) {
  auto changes              = LogicTerm(true);
  const auto [xCtrl, xTrgt] = tvars->twoQubitXChange(pos, ctrl, trgt);
  const auto [zCtrl, zTrgt] = tvars->twoQubitZChange(pos, ctrl, trgt);

  changes = changes && (tvars->x[pos + 1][ctrl] == xCtrl);
  changes = changes && (tvars->x[pos + 1][trgt] == xTrgt);
  changes = changes && (tvars->z[pos + 1][ctrl] == zCtrl);
  changes = changes && (tvars->z[pos + 1][trgt] == zTrgt);

  rChanges =
      rChanges ^ LogicTerm::ite(vars.gC[pos][ctrl][trgt],
                                tvars->twoQubitRChange(pos, ctrl, trgt),
                                LogicTerm(0, static_cast<std::int16_t>(S)));

  return changes;
}

void STQGatesEncoder::assertSingleQubitGateOrderConstraints(
    const std::size_t pos, const std::size_t qubit) {
  // nothing to assert at the end
  if (pos == T - 1U) {
    return;
  }

  // gate variables of the current and the next time step
  const auto& gSNow  = vars.gS[pos];
  const auto& gSNext = vars.gS[pos + 1];

  // once no gate is applied to the considered qubit, no single qubit gate can
  // be applied to it in the next time step.
  auto noSingleQubitGate = LogicTerm(true);
  for (const auto gate : SINGLE_QUBIT_GATES) {
    if (gate == qc::OpType::None) {
      continue;
    }
    noSingleQubitGate = noSingleQubitGate && !gSNext[gateToIndex(gate)][qubit];
  }
  lb->assertFormula(LogicTerm::implies(
      gSNow[gateToIndex(qc::OpType::None)][qubit], noSingleQubitGate));
}

void STQGatesEncoder::assertTwoQubitGateOrderConstraints(
    const std::size_t pos, const std::size_t ctrl, const std::size_t trgt) {
  // nothing to assert at the end
  if (pos == T - 1U) {
    return;
  }

  // gate variables of the current and the next time step
  const auto& current = vars.gC[pos][ctrl][trgt];
  const auto& gSNow   = vars.gS[pos];
  const auto& gCNext  = vars.gC[pos + 1];

  // two identical CNOTs may not be applied in a row because they would cancel.
  lb->assertFormula(LogicTerm::implies(current, !gCNext[ctrl][trgt]));

  // if no gate is applied to both qubits, no CNOT on them can be applied in the
  // next time step.
  const auto noneIndex     = gateToIndex(qc::OpType::None);
  const auto noGate        = gSNow[noneIndex][ctrl] && gSNow[noneIndex][trgt];
  const auto noFurtherCnot = !gCNext[ctrl][trgt] && !gCNext[trgt][ctrl];
  lb->assertFormula(LogicTerm::implies(noGate, noFurtherCnot));
}

} // namespace cs::encoding
