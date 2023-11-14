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
      if (singleQubitGates.paulis().size() >
          1) { // if no Paulis are present we can relax the exactlyone
               // constraint
        assertExactlyOne(gateVariables);

      } else {
        auto atLeastOne = LogicTerm(false);
        for (const auto& var : gateVariables) {
          atLeastOne = atLeastOne || var;
        }
        lb->assertFormula(atLeastOne);
      }
    }
  }
}

void encoding::MultiGateEncoder::assertGateConstraints() {
  DEBUG() << "Asserting gate constraints";
  xorHelpers = logicbase::LogicMatrix{T};
  for (std::size_t t = 0U; t < T; ++t) {
    TRACE() << "Asserting gate constraints at time " << t;
    if (!ignoreRChanges) {
      rChanges = tvars->r[t];
      splitXorR(tvars->r[t], t);
    }
    assertSingleQubitGateConstraints(t);
    assertTwoQubitGateConstraints(t);
    TRACE() << "Asserting r changes at time " << t;

    if (!ignoreRChanges) {
      lb->assertFormula(tvars->r[t + 1] == xorHelpers[t].back());
    }
  }
}

void MultiGateEncoder::assertRConstraints(const std::size_t pos,
                                          const std::size_t qubit) {
  for (const auto gate : singleQubitGates) {
    const auto& change =
        LogicTerm::ite(vars.gS[pos][singleQubitGates.gateToIndex(gate)][qubit],
                       tvars->singleQubitRChange(pos, qubit, gate),
                       LogicTerm(0, static_cast<std::int16_t>(S)));
    splitXorR(change, pos);
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

  return changes;
}

[[nodiscard]] logicbase::LogicTerm
MultiGateEncoder::createTwoQubitRConstraint(std::size_t pos, std::size_t ctrl,
                                            std::size_t trgt) {
  const auto& newRChanges = LogicTerm::ite(
      vars.gC[pos][ctrl][trgt], tvars->twoQubitRChange(pos, ctrl, trgt),
      LogicTerm(0, static_cast<std::int16_t>(S)));
  splitXorR(newRChanges, pos);
  return logicbase::LogicTerm(true);
}

void MultiGateEncoder::assertSingleQubitGateOrderConstraints(
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
  for (const auto gate : singleQubitGates) {
    if (gate == qc::OpType::None) {
      continue;
    }
    noSingleQubitGate =
        noSingleQubitGate && !gSNext[singleQubitGates.gateToIndex(gate)][qubit];
  }
  lb->assertFormula(LogicTerm::implies(
      gSNow[singleQubitGates.gateToIndex(qc::OpType::None)][qubit],
      noSingleQubitGate));
}

void MultiGateEncoder::assertTwoQubitGateOrderConstraints(
    const std::size_t pos, const std::size_t ctrl, const std::size_t trgt) {
  // nothing to assert at the end
  if (pos == T - 1U) {
    return;
  }

  // gate variables of the current and the next time step
  const auto& gSNow  = vars.gS[pos];
  const auto& gCNext = vars.gC[pos + 1];

  // two identical CNOTs may not be applied in a row because they would cancel.
  lb->assertFormula(
      LogicTerm::implies(vars.gC[pos][ctrl][trgt], !gCNext[ctrl][trgt]));
  lb->assertFormula(
      LogicTerm::implies(vars.gC[pos][trgt][ctrl], !gCNext[trgt][ctrl]));

  // no gate on both qubits => no CNOT on them in the next time step.
  // hadamards on both qubits => no CNOT on them in the next time step (CNOT can
  // be conjugated with Hadamards) No Combination of Paulis on both Qubits
  // before a CNOT These gates can be just pushed through to the other side
  auto       noneIndex     = singleQubitGates.gateToIndex(qc::OpType::None);
  const auto noGate        = gSNow[noneIndex][ctrl] && gSNow[noneIndex][trgt];
  const auto noFurtherCnot = !gCNext[ctrl][trgt] && !gCNext[trgt][ctrl];

  auto       hIndex = singleQubitGates.gateToIndex(qc::OpType::H);
  auto       xIndex = singleQubitGates.gateToIndex(qc::OpType::X);
  auto       zIndex = singleQubitGates.gateToIndex(qc::OpType::Z);
  auto       yIndex = singleQubitGates.gateToIndex(qc::OpType::Y);
  const auto hh     = gSNow[hIndex][ctrl] && gSNow[hIndex][trgt];
  const auto gateBeforeCtrl =
      gSNow[zIndex][ctrl] || gSNow[xIndex][ctrl] || gSNow[yIndex][ctrl];
  const auto gateBeforeTarget =
      gSNow[zIndex][trgt] || gSNow[xIndex][trgt] || gSNow[yIndex][ctrl];
  lb->assertFormula(LogicTerm::implies(
      noGate || hh || (gateBeforeCtrl && gateBeforeTarget), noFurtherCnot));
}

void MultiGateEncoder::splitXorR(const logicbase::LogicTerm& changes,
                                 std::size_t                 pos) {
  auto&             xorHelper = xorHelpers[pos];
  const std::string hName =
      "h_" + std::to_string(pos) + "_" + std::to_string(xorHelper.size());
  DEBUG() << "Creating helper variable for RChange XOR " << hName;
  const auto n = static_cast<std::int16_t>(S);
  xorHelper.emplace_back(lb->makeVariable(hName, CType::BITVECTOR, n));
  if (xorHelper.size() == 1) {
    lb->assertFormula(xorHelper.back() == changes);
  } else {
    lb->assertFormula(xorHelper.back() ==
                      (xorHelper[xorHelpers[pos].size() - 2] ^ changes));
  }
}

} // namespace cs::encoding
