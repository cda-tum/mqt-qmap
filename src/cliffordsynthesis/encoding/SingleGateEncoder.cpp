//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/encoding/SingleGateEncoder.hpp"

#include "utils/logging.hpp"

namespace cs::encoding {

using namespace logicbase;

void SingleGateEncoder::assertConsistency() const {
  DEBUG() << "Asserting gate consistency";
  for (std::size_t t = 0U; t < T; ++t) {
    LogicVector gateVariables{};

    for (std::size_t q = 0U; q < N; ++q) {
      vars.collectSingleQubitGateVariables(t, q, gateVariables);
      vars.collectTwoQubitGateVariables(t, q, true, gateVariables);
    }
    IF_PLOG(plog::verbose) {
      TRACE() << "Gate variables at time " << t;
      for (const auto& var : gateVariables) {
        TRACE() << var.getName();
      }
    }
    assertExactlyOne(gateVariables);
  }
}

void SingleGateEncoder::assertGateConstraints() {
  DEBUG() << "Asserting gate constraints";
  for (std::size_t t = 0U; t < T; ++t) {
    TRACE() << "Asserting gate constraints at time " << t;
    assertSingleQubitGateConstraints(t);
    assertTwoQubitGateConstraints(t);
  }
}

void SingleGateEncoder::assertSingleQubitGateConstraints(std::size_t pos) {
  const auto& singleQubitGates = vars.gS[pos];
  for (std::size_t q = 0U; q < N; ++q) {
    for (const auto gate : SINGLE_QUBIT_GATES) {
      const auto changes = createSingleQubitGateConstraint(pos, q, gate);

      DEBUG() << "Asserting " << toString(gate) << " on " << q;

      lb->assertFormula(
          LogicTerm::implies(singleQubitGates[gateToIndex(gate)][q], changes));
    }
  }
}

LogicTerm SingleGateEncoder::createSingleQubitGateConstraint(
    const std::size_t pos, const std::size_t qubit, const qc::OpType gate) {
  auto changes = LogicTerm(true);

  changes = changes && (tvars->x[pos + 1][qubit] ==
                        tvars->singleQubitXChange(pos, qubit, gate));
  changes = changes && (tvars->z[pos + 1][qubit] ==
                        tvars->singleQubitZChange(pos, qubit, gate));
  changes = changes &&
            (tvars->r[pos + 1] ==
             (tvars->r[pos] ^ tvars->singleQubitRChange(pos, qubit, gate)));

  return changes && createNoChange(pos, qubit, std::nullopt);
}

void SingleGateEncoder::assertTwoQubitGateConstraints(const std::size_t pos) {
  const auto& twoQubitGates = vars.gC[pos];
  for (std::size_t ctrl = 0U; ctrl < N; ++ctrl) {
    for (std::size_t trgt = 0U; trgt < N; ++trgt) {
      if (ctrl == trgt) {
        continue;
      }
      const auto changes = createTwoQubitGateConstraint(pos, ctrl, trgt);

      DEBUG() << "Asserting CNOT on " << ctrl << " and " << trgt;

      lb->assertFormula(LogicTerm::implies(twoQubitGates[ctrl][trgt], changes));
    }
  }
}

LogicTerm SingleGateEncoder::createTwoQubitGateConstraint(
    const std::size_t pos, const std::size_t ctrl, const std::size_t trgt) {
  auto changes              = LogicTerm(true);
  const auto [xCtrl, xTrgt] = tvars->twoQubitXChange(pos, ctrl, trgt);
  const auto [zCtrl, zTrgt] = tvars->twoQubitZChange(pos, ctrl, trgt);

  changes = changes && (tvars->x[pos + 1][ctrl] == xCtrl);
  changes = changes && (tvars->x[pos + 1][trgt] == xTrgt);
  changes = changes && (tvars->z[pos + 1][ctrl] == zCtrl);
  changes = changes && (tvars->z[pos + 1][trgt] == zTrgt);
  changes =
      changes && (tvars->r[pos + 1] ==
                  (tvars->r[pos] ^ tvars->twoQubitRChange(pos, ctrl, trgt)));

  return changes && createNoChange(pos, ctrl, trgt);
}

LogicTerm SingleGateEncoder::createNoChange(
    const std::size_t pos, const std::size_t except,
    const std::optional<std::size_t> except2) const {
  auto changes = LogicTerm(true);
  for (std::size_t q = 0U; q < N; ++q) {
    if (q == except) {
      continue;
    }
    if (except2.has_value() && q == except2.value()) {
      continue;
    }

    changes = changes && (tvars->x[pos + 1][q] == tvars->x[pos][q]);
    changes = changes && (tvars->z[pos + 1][q] == tvars->z[pos][q]);
  }
  return changes;
}

void SingleGateEncoder::assertSingleQubitGateOrderConstraints(
    const std::size_t pos, const std::size_t qubit) {
  // nothing to assert at the end
  if (pos == T - 1U) {
    return;
  }

  // gate variables of the current and the next time step
  const auto& gSNow  = vars.gS[pos];
  const auto& gSNext = vars.gS[pos + 1];

  // collect variables of single-qubit gates that could be applied to `qubit`
  auto singleQubitGate = LogicTerm(false);
  for (const auto& gate : SINGLE_QUBIT_GATES) {
    if (gate == qc::OpType::None) {
      continue;
    }
    // any single-qubit gate on qubit q
    singleQubitGate = singleQubitGate || gSNow[gateToIndex(gate)][qubit];
  }

  // collect gate variables of the next timestep that should not be applied.
  auto disallowed = LogicTerm(true);

  // no single-qubit gate on a lower qubit
  for (std::size_t lower = 0U; lower < qubit; ++lower) {
    for (const auto& gate : SINGLE_QUBIT_GATES) {
      if (gate == qc::OpType::None) {
        continue;
      }
      disallowed = disallowed && !gSNext[gateToIndex(gate)][lower];
    }
  }

  lb->assertFormula(LogicTerm::implies(singleQubitGate, disallowed));

  // once no gate is applied, no other gate can be applied, i.e., in the next
  // timestep any of the `None` gate variables must be selected.
  const auto noneIndex = gateToIndex(qc::OpType::None);
  auto       noGate    = LogicTerm(false);
  for (std::size_t q = 0U; q < N; ++q) {
    noGate = noGate || gSNext[noneIndex][q];
  }
  lb->assertFormula(LogicTerm::implies(gSNow[noneIndex][qubit], noGate));
}

void SingleGateEncoder::assertTwoQubitGateOrderConstraints(
    const std::size_t pos, const std::size_t ctrl, const std::size_t trgt) {
  // nothing to assert at the end
  if (pos == T - 1U) {
    return;
  }

  // gate variables of the current and the next time step
  const auto& current = vars.gC[pos][ctrl][trgt];
  const auto& gCNext  = vars.gC[pos + 1];
  const auto& gSNext  = vars.gS[pos + 1];

  // two identical CNOTs may not be applied in a row because they would cancel.
  auto disallowed = !gCNext[ctrl][trgt];

  // any single-qubit gate independent of the CNOT should be applied prior.
  for (std::size_t q = 0U; q < N; ++q) {
    if (q == ctrl || q == trgt) {
      continue;
    }
    for (const auto& gate : SINGLE_QUBIT_GATES) {
      if (gate == qc::OpType::None) {
        continue;
      }
      disallowed = disallowed && !gSNext[gateToIndex(gate)][q];
    }
  }

  // no X gate may be placed on the target qubit since it would commute.
  disallowed = disallowed && !gSNext[gateToIndex(qc::OpType::X)][trgt];

  // no diagonal gate may be placed on the control qubit since it would commute.
  disallowed = disallowed && !gSNext[gateToIndex(qc::OpType::Z)][ctrl] &&
               !gSNext[gateToIndex(qc::OpType::S)][ctrl] &&
               !gSNext[gateToIndex(qc::OpType::Sdag)][ctrl];

  // no CNOT with the same control and a lower target qubit may be placed.
  for (std::size_t t = 0U; t < trgt; ++t) {
    if (t == ctrl) {
      continue;
    }
    disallowed = disallowed && !gCNext[ctrl][t];
  }

  // no CNOT with a lower control different from trgt may be placed.
  for (std::size_t c = 0U; c < ctrl; ++c) {
    for (std::size_t t = 0U; t < N; ++t) {
      if (c == t || c == trgt) {
        continue;
      }
      disallowed = disallowed && !gCNext[c][t];
    }
  }

  lb->assertFormula(LogicTerm::implies(current, disallowed));
}

} // namespace cs::encoding
