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

} // namespace cs::encoding
