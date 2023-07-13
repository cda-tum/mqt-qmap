//
// Created by Velsh Aleksei on 16.06.23.
//

#include "cliffordsynthesis/encoding/STQGatesEncoder.hpp"
#include "LogicTerm/LogicTerm.hpp"
#include "utils/logging.hpp"

namespace cs::encoding {

using namespace logicbase;

void encoding::STQGatesEncoder::collectTwoQubitGateVariables(
    const std::size_t pos, const std::size_t qubit, const bool target,
    LogicVector& variables) const {
  const auto& twoQubitGates = vars.gC[pos];
  const auto  n             = twoQubitGates.size();
  for (std::size_t q = 0; q < n; ++q) {
    if (q == qubit) {
      if (!target) {
        variables.emplace_back(twoQubitGates[qubit].back());
      }
      continue;
    }
    if (target) {
      variables.emplace_back(twoQubitGates[q][qubit]);
    } else {
      variables.emplace_back(twoQubitGates[qubit][q]);
    }
  }
}

void encoding::STQGatesEncoder::assertConsistency() const {
  DEBUG() << "Asserting gate consistency";
  for (std::size_t t = 0U; t < T/2U; ++t) {
    // asserting only a single gate is applied on each qubit.
    for (std::size_t q = 0U; q < N; ++q) {
      LogicVector singleQubitGateVariables{};
      LogicVector twoQubitGateVariables{};
      vars.collectSingleQubitGateVariables(t, q, singleQubitGateVariables);
      //TODO: may need to return the vars back
      collectTwoQubitGateVariables(t, q, true, twoQubitGateVariables);
      collectTwoQubitGateVariables(t, q, false, twoQubitGateVariables);

      IF_PLOG(plog::verbose) {
        TRACE() << "Single Qubit Gate variables at time " << t << " and qubit " << q;
        for (const auto& var : singleQubitGateVariables) {
          TRACE() << var.getName();
        }

        TRACE() << "Two Qubit Gate variables at time " << t << " and qubit " << q;
        for (const auto& var : twoQubitGateVariables) {
          TRACE() << var.getName();
        }
      }

      assertExactlyOne(singleQubitGateVariables);
      assertExactlyOne(twoQubitGateVariables);
    }
  }
}

void encoding::STQGatesEncoder::assertGateConstraints() {
  DEBUG() << "Asserting gate constraints";
  for (std::size_t t = 0U; t < T; ++t) {
    TRACE() << "Asserting gate constraints at time " << t;

    const std::size_t pos = t < L  ? t : t - L;
    t % 2U == 0
        ? assertSingleQubitGateConstraints(pos)
        : assertTwoQubitGateConstraints(pos);

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
        const auto changes = createIdentityConstraintOnTQG(pos, ctrl, trgt);
        lb->assertFormula(LogicTerm::implies(twoQubitGates[ctrl].back(), changes));
        DEBUG() << "Asserting Identity gate on " << ctrl << " and " << trgt;
      } else {
        const auto changes = createTwoQubitGateConstraint(pos, ctrl, trgt);
        lb->assertFormula(LogicTerm::implies(twoQubitGates[ctrl][trgt], changes));
        DEBUG() << "Asserting CNOT on " << ctrl << " and " << trgt;
      }
    }
  }
}
//TODO: to implement
LogicTerm encoding::STQGatesEncoder::createIdentityConstraintOnTQG(
    std::size_t pos, std::size_t ctrl, std::size_t trgt) {
  auto changes              = LogicTerm(true);

  return changes;
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

void STQGatesEncoder::extractCircuitFromModel(Results& res, Model& model) {
  std::size_t nSingleQubitGates = 0U;
  std::size_t nTwoQubitGates    = 0U;

  qc::QuantumComputation qc(N);
  for (std::size_t t = 0; t < T; ++t) {
    const std::size_t pos = t < L  ? t : t - L;
    t % 2U == 0
        ? extractSingleQubitGatesFromModel(t, model, qc, nSingleQubitGates)
        : extractTwoQubitGatesFromModel(t, model, qc, nTwoQubitGates);
  }

  res.setSingleQubitGates(nSingleQubitGates);
  res.setTwoQubitGates(nTwoQubitGates);
  res.setDepth(qc.getDepth());
  res.setResultCircuit(qc);
}

// mock-functions
void STQGatesEncoder::assertSingleQubitGateOrderConstraints(
    const std::size_t pos, const std::size_t qubit) {
  std::ostringstream posToStr;
  std::ostringstream qubitToStr;

  posToStr << pos;
  qubitToStr << qubit;
  throw std::runtime_error("Function assertSingleQubitGateOrderConstraints"
                           "with pos: " + posToStr.str() + " and qubit: " +
                           qubitToStr.str() + " is not provided for the STQGatesEncoder");
}

// mock-functions
void STQGatesEncoder::assertTwoQubitGateOrderConstraints(
    const std::size_t pos, const std::size_t ctrl, const std::size_t trgt) {
  std::ostringstream posToStr;
  std::ostringstream ctrlToStr;
  std::ostringstream trgtToStr;


  posToStr << pos;
  ctrlToStr << ctrl;
  trgtToStr << trgt;

  throw std::runtime_error("Function assertTwoQubitGateOrderConstraints"
                           "with pos: " + posToStr.str() + " , ctrl: " +
                           ctrlToStr.str() + " and target: " + trgtToStr.str() +
                            "is not provided for the STQGatesEncoder");
}
} // namespace cs::encoding
