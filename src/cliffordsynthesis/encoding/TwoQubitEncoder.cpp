//
// Created by Velsh Aleksei on 16.06.23.
//

#include "cliffordsynthesis/encoding/TwoQubitEncoder.hpp"

#include "LogicTerm/LogicTerm.hpp"
#include "operations/OpType.hpp"
#include "utils/logging.hpp"

namespace cs::encoding {

using namespace logicbase;

void TwoQubitEncoder::createSingleQubitGateVariables() {
  DEBUG() << "Creating single-qubit gate variables.";
  vars.gS.reserve(T / 2);
  for (std::size_t t = 0U; t < T; t += 2U) {
    auto& timeStep = vars.gS.emplace_back();
    timeStep.reserve(SINGLE_QUBIT_GATES.size());
    for (const auto gate : SINGLE_QUBIT_GATES) {
      auto& g = timeStep.emplace_back();
      g.reserve(N);
      for (std::size_t q = 0U; q < N; ++q) {
        const std::string gName = "g_" + std::to_string(t) + "_" +
                                  toString(gate) + "_" + std::to_string(q);
        TRACE() << "Creating variable " << gName;
        g.emplace_back(lb->makeVariable(gName));
      }
    }
  }
  createPauliGateVariables();
}

void TwoQubitEncoder::createTwoQubitGateVariables() {
  DEBUG() << "Creating two-qubit gate variables.";
  vars.gC.reserve(T / 2);
  for (std::size_t t = 1U; t < T; t += 2U) {
    auto& timeStep = vars.gC.emplace_back();
    timeStep.reserve(N);
    for (std::size_t ctrl = 0U; ctrl < N; ++ctrl) {
      auto& control = timeStep.emplace_back();
      control.reserve(N);
      for (std::size_t trgt = 0U; trgt < N; ++trgt) {
        const std::string gName = "g_" + std::to_string(t) + "_cx_" +
                                  std::to_string(ctrl) + "_" +
                                  std::to_string(trgt);
        TRACE() << "Creating variable " << gName;
        control.emplace_back(lb->makeVariable(gName));
      }
    }
  }
}

void TwoQubitEncoder::createPauliGateVariables() {
  DEBUG() << "Creating pauli-qubit gate variables.";
  gP.reserve(PAULI_GATES.size());
  for (const auto gate : PAULI_GATES) {
    auto& g = gP.emplace_back();
    g.reserve(N);
    for (std::size_t q = 0U; q < N; ++q) {
      const std::string gName = "g_" + std::to_string(T) + "_" +
                                toString(gate) + "_" + std::to_string(q);
      TRACE() << "Creating variable " << gName;
      g.emplace_back(lb->makeVariable(gName));
    }
  }
}

void encoding::TwoQubitEncoder::assertConsistency() const {
  DEBUG() << "Asserting gate consistency";
  for (std::size_t t = 0U; t < T; t += 2U) {
    // asserting only a single gate is applied on each qubit.
    for (std::size_t q = 0U; q < N; ++q) {
      LogicVector singleQubitGateVariables{};
      LogicVector twoQubitGateVariables{};

      vars.collectSingleQubitGateVariables(t / 2, q, singleQubitGateVariables);
      assertExactlyOne(singleQubitGateVariables);
      IF_PLOG(plog::verbose) {
        TRACE() << "Single Qubit Gate variables at time " << t << " and qubit "
                << q;
        for (const auto& var : singleQubitGateVariables) {
          TRACE() << var.getName();
        }
      }

      vars.collectTwoQubitGateVariables(t / 2, q, true, twoQubitGateVariables);
      vars.collectTwoQubitGateVariables(t / 2, q, false, twoQubitGateVariables);

      // add a variable that will be treated as identity
      twoQubitGateVariables.emplace_back(vars.gC[t / 2][q][q]);

      assertExactlyOne(twoQubitGateVariables);
      IF_PLOG(plog::verbose) {
        std::cout << twoQubitGateVariables.size() << std::endl;

        TRACE() << "Two Qubit Gate variables at time " << t + 1 << " and qubit "
                << q;
        for (const auto& var : twoQubitGateVariables) {
          TRACE() << var.getName();
        }
      }
    }
  }

  LogicVector pauliGateVariables{};
  for (std::size_t q = 0U; q < N; ++q) {
    collectPauliGateVariables(q, pauliGateVariables);
    assertExactlyOne(pauliGateVariables);
    IF_PLOG(plog::verbose) {
      TRACE() << "Pauli Qubit Gate variables at time " << T << " and qubit "
              << q;
      for (const auto& var : pauliGateVariables) {
        TRACE() << var.getName();
      }
    }
  }
}

void encoding::TwoQubitEncoder::assertGateConstraints() {
  DEBUG() << "Asserting gate constraints";
  xorHelpers = logicbase::LogicMatrix{T + 1};
  for (std::size_t t = 0U; t < T; ++t) {
    TRACE() << "Asserting gate constraints at time " << t;
    splitXorR(tvars->r[t], t);
    if (t % 2 == 0) {
      assertSingleQubitGateConstraints(t);
    } else {
      assertTwoQubitGateConstraints(t);
    }
    TRACE() << "Asserting r changes at time " << t;
    lb->assertFormula(tvars->r[t + 1] == xorHelpers[t].back());
  }

  TRACE() << "Asserting pauli gate constraints at time " << T;
  splitXorR(tvars->r[T], T);

  assertPauliGateConstraints(T);

  TRACE() << "Asserting r changes at time " << T;
  lb->assertFormula(tvars->r[T + 1] == xorHelpers[T].back());
}

void encoding::TwoQubitEncoder::assertSingleQubitGateConstraints(
    const std::size_t pos) {
  for (std::size_t q = 0U; q < N; ++q) {
    assertZConstraints(pos, q);
    assertXConstraints(pos, q);
    assertRConstraints(pos, q);
  }
}

void TwoQubitEncoder::assertRConstraints(const std::size_t pos,
                                         const std::size_t qubit) {
  for (const auto gate : SINGLE_QUBIT_GATES) {
    const auto& change =
        LogicTerm::ite(vars.gS[pos / 2][gateToIndex(gate)][qubit],
                       tvars->singleQubitRChange(pos, qubit, gate),
                       LogicTerm(0, static_cast<std::int16_t>(S)));
    splitXorR(change, pos);
  }
}

void encoding::TwoQubitEncoder::assertTwoQubitGateConstraints(
    const std::size_t pos) {
  const auto& twoQubitGates = vars.gC[pos / 2];
  for (std::size_t ctrl = 0U; ctrl < N; ++ctrl) {
    for (std::size_t trgt = 0U; trgt < N; ++trgt) {
      if (ctrl == trgt) {
        // create identity constraint on TQG
        auto changes = tvars->x[pos + 1][ctrl] == tvars->x[pos][ctrl];
        changes      = changes && (tvars->z[pos + 1][ctrl] ==
                              tvars->z[pos][ctrl]); // && here is overloaded

        lb->assertFormula(
            LogicTerm::implies(twoQubitGates[ctrl][ctrl], changes));
        splitXorR(tvars->singleQubitRChange(pos, ctrl, qc::OpType::None), pos);
      } else {
        const auto changes = createTwoQubitGateConstraint(pos, ctrl, trgt);
        lb->assertFormula(
            LogicTerm::implies(twoQubitGates[ctrl][trgt], changes));
        DEBUG() << "Asserting CNOT on " << ctrl << " and " << trgt;
      }
    }
  }
}

void encoding::TwoQubitEncoder::assertPauliGateConstraints(
    const std::size_t pos) {
  for (std::size_t q = 0U; q < N; ++q) {
    for (const auto gate : PAULI_GATES) {
      const auto& change = LogicTerm::ite(
          gP[gateToIndex(gate)][q], tvars->singleQubitRChange(pos, q, gate),
          LogicTerm(0, static_cast<std::int16_t>(S)));

      splitXorR(change, pos);
    }
    lb->assertFormula(tvars->x[pos][q] == tvars->x[pos + 1][q]);
    lb->assertFormula(tvars->z[pos][q] == tvars->z[pos + 1][q]);
  }
}

void TwoQubitEncoder::assertGatesImplyTransform(
    const std::size_t pos, const std::size_t qubit,
    const std::vector<TransformationFamily>& transformations) {
  const auto& singleQubitGates = vars.gS[pos / 2];
  for (const auto& [transformation, gates] : transformations) {
    auto gateOr = LogicTerm(false);
    for (const auto& gate : gates) {
      gateOr = gateOr || singleQubitGates[gateToIndex(gate)][qubit];
    }
    lb->assertFormula(LogicTerm::implies(gateOr, transformation));
  }
}

void TwoQubitEncoder::splitXorR(const logicbase::LogicTerm& changes,
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

void encoding::TwoQubitEncoder::collectPauliGateVariables(
    const std::size_t qubit, LogicVector& variables) const {
  for (const auto& gate : gP) {
    variables.emplace_back(gate[qubit]);
  }
}

LogicTerm encoding::TwoQubitEncoder::createTwoQubitGateConstraint(
    std::size_t pos, std::size_t ctrl, std::size_t trgt) {
  auto changes              = LogicTerm(true);
  const auto [xCtrl, xTrgt] = tvars->twoQubitXChange(pos, ctrl, trgt);
  const auto [zCtrl, zTrgt] = tvars->twoQubitZChange(pos, ctrl, trgt);

  changes = changes && (tvars->x[pos + 1][ctrl] == xCtrl);
  changes = changes && (tvars->x[pos + 1][trgt] == xTrgt);
  changes = changes && (tvars->z[pos + 1][ctrl] == zCtrl);
  changes = changes && (tvars->z[pos + 1][trgt] == zTrgt);

  const auto& newRChanges = LogicTerm::ite(
      vars.gC[pos / 2][ctrl][trgt], tvars->twoQubitRChange(pos, ctrl, trgt),
      LogicTerm(0, static_cast<std::int16_t>(S)));
  splitXorR(newRChanges, pos);

  return changes;
}

void TwoQubitEncoder::extractCircuitFromModel(Results& res, Model& model) {
  std::size_t nSingleQubitGates = 0U;
  std::size_t nTwoQubitGates    = 0U;

  qc::QuantumComputation qc(N);
  for (std::size_t t = 0; t < T / 2; ++t) {
    extractSingleQubitGatesFromModel(t, model, qc, nSingleQubitGates);
    extractTwoQubitGatesFromModel(t, model, qc, nTwoQubitGates);
  }

  res.setSingleQubitGates(nSingleQubitGates);
  res.setTwoQubitGates(nTwoQubitGates);
  res.setDepth(qc.getDepth());
  res.setTwoQubitDepth(qc.getTwoQubitDepth());
  res.setResultCircuit(qc);
}

// mock-functions
void TwoQubitEncoder::encodeSymmetryBreakingConstraints() {
  DEBUG() << "Encoding symmetry breaking constraints IS NOT supported for the "
             "two qubit encoding.";
}

void TwoQubitEncoder::assertSingleQubitGateOrderConstraints(
    const std::size_t pos, const std::size_t qubit) {
  std::ostringstream posToStr;
  std::ostringstream qubitToStr;

  posToStr << pos;
  qubitToStr << qubit;
}

// mock-functions
void TwoQubitEncoder::assertTwoQubitGateOrderConstraints(
    const std::size_t pos, const std::size_t ctrl, const std::size_t trgt) {
  std::ostringstream posToStr;
  std::ostringstream ctrlToStr;
  std::ostringstream trgtToStr;

  posToStr << pos;
  ctrlToStr << ctrl;
  trgtToStr << trgt;
}
} // namespace cs::encoding
