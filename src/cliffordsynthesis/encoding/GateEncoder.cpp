//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/encoding/GateEncoder.hpp"

#include "Encodings/Encodings.hpp"
#include "utils/logging.hpp"

namespace cs::encoding {

using namespace logicbase;

void GateEncoder::createSingleQubitGateVariables() {
  DEBUG() << "Creating single-qubit gate variables.";
  vars.gS.reserve(T);
  for (std::size_t t = 0U; t < T; ++t) {
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
}

void GateEncoder::createTwoQubitGateVariables() {
  DEBUG() << "Creating two-qubit gate variables.";
  vars.gC.reserve(T);
  for (std::size_t t = 0U; t < T; ++t) {
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

void GateEncoder::Variables::collectSingleQubitGateVariables(
    const std::size_t pos, const std::size_t qubit,
    LogicVector& variables) const {
  const auto& singleQubitGates = gS[pos];
  for (const auto& gate : singleQubitGates) {
    variables.emplace_back(gate[qubit]);
  }
}

void GateEncoder::Variables::collectTwoQubitGateVariables(
    const std::size_t pos, const std::size_t qubit, const bool target,
    LogicVector& variables) const {
  const auto& twoQubitGates = gC[pos];
  const auto  n             = twoQubitGates.size();
  for (std::size_t q = 0; q < n; ++q) {
    if (q == qubit) {
      continue;
    }
    if (target) {
      variables.emplace_back(twoQubitGates[q][qubit]);
    } else {
      variables.emplace_back(twoQubitGates[qubit][q]);
    }
  }
}

void GateEncoder::assertExactlyOne(const LogicVector& variables) const {
  const auto variableGrouping = encodings::groupVars(variables, 3U);
  lb->assertFormula(encodings::exactlyOneCmdr(variableGrouping,
                                              LogicTerm::noneTerm(), lb.get()));
}

std::vector<GateEncoder::TransformationFamily>
GateEncoder::collectGateTransformations(
    const std::size_t pos, const std::size_t qubit,
    const GateToTransformation& gateToTransformation) {
  std::vector<TransformationFamily> transformations;

  for (const auto& gate : SINGLE_QUBIT_GATES) {
    const auto& transformation = gateToTransformation(pos, qubit, gate);
    const auto& it             = std::find_if(
        transformations.begin(), transformations.end(), [&](const auto& entry) {
          return entry.first.deepEquals(transformation);
        });
    if (it != transformations.end()) {
      it->second.emplace_back(gate);
    } else {
      transformations.emplace_back(transformation,
                                   std::vector<qc::OpType>{gate});
    }
  }
  return transformations;
}

void GateEncoder::assertGatesImplyTransform(
    const std::size_t pos, const std::size_t qubit,
    const std::vector<TransformationFamily>& transformations) {
  const auto& singleQubitGates = vars.gS[pos];
  for (const auto& [transformation, gates] : transformations) {
    auto gateOr = LogicTerm(false);
    for (const auto& gate : gates) {
      gateOr = gateOr || singleQubitGates[gateToIndex(gate)][qubit];
    }
    lb->assertFormula(LogicTerm::implies(gateOr, transformation));
  }
}

void GateEncoder::assertZConstraints(const std::size_t pos,
                                     const std::size_t qubit) {
  const auto& gatesToZTransformations = [this](const auto& p1, const auto& p2,
                                               const auto& p3) {
    return tvars->singleQubitZChange(p1, p2, p3);
  };
  auto gateTransformations =
      collectGateTransformations(pos, qubit, gatesToZTransformations);
  for (auto& [transformation, _] : gateTransformations) {
    transformation = tvars->z[pos + 1][qubit] == transformation;
  }
  assertGatesImplyTransform(pos, qubit, gateTransformations);
}

void GateEncoder::assertXConstraints(const std::size_t pos,
                                     const std::size_t qubit) {
  const auto& gatesToXTransformations = [this](const auto& p1, const auto& p2,
                                               const auto& p3) {
    return tvars->singleQubitXChange(p1, p2, p3);
  };
  auto gateTransformations =
      collectGateTransformations(pos, qubit, gatesToXTransformations);
  for (auto& [transformation, _] : gateTransformations) {
    transformation = tvars->x[pos + 1][qubit] == transformation;
  }
  assertGatesImplyTransform(pos, qubit, gateTransformations);
}

void GateEncoder::assertRConstraints(const std::size_t pos,
                                     const std::size_t qubit) {
  const auto& gatesToRTransformations = [this](const auto& p1, const auto& p2,
                                               const auto& p3) {
    return tvars->singleQubitRChange(p1, p2, p3);
  };
  auto gateTransformations =
      collectGateTransformations(pos, qubit, gatesToRTransformations);
  for (auto& [transformation, _] : gateTransformations) {
    transformation = tvars->r[pos + 1] == (tvars->r[pos] ^ transformation);
  }
  assertGatesImplyTransform(pos, qubit, gateTransformations);
}

void GateEncoder::extractCircuitFromModel(Results& res, Model& model) {
  std::size_t nSingleQubitGates = 0U;
  std::size_t nTwoQubitGates    = 0U;

  qc::QuantumComputation qc(N);
  for (std::size_t t = 0; t < T; ++t) {
    extractSingleQubitGatesFromModel(t, model, qc, nSingleQubitGates);
    extractTwoQubitGatesFromModel(t, model, qc, nTwoQubitGates);
  }

  res.setSingleQubitGates(nSingleQubitGates);
  res.setTwoQubitGates(nTwoQubitGates);
  res.setDepth(qc.getDepth());
  res.setResultCircuit(qc);
}

void GateEncoder::extractSingleQubitGatesFromModel(
    const std::size_t pos, Model& model, qc::QuantumComputation& qc,
    std::size_t& nSingleQubitGates) {
  const auto& singleQubitGates = vars.gS[pos];
  for (std::size_t q = 0U; q < N; ++q) {
    for (const auto gate : SINGLE_QUBIT_GATES) {
      if (gate == qc::OpType::None) {
        continue;
      }
      if (model.getBoolValue(singleQubitGates[gateToIndex(gate)][q],
                             lb.get())) {
        qc.emplace_back<qc::StandardOperation>(N, q, gate);
        ++nSingleQubitGates;
        DEBUG() << toString(gate) << "(" << q << ")";
      }
    }
  }
}

void GateEncoder::extractTwoQubitGatesFromModel(const std::size_t       pos,
                                                Model&                  model,
                                                qc::QuantumComputation& qc,
                                                size_t& nTwoQubitGates) {
  const auto& twoQubitGates = vars.gC[pos];
  for (std::size_t ctrl = 0U; ctrl < N; ++ctrl) {
    for (std::size_t trgt = 0U; trgt < N; ++trgt) {
      if (ctrl == trgt) {
        continue;
      }
      const auto control =
          qc::Control{static_cast<qc::Qubit>(ctrl), qc::Control::Type::Pos};
      if (model.getBoolValue(twoQubitGates[ctrl][trgt], lb.get())) {
        qc.emplace_back<qc::StandardOperation>(N, control, trgt, qc::OpType::X);
        ++nTwoQubitGates;
        DEBUG() << "CX(" << ctrl << ", " << trgt << ")";
      }
    }
  }
}

void GateEncoder::encodeSymmetryBreakingConstraints() {
  DEBUG() << "Encoding symmetry breaking constraints.";
  for (std::size_t t = 0U; t < T; ++t) {
    assertSingleQubitGateSymmetryBreakingConstraints(t);
    assertTwoQubitGateSymmetryBreakingConstraints(t);
  }
}

void GateEncoder::assertSingleQubitGateCancellationConstraints(
    const std::size_t pos, const std::size_t qubit) {
  // nothing to assert for the last timestep
  if (pos == T - 1U) {
    return;
  }

  // gate variables of the current and the next timestep
  const auto& gSNow  = vars.gS[pos];
  const auto& gSNext = vars.gS[pos + 1U];

  for (const auto gate : SINGLE_QUBIT_GATES) {
    if (gate == qc::OpType::None) {
      continue;
    }
    const auto gateIndex = gateToIndex(gate);
    switch (gate) {
    case qc::OpType::X:
    case qc::OpType::Y:
    case qc::OpType::H:
      // self-inverse gates
      lb->assertFormula(LogicTerm::implies(gSNow[gateIndex][qubit],
                                           !gSNext[gateIndex][qubit]));
      break;
    case qc::OpType::Z:
      lb->assertFormula(LogicTerm::implies(gSNow[gateIndex][qubit],
                                           !gSNext[gateIndex][qubit]));
      lb->assertFormula(LogicTerm::implies(
          gSNow[gateIndex][qubit], !gSNext[gateToIndex(qc::OpType::X)][qubit]));
      break;
    case qc::OpType::S:
    case qc::OpType::Sdag: {
      // Sdag * S = S * Sdag = I
      // Sdag * Sdag = S * S = Z
      auto disallowed = !gSNext[gateToIndex(qc::OpType::S)][qubit] &&
                        !gSNext[gateToIndex(qc::OpType::Sdag)][qubit];
      // Z and Sdag commute. Z should always precede S and Sdag
      disallowed = disallowed && !gSNext[gateToIndex(qc::OpType::Z)][qubit];
      lb->assertFormula(
          LogicTerm::implies(gSNow[gateIndex][qubit], disallowed));
      break;
    }
    default:
      break;
    }
  }

  if (pos >= T - 2U) {
    return;
  }
  const auto& gSNextNext = vars.gS[pos + 2U];

  for (const auto gate : SINGLE_QUBIT_GATES) {
    if (gate == qc::OpType::None) {
      continue;
    }
    const auto gateIndex = gateToIndex(gate);
    switch (gate) {
    case qc::OpType::X:
      lb->assertFormula(
          LogicTerm::implies(gSNow[gateIndex][qubit],
                             !(gSNext[gateToIndex(qc::OpType::Z)][qubit] &&
                               gSNextNext[gateToIndex(qc::OpType::Y)][qubit])));
      break;
    case qc::OpType::Y:
      lb->assertFormula(
          LogicTerm::implies(gSNow[gateIndex][qubit],
                             !(gSNext[gateToIndex(qc::OpType::Z)][qubit] &&
                               gSNextNext[gateToIndex(qc::OpType::X)][qubit])));
      break;
    case qc::OpType::Z:
      lb->assertFormula(
          LogicTerm::implies(gSNow[gateIndex][qubit],
                             !(gSNext[gateToIndex(qc::OpType::Y)][qubit] &&
                               gSNextNext[gateToIndex(qc::OpType::X)][qubit])));
      lb->assertFormula(
          LogicTerm::implies(gSNow[gateIndex][qubit],
                             !(gSNext[gateToIndex(qc::OpType::S)][qubit] &&
                               gSNextNext[gateToIndex(qc::OpType::S)][qubit])));
      break;
    case qc::OpType::H:
      lb->assertFormula(
          LogicTerm::implies(gSNow[gateIndex][qubit],
                             !(gSNext[gateToIndex(qc::OpType::X)][qubit] &&
                               gSNextNext[gateToIndex(qc::OpType::H)][qubit])));
      lb->assertFormula(
          LogicTerm::implies(gSNow[gateIndex][qubit],
                             !(gSNext[gateToIndex(qc::OpType::Z)][qubit] &&
                               gSNextNext[gateToIndex(qc::OpType::H)][qubit])));
    default:
      break;
    }
  }
}

void GateEncoder::assertSingleQubitGateSymmetryBreakingConstraints(
    const std::size_t pos) {
  for (std::size_t q = 0U; q < N; ++q) {
    assertSingleQubitGateOrderConstraints(pos, q);
    assertSingleQubitGateCancellationConstraints(pos, q);
  }
}

void GateEncoder::assertTwoQubitGateSymmetryBreakingConstraints(
    const std::size_t pos) {
  for (std::size_t ctrl = 1U; ctrl < N; ++ctrl) {
    for (std::size_t trgt = 0U; trgt < ctrl; ++trgt) {
      assertTwoQubitGateOrderConstraints(pos, ctrl, trgt);
    }
  }
}

} // namespace cs::encoding
