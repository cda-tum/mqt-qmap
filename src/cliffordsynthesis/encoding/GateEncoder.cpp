//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/encoding/GateEncoder.hpp"

#include "Definitions.hpp"
#include "Encodings/Encodings.hpp"
#include "LogicTerm/LogicTerm.hpp"
#include "QuantumComputation.hpp"
#include "cliffordsynthesis/GateSet.hpp"
#include "cliffordsynthesis/Tableau.hpp"
#include "cliffordsynthesis/Utils.hpp"
#include "operations/OpType.hpp"
#include "utils/logging.hpp"

#include <vector>

namespace cs::encoding {

using namespace logicbase;

void GateEncoder::createSingleQubitGateVariables() {
  DEBUG() << "Creating single-qubit gate variables.";
  vars.gS.reserve(T);
  for (std::size_t t = 0U; t < T; ++t) {
    auto& timeStep = vars.gS.emplace_back();
    timeStep.reserve(singleQubitGates.size());
    for (const auto gate : singleQubitGates) {
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
  const auto& singleQubitGateVariables = gS[pos];
  for (const auto& gate : singleQubitGateVariables) {
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

LogicTerm
GateEncoder::createExactlyOne(const logicbase::LogicVector& variables) const {
  const auto variableGrouping = encodings::groupVars(variables, 3U);
  return encodings::exactlyOneCmdr(variableGrouping, LogicTerm::noneTerm(),
                                   lb.get());
}
void GateEncoder::assertExactlyOne(const LogicVector& variables) const {
  lb->assertFormula(createExactlyOne(variables));
}

std::vector<GateEncoder::TransformationFamily>
GateEncoder::collectGateTransformations(
    const std::size_t pos, const std::size_t qubit,
    const GateToTransformation& gateToTransformation) {
  std::vector<TransformationFamily> transformations;

  for (const auto& gate : singleQubitGates) {
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
  const auto& singleQubitGateVars = vars.gS[pos];
  for (const auto& [transformation, gates] : transformations) {
    auto gateOr = LogicTerm(false);
    for (const auto& gate : gates) {
      gateOr = gateOr ||
               singleQubitGateVars[singleQubitGates.gateToIndex(gate)][qubit];
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

qc::QuantumComputation GateEncoder::extractCircuitFromModel(Results& res,
                                                            Model&   model) {
  std::size_t nSingleQubitGates = 0U;
  std::size_t nTwoQubitGates    = 0U;

  Tableau                intermediateTableau = init;
  qc::QuantumComputation qc(N);
  for (std::size_t t = 0; t < T; ++t) {
    std::vector<bool> hasGate(N, false);
    extractSingleQubitGatesFromModel(t, model, qc, nSingleQubitGates, hasGate,
                                     intermediateTableau);
    extractTwoQubitGatesFromModel(t, model, qc, nTwoQubitGates, hasGate,
                                  intermediateTableau);
  }

  res.setSingleQubitGates(nSingleQubitGates);
  res.setTwoQubitGates(nTwoQubitGates);
  res.setDepth(qc.getDepth());
  res.setResultCircuit(qc);
  return qc;
}

void GateEncoder::assertSingleQubitGateConstraints(std::size_t pos) {
  for (std::size_t q = 0U; q < N; ++q) {
    DEBUG() << "Asserting gates on " << q;
    assertZConstraints(pos, q);
    assertXConstraints(pos, q);
    if (!ignoreRChanges) {
      assertRConstraints(pos, q);
    }
  }
}

void GateEncoder::assertTwoQubitGateConstraints(std::size_t pos) {
  const auto& twoQubitGates = vars.gC[pos];
  for (std::size_t ctrl = 0U; ctrl < N; ++ctrl) {
    for (std::size_t trgt = 0U; trgt < N; ++trgt) {
      if (ctrl == trgt) {
        continue;
      }
      auto changes = createTwoQubitGateConstraint(pos, ctrl, trgt);
      if (!ignoreRChanges) {
        changes = changes && createTwoQubitRConstraint(pos, ctrl, trgt);
      }
      lb->assertFormula(LogicTerm::implies(twoQubitGates[ctrl][trgt], changes));
      DEBUG() << "Asserting CNOT on " << ctrl << " and " << trgt;
    }
  }
}

void GateEncoder::extractSingleQubitGatesFromModel(
    const std::size_t pos, Model& model, qc::QuantumComputation& qc,
    std::size_t& nSingleQubitGates, std::vector<bool>& hasGate,
    Tableau& intermediateTableau) {
  const auto& singleQubitGateVars = vars.gS[pos];
  auto        validPaulis         = singleQubitGates.paulis();
  for (std::size_t q = 0U; q < N; ++q) {
    if (hasGate[q]) {
      continue;
    }
    for (const auto gate : singleQubitGates) {
      if (gate == qc::OpType::None) {
        continue;
      }

      if (model.getBoolValue(
              singleQubitGateVars[singleQubitGates.gateToIndex(gate)][q],
              lb.get())) {
        auto temp = intermediateTableau;
        temp.applySingleQGate(gate, q);
        if ((ignoreRChanges && intermediateTableau.equalUpToPhase(temp)) ||
            intermediateTableau == temp) {
          continue;
        }
        intermediateTableau.applySingleQGate(gate, q);
        qc.emplace_back<qc::StandardOperation>(N, q, gate);
        hasGate[q] = true;
        ++nSingleQubitGates;
        DEBUG() << toString(gate) << "(" << q << ")";
        break;
      }
    }
  }
}

void GateEncoder::extractTwoQubitGatesFromModel(const std::size_t       pos,
                                                Model&                  model,
                                                qc::QuantumComputation& qc,
                                                size_t& nTwoQubitGates,
                                                std::vector<bool>& hasGate,
                                                Tableau& intermediateTableau) {
  const auto& twoQubitGates = vars.gC[pos];
  for (std::size_t ctrl = 0U; ctrl < N; ++ctrl) {
    if (hasGate[ctrl]) {
      continue;
    }
    for (std::size_t trgt = 0U; trgt < N; ++trgt) {
      if (hasGate[trgt] || ctrl == trgt) {
        continue;
      }
      const auto control =
          qc::Control{static_cast<qc::Qubit>(ctrl), qc::Control::Type::Pos};
      if (model.getBoolValue(twoQubitGates[ctrl][trgt], lb.get())) {
        auto temp = intermediateTableau;
        temp.applyCX(ctrl, trgt);
        DEBUG() << temp.toString();
        DEBUG() << intermediateTableau.toString();
        if ((ignoreRChanges && intermediateTableau.equalUpToPhase(temp)) ||
            intermediateTableau == temp) { // gate doesn't do anything
          continue;
        }
        intermediateTableau.applyCX(ctrl, trgt);
        qc.emplace_back<qc::StandardOperation>(N, control, trgt, qc::OpType::X);
        hasGate[ctrl] = true;
        hasGate[trgt] = true;
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

  // Any Pauli must not be followed by another Pauli since -iXYZ = I.
  std::vector<qc::OpType> paulis{};
  if (singleQubitGates.containsX()) {
    paulis.emplace_back(qc::OpType::X);
  }
  if (singleQubitGates.containsY()) {
    paulis.emplace_back(qc::OpType::Y);
  }
  if (singleQubitGates.containsZ()) {
    paulis.emplace_back(qc::OpType::Z);
  }
  const bool containsPaulis = singleQubitGates.containsX() ||
                              singleQubitGates.containsY() ||
                              singleQubitGates.containsZ();
  if (containsPaulis) {
    auto gates      = gSNow[singleQubitGates.gateToIndex(paulis[0])][qubit];
    auto disallowed = !gSNext[singleQubitGates.gateToIndex(paulis[0])][qubit];
    for (std::size_t i = 1U; i < paulis.size(); ++i) {
      gates = gates || gSNow[singleQubitGates.gateToIndex(paulis[i])][qubit];
      disallowed =
          disallowed && !gSNext[singleQubitGates.gateToIndex(paulis[i])][qubit];
    }

    if (singleQubitGates.containsH()) {
      // -(X|Y|Z)-H- ~= -H-(Z|Y|X)-
      auto gateIndex = singleQubitGates.gateToIndex(qc::OpType::H);
      disallowed     = disallowed && !gSNext[gateIndex][qubit];
    }

    if (singleQubitGates.containsS() && singleQubitGates.containsSdg()) {
      const auto gateIndexS   = singleQubitGates.gateToIndex(qc::OpType::S);
      const auto gateIndexSdg = singleQubitGates.gateToIndex(qc::OpType::Sdg);

      // -X-(S|Sd)- ~= -(Sd|S)-X-
      // -Y-(S|Sd)- ~= -(Sd|S)-Y-
      // -Z-(S|Sd)-  = -(S|Sd)-Z-
      disallowed = disallowed && !gSNext[gateIndexS][qubit] &&
                   !gSNext[gateIndexSdg][qubit];
    }

    lb->assertFormula(LogicTerm::implies(gates, disallowed));
  }

  // H is self-inverse
  if (singleQubitGates.containsH()) {
    auto gateIndex = singleQubitGates.gateToIndex(qc::OpType::H);
    lb->assertFormula(
        LogicTerm::implies(gSNow[gateIndex][qubit], !gSNext[gateIndex][qubit]));
  }

  if (singleQubitGates.containsS()) {
    auto gateIndexS = singleQubitGates.gateToIndex(qc::OpType::S);

    if (singleQubitGates.containsZ()) {
      auto gateIndexZ = singleQubitGates.gateToIndex(qc::OpType::Z);

      // -S-S- = -Z-
      auto gates      = gSNow[gateIndexS][qubit];
      auto disallowed = !gSNext[gateIndexS][qubit];

      if (singleQubitGates.containsSdg()) {
        auto gateIndexSdag = singleQubitGates.gateToIndex(qc::OpType::Sdg);

        // -Sd-Sd- = -Z-
        // -Sd-S-  = -I-
        // -Sd-Z-  = -S-
        // -S-Sd-  = -I-
        // -S-Z-   = -Sd-
        gates      = gates || gSNow[gateIndexSdag][qubit];
        disallowed = disallowed && !gSNext[gateIndexSdag][qubit] &&
                     !gSNext[gateIndexZ][qubit];
      }

      lb->assertFormula(LogicTerm::implies(gates, disallowed));
    } else {
      if (singleQubitGates.containsSdg()) {
        auto gateIndexSdag = singleQubitGates.gateToIndex(qc::OpType::Sdg);

        // -S-Sd- = -I-
        // -Sd-S- = -I-
        lb->assertFormula(LogicTerm::implies(gSNow[gateIndexS][qubit],
                                             !gSNext[gateIndexSdag][qubit]));
        lb->assertFormula(LogicTerm::implies(gSNow[gateIndexSdag][qubit],
                                             !gSNext[gateIndexS][qubit]));
      }
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
