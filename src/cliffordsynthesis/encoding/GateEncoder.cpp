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
  const auto variableGrouping =
      encodings::groupVars(variables, variables.size() / 2U);
  lb->assertFormula(encodings::exactlyOneCmdr(variableGrouping,
                                              LogicTerm::noneTerm(), lb.get()));
}

void GateEncoder::extractCircuitFromModel(Results& res, Model& model) {
  std::size_t nSingleQubitGates = 0U;
  std::size_t nTwoQubitGates    = 0U;
  std::size_t depth             = 0U;

  qc::QuantumComputation qc(N);
  for (std::size_t t = 0; t < T; ++t) {
    const auto gateCount = nSingleQubitGates + nTwoQubitGates;
    DEBUG() << "Timestep: " << t << ", gate count: " << gateCount
            << ", depth: " << depth;

    extractSingleQubitGatesFromModel(t, model, qc, nSingleQubitGates);
    extractTwoQubitGatesFromModel(t, model, qc, nTwoQubitGates);

    if ((nSingleQubitGates + nTwoQubitGates) > gateCount) {
      ++depth;
    }
  }

  res.setSingleQubitGates(nSingleQubitGates);
  res.setTwoQubitGates(nTwoQubitGates);
  res.setDepth(depth);
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
          dd::Control{static_cast<dd::Qubit>(ctrl), dd::Control::Type::pos};
      if (model.getBoolValue(twoQubitGates[ctrl][trgt], lb.get())) {
        qc.emplace_back<qc::StandardOperation>(N, control, trgt, qc::OpType::X);
        ++nTwoQubitGates;
        DEBUG() << "CX(" << ctrl << ", " << trgt << ")";
      }
    }
  }
}

} // namespace cs::encoding