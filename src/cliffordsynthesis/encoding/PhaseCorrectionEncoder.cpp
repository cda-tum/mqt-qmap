#include "cliffordsynthesis/encoding/PhaseCorrectionEncoder.hpp"

#include "LogicTerm/Logic.hpp"
#include "LogicTerm/LogicTerm.hpp"
#include "cliffordsynthesis/GateSet.hpp"
#include "cliffordsynthesis/Utils.hpp"
#include "utils/logging.hpp"

#include <bitset>
#include <cstddef>
#include <string>
#include <vector>

namespace cs::encoding {
using namespace logicbase;

std::vector<qc::OpType> PhaseCorrectionEncoder::phaseCorrection() {
  createVariables();
  encodePauliConstraints();
  lb->solve();
  return extractResult();
}

void PhaseCorrectionEncoder::createVariables() {
  DEBUG() << "Creating phase correction gate variables.";
  for (std::size_t i = 0U; i < N; ++i) {
    paulis[i].reserve(4);
    for (const auto& pauli : PAULIS) {
      const std::string gName =
          "g_" + std::to_string(i) + "_" + toString(pauli);
      TRACE() << "Creating pauli variable " << gName;
      paulis[i].emplace_back(lb->makeVariable(gName));
    }
  }

  std::bitset<64> init = uncorrected.getBVFrom(2 * N);
  initialPhase         = vectorFromBitset(init);
  std::bitset<64> tar  = target.getBVFrom(2 * N);
  targetPhase          = vectorFromBitset(tar);
}

void PhaseCorrectionEncoder::encodePauliConstraints() {
  DEBUG() << "Encoding pauli constraints.";
  splitXorR(initialPhase);
  for (std::size_t q = 0U; q < N; ++q) {
    const auto& xCol = vectorFromBitset(uncorrected.getBVFrom(q));
    const auto& zCol = vectorFromBitset(uncorrected.getBVFrom(N + q));

    LogicVector iChange{};
    LogicVector xChange{};
    LogicVector zChange{};
    LogicVector yChange{};
    iChange.reserve(S);
    xChange.reserve(S);
    zChange.reserve(S);
    yChange.reserve(S);

    for (std::size_t row = 0U; row < S; ++row) {
      DEBUG() << "Encoding identity constraint.";
      iChange.emplace_back(logicbase::LogicTerm::ite(
          paulis[q][0], LogicTerm(false), LogicTerm(false)));
      DEBUG() << "Encoding X constraint.";
      xChange.emplace_back(
          logicbase::LogicTerm::ite(paulis[q][1], zCol[row], LogicTerm(false)));
      DEBUG() << "Encoding Y constraint.";
      yChange.emplace_back(logicbase::LogicTerm::ite(
          paulis[q][2], !(zCol[row] == xCol[row]), LogicTerm(false)));
      DEBUG() << "Encoding Z constraint.";
      zChange.emplace_back(
          logicbase::LogicTerm::ite(paulis[q][3], xCol[row], LogicTerm(false)));
    }
    splitXorR(iChange);
    splitXorR(xChange);
    splitXorR(yChange);
    splitXorR(zChange);

    // at least one gate per qubit
    lb->assertFormula(paulis[q][0] || paulis[q][1] || paulis[q][2] ||
                      paulis[q][3]);
  }
  DEBUG() << "Assert Phase.";
  for (std::size_t row = 0U; row < S; ++row) {
    lb->assertFormula(targetPhase[row] == xorHelpers.back()[row]);
  }
  DEBUG() << "Finished Phase Encoding.";
}

void PhaseCorrectionEncoder::splitXorR(const LogicVector& changes) {
  auto xorHelper = LogicVector{};
  xorHelper.reserve(S);
  for (std::size_t row = 0U; row < S; ++row) {
    const std::string hName =
        "h_" + std::to_string(row) + "_" + std::to_string(xorHelpers.size());
    DEBUG() << "Creating helper variable for RChange XOR " << hName;
    xorHelper.emplace_back(lb->makeVariable(hName));
  }
  xorHelpers.emplace_back(xorHelper);
  if (xorHelpers.size() == 1) {
    for (std::size_t row = 0U; row < S; ++row) {
      lb->assertFormula(xorHelpers[0][row] == changes[row]);
    }
  } else {
    for (std::size_t row = 0U; row < S; ++row) {
      lb->assertFormula(
          xorHelpers.back()[row] ==
          (changes[row] != xorHelpers[xorHelpers.size() - 2][row]));
    }
  }
}

LogicVector
PhaseCorrectionEncoder::vectorFromBitset(const std::bitset<64>& bs) const {
  LogicVector result{};
  result.reserve(S);
  for (std::size_t i = 0U; i < S; ++i) {
    result.emplace_back(LogicTerm(bs[i]));
  }
  return result;
}

std::vector<qc::OpType> PhaseCorrectionEncoder::extractResult() {
  const auto&             model = lb->getModel();
  std::vector<qc::OpType> result{};
  result.reserve(N);
  for (std::size_t q = 0U; q < N; ++q) {
    GateSet pauliGates{};
    pauliGates.removePaulis();
    if (model->getBoolValue(paulis[q][0], lb.get())) {
      pauliGates.emplace_back(qc::OpType::None);
      DEBUG() << "Identity gate for qubit " << q;
    }
    if (model->getBoolValue(paulis[q][1], lb.get())) {
      pauliGates.emplace_back(qc::OpType::X);
      DEBUG() << "X gate for qubit " << q;
    }
    if (model->getBoolValue(paulis[q][2], lb.get())) {
      pauliGates.emplace_back(qc::OpType::Y);
      DEBUG() << "Y gate for qubit " << q;
    }
    if (model->getBoolValue(paulis[q][3], lb.get())) {
      pauliGates.emplace_back(qc::OpType::Z);
      DEBUG() << "Z gate for qubit " << q;
    }
    result.emplace_back(multiplyPaulis(pauliGates));
  }
  return result;
}
} // namespace cs::encoding
