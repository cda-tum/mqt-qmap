//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//
#include "hybridmap/HybridSynthesisMapper.hpp"

#include "Definitions.hpp"
#include "QuantumComputation.hpp"
#include "hybridmap/HybridNeutralAtomMapper.hpp"
#include "hybridmap/NeutralAtomDefinitions.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <utility>
#include <vector>

namespace na {

size_t HybridSynthesisMapper::evaluateSynthesisSteps(qcs& synthesisSteps,
                                                     bool directlyMap) {
  std::vector<std::pair<qc::QuantumComputation, qc::fp>> costs;
  for (auto& qc : synthesisSteps) {
    costs.emplace_back(qc, this->evaluateSynthesisStep(qc));
  }
  const auto bestQc = std::min_element(
      costs.begin(), costs.end(),
      [](const auto& a, const auto& b) { return a.second < b.second; });
  if (directlyMap) {
    this->appendWithoutMapping(bestQc->first);
  }
  return static_cast<size_t>(std::distance(costs.begin(), bestQc));
}

qc::fp
HybridSynthesisMapper::evaluateSynthesisStep(qc::QuantumComputation& qc) {
  NeutralAtomMapper tempMapper(arch, parameters);
  tempMapper.loadHwQubits(hardwareQubits);
  const auto mappedQc = tempMapper.map(qc, mapping);
  const auto results  = tempMapper.schedule();
  return results.totalFidelities;
}

void HybridSynthesisMapper::appendWithoutMapping(
    const qc::QuantumComputation& qc) {
  for (const auto& op : qc) {
    this->synthesizedQc.emplace_back(op->clone());
    this->mapGate(op.get());
  }
}

AdjacencyMatrix HybridSynthesisMapper::getCircuitAdjacencyMatrix() const {
  auto            numCircQubits = mappedQc.getNqubits();
  AdjacencyMatrix adjMatrix(numCircQubits);

  for (uint32_t i = 0; i < numCircQubits; ++i) {
    for (uint32_t j = 0; j < i; ++j) {
      auto mappedI = this->mapping.getHwQubit(i);
      auto mappedJ = this->mapping.getHwQubit(j);
      if (this->arch->getSwapDistance(mappedI, mappedJ) == 0) {
        adjMatrix(i, j) = 1;
      } else {
        adjMatrix(i, j) = 0;
      }
    }
  }
  return adjMatrix;
}

} // namespace na
