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
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

namespace na {

size_t HybridSynthesisMapper::evaluateSynthesisSteps(qcs& synthesisSteps,
                                                     bool alsoMap) {
  std::vector<std::pair<qc::QuantumComputation, qc::fp>> costs;
  size_t                                                 qcIndex = 0;
  for (auto& qc : synthesisSteps) {
    if (this->parameters->verbose) {
      std::cout << "Evaluating synthesis step number " << qcIndex++ << "\n";
    }
    costs.emplace_back(qc, this->evaluateSynthesisStep(qc));
    if (this->parameters->verbose) {
      std::cout << "Fidelity: " << costs.back().second << "\n";
    }
    qcIndex++;
  }
  const auto bestQc = std::max_element(
      costs.begin(), costs.end(),
      [](const auto& a, const auto& b) { return a.second < b.second; });
  if (alsoMap) {
    this->appendWithMapping(bestQc->first);
  }
  return static_cast<size_t>(std::distance(costs.begin(), bestQc));
}

qc::fp
HybridSynthesisMapper::evaluateSynthesisStep(qc::QuantumComputation& qc) {
  NeutralAtomMapper tempMapper;
  tempMapper.copyStateFrom(*this);
  auto mappedQc = tempMapper.map(qc, mapping);
  tempMapper.convertToAod();
  const auto results = tempMapper.schedule();
  return results.totalFidelities;
}

void HybridSynthesisMapper::appendWithoutMapping(
    const qc::QuantumComputation& qc) {
  for (const auto& op : qc) {
    this->synthesizedQc.emplace_back(op->clone());
    this->mapGate(op.get());
  }
}

void HybridSynthesisMapper::appendWithMapping(qc::QuantumComputation& qc) {
  if (mappedQc.empty()) {
    initMapping(qc.getNqubits());
  }
  mapAppend(qc, this->mapping);
  for (const auto& op : qc) {
    this->synthesizedQc.emplace_back(op->clone());
  }
}

AdjacencyMatrix HybridSynthesisMapper::getCircuitAdjacencyMatrix() const {
  auto            numCircQubits = synthesizedQc.getNqubits();
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