//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//
#include "hybridmap/HybridSynthesisMapper.hpp"

#include "CircuitOptimizer.hpp"
#include "QuantumComputation.hpp"
#include "hybridmap/NeutralAtomDefinitions.hpp"
#include "hybridmap/NeutralAtomUtils.hpp"

#include <cstdint>

namespace na {

qc::QuantumComputation
HybridSynthesisMapper::completelyRemap(InitialMapping initialMapping) {
  this->map(unmappedQc, initialMapping);
  return this->mappedQc;
}

uint32_t HybridSynthesisMapper::evaluateSynthesisSteps(qcs& synthesisSteps,
                                                       bool directlyMap) {
  return 0;
}
void HybridSynthesisMapper::directlyMap(const qc::QuantumComputation& qc) {
  for (const auto& op : qc) {
    this->mapGate(op.get());
  }
}

AdjacencyMatrix HybridSynthesisMapper::getAdjacencyMatrix() const {
  auto            numCircQubits = mappedQc.getNqubits();
  AdjacencyMatrix adjMatrix(numCircQubits);

  for (uint32_t i = 0; i < numCircQubits; ++i) {
    for (uint32_t j = 0; j < i; ++j) {
      auto mappedI = this->mapping.getHwQubit(i);
      auto mappedJ = this->mapping.getHwQubit(j);
      if (this->arch.getSwapDistance(mappedI, mappedJ) == 0) {
        adjMatrix(i, j) = 1;
      } else {
        adjMatrix(i, j) = 0;
      }
    }
  }
  return adjMatrix;
}

} // namespace na
