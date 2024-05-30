//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//
#include "hybridmap/HybridSynthesisMapper.hpp"

#include "CircuitOptimizer.hpp"
#include "operations/OpType.hpp"

#include <cstdint>

namespace qc {

QuantumComputation
HybridSynthesisMapper::completelyRemap(InitialMapping initialMapping) {
  // copy mapped circuit, removing SWAPs and MOVEs
  QuantumComputation newCirc = this->mappedQc;
  CircuitOptimizer::removeOpTypes(newCirc, {OpType::SWAP, OpType::Move});
  this->map(newCirc, initialMapping, false);
  return this->mappedQc;
}

uint32_t HybridSynthesisMapper::evaluateSynthesisSteps(qcs& synthesisSteps,
                                                       bool directlyMap) {
  return 0;
}
void HybridSynthesisMapper::directlyMap(const QuantumComputation& qc) {}

} // namespace qc
