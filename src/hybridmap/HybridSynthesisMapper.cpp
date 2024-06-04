//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//
#include "hybridmap/HybridSynthesisMapper.hpp"

#include "CircuitOptimizer.hpp"
#include "QuantumComputation.hpp"
#include "hybridmap/NeutralAtomUtils.hpp"
#include "operations/OpType.hpp"

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
void HybridSynthesisMapper::directlyMap(const qc::QuantumComputation& qc) {}

} // namespace na
