//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "hybridmap/Mapping.hpp"

#include "hybridmap/NeutralAtomDefinitions.hpp"

#include <stdexcept>

namespace na {
void Mapping::applySwap(Swap swap) {
  auto q1 = swap.first;
  auto q2 = swap.second;
  if (this->isMapped(q1) && this->isMapped(q2)) {
    auto circQ1 = this->getCircQubit(q1);
    auto circQ2 = this->getCircQubit(q2);
    this->setCircuitQubit(circQ2, q1);
    this->setCircuitQubit(circQ1, q2);
  } else if (this->isMapped(q1) && !this->isMapped(q2)) {
    this->setCircuitQubit(this->getCircQubit(q1), q2);
  } else if (this->isMapped(q2) && !this->isMapped(q1)) {
    this->setCircuitQubit(this->getCircQubit(q2), q1);
  } else {
    throw std::runtime_error("Cannot swap unmapped qubits");
  }
}

} // namespace na
