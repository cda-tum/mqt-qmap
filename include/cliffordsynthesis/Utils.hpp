//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "cliffordsynthesis/GateSet.hpp"
#include "operations/OpType.hpp"

namespace cs {

const GateSet PAULIS{
    qc::OpType::None,
    qc::OpType::X,
    qc::OpType::Y,
    qc::OpType::Z,
};

inline bool isPauli(const qc::OpType& gate) {
  return PAULIS.containsGate(gate);
}

qc::OpType multiplyPaulis(const GateSet& paulis);

} // namespace cs
