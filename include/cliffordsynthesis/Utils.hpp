#pragma once

#include "cliffordsynthesis/GateSet.hpp"
#include "operations/OpType.hpp"

namespace cs {
    inline bool isPauli(const qc::OpType& gate) {
    return gate == qc::OpType::X || gate == qc::OpType::Y ||
           gate == qc::OpType::Z || gate == qc::OpType::I;
  }

  qc::OpType multiplyPaulis(const GateSet& paulis);

} // namespace cs;
