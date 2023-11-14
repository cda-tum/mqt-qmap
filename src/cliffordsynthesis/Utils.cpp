#include "cliffordsynthesis/Utils.hpp"

#include "operations/OpType.hpp"

namespace cs {
qc::OpType multiplyPaulis(const GateSet& paulis) {
  auto pauliCopy = paulis;
  pauliCopy.removeGate(qc::OpType::None);
  if (paulis.empty()) {
    return qc::OpType::None;
  }

  if (pauliCopy.size() == 1) {
    return pauliCopy.at(0);
  }

  if (pauliCopy.size() == 2) {
    if (pauliCopy.containsX()) {
      if (pauliCopy.containsY()) {
        return qc::OpType::Z;
      }
      return qc::OpType::Y;
    }
    return qc::OpType::X;
  }

  if (pauliCopy.size() == 3) {
    return qc::OpType::None;
  }
  return qc::OpType::None;
};
} // namespace cs
