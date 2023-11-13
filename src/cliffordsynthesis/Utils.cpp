#include "cliffordsynthesis/Utils.hpp"
#include "operations/OpType.hpp"

namespace cs {
qc::OpType multiplyPaulis(const GateSet& paulis                        ) {
  if(paulis.empty()) {
    return qc::OpType::None;
  }

  qc::OpType combinedPauli = qc::OpType::None;
  if(paulis.size() == 1) {
    combinedPauli = paulis.at(0);
  } else if(paulis.size() == 2) {
    if(paulis.containsX()) {
      if(paulis.containsY()) {
        combinedPauli = qc::OpType::Z;
      } else {
        combinedPauli = qc::OpType::Y;
      }
    } else {
      combinedPauli = qc::OpType::X;
    }
  } else if(paulis.size() == 3) {
    if(!paulis.containsX()) {
      combinedPauli = qc::OpType::X;
    } else if(!paulis.containsY()) {
      combinedPauli = qc::OpType::Y;
    } else {
      combinedPauli = qc::OpType::Z;
    }
  }
     return combinedPauli;
};
}  // namespace cs
