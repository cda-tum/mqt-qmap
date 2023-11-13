//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "operations/OpType.hpp"

#include <vector>

namespace cs {
class GateSet : public std::vector<qc::OpType> {
  static constexpr std::array<qc::OpType, 9> SINGLE_QUBIT_CLIFFORDS = {
      qc::OpType::I,   qc::OpType::H,  qc::OpType::X,
      qc::OpType::Y,   qc::OpType::Z,  qc::OpType::S,
      qc::OpType::Sdg, qc::OpType::SX, qc::OpType::SXdg};

public:
  using std::vector<qc::OpType>::vector;

  template <qc::OpType Gate> [[nodiscard]] bool containsGate() const {
    for (const auto& g : // NOLINT(readability-use-anyofallof)
         *this) {
      if (g == Gate) {
        return true;
      }
    }
    return false;
  }
  [[nodiscard]] bool containsX() const { return containsGate<qc::OpType::X>(); }
  [[nodiscard]] bool containsY() const { return containsGate<qc::OpType::Y>(); }
  [[nodiscard]] bool containsZ() const { return containsGate<qc::OpType::Z>(); }
  [[nodiscard]] bool containsH() const { return containsGate<qc::OpType::H>(); }
  [[nodiscard]] bool containsS() const { return containsGate<qc::OpType::S>(); }
  [[nodiscard]] bool containsSdg() const {
    return containsGate<qc::OpType::Sdg>();
  }
  [[nodiscard]] bool containsSX() const {
    return containsGate<qc::OpType::SX>();
  }
  [[nodiscard]] bool containsSXdg() const {
    return containsGate<qc::OpType::SXdg>();
  }
  [[nodiscard]] std::size_t gateToIndex(const qc::OpType type) const {
    for (std::size_t i = 0; i < this->size(); ++i) {
      if (this->at(i) == type) {
        return i;
      }
    }
    return 0;
  }

  [[nodiscard]] bool isValidGateSet() const {
    return std::all_of(this->begin(), this->end(), [](const auto& g) {
      return std::find(SINGLE_QUBIT_CLIFFORDS.begin(),
                       SINGLE_QUBIT_CLIFFORDS.end(),
                       g) != SINGLE_QUBIT_CLIFFORDS.end();
    });
  }

  // Any single-qubit Clifford gate can be obtained from a product of PI/2
  // rotations around different axes.
  [[nodiscard]] bool isComplete() const {
    if (containsS() || containsSdg()) {
      if (containsSX() || containsSXdg() || containsH()) {
        return true;
      }
    }

    if (containsSX() || containsSXdg()) {
      if (containsH()) {
        return true;
      }
    }
    return false;
  }

  [[nodiscard]] std::string toString() const {
    std::string result = "{";
    for (const auto& g : *this) {
      result += qc::toString(g) + ", ";
    }
    result += "}";
    return result;
  }
};
} // namespace cs
