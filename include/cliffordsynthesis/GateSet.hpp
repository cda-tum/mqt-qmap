//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "operations/OpType.hpp"

#include <initializer_list>
#include <vector>

namespace cs {
class GateSet {
  static constexpr std::array<qc::OpType, 10> SINGLE_QUBIT_CLIFFORDS = {
      qc::OpType::I,    qc::OpType::H,   qc::OpType::X,   qc::OpType::Y,
      qc::OpType::Z,    qc::OpType::S,   qc::OpType::Sdg, qc::OpType::SX,
      qc::OpType::SXdg, qc::OpType::None};

  using iterator       = typename std::vector<qc::OpType>::iterator;
  using const_iterator = typename std::vector<qc::OpType>::const_iterator;

private:
  // Check if None is already contained and if not, append it
  void appendNone() {
    if (!containsGate(qc::OpType::None)) {
      gates.push_back(qc::OpType::None);
    }
  }
  std::vector<qc::OpType> gates{};
public:

  
  GateSet() { appendNone(); };

  explicit GateSet(std::vector<qc::OpType> gateSet)
      : gates(std::move(gateSet)) {
    appendNone();
  };

  GateSet(const GateSet& gateSet) : gates(gateSet.gates) { appendNone(); };

  GateSet(GateSet&& gateSet) noexcept : gates(std::move(gateSet.gates)) {
    appendNone();
  };

  GateSet(std::initializer_list<qc::OpType> gateSet) : gates(gateSet) {
    appendNone();
  }

  GateSet& operator=(const GateSet& other) {
    gates = other.gates;
    appendNone();
    return *this;
  }

  GateSet& operator=(GateSet&& other) noexcept {
    gates = std::move(other.gates);
    appendNone();
    return *this;
  }

  GateSet& operator=(std::initializer_list<qc::OpType> other) {
    gates = other;
    appendNone();
    return *this;
  }

  void removePaulis() {
    gates.erase(std::remove_if(gates.begin(), gates.end(),
                               [](qc::OpType gate) {
                                 return gate == qc::OpType::X ||
                                        gate == qc::OpType::Y ||
                                        gate == qc::OpType::Z;
                               }),
                gates.end());
  }
  [[nodiscard]] bool containsGate(qc::OpType gate) const {
    for (const auto& g : // NOLINT(readability-use-anyofallof)
         gates) {
      if (g == gate) {
        return true;
      }
    }
    return false;
  }
  [[nodiscard]] bool containsX() const { return containsGate(qc::OpType::X); }
  [[nodiscard]] bool containsY() const { return containsGate(qc::OpType::Y); }
  [[nodiscard]] bool containsZ() const { return containsGate(qc::OpType::Z); }
  [[nodiscard]] bool containsH() const { return containsGate(qc::OpType::H); }
  [[nodiscard]] bool containsS() const { return containsGate(qc::OpType::S); }
  [[nodiscard]] bool containsSdg() const {
    return containsGate(qc::OpType::Sdg);
  }
  [[nodiscard]] bool containsSX() const { return containsGate(qc::OpType::SX); }
  [[nodiscard]] bool containsSXdg() const {
    return containsGate(qc::OpType::SXdg);
  }
  [[nodiscard]] std::size_t gateToIndex(const qc::OpType type) const {
    for (std::size_t i = 0; i < gates.size(); ++i) {
      if (gates.at(i) == type) {
        return i;
      }
    }
    return 0;
  }

  [[nodiscard]] GateSet paulis() const {
    std::vector<qc::OpType> result;
    for (const auto& g : gates) {
      if (g == qc::OpType::X || g == qc::OpType::Y || g == qc::OpType::Z || g == qc::OpType::I) {
        result.push_back(g);
      }
    }
    return GateSet(result);
  }
  
  [[nodiscard]] bool isValidGateSet() const {
    return std::all_of(gates.begin(), gates.end(), [](const auto& g) {
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
    for (const auto& g : gates) {
      result += qc::toString(g) + ", ";
    }
    result += "}";
    return result;
  }

  /**
   * Pass-Through
   */

  // Iterators (pass-through)
  auto               begin() noexcept { return gates.begin(); }
  [[nodiscard]] auto begin() const noexcept { return gates.begin(); }
  [[nodiscard]] auto cbegin() const noexcept { return gates.cbegin(); }
  auto               end() noexcept { return gates.end(); }
  [[nodiscard]] auto end() const noexcept { return gates.end(); }
  [[nodiscard]] auto cend() const noexcept { return gates.cend(); }
  auto               rbegin() noexcept { return gates.rbegin(); }
  [[nodiscard]] auto rbegin() const noexcept { return gates.rbegin(); }
  [[nodiscard]] auto crbegin() const noexcept { return gates.crbegin(); }
  auto               rend() noexcept { return gates.rend(); }
  [[nodiscard]] auto rend() const noexcept { return gates.rend(); }
  [[nodiscard]] auto crend() const noexcept { return gates.crend(); }

  // Capacity (pass-through)
  [[nodiscard]] bool        empty() const noexcept { return gates.empty(); }
  [[nodiscard]] std::size_t size() const noexcept { return gates.size(); }
  // NOLINTNEXTLINE(readability-identifier-naming)
  [[nodiscard]] std::size_t max_size() const noexcept {
    return gates.max_size();
  }
  [[nodiscard]] std::size_t capacity() const noexcept {
    return gates.capacity();
  }

  void reserve(const std::size_t newCap) { gates.reserve(newCap); }
  // NOLINTNEXTLINE(readability-identifier-naming)
  void shrink_to_fit() { gates.shrink_to_fit(); }

  // Modifiers (pass-through)
  void clear() noexcept { gates.clear(); }
  // NOLINTNEXTLINE(readability-identifier-naming)
  void     pop_back() { return gates.pop_back(); }
  void     resize(std::size_t count) { gates.resize(count); }
  iterator erase(const_iterator pos) { return gates.erase(pos); }
  iterator erase(const_iterator first, const_iterator last) {
    return gates.erase(first, last);
  }

  // NOLINTNEXTLINE(readability-identifier-naming)
  void push_back(const qc::OpType& gate) {
    if (!containsGate(gate)) {
      gates.push_back(gate);
    }
  }

  // NOLINTNEXTLINE(readability-identifier-naming)
  template <class T, class... Args> void emplace_back(Args&&... args) {
    gates.emplace_back(args...);
  }

  // NOLINTNEXTLINE(readability-identifier-naming)
  void emplace_back(qc::OpType& gate) {
    if (!containsGate(gate)) {
      gates.emplace_back(gate);
    }
  }

  // NOLINTNEXTLINE(readability-identifier-naming)
  void emplace_back(qc::OpType&& gate) {
    if (!containsGate(gate)) {
      gates.emplace_back(gate);
    }
  }

  [[nodiscard]] const auto& at(const std::size_t i) const {
    return gates.at(i);
  }
  [[nodiscard]] auto&       at(const std::size_t i) { return gates.at(i); }
  [[nodiscard]] const auto& front() const { return gates.front(); }
  [[nodiscard]] const auto& back() const { return gates.back(); }
};


} // namespace cs
