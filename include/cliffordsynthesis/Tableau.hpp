//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "QuantumComputation.hpp"
#include "utils/logging.hpp"

#include <cstdint>
#include <fstream>
#include <limits>
#include <ostream>
#include <utility>
#include <vector>

namespace cs {
class Tableau {
  using EntryType   = std::uint8_t;
  using RowType     = std::vector<EntryType>;
  using TableauType = std::vector<RowType>;
  std::size_t nQubits{};
  TableauType tableau;

public:
  Tableau() = default;
  explicit Tableau(const qc::QuantumComputation& qc, std::size_t begin = 0,
                   std::size_t end = std::numeric_limits<std::size_t>::max());
  explicit Tableau(const std::size_t nQubits) : nQubits(nQubits) {
    createDiagonalTableau(nQubits);
  }
  explicit Tableau(const std::string& description) {
    fromString(description);
    nQubits = tableau.size();
  }

  [[nodiscard]] RowType operator[](const std::size_t index) {
    return tableau[index];
  }
  [[nodiscard]] RowType operator[](const std::size_t index) const {
    return tableau[index];
  }

  [[nodiscard]] RowType at(const std::size_t index) {
    return tableau.at(index);
  }

  [[nodiscard]] auto getQubitCount() const { return nQubits; }

  void dump(const std::string& filename) const;

  void dump(std::ostream& of) const;

  void import(const std::string& filename);
  void import(std::istream& is);

  template <std::size_t N>
  void populateTableauFrom(const std::bitset<N> bv, const std::size_t nQ,
                           const std::size_t column) {
    assert(column <= 2 * nQ);
    assert(nQ <= nQubits);
    assert(nQ <= N);
    for (std::size_t i = 0U; i < nQ; ++i) {
      if (bv[i]) {
        tableau[i][column] = 1;
      } else {
        tableau[i][column] = 0;
      }
    }
  }
  void populateTableauFrom(const std::uint64_t bv, const std::size_t nQ,
                           const std::size_t column) {
    populateTableauFrom(std::bitset<64>(bv), nQ, column);
  }

  void applyGate(const qc::Operation* gate);
  void applyH(std::size_t target);
  void applyS(std::size_t target);
  void applySdag(std::size_t target);
  void applySx(std::size_t target);
  void applySxdag(std::size_t target);
  void applyX(std::size_t target);
  void applyY(std::size_t target);
  void applyZ(std::size_t target);
  void applyCX(std::size_t control, std::size_t target);
  void applyCY(std::size_t control, std::size_t target);
  void applyCZ(std::size_t control, std::size_t target);
  void applySwap(std::size_t q1, std::size_t q2);

  friend bool operator==(const Tableau& lhs, const Tableau& rhs) {
    return lhs.tableau == rhs.tableau;
  }
  friend bool operator!=(const Tableau& lhs, const Tableau& rhs) {
    return !(lhs == rhs);
  }

  void createDiagonalTableau(std::size_t nQ);

  friend std::ostream& operator<<(std::ostream& os, const Tableau& dt) {
    os << dt.toString();
    return os;
  }
  friend std::istream& operator>>(std::istream& is, Tableau& dt) {
    if (is.good()) {
      std::stringstream ss;
      ss << is.rdbuf();
      dt.fromString(ss.str());
    }
    return is;
  }

  [[nodiscard]] std::string toString() const;
  void                      fromString(const std::string& str);

  template <std::size_t N>
  [[nodiscard]] std::bitset<N> getBVFrom(const std::size_t column) const {
    assert(column <= 2 * nQubits);
    assert(nQubits <= N);
    std::bitset<N> bv;
    for (std::size_t i = 0U; i < nQubits; ++i) {
      if (tableau[i][column] == 1U) {
        bv[i] = 1;
      }
    }
    return bv;
  }
  [[nodiscard]] std::uint64_t getBVFrom(const std::size_t column) const {
    return getBVFrom<64>(column).to_ullong();
  }
};
} // namespace cs
