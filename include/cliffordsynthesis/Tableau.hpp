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
  TableauType tableau;

public:
  [[nodiscard]] Tableau() = default;
  [[nodiscard]] explicit Tableau(TableauType tableau)
      : tableau(std::move(tableau)) {}

  [[nodiscard]] explicit Tableau(
      const qc::QuantumComputation& qc, std::size_t begin = 0,
      std::size_t end = std::numeric_limits<std::size_t>::max());
  [[nodiscard]] explicit Tableau(std::size_t nQubits);
  [[nodiscard]] explicit Tableau(const std::string& description);

  [[nodiscard]] RowType operator[](const std::size_t index) {
    return tableau[index];
  }
  [[nodiscard]] RowType operator[](const std::size_t index) const {
    return tableau[index];
  }

  [[nodiscard]] RowType at(const std::size_t index) {
    return tableau.at(index);
  }

  [[nodiscard]] auto getQubitCount() const { return tableau.size(); }

  void resize(const std::size_t size) { tableau.resize(size); }

  void clear() { tableau.clear(); }

  [[nodiscard]] bool empty() const { return tableau.empty(); }

  [[nodiscard]] auto back() const { return tableau.back(); }

  [[nodiscard]] auto begin() const { return tableau.begin(); }
  [[nodiscard]] auto end() const { return tableau.end(); }

  void dump(const std::string& filename) const;

  void dump(std::ostream& of) const;

  void import(const std::string& filename);
  void import(std::istream& is);

  void init(std::size_t nQubits);

  void populateTableauFrom(std::uint64_t bv, std::size_t nQubits,
                           std::size_t column);

  void applyGate(const qc::Operation* gate);

  [[nodiscard]] bool operator==(const Tableau& other) const;
  [[nodiscard]] bool operator!=(const Tableau& other) const {
    return !(other == *this);
  }

  [[nodiscard]] static Tableau getDiagonalTableau(std::size_t nQubits);

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

  [[nodiscard]] std::uint64_t getBVFrom(std::size_t column) const;

private:
  void applyGateH(std::uint16_t target, std::size_t nqubits);
  void applyGateS(std::uint16_t target, std::size_t nqubits);
  void applyGateSdag(std::uint16_t target, std::size_t nqubits);
  void applyGateX(std::uint16_t target, std::size_t nqubits);
  void applyGateY(std::uint16_t target, std::size_t nqubits);
  void applyGateZ(std::uint16_t target, std::size_t nqubits);
  void applyGateCX(std::uint16_t control, std::uint16_t target,
                   std::size_t nqubits);
};
} // namespace cs
