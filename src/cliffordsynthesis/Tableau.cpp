//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/Tableau.hpp"

#include "utils.hpp"

namespace cs {
void Tableau::dump(const std::string& filename) const {
  auto of = std::ofstream(filename);
  if (!of.good()) {
    util::fatal("Error opening file " + filename);
  }
  dump(of);
}

void Tableau::dump(std::ostream& of) const { of << *this; }

void Tableau::import(const std::string& filename) {
  auto is = std::ifstream(filename);
  if (!is.good()) {
    util::fatal("Error opening file " + filename);
  }
  import(is);
}

void Tableau::import(std::istream& is) {
  tableau.clear();

  std::string              line;
  std::vector<std::string> data{};
  char                     delimiter = '|';

  while (std::getline(is, line)) {
    if (line.find('|', 0) == std::string::npos) {
      delimiter = ';';
    }
    tableau.emplace_back();
    ::parse_line(line, delimiter, {'\"'}, {'\\', '\r', '\n', '\t'}, data);
    for (const auto& datum : data) {
      if (datum.empty()) {
        continue;
      }
      tableau.back().emplace_back(static_cast<EntryType>(std::stoul(datum)));
    }
  }
}

void Tableau::populateTableauFrom(const std::uint64_t bv,
                                  const std::size_t   nQubits,
                                  const std::size_t   column) {
  for (std::size_t j = 0; j < nQubits; ++j) {
    if ((bv & (1U << j)) != 0U) {
      tableau[j][column] = 1;
    }
  }
}

void Tableau::applyGate(const qc::Operation* const gate) {
  const auto nqubits = getQubitCount();
  switch (gate->getType()) {
  case qc::OpType::H: { // HADAMARD
    if (gate->isControlled()) {
      util::fatal("Expected single-qubit gate");
    }
    const auto a = gate->getTargets().at(0U);
    applyGateH(a, nqubits);
  } break;
  case qc::OpType::S: { // PHASE
    if (gate->isControlled()) {
      util::fatal("Expected single-qubit gate");
    }
    const auto a = gate->getTargets().at(0U);
    applyGateS(a, nqubits);
  } break;
  case qc::OpType::X: { // CNOT and X
    if (gate->getNcontrols() != 1U) {
      const auto a = gate->getTargets().at(0U);
      applyGateX(a, nqubits);
    } else {
      const auto a = (*gate->getControls().begin()).qubit;
      const auto b = gate->getTargets().at(0);
      if (a == b) {
        util::fatal("Invalid CX with same control and target.");
      }
      applyGateCX(a, b, nqubits);
    }
  } break;
  case qc::OpType::Sdag: {
    if (gate->isControlled()) {
      util::fatal("Expected single-qubit gate");
    }
    const auto a = gate->getTargets().at(0U);
    applyGateSdag(a, nqubits);
  } break;
  case qc::OpType::Z: { // Z = S * S
    if (!gate->isControlled()) {
      const auto a = gate->getTargets().at(0U);
      applyGateZ(a, nqubits);
    } else { // CZ = H(1) x CX(0,1) x H(1)
      const auto a = (*gate->getControls().begin()).qubit;
      const auto b = gate->getTargets().at(0);
      if (a == b) {
        util::fatal("Invalid CZ with same control and target.");
      }
      applyGateH(b, nqubits);
      applyGateCX(a, b, nqubits);
      applyGateH(b, nqubits);
    }
  } break;
  case qc::OpType::Y: {
    if (!gate->isControlled()) {
      const auto a = gate->getTargets().at(0U);
      applyGateY(a, nqubits);
    } else { // CY = Sdag(1) * CX(0,1) * S(1)
      const auto a = (*gate->getControls().begin()).qubit;
      const auto b = gate->getTargets().at(0);
      if (a == b) {
        util::fatal("Invalid CY with same control and target.");
      }
      applyGateSdag(b, nqubits);
      applyGateCX(a, b, nqubits);
      applyGateS(b, nqubits);
    }
  } break;
  case qc::OpType::SWAP: {
    const auto a = gate->getTargets().at(0);
    const auto b = gate->getTargets().at(1);
    if (a == b) {
      util::fatal("Invalid SWAP with same control and target.");
    }
    applyGateCX(a, b, nqubits);
    applyGateCX(b, a, nqubits);
    applyGateCX(a, b, nqubits);
  } break;
  default:
    util::fatal("Unsupported gate encountered: " +
                std::to_string(gate->getType()));
    break;
  }
}

Tableau Tableau::getDiagonalTableau(const std::size_t nQubits) {
  TableauType result{};
  result.resize(nQubits);
  for (std::size_t i = 0U; i < nQubits; ++i) {
    result[i].resize((2U * nQubits) + 1U);
    for (std::size_t j = 0U; j < (2 * nQubits); ++j) {
      if (i == (j - nQubits)) {
        result[i][j] = 1;
      }
    }
    result[i][2U * nQubits] = 0;
  }

  return Tableau(result);
}

std::uint64_t Tableau::getBVFrom(const std::size_t column) const {
  std::uint64_t result  = 0UL;
  const auto    nQubits = getQubitCount();
  for (std::size_t j = 0U; j < nQubits; ++j) {
    if (tableau[j][column] == 1) {
      result |= (1U << j);
    }
  }
  return result;
}
void Tableau::init(const std::size_t nQubits) {
  tableau.clear();
  tableau.resize(nQubits);
  this->tableau = Tableau::getDiagonalTableau(nQubits).tableau;
}

bool Tableau::operator==(const Tableau& other) const {
  return tableau == other.tableau;
}

std::string Tableau::toString() const {
  std::stringstream ss;

  if (empty()) {
    util::debug("Tableau is empty.");
    return "";
  }
  for (const auto& row : tableau) {
    if (row.size() != back().size()) {
      util::fatal("Tableau is not rectangular.");
      return "";
    }
    for (const auto& s : row) {
      ss << std::to_string(s) << ';';
    }
    ss << "\n";
  }
  return ss.str();
}
void Tableau::fromString(const std::string& str) {
  std::stringstream ss(str);
  std::string       line;
  std::getline(ss, line);
  if (line.empty()) {
    return;
  }
  const auto  rStabilizer = std::regex("([\\+-])([IYZX]+)");
  std::smatch m;
  if (std::regex_search(line, rStabilizer)) {
    // string is a list of stabilizers
    auto iter = line.cbegin();
    while (std::regex_search(iter, line.cend(), m, rStabilizer)) {
      std::string s = m.str(0U);
      RowType     row;

      for (const auto c : s) {
        if (c == 'I' || c == 'Z') {
          row.push_back(0);
        } else if (c == 'X' || c == 'Y') {
          row.push_back(1);
        }
      }
      for (const auto c : s) {
        if (c == 'I' || c == 'X') {
          row.push_back(0);
        } else if (c == 'Y' || c == 'Z') {
          row.push_back(1);
        }
      }
      if (s[0U] == '-') {
        row.push_back(1);
      } else {
        row.push_back(0);
      }
      tableau.push_back(row);
      iter = m[0].second;
    }
  } else {
    // assume string is a semicolon separated binary matrix
    ss = std::stringstream(str);
    import(ss);
  }
}
void Tableau::applyGateH(const std::uint16_t target,
                         const std::size_t   nqubits) {
  for (auto i = 0U; i < nqubits; ++i) {
    tableau[i][2U * nqubits] ^=
        (tableau[i][target] & tableau[i][target + nqubits]);
    std::swap(tableau[i][target], tableau[i][target + nqubits]);
  }
}
void Tableau::applyGateS(const std::uint16_t target,
                         const std::size_t   nqubits) {
  for (auto i = 0U; i < nqubits; ++i) {
    tableau[i][2U * nqubits] ^=
        tableau[i][target] & tableau[i][target + nqubits];
    tableau[i][target + nqubits] ^= tableau[i][target];
  }
}
void Tableau::applyGateCX(const std::uint16_t control,
                          const std::uint16_t target,
                          const std::size_t   nqubits) {
  for (auto i = 0U; i < nqubits; ++i) {
    const auto xa = tableau[i][target];
    const auto za = tableau[i][target + nqubits];
    const auto xb = tableau[i][control];
    const auto zb = tableau[i][control + nqubits];
    tableau[i][2 * nqubits] ^= (xa & zb) & ((xb ^ za) ^ 1);
    tableau[i][control + nqubits] = za ^ zb;
    tableau[i][target]            = xb ^ xa;
  }
}
void Tableau::applyGateSdag(const std::uint16_t target,
                            const std::size_t   nqubits) {
  // Sdag = S * S * S
  applyGateS(target, nqubits);
  applyGateS(target, nqubits);
  applyGateS(target, nqubits);
}
void Tableau::applyGateX(const std::uint16_t target,
                         const std::size_t   nqubits) {
  // X = H * Z * H
  applyGateH(target, nqubits);
  applyGateZ(target, nqubits);
  applyGateH(target, nqubits);
}
void Tableau::applyGateY(const std::uint16_t target,
                         const std::size_t   nqubits) {
  // Y = X * Z
  applyGateX(target, nqubits);
  applyGateZ(target, nqubits);
}
void Tableau::applyGateZ(const std::uint16_t target,
                         const std::size_t   nqubits) {
  // Z = S * S
  applyGateS(target, nqubits);
  applyGateS(target, nqubits);
}

Tableau::Tableau(const qc::QuantumComputation& qc, const std::size_t begin,
                 const std::size_t end) {
  init(qc.getNqubits());
  std::size_t currentG = 0;
  for (const auto& gate : qc) {
    if (gate->getType() == qc::OpType::Compound) {
      const auto* const compOp =
          dynamic_cast<const qc::CompoundOperation* const>(gate.get());
      auto cit = compOp->begin();
      while (cit != compOp->end()) {
        if (currentG >= begin && (currentG < end)) {
          applyGate((*cit).get());
        }
        ++cit;
        ++currentG;
      }
    } else {
      if (currentG >= begin && (currentG < end)) {
        applyGate(gate.get());
      }
      ++currentG;
    }
    if (currentG >= end) {
      break;
    }
  }
}
Tableau::Tableau(const std::size_t nQubits) {
  tableau.clear();
  tableau.resize(nQubits);
  this->tableau = Tableau::getDiagonalTableau(nQubits).tableau;
}

Tableau::Tableau(const std::string& description) { fromString(description); }
} // namespace cs
