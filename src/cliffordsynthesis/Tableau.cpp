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
    ::parseLine(line, delimiter, {'\"'}, {'\\', '\r', '\n', '\t'}, data);
    for (const auto& datum : data) {
      if (datum.empty()) {
        continue;
      }
      tableau.back().emplace_back(static_cast<EntryType>(std::stoul(datum)));
    }
  }
}

void Tableau::applyGate(const qc::Operation* const gate) {
  if (gate->getNcontrols() > 1U) {
    util::fatal("Tableau::applyGate: Only operations with up to one control "
                "are supported.");
  }

  const auto target = static_cast<std::size_t>(gate->getTargets().at(0U));

  // non-controlled gates
  if (!gate->isControlled()) {
    switch (gate->getType()) {
    case qc::OpType::H:
      applyH(target);
      break;
    case qc::OpType::S:
      applyS(target);
      break;
    case qc::OpType::Sdag:
      applySdag(target);
      break;
    case qc::OpType::SX:
      applySx(target);
      break;
    case qc::OpType::SXdag:
      applySxdag(target);
      break;
    case qc::OpType::X:
      applyX(target);
      break;
    case qc::OpType::Y:
      applyY(target);
      break;
    case qc::OpType::Z:
      applyZ(target);
      break;
    case qc::OpType::SWAP: {
      const auto target2 = static_cast<std::size_t>(gate->getTargets().at(1U));
      applySwap(target, target2);
      break;
    }
    default:
      // unsupported non-controlled gate type
      util::fatal("Tableau::applyGate: Unsupported non-controlled gate type " +
                  qc::toString(gate->getType()));
    }
  } else {
    const auto control =
        static_cast<std::size_t>((*gate->getControls().begin()).qubit);
    switch (gate->getType()) {
    case qc::OpType::X:
      applyCX(control, target);
      break;
    case qc::OpType::Y:
      applyCY(control, target);
      break;
    case qc::OpType::Z:
      applyCZ(control, target);
      break;
    default:
      // unsupported controlled gate type
      util::fatal("Tableau::applyGate: Unsupported controlled gate type " +
                  qc::toString(gate->getType()));
    }
  }
}

void Tableau::createDiagonalTableau(const std::size_t nQ,
                                    const bool        useFullsizeTableau) {
  nQubits = nQ;
  tableau.clear();
  if (useFullsizeTableau) {
    tableau.resize(2U * nQubits);
  } else {
    tableau.resize(nQubits);
  }
  for (std::size_t i = 0U; i < getTableauSize(); ++i) {
    tableau[i].resize((2U * nQubits) + 1U);
    if (useFullsizeTableau) {
      for (std::size_t j = 0; j < (2U * nQubits); ++j) {
        tableau[i][j] = 0;
        if (i == j) {
          tableau[i][j] = 1;
        }
      }
    } else {
      for (std::size_t j = nQubits; j < (2U * nQubits); ++j) {
        tableau[i][j] = 0;
        if (i == (j - nQubits)) {
          tableau[i][j] = 1;
        }
      }
    }
  }
}

std::string Tableau::toString() const {
  std::stringstream ss;
  for (const auto& row : tableau) {
    if (row.size() != tableau.back().size()) {
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
  const auto  rStabilizer = std::regex("([\\+-]?)([IYZX]+)");
  std::smatch m;
  if (std::regex_search(line, rStabilizer)) {
    // string is a list of stabilizers
    loadStabilizerDestabilizerString(str);
  } else {
    // assume string is a semicolon separated binary matrix
    ss = std::stringstream(str);
    import(ss);
  }
}

void Tableau::applyH(const std::size_t target) {
  assert(target < nQubits);
  for (std::size_t i = 0U; i < getTableauSize(); ++i) {
    tableau[i][2U * nQubits] ^= static_cast<EntryType>(
        tableau[i][target] & tableau[i][target + nQubits]);
    std::swap(tableau[i][target], tableau[i][target + nQubits]);
  }
}

void Tableau::applyS(const std::size_t target) {
  assert(target < nQubits);
  for (std::size_t i = 0U; i < getTableauSize(); ++i) {
    tableau[i][2U * nQubits] ^= static_cast<EntryType>(
        tableau[i][target] & tableau[i][target + nQubits]);
    tableau[i][target + nQubits] ^= tableau[i][target];
  }
}

// Sdag = S * S * S
void Tableau::applySdag(const std::size_t target) {
  assert(target < nQubits);
  applyS(target);
  applyS(target);
  applyS(target);
}

// Sx = Sdag * H * Sdag
void Tableau::applySx(const std::size_t target) {
  assert(target < nQubits);
  applySdag(target);
  applyH(target);
  applySdag(target);
}

// Sxdag = S * H * S
void Tableau::applySxdag(const std::size_t target) {
  assert(target < nQubits);
  applyS(target);
  applyH(target);
  applyS(target);
}

// X = H * Z * H
void Tableau::applyX(const std::size_t target) {
  assert(target < nQubits);
  applyH(target);
  applyZ(target);
  applyH(target);
}

// Y = X * Z
void Tableau::applyY(const std::size_t target) {
  assert(target < nQubits);
  applyX(target);
  applyZ(target);
}

// Z = S * S
void Tableau::applyZ(const std::size_t target) {
  assert(target < nQubits);
  applyS(target);
  applyS(target);
}

void Tableau::applyCX(const std::size_t control, const std::size_t target) {
  assert(control < nQubits);
  assert(target < nQubits);
  assert(control != target);
  for (auto i = 0U; i < getTableauSize(); ++i) {
    const auto xa = tableau[i][target];
    const auto za = tableau[i][target + nQubits];
    const auto xb = tableau[i][control];
    const auto zb = tableau[i][control + nQubits];
    tableau[i][2 * nQubits] ^= (xa & zb) & ((xb ^ za) ^ 1);
    tableau[i][control + nQubits] = za ^ zb;
    tableau[i][target]            = xb ^ xa;
  }
}

void Tableau::applyCY(const std::size_t control, const std::size_t target) {
  assert(control < nQubits);
  assert(target < nQubits);
  assert(control != target);
  applySdag(target);
  applyCX(control, target);
  applyS(target);
}

void Tableau::applyCZ(const std::size_t control, const std::size_t target) {
  assert(control < nQubits);
  assert(target < nQubits);
  assert(control != target);
  applyH(target);
  applyCX(control, target);
  applyH(target);
}

void Tableau::applySwap(const std::size_t q1, const std::size_t q2) {
  assert(q1 < nQubits);
  assert(q2 < nQubits);
  assert(q1 != q2);
  applyCX(q1, q2);
  applyCX(q2, q1);
  applyCX(q1, q2);
}

Tableau::Tableau(const qc::QuantumComputation& qc, const std::size_t begin,
                 const std::size_t end, const bool useFullsizeTableau)
    : Tableau(qc.getNqubits(), useFullsizeTableau) {
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
void Tableau::fromStabilizersDestabilizers(const std::string& stabilizers,
                                           const std::string& destabilizers) {
  loadStabilizerDestabilizerString(destabilizers);
  loadStabilizerDestabilizerString(stabilizers);
}
void Tableau::loadStabilizerDestabilizerString(const std::string& string) {
  std::stringstream ss(string);
  std::string       line;
  std::getline(ss, line);
  if (line.empty()) {
    return;
  }
  const auto  rStabilizer = std::regex("([\\+-]?)([IYZX]+)");
  std::smatch m;
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
}
} // namespace cs
