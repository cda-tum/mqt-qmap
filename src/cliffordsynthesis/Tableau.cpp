/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
*/

#include "cliffordsynthesis/Tableau.hpp"

#include "utils.hpp"

void Tableau::dump(const std::string& filename) const {
    auto of = std::ofstream(filename);
    if (!of.good()) {
        FATAL() << "Error opening file " << filename;
    }
    dump(of);
}

void Tableau::dump(std::ostream& of) const {
    of << *this;
}

void Tableau::import(const std::string& filename) {
    auto is = std::ifstream(filename);
    if (!is.good()) {
        FATAL() << "Error opening file " << filename;
    }
    import(is);
}

void Tableau::import(std::istream& is) {
    tableau.clear();

    std::string              token;
    std::string              line;
    std::vector<std::string> data{};
    char                     delimiter = '|';

    while (std::getline(is, line)) {
        if (line.find('|', 0) == std::string::npos) {
            delimiter = ';';
        }
        tableau.emplace_back();
        parse_line(line, delimiter, {'\"'}, {'\\', '\r', '\n', '\t'}, data);
        for (const auto& datum: data) {
            if (datum.empty()) {
                continue;
            }
            tableau.back().emplace_back(static_cast<EntryType>(std::stoul(datum)));
        }
    }
}

void Tableau::populateTableauFrom(std::uint64_t bv, std::size_t nQubits,
                                  int column) {
    for (std::size_t j = 0; j < nQubits; ++j) {
        if ((bv & (1U << j)) != 0U) {
            tableau[j][column] = 1;
        }
    }
}

void Tableau::applyGate(const std::unique_ptr<qc::Operation>& gate) {
    auto nqubits = getQubitCount();
    switch (gate->getType()) {
        case qc::OpType::H: // HADAMARD
        {
            if (gate->isControlled()) {
                util::fatal("Expected single-qubit gate");
            }
            const auto a = gate->getTargets().at(0U);
            applyGateH(a, nqubits);
        } break;
        case qc::OpType::S: // PHASE
        {
            if (gate->isControlled()) {
                util::fatal("Expected single-qubit gate");
            }
            const auto a = gate->getTargets().at(0U);
            applyGateS(a, nqubits);
        } break;
        case qc::OpType::X: // CNOT
        {
            if (gate->getNcontrols() != 1U) { // NOT = H x S x S x H
                const auto a = gate->getTargets().at(0U);
                applyGateH(a, nqubits);
                applyGateS(a, nqubits);
                applyGateS(a, nqubits);
                applyGateH(a, nqubits);
            } else {
                const auto a = (*gate->getControls().begin()).qubit;
                const auto b = gate->getTargets().at(0);
                if (a == b) {
                    util::fatal("Invalid CX with same control and target.");
                }
                applyGateCX(a, b, nqubits);
            }
        } break;
        case qc::OpType::Sdag: { // Sdag  = S x S x S
            if (gate->isControlled()) {
                util::fatal("Expected single-qubit gate");
            }
            const auto a = gate->getTargets().at(0U);
            applyGateS(a, nqubits);
            applyGateS(a, nqubits);
            applyGateS(a, nqubits);
        } break;
        case qc::OpType::Z: { // Z = S x S
            if (!gate->isControlled()) {
                const auto a = gate->getTargets().at(0U);
                applyGateS(a, nqubits);
                applyGateS(a, nqubits);
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
        case qc::OpType::Y: { // Y = H x S x S x H x S x S
            if (!gate->isControlled()) {
                const auto a = gate->getTargets().at(0U);
                applyGateH(a, nqubits);
                applyGateS(a, nqubits);
                applyGateS(a, nqubits);
                applyGateH(a, nqubits);
                applyGateS(a, nqubits);
                applyGateS(a, nqubits);
            } else { // CY = Sdag(1) x CX(0,1) x S(1)
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
            const auto a = (*gate->getControls().begin()).qubit;
            const auto b = gate->getTargets().at(0);
            if (a == b) {
                util::fatal("Invalid SWAP with same control and target.");
            }
            applyGateCX(a, b, nqubits);
            applyGateCX(b, a, nqubits);
            applyGateCX(a, b, nqubits);
        } break;
        default:
            util::fatal("Unsupported gate encountered: " + std::to_string(gate->getType()));
            break;
    }
}

Tableau Tableau::getDiagonalTableau(std::size_t nQubits) {
    TableauType result{};
    result.resize(nQubits);
    for (std::size_t i = 0; i < nQubits; i++) {
        result[i].resize(2U * nQubits + 1U);
        for (std::size_t j = 0; j < 2 * nQubits; j++) {
            if (i == j - nQubits) {
                result[i][j] = 1;
            }
        }
        result[i][2U * nQubits] = 0;
    }

    return Tableau(result);
}

std::uint64_t Tableau::getBVFrom(int column) const {
    std::uint64_t result = 0UL;
    for (std::size_t j = 0; j < getQubitCount(); ++j) {
        if (tableau[j][column] == 1) {
            result |= (1U << j);
        }
    }
    return result;
}
void Tableau::init(std::size_t nQubits) {
    tableau.clear();
    tableau.resize(nQubits);
    this->tableau = Tableau::getDiagonalTableau(nQubits).tableau;
}

bool Tableau::operator==(const Tableau& other) const {
    if (tableau.size() != other.tableau.size()) {
        return false;
    }
    for (std::size_t i = 0; i < getQubitCount(); ++i) {
        const auto& row1 = tableau[i];
        const auto& row2 = other.tableau[i];
        if (row1.size() != row2.size()) {
            return false;
        }
        for (std::size_t j = 0; j < 2 * getQubitCount() + 1; ++j) {
            const auto& col1 = row1[j];
            const auto& col2 = row2[j];
            if (col1 != col2) {
                return false;
            }
        }
    }
    return true;
}

std::string Tableau::toString() const {
    std::stringstream ss;

    if (empty()) {
        DEBUG() << "Empty tableau";
        return "";
    }
    for (const auto& row: tableau) {
        if (row.size() != back().size()) {
            FATAL() << "Tableau is not rectangular";
            return "";
        }
        for (const auto& s: row) {
            ss << std::to_string(s) << ';';
        }
        ss << std::endl;
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
    auto        rStabilizer = std::regex("([\\+-])([IYZX]+)");
    std::smatch m;
    if (std::regex_search(line, rStabilizer)) {
        // string is a list of stabilizers
        auto iter = line.cbegin();
        while (std::regex_search(iter, line.cend(), m, rStabilizer)) {
            std::string s = m.str(0U);
            RowType     row;

            for (const auto c: s) {
                if (c == 'I' || c == 'Z') {
                    row.push_back(0);
                } else if (c == 'X' || c == 'Y') {
                    row.push_back(1);
                }
            }
            for (const auto c: s) {
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
        import(ss);
    }
}
void Tableau::applyGateH(dd::Qubit target, std::size_t nqubits) {
    for (auto i = 0U; i < nqubits; i++) {
        tableau[i][2U * nqubits] ^= (tableau[i][target] & tableau[i][target + nqubits]);
        std::swap(tableau[i][target], tableau[i][target + nqubits]);
    }
}
void Tableau::applyGateS(dd::Qubit target, std::size_t nqubits) {
    for (auto i = 0U; i < nqubits; i++) {
        tableau[i][2U * nqubits] ^= tableau[i][target] & tableau[i][target + nqubits];
        tableau[i][target + nqubits] ^= tableau[i][target];
    }
}
void Tableau::applyGateCX(dd::Qubit control, dd::Qubit target, std::size_t nqubits) {
    for (auto i = 0U; i < nqubits; i++) {
        const auto xa = tableau[i][target];
        const auto za = tableau[i][target + nqubits];
        const auto xb = tableau[i][control];
        const auto zb = tableau[i][control + nqubits];
        tableau[i][2 * nqubits] ^= (xa & zb) & ((xb ^ za) ^ 1);
        tableau[i][target + nqubits] = za ^ zb;
        tableau[i][control]          = xb ^ xa;
    }
}
void Tableau::applyGateSdag(dd::Qubit target, std::size_t nqubits) {
    applyGateS(target, nqubits);
    applyGateS(target, nqubits);
    applyGateS(target, nqubits);
}
Tableau::Tableau(const qc::QuantumComputation& qc, std::size_t begin, std::size_t end) {
    init(qc.getNqubits());
    std::size_t currentG = 0;
    for (const auto& gate: qc) {
        if (currentG >= begin && (currentG < end)) {
            if (gate->getType() == qc::OpType::Compound) {
                auto* compOp = dynamic_cast<qc::CompoundOperation*>(gate.get());
                auto  cit    = compOp->begin();
                while (cit != compOp->end() && currentG >= begin &&
                       (currentG < end)) {
                    applyGate((*cit));
                    ++cit;
                    ++currentG;
                }
            } else {
                applyGate(gate);
                ++currentG;
            }
        }
    }
}
Tableau::Tableau(std::size_t nQubits) {
    tableau.clear();
    tableau.resize(nQubits);
    this->tableau = Tableau::getDiagonalTableau(nQubits).tableau;
}
