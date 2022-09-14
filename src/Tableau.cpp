/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
*/

#include "Tableau.hpp"

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

    std::size_t              nQubits = 0;
    std::string              token;
    std::string              line;
    std::vector<std::string> data{};
    char                     delimiter = '|';
    // Try to find out size by reading first line
    if (std::getline(is, line)) {
        if (line.find('|', 0) == std::string::npos)
            delimiter = ',';
        parse_line(line, delimiter, {'\"'}, {'\\'}, data);
        nQubits = static_cast<std::size_t>(std::stoul(data.at(0)));
    }

    tableau.reserve(nQubits);

    while (std::getline(is, line)) {
        tableau.emplace_back();
        tableau.back().reserve(2 * nQubits + 1);
        parse_line(line, delimiter, {'\"'}, {'\\', '\r', '\n', '\t'}, data);
        bool skipFirst = true;
        for (const auto& datum: data) {
            if (skipFirst) {
                skipFirst = false;
                continue;
            }
            if (datum.empty())
                continue;
            tableau.back().emplace_back(static_cast<EntryType>(std::stoul(datum)));
        }
    }
}

void Tableau::populateTableauFrom(unsigned long bv, int nQubits,
                                  int column) {
    for (int j = 0; j < nQubits; ++j) {
        if ((bv & (1U << j)) != 0U) {
            tableau[j][column] = 1;
        }
    }
}

void Tableau::generateTableau(Tableau& tableau, const qc::QuantumComputation& circuit, int begin, int end) {
    initTableau(tableau, circuit.getNqubitsWithoutAncillae());
    int current_g = 0;
    for (const auto& gate: circuit) {
        if (current_g >= begin && (current_g < end || end < 0)) {
            if (gate->getType() == qc::OpType::Compound) {
                auto compOp = dynamic_cast<qc::CompoundOperation*>(gate.get());
                auto cit    = compOp->begin();
                while (cit != compOp->end() && current_g >= begin &&
                       (current_g < end || end < 0)) {
                    tableau.applyGate((*cit));
                    ++cit;
                    ++current_g;
                }
            } else {
                tableau.applyGate(gate);
                ++current_g;
            }
        }
    }
}

void Tableau::initTableau(Tableau& tableau, std::size_t nqubits) {
    tableau.init(nqubits);
}

int Tableau::applyGate(const std::unique_ptr<qc::Operation>& gate) {
    auto nqubits = getQubitCount();
    switch (gate->getType()) {
        case qc::OpType::H: // HADAMARD
        {
            if (gate->isControlled()) {
                util::fatal("Expected single-qubit gate");
            }
            const auto a = gate->getTargets().at(0U);
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= (tableau[i][a] & tableau[i][a + nqubits]);
                std::swap(tableau[i][a], tableau[i][a + nqubits]);
            }
            return 1U;
        }
        case qc::OpType::S: // PHASE
        {
            if (gate->isControlled()) {
                util::fatal("Expected single-qubit gate");
            }
            const auto a = gate->getTargets().at(0U);
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= tableau[i][a] & tableau[i][a + nqubits];
                tableau[i][a + nqubits] ^= tableau[i][a];
            }
            return 1U;
        }
        case qc::OpType::X: // CNOT
        {
            if (gate->getNcontrols() != 1U) { // NOT = H x S x S x H
                const auto a = gate->getTargets().at(0U);
                for (auto i = 0U; i < nqubits; i++) {
                    tableau[i][2U * nqubits] ^= tableau[i][a] & tableau[i][a + nqubits];
                    std::swap(tableau[i][a], tableau[i][a + nqubits]);
                }
                for (auto i = 0U; i < nqubits; i++) {
                    tableau[i][2U * nqubits] ^= tableau[i][a] & tableau[i][a + nqubits];
                    tableau[i][a + nqubits] ^= tableau[i][a];
                }
                for (auto i = 0U; i < nqubits; i++) {
                    tableau[i][2U * nqubits] ^= tableau[i][a] & tableau[i][a + nqubits];
                    tableau[i][a + nqubits] ^= tableau[i][a];
                }
                for (auto i = 0U; i < nqubits; i++) {
                    tableau[i][2U * nqubits] ^= tableau[i][a] & tableau[i][a + nqubits];
                    std::swap(tableau[i][a], tableau[i][a + nqubits]);
                }
                return 4U;
            } else {
                const auto a = (*gate->getControls().begin()).qubit;
                const auto b = gate->getTargets().at(0);
                if (a == b) {
                    util::fatal("Invalid CNOT with same control and target.");
                }
                for (auto i = 0U; i < nqubits; i++) {
                    const auto xa = tableau[i][a];
                    const auto za = tableau[i][a + nqubits];
                    const auto xb = tableau[i][b];
                    const auto zb = tableau[i][b + nqubits];
                    tableau[i][2 * nqubits] ^= (xa & zb) & ((xb ^ za) ^ 1);
                    tableau[i][a + nqubits] = za ^ zb;
                    tableau[i][b]           = xb ^ xa;
                }
                return 1U;
            }
        }

        case qc::OpType::Sdag: { // Sdag  = S x S x S
            if (gate->isControlled()) {
                util::fatal("Expected single-qubit gate");
            }
            const auto a = gate->getTargets().at(0U);
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= tableau[i][a] & tableau[i][a + nqubits];
                tableau[i][a + nqubits] ^= tableau[i][a];
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= tableau[i][a] & tableau[i][a + nqubits];
                tableau[i][a + nqubits] ^= tableau[i][a];
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= tableau[i][a] & tableau[i][a + nqubits];
                tableau[i][a + nqubits] ^= tableau[i][a];
            }
            return 3U;
        }
        case qc::OpType::Z: { // Z = S x S
            if (gate->isControlled()) {
                util::fatal("Expected single-qubit gate");
            }
            const auto a = gate->getTargets().at(0U);
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= tableau[i][a] & tableau[i][a + nqubits];
                tableau[i][a + nqubits] ^= tableau[i][a];
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= tableau[i][a] & tableau[i][a + nqubits];
                tableau[i][a + nqubits] ^= tableau[i][a];
            }
            return 2U;
        }
        case qc::OpType::Y: { // Y = H x S x S x H x S x S
            if (gate->isControlled()) {
                util::fatal("Expected single-qubit gate");
            }
            const auto a = gate->getTargets().at(0U);
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= tableau[i][a] & tableau[i][a + nqubits];
                std::swap(tableau[i][a], tableau[i][a + nqubits]);
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= tableau[i][a] & tableau[i][a + nqubits];
                tableau[i][a + nqubits] ^= tableau[i][a];
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= tableau[i][a] & tableau[i][a + nqubits];
                tableau[i][a + nqubits] ^= tableau[i][a];
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= tableau[i][a] & tableau[i][a + nqubits];
                std::swap(tableau[i][a], tableau[i][a + nqubits]);
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= tableau[i][a] & tableau[i][a + nqubits];
                tableau[i][a + nqubits] ^= tableau[i][a];
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= tableau[i][a] & tableau[i][a + nqubits];
                tableau[i][a + nqubits] ^= tableau[i][a];
            }
            return 6U;
        }
        default:
            util::fatal("Unsupported gate encountered: " + std::to_string(gate->getType()));
            break;
    }
    return 0U;
}

Tableau Tableau::getDiagonalTableau(int nQubits) {
    TableauType result{};
    result.resize(nQubits);
    for (auto i = 0; i < nQubits; i++) {
        result[i].resize(2U * nQubits + 1U);
        for (auto j = 0; j < 2 * nQubits; j++) {
            if (i == j - nQubits) {
                result[i][j] = 1;
            } else {
                result[i][j] = 0;
            }
        }
        result[i][2U * nQubits] = 0;
    }

    return Tableau(result);
}

double Tableau::tableauDistance(Tableau other, int nQubits) {
    double result = 0.0;
    if (tableau.size() != other.tableau.size()) {
        result = std::numeric_limits<double>::max();
    } else {
        for (int i = 0; i < nQubits; ++i) {
            auto first  = std::find_if(tableau[i].begin(), tableau[i].end(),
                                       [](const auto x) { return x == 1; });
            auto last   = std::find_if(tableau[i].rbegin(), tableau[i].rend(),
                                       [](const auto x) { return x == 1; });
            auto first2 = std::find_if(other.tableau[i].begin(), other.tableau[i].end(),
                                       [](const auto x) { return x == 1; });
            auto last2  = std::find_if(other.tableau[i].rbegin(), other.tableau[i].rend(),
                                       [](const auto x) { return x == 1; });
            auto d1     = std::distance(tableau[i].begin(), first);
            auto d2     = std::distance(tableau[i].rbegin(), last);
            auto d3     = std::distance(other.tableau[i].begin(), first2);
            auto d4     = std::distance(other.tableau[i].rbegin(), last2);
            result += static_cast<double>(std::abs(d1 - d3)) / 2.0 + static_cast<double>(std::abs(d2 - d4)) / 2.0;
        }
    }
    return result;
}

Tableau Tableau::embedTableau(int nQubits) {
    TableauType      result{};
    auto             diagonal = getDiagonalTableau(nQubits);
    std::vector<int> indices{};
    auto             m = tableau.size();

    for (unsigned long i = 0U; i < static_cast<unsigned long>(nQubits); i++) {
        indices.push_back(i < m ? 0 : 1);
    }
    do {
        TableauType intermediate_result{};
        int         i = 0;
        intermediate_result.resize(nQubits);
        for (auto k = 0; k < nQubits; k++) {
            intermediate_result[k].resize(2 * nQubits + 1);
            int n = 0;
            for (auto j = 0; j < 2 * nQubits; j++) {
                if (indices[k] == 1 || (j < nQubits && indices[j] == 1) ||
                    (j >= nQubits && indices[j - nQubits] == 1)) {
                    intermediate_result[k][j] = diagonal[k][j];
                } else {
                    intermediate_result[k][j] = tableau[i][n];
                    n++;
                }
            }
            if (indices[k] == 1) {
                intermediate_result[k][2 * nQubits] = 0;
            } else {
                intermediate_result[k][2 * nQubits] = tableau[i][2 * getQubitCount()];
                i++;
            }
        }
        if (Tableau::tableauDistance(diagonal.tableau, intermediate_result, nQubits) <
            Tableau::tableauDistance(diagonal.tableau, result, nQubits)) {
            result = intermediate_result;
        }
    } while (std::next_permutation(indices.begin(), indices.end()));
    return Tableau(result);
}
double Tableau::tableauDistance(const TableauType& tableau1, const TableauType& tableau2, std::size_t nQubits) {
    double result = 0.0;
    if (tableau1.size() != tableau2.size()) {
        result = std::numeric_limits<double>::max();
    } else {
        for (std::size_t i = 0U; i < nQubits; ++i) {
            auto first  = std::find_if(tableau1[i].begin(), tableau1[i].end(),
                                       [](const auto x) { return x == 1; });
            auto last   = std::find_if(tableau1[i].rbegin(), tableau1[i].rend(),
                                       [](const auto x) { return x == 1; });
            auto first2 = std::find_if(tableau2[i].begin(), tableau2[i].end(),
                                       [](const auto x) { return x == 1; });
            auto last2  = std::find_if(tableau2[i].rbegin(), tableau2[i].rend(),
                                       [](const auto x) { return x == 1; });
            auto d1     = std::distance(tableau1[i].begin(), first);
            auto d2     = std::distance(tableau1[i].rbegin(), last);
            auto d3     = std::distance(tableau2[i].begin(), first2);
            auto d4     = std::distance(tableau2[i].rbegin(), last2);
            result += static_cast<double>(std::abs(d1 - d3)) / 2.0 + static_cast<double>(std::abs(d2 - d4)) / 2.0;
        }
    }
    return result;
}

unsigned long Tableau::getBVFrom(int column) const {
    unsigned long result = 0UL;
    for (std::size_t j = 0; j < getQubitCount(); ++j) {
        if (tableau[j][column] == 1) {
            result |= (1U << j);
        }
    }
    return result;
}
void Tableau::init(std::size_t nQubits) {
    tableau.resize(nQubits);
    for (auto i = 0U; i < nQubits; i++) {
        tableau[i].resize(2U * nQubits + 1U);
        for (auto j = 0U; j < 2U * nQubits; j++) {
            if (i == j - nQubits) {
                tableau[i][j] = 1;
            } else {
                tableau[i][j] = 0;
            }
        }
        tableau[i][2U * nQubits] = 0;
    }
}

std::string Tableau::getStrRepresentation() const {
    std::stringstream out;
    out << *this;
    return out.str();
}
void Tableau::importString(const std::string& tableauRepr) {
    std::stringstream in(tableauRepr);
    in >> *this;
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

std::string Tableau::getRepresentation() const {
    std::stringstream result;
    result << *this;
    return result.str();
}
