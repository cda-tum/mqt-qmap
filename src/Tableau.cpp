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

    size_t                   nQubits = 0;
    std::string              token;
    std::string              line;
    std::vector<std::string> data{};
    char                     delimiter = '|';
    // Try to find out size by reading first line
    if (std::getline(is, line)) {
        if (line.find('|', 0) == std::string::npos)
            delimiter = ',';
        parse_line(line, delimiter, {'\"'}, {'\\'}, data);
        nQubits = static_cast<size_t>(std::stoul(data.at(0)));
    }

    tableau.reserve(nQubits);

    while (std::getline(is, line)) {
        if (line.find('-', 0) != std::string::npos)
            continue;
        tableau.emplace_back();
        tableau.back().reserve(2 * nQubits + 1);
        parse_line(line, delimiter, {'\"'}, {'\\'}, data);
        bool skipFirst = true;
        for (const auto& datum: data) {
            if (skipFirst) {
                skipFirst = false;
                continue;
            }
            if (datum == "")
                continue;
            tableau.back().emplace_back(static_cast<int32_t>(std::stoul(datum)));
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

void Tableau::generateTableau(Tableau& tableau, qc::QuantumComputation& circuit, int begin, int end) {
    initTableau(tableau, circuit.getNqubitsWithoutAncillae());
    int current_g = 0;
    for (auto& gate: circuit) {
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

void Tableau::initTableau(Tableau& tableau, size_t nqubits) {
    tableau.init(nqubits);
}

int Tableau::applyGate(std::unique_ptr<qc::Operation>& gate) {
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
                tableau[i][2U * nqubits] ^= static_cast<int32_t>(tableau[i][a] & tableau[i][a + nqubits]);
                tableau[i][a + nqubits] ^= static_cast<int32_t>(tableau[i][a]);
            }
            return 1U;
        }
        case qc::OpType::X: // CNOT
        {
            if (gate->getNcontrols() != 1U) { // NOT = H x S x S x H
                const auto a = gate->getTargets().at(0U);
                for (auto i = 0U; i < nqubits; i++) {
                    tableau[i][2U * nqubits] ^= static_cast<int32_t>(tableau[i][a] & tableau[i][a + nqubits]);
                    std::swap(tableau[i][a], tableau[i][a + nqubits]);
                }
                for (auto i = 0U; i < nqubits; i++) {
                    tableau[i][2U * nqubits] ^= static_cast<int32_t>(tableau[i][a] & tableau[i][a + nqubits]);
                    tableau[i][a + nqubits] ^= static_cast<int32_t>(tableau[i][a]);
                }
                for (auto i = 0U; i < nqubits; i++) {
                    tableau[i][2U * nqubits] ^= static_cast<int32_t>(tableau[i][a] & tableau[i][a + nqubits]);
                    tableau[i][a + nqubits] ^= static_cast<int32_t>(tableau[i][a]);
                }
                for (auto i = 0U; i < nqubits; i++) {
                    tableau[i][2U * nqubits] ^= static_cast<int32_t>(tableau[i][a] & tableau[i][a + nqubits]);
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
                    const auto xa = static_cast<int32_t>(tableau[i][a]);
                    const auto za = static_cast<int32_t>(tableau[i][a + nqubits]);
                    const auto xb = static_cast<int32_t>(tableau[i][b]);
                    const auto zb = static_cast<int32_t>(tableau[i][b + nqubits]);
                    tableau[i][2 * nqubits] ^= static_cast<int32_t>((xa & zb) & ((xb ^ za) ^ 1));
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
                tableau[i][2U * nqubits] ^= static_cast<int32_t>(tableau[i][a] & tableau[i][a + nqubits]);
                tableau[i][a + nqubits] ^= static_cast<int32_t>(tableau[i][a]);
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= static_cast<int32_t>(tableau[i][a] & tableau[i][a + nqubits]);
                tableau[i][a + nqubits] ^= static_cast<int32_t>(tableau[i][a]);
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= static_cast<int32_t>(tableau[i][a] & tableau[i][a + nqubits]);
                tableau[i][a + nqubits] ^= static_cast<int32_t>(tableau[i][a]);
            }
            return 3U;
        }
        case qc::OpType::Z: { // Z = S x S
            if (gate->isControlled()) {
                util::fatal("Expected single-qubit gate");
            }
            const auto a = gate->getTargets().at(0U);
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= static_cast<int32_t>(tableau[i][a] & tableau[i][a + nqubits]);
                tableau[i][a + nqubits] ^= static_cast<int32_t>(tableau[i][a]);
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= static_cast<int32_t>(tableau[i][a] & tableau[i][a + nqubits]);
                tableau[i][a + nqubits] ^= static_cast<int32_t>(tableau[i][a]);
            }
            return 2U;
        }
        case qc::OpType::Y: { // Y = H x S x S x H x S x S
            if (gate->isControlled()) {
                util::fatal("Expected single-qubit gate");
            }
            const auto a = gate->getTargets().at(0U);
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= static_cast<int32_t>(tableau[i][a] & tableau[i][a + nqubits]);
                std::swap(tableau[i][a], tableau[i][a + nqubits]);
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= static_cast<int32_t>(tableau[i][a] & tableau[i][a + nqubits]);
                tableau[i][a + nqubits] ^= static_cast<int32_t>(tableau[i][a]);
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= static_cast<int32_t>(tableau[i][a] & tableau[i][a + nqubits]);
                tableau[i][a + nqubits] ^= static_cast<int32_t>(tableau[i][a]);
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= static_cast<int32_t>(tableau[i][a] & tableau[i][a + nqubits]);
                std::swap(tableau[i][a], tableau[i][a + nqubits]);
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= static_cast<int32_t>(tableau[i][a] & tableau[i][a + nqubits]);
                tableau[i][a + nqubits] ^= static_cast<int32_t>(tableau[i][a]);
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= static_cast<int32_t>(tableau[i][a] & tableau[i][a + nqubits]);
                tableau[i][a + nqubits] ^= static_cast<int32_t>(tableau[i][a]);
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
    innerTableau result{};
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
                                       [](int32_t x) { return x == 1; });
            auto last   = std::find_if(tableau[i].rbegin(), tableau[i].rend(),
                                       [](int32_t x) { return x == 1; });
            auto first2 = std::find_if(other.tableau[i].begin(), other.tableau[i].end(),
                                       [](int32_t x) { return x == 1; });
            auto last2  = std::find_if(other.tableau[i].rbegin(), other.tableau[i].rend(),
                                       [](int32_t x) { return x == 1; });
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
    innerTableau     result{};
    auto             diagonal = getDiagonalTableau(nQubits);
    std::vector<int> indices{};
    auto             m = tableau.size();

    for (unsigned long i = 0U; i < static_cast<unsigned long>(nQubits); i++) {
        indices.push_back(i < m ? 0 : 1);
    }
    do {
        innerTableau intermediate_result{};
        int          i = 0;
        intermediate_result.resize(nQubits);
        for (auto k = 0; k < nQubits; k++) {
            intermediate_result[k].resize(2 * nQubits + 1);
            int n = 0;
            for (auto j = 0; j < 2 * nQubits; j++) {
                DEBUG() << "i = " << i << " n = " << n << " j = " << j << " k = " << k
                        << std::endl;
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
                intermediate_result[k][2 * nQubits] = tableau[i][2 * nQubits];
                i++;
            }
            // intermediate_result[2 * nQubits] = tableau[2 * nQubits];
        }
        if (Tableau::tableauDistance(diagonal.tableau, intermediate_result, nQubits) <
            Tableau::tableauDistance(diagonal.tableau, result, nQubits)) {
            result = intermediate_result;
        }
    } while (std::next_permutation(indices.begin(), indices.end()));
    return Tableau(result);
}
double Tableau::tableauDistance(innerTableau tableau1, innerTableau tableau2, int nQubits) {
    double result = 0.0;
    if (tableau1.size() != tableau2.size()) {
        result = std::numeric_limits<double>::max();
    } else {
        for (int i = 0; i < nQubits; ++i) {
            auto first  = std::find_if(tableau1[i].begin(), tableau1[i].end(),
                                       [](int32_t x) { return x == 1; });
            auto last   = std::find_if(tableau1[i].rbegin(), tableau1[i].rend(),
                                       [](int32_t x) { return x == 1; });
            auto first2 = std::find_if(tableau2[i].begin(), tableau2[i].end(),
                                       [](int32_t x) { return x == 1; });
            auto last2  = std::find_if(tableau2[i].rbegin(), tableau2[i].rend(),
                                       [](int32_t x) { return x == 1; });
            auto d1     = std::distance(tableau1[i].begin(), first);
            auto d2     = std::distance(tableau1[i].rbegin(), last);
            auto d3     = std::distance(tableau2[i].begin(), first2);
            auto d4     = std::distance(tableau2[i].rbegin(), last2);
            result += static_cast<double>(std::abs(d1 - d3)) / 2.0 + static_cast<double>(std::abs(d2 - d4)) / 2.0;
        }
    }
    // DEBUG() << "Tableau distance: " << result << std::endl;
    return result;
}

unsigned long Tableau::getBVFrom(int column) const {
    unsigned long result = 0UL;
    for (size_t j = 0; j < getQubitCount(); ++j) {
        if (tableau[j][column] == 1) {
            result |= (1U << j);
        }
    }
    return result;
}
void Tableau::init(size_t nQubits) {
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

std::ostream& operator<<(std::ostream& os, const Tableau& dt) {
    size_t nQubits = dt.getQubitCount();
    if (dt.empty()) {
        DEBUG() << "Empty tableau";
        return os;
    }
    os << nQubits << '|';
    for (std::size_t i = 1; i < dt.back().size(); ++i) {
        os << i << '|';
    }
    //    os << std::string(nQubits * 2, '-') << std::endl;
    os << std::endl;
    auto i = 1;
    for (const auto& row: dt) {
        if (row.size() != dt.back().size()) {
            FATAL() << "Tableau is not rectangular";
            return os;
        }
        os << i++ << "|";
        for (const auto& s: row)
            os << s << '|';
        os << std::endl;
        //        os << std::string(nQubits * 2, '-') << std::endl;
    }
    return os;
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
std::istream& operator>>(std::istream& is, Tableau& dt) {
    std::string line;
    std::getline(is, line);
    if (line.empty()) {
        return is;
    }
    std::regex  r_stabilizer = std::regex("([\\+-])([IYZX]+)");
    std::smatch m;
    if (line.find("Destabilizer") != std::string::npos) {
        std::string::const_iterator iter = line.cbegin();
        while (std::regex_search(iter, line.cend(), m, r_stabilizer)) {
            std::string          s = m.str(0);
            std::vector<int32_t> row;

            for (auto c: s) {
                if (c == 'I')
                    row.push_back(0);
                else if (c == 'X')
                    row.push_back(0);
                else if (c == 'Y')
                    row.push_back(1);
                else if (c == 'Z')
                    row.push_back(1);
            }
            for (auto c: s) {
                if (c == 'I')
                    row.push_back(0);
                else if (c == 'X')
                    row.push_back(1);
                else if (c == 'Y')
                    row.push_back(1);
                else if (c == 'Z')
                    row.push_back(0);
            }
            if (s[0] == '-') {
                row.push_back(1);
            } else {
                row.push_back(0);
            }
            //            std::cout << "Destabilizer:" << s << std::endl;
            dt.tableau.push_back(row);
            iter = m[0].second;
        }
    } else {
        std::string::const_iterator iter = line.cbegin();
        while (std::regex_search(iter, line.cend(), m, r_stabilizer)) {
            std::string          s = m.str(0);
            std::vector<int32_t> row;

            for (auto c: s) {
                if (c == 'I')
                    row.push_back(0);
                else if (c == 'X')
                    row.push_back(1);
                else if (c == 'Y')
                    row.push_back(1);
                else if (c == 'Z')
                    row.push_back(0);
            }
            for (auto c: s) {
                if (c == 'I')
                    row.push_back(0);
                else if (c == 'X')
                    row.push_back(0);
                else if (c == 'Y')
                    row.push_back(1);
                else if (c == 'Z')
                    row.push_back(1);
            }
            if (s[0] == '-') {
                row.push_back(1);
            } else {
                row.push_back(0);
            }
            //            std::cout << "Stabilizer:" << s << std::endl;
            dt.tableau.push_back(row);
            iter = m[0].second;
        }
    }
    std::cout << "Tableau: " << dt << std::endl;
    return is;
}
