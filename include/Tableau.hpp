/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
*/

#ifndef QMAP_TABLEAU_HPP
#define QMAP_TABLEAU_HPP

#include "QuantumComputation.hpp"
#include "utils/logging.hpp"

#include <cstdint>
#include <fstream>
#include <limits>
#include <ostream>
#include <utility>
#include <vector>

class Tableau {
    using EntryType   = std::int32_t;
    using RowType     = std::vector<EntryType>;
    using TableauType = std::vector<RowType>;
    TableauType tableau;

public:
    [[nodiscard]] Tableau() = default;
    [[nodiscard]] explicit Tableau(TableauType tableau):
        tableau(std::move(tableau)) {}

    [[nodiscard]] RowType operator[](std::size_t index) {
        return tableau[index];
    }
    [[nodiscard]] RowType operator[](std::size_t index) const {
        return tableau[index];
    }

    [[nodiscard]] RowType at(std::size_t index) {
        return tableau.at(index);
    }

    [[nodiscard]] auto getQubitCount() const {
        return tableau.size();
    }

    inline void resize(std::size_t size) {
        tableau.resize(size);
    }

    inline void clear() {
        tableau.clear();
    }

    [[nodiscard]] bool empty() const {
        return tableau.empty();
    }

    [[nodiscard]] auto back() const {
        return tableau.back();
    }

    [[nodiscard]] auto begin() const {
        return tableau.begin();
    }
    [[nodiscard]] auto end() const {
        return tableau.end();
    }

    void dump(const std::string& filename) const;

    void dump(std::ostream& of) const;

    void importString(const std::string& tableauRepr);

    void import(const std::string& filename);
    void import(std::istream& is);

    [[nodiscard]] std::string getRepresentation() const;

    void init(std::size_t nQubits);

    void populateTableauFrom(unsigned long bv, int nQubits,
                             int column);

    static void generateTableau(Tableau& tableau, const qc::QuantumComputation& circ, int begin = 0, int end = -1);
    static void initTableau(Tableau& tableau, std::size_t nqubits);

    [[nodiscard]] int applyGate(const std::unique_ptr<qc::Operation>& gate);

    [[nodiscard]] bool operator==(const Tableau& other) const;

    [[nodiscard]] static Tableau getDiagonalTableau(int nQubits);
    [[nodiscard]] double         tableauDistance(Tableau other, int nQubits);
    [[nodiscard]] Tableau        embedTableau(int nQubits);
    friend std::ostream&         operator<<(std::ostream& os, const Tableau& dt) {
                std::size_t nQubits = dt.getQubitCount();
                if (dt.empty()) {
                    DEBUG() << "Empty tableau";
                    return os;
        }
                os << nQubits << '|';
                for (std::size_t i = 1U; i < dt.back().size(); ++i) {
                    os << i << '|';
        }
                os << "R|";
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
        }
                return os;
    }
    friend std::istream& operator>>(std::istream& is, Tableau& dt) {
        std::string line;
        std::getline(is, line);
        if (line.empty()) {
            return is;
        }
        auto        r_stabilizer = std::regex("([\\+-])([IYZX]+)");
        std::smatch m;
        if (line.find("Destabilizer") != std::string::npos) {
            auto iter = line.cbegin();
            while (std::regex_search(iter, line.cend(), m, r_stabilizer)) {
                std::string s = m.str(0U);
                RowType     row;

                for (const auto c: s) {
                    if (c == 'I' || c == 'X') {
                        row.push_back(0);
                    } else if (c == 'Y' || c == 'Z') {
                        row.push_back(1);
                    }
                }
                for (auto c: s) {
                    if (c == 'I' || c == 'Z') {
                        row.push_back(0);
                    } else if (c == 'X' || c == 'Y') {
                        row.push_back(1);
                    }
                }
                if (s[0U] == '-') {
                    row.push_back(1);
                } else {
                    row.push_back(0);
                }
                dt.tableau.push_back(row);
                iter = m[0].second;
            }
        } else {
            auto iter = line.cbegin();
            while (std::regex_search(iter, line.cend(), m, r_stabilizer)) {
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
                dt.tableau.push_back(row);
                iter = m[0].second;
            }
        }
        return is;
    }

    [[nodiscard]] static double tableauDistance(const TableauType& tableau1, const TableauType& tableau2, std::size_t nQubits);

    [[nodiscard]] unsigned long getBVFrom(int column) const;

    [[nodiscard]] std::string getStrRepresentation() const;
};
#endif //QMAP_TABLEAU_HPP
