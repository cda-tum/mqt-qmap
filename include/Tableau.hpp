/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
*/

#ifndef QMAP_TABLEAU_HPP
#define QMAP_TABLEAU_HPP

#include "QuantumComputation.hpp"
#include "utils/logging.hpp"

#include <fstream>
#include <limits>
#include <ostream>
#include <vector>

using innerTableau = std::vector<std::vector<int32_t>>;

class Tableau {
private:
    innerTableau tableau;

public:
    [[nodiscard]] Tableau() = default;
    [[nodiscard]] explicit Tableau(innerTableau& inner):
        tableau(inner) {}

    Tableau(const Tableau& other) {
        this->tableau = other.tableau;
    }
    Tableau(Tableau& other) {
        this->tableau = other.tableau;
    }

    Tableau operator=(Tableau other) {
        tableau = other.tableau;
        return *this;
    }

    std::vector<int32_t> operator[](std::size_t index) {
        return tableau[index];
    }

    std::vector<int32_t> at(std::size_t index) {
        return tableau.at(index);
    }

    [[nodiscard]] size_t getQubitCount() const {
        return tableau.size();
    }

    void resize(std::size_t size) {
        tableau.resize(size);
    }

    void clear() {
        tableau.clear();
    }

    [[nodiscard]] bool empty() const {
        return tableau.empty();
    }

    [[nodiscard]] std::vector<int32_t> back() const {
        return tableau.back();
    }

    [[nodiscard]] std::vector<std::vector<int32_t>>::const_iterator begin() const {
        return tableau.cbegin();
    }
    [[nodiscard]] std::vector<std::vector<int32_t>>::const_iterator end() const {
        return tableau.cend();
    }

    void dump(const std::string& filename) const;

    void dump(std::ostream& of) const;

    void import(const std::string& filename);
    void import(std::istream& is);

    [[nodiscard]] std::string getRepresentation() const {
        std::stringstream result;
        result << *this;
        return result.str();
    }

    void init(size_t nQubits);

    void populateTableauFrom(unsigned long bv, int nQubits,
                             int column);

    static void generateTableau(Tableau& tableau, qc::QuantumComputation& circ, int begin = 0, int end = -1);
    static void initTableau(Tableau& tableau, size_t nqubits);

    int applyGate(std::unique_ptr<qc::Operation>& gate);

    bool operator==(const Tableau& other) const {
        if (tableau.size() != other.tableau.size()) {
            return false;
        }
        for (size_t i = 0; i < getQubitCount(); ++i) {
            const auto& row1 = tableau[i];
            const auto& row2 = other.tableau[i];
            if (row1.size() != row2.size()) {
                return false;
            }
            for (size_t j = 0; j < 2 * getQubitCount() + 1; ++j) {
                const auto& col1 = row1[j];
                const auto& col2 = row2[j];
                if (col1 != col2) {
                    return false;
                }
            }
        }
        return true;
    }

    static Tableau       getDiagonalTableau(int nQubits);
    double               tableauDistance(Tableau other, int nQubits);
    Tableau              embedTableau(int nQubits);
    friend std::ostream& operator<<(std::ostream& os, const Tableau& dt);

    static double tableauDistance(innerTableau tableau1, innerTableau tableau2, int nQubits);

    unsigned long getBVFrom(int column) const;
};
#endif //QMAP_TABLEAU_HPP
