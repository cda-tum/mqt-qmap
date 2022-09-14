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

    [[nodiscard]] Tableau(const Tableau& other) {
        this->tableau = other.tableau;
    }
    [[nodiscard]] Tableau(Tableau& other) {
        this->tableau = other.tableau;
    }

    [[nodiscard]] Tableau& operator=(const Tableau& other) {
        tableau = other.tableau;
        return *this;
    }

    [[nodiscard]] std::vector<int32_t> operator[](std::size_t index) {
        return tableau[index];
    }
    [[nodiscard]] std::vector<int32_t> operator[](std::size_t index) const {
        return tableau[index];
    }

    [[nodiscard]] std::vector<int32_t> at(std::size_t index) {
        return tableau.at(index);
    }

    [[nodiscard]] size_t getQubitCount() const {
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

    void importString(const std::string& tableauRepr);

    void import(const std::string& filename);
    void import(std::istream& is);

    [[nodiscard]] std::string getRepresentation() const;

    void init(size_t nQubits);

    void populateTableauFrom(unsigned long bv, int nQubits,
                             int column);

    static void generateTableau(Tableau& tableau, qc::QuantumComputation& circ, int begin = 0, int end = -1);
    static void initTableau(Tableau& tableau, size_t nqubits);

    [[nodiscard]] int applyGate(std::unique_ptr<qc::Operation>& gate);

    [[nodiscard]] bool operator==(const Tableau& other) const;

    [[nodiscard]] static Tableau getDiagonalTableau(int nQubits);
    [[nodiscard]] double         tableauDistance(Tableau other, int nQubits);
    [[nodiscard]] Tableau        embedTableau(int nQubits);
    friend std::ostream&         operator<<(std::ostream& os, const Tableau& dt);
    friend std::istream&         operator>>(std::istream& is, Tableau& dt);

    [[nodiscard]] static double tableauDistance(innerTableau tableau1, innerTableau tableau2, int nQubits);

    [[nodiscard]] unsigned long getBVFrom(int column) const;

    [[nodiscard]] std::string getStrRepresentation() const;
};
#endif //QMAP_TABLEAU_HPP
