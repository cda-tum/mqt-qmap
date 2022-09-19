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

    void import(const std::string& filename);
    void import(std::istream& is);

    void init(std::size_t nQubits);

    void populateTableauFrom(unsigned long bv, std::size_t nQubits,
                             int column);

    static void generateTableau(Tableau& tableau, const qc::QuantumComputation& circ, std::size_t begin = 0, std::size_t end = -1);
    static void initTableau(Tableau& tableau, std::size_t nqubits);

    [[nodiscard]] int applyGate(const std::unique_ptr<qc::Operation>& gate);

    [[nodiscard]] bool operator==(const Tableau& other) const;
    [[nodiscard]] bool operator!=(const Tableau& other) const {
        return !(other == *this);
    }

    [[nodiscard]] static Tableau getDiagonalTableau(std::size_t nQubits);
    [[nodiscard]] double         tableauDistance(const Tableau& other, std::size_t nQubits);
    [[nodiscard]] Tableau        embedTableau(std::size_t nQubits);

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
    void fromString(const std::string& str);

    [[nodiscard]] static double tableauDistance(const TableauType& tableau1, const TableauType& tableau2, std::size_t nQubits);

    [[nodiscard]] unsigned long getBVFrom(int column) const;

private:
    void applyGateH(dd::Qubit target, std::size_t nqubits);
    void applyGateS(dd::Qubit target, std::size_t nqubits);
    void applyGateSdag(dd::Qubit target, std::size_t nqubits);
    void applyGateCX(dd::Qubit control, dd::Qubit target, std::size_t nqubits);

};
#endif //QMAP_TABLEAU_HPP
