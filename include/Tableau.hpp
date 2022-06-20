//
// Created by Sarah on 20.06.2022.
//

#ifndef QMAP_TABLEAU_HPP
#define QMAP_TABLEAU_HPP

#include "utils/logging.hpp"

#include <ostream>
#include <vector>

using _tableau = std::vector<std::vector<short>>;

class Tableau {
private:
    _tableau tableau;

public:
    [[nodiscard]] Tableau() {}
    [[nodiscard]] Tableau(_tableau tableau1):
        tableau(tableau1) {}

    std::vector<short> operator[](std::size_t index) {
        return tableau.at(index);
    }

    std::vector<short> at(std::size_t index) {
        return tableau.at(index);
    }

    [[nodiscard]] bool empty() const {
        return tableau.empty();
    }

    [[nodiscard]] std::vector<short> back() const {
        return tableau.back();
    }

    [[nodiscard]] std::vector<std::vector<short>>::const_iterator begin() const {
        return tableau.cbegin();
    }
    [[nodiscard]] std::vector<std::vector<short>>::const_iterator end() const {
        return tableau.cend();
    }

    [[nodiscard]] std::string getRepresentation() const {
        std::stringstream result;
        result << *this;
        return result.str();
    }

    friend std::ostream& operator<<(std::ostream& os, const Tableau& dt);
};
std::ostream& operator<<(std::ostream& os, const Tableau& dt) {
    if (dt.empty()) {
        DEBUG() << "Empty tableau";
        return os;
    }
    for (std::size_t i = 0; i < dt.back().size(); ++i) {
        os << i << "\t";
    }
    os << std::endl;
    auto i = 1;
    for (const auto& row: dt) {
        if (row.size() != dt.back().size()) {
            FATAL() << "Tableau is not rectangular";
        }
        os << i++ << "\t";
        for (const auto& s: row)
            os << s << '\t';
        os << std::endl;
    };
}
#endif //QMAP_TABLEAU_HPP
