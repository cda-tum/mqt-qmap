/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
*/

#include "Tableau.hpp"

void Tableau::populateTableauFrom(unsigned long bv, int nQubits,
                         int                              column) {
    for (int j = 0; j < nQubits; ++j) {
        if ((bv & (1U << j)) != 0U) {
            tableau[j][column] = 1;
        }
    }
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

 double Tableau::tableauDistance(Tableau other,int nQubits) {
    double result = 0.0;
    if (tableau.size() != other.tableau.size()) {
        result = std::numeric_limits<double>::max();
    } else {
        for (int i = 0; i < nQubits; ++i) {
            auto first  = std::find_if(tableau[i].begin(), tableau[i].end(),
                                       [](short x) { return x == 1; });
            auto last   = std::find_if(tableau[i].rbegin(), tableau[i].rend(),
                                       [](short x) { return x == 1; });
            auto first2 = std::find_if(other.tableau[i].begin(), other.tableau[i].end(),
                                       [](short x) { return x == 1; });
            auto last2  = std::find_if(other.tableau[i].rbegin(), other.tableau[i].rend(),
                                       [](short x) { return x == 1; });
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
    innerTableau                    result{};
    auto                            diagonal = getDiagonalTableau(nQubits);
    std::vector<int>                indices{};
    auto                            m = tableau.size();

    for (unsigned long i = 0U; i < static_cast<unsigned long>(nQubits); i++) {
        indices.push_back(i < m ? 0 : 1);
    }
    do {
        innerTableau                    intermediate_result{};
        int                             i = 0;
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
                                       [](short x) { return x == 1; });
            auto last   = std::find_if(tableau1[i].rbegin(), tableau1[i].rend(),
                                       [](short x) { return x == 1; });
            auto first2 = std::find_if(tableau2[i].begin(), tableau2[i].end(),
                                       [](short x) { return x == 1; });
            auto last2  = std::find_if(tableau2[i].rbegin(), tableau2[i].rend(),
                                       [](short x) { return x == 1; });
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
     }
     return os;
 }