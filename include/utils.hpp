/*
 * This file is part of the MQT QMAP library which is released under the MIT license.
 * See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
 */

#ifndef QMAP_UTILS_HPP
#define QMAP_UTILS_HPP

#include "QuantumComputation.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <operations/Operation.hpp>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using Matrix      = std::vector<std::vector<double>>;
using Edge        = std::pair<unsigned short, unsigned short>;
using CouplingMap = std::set<Edge>;

struct Exchange {
    Exchange(unsigned short first, unsigned short second, qc::OpType op):
        first(first), second(second), middle_ancilla(std::numeric_limits<decltype(middle_ancilla)>::max()), op(op) {}
    Exchange(unsigned short first, unsigned short second, unsigned short middle_anc, qc::OpType op):
        first(first), second(second), middle_ancilla(middle_anc), op(op) {}
    unsigned short first;
    unsigned short second;
    unsigned short middle_ancilla;
    qc::OpType     op;
};

class QMAPException: public std::runtime_error {
    std::string msg;

public:
    explicit QMAPException(std::string msg):
        std::runtime_error("QMAP Exception"), msg(std::move(msg)) {}

    [[nodiscard]] const char* what() const noexcept override {
        return msg.c_str();
    }
};

/// Computes n! recursively
/// \param n integer to compute factorial of
/// \return n!
static inline unsigned long factorial(unsigned long n) {
    if (n == 1)
        return 1;
    else
        return n * factorial(n - 1);
}

class Dijkstra {
    static constexpr bool VERBOSE = false;

public:
    struct Node {
        bool   contains_correct_edge = false;
        bool   visited               = false;
        int    pos                   = -1;
        double cost                  = -1.;
    };

    static void build_table(unsigned short n, const std::set<Edge>& graph, Matrix& distanceTable, const std::function<double(const Node&)>& cost);

protected:
    static void dijkstra(const CouplingMap& couplingMap, std::vector<Node>& nodes, unsigned short start);
};

inline bool operator<(const Dijkstra::Node& x, const Dijkstra::Node& y) {
    if (x.cost != y.cost) {
        return x.cost < y.cost;
    }
    return x.contains_correct_edge && !y.contains_correct_edge;
}

/// Iterating routine through all combinations
/// \tparam Iterator iterator type
/// \param first iterator to beginning
/// \param k current iterator
/// \param last iterator to end
/// \return true if another combination was found
template<typename Iterator>
bool next_combination(Iterator first, Iterator k, Iterator last) {
    /* Credits: Thomas Draper */
    if ((first == last) || (first == k) || (last == k))
        return false;
    Iterator itr1 = first;
    Iterator itr2 = last;
    ++itr1;
    if (last == itr1)
        return false;
    itr1 = last;
    --itr1;
    itr1 = k;
    --itr2;
    while (first != itr1) {
        if (*--itr1 < *itr2) {
            Iterator j = k;
            while (!(*itr1 < *j)) ++j;
            std::iter_swap(itr1, j);
            ++itr1;
            ++j;
            itr2 = k;
            std::rotate(itr1, j, last);
            while (last != j) {
                ++j;
                ++itr2;
            }
            std::rotate(k, itr2, last);
            return true;
        }
    }
    std::rotate(first, k, last);
    return false;
}

/// Create a string representation of a given permutation
/// \param pi permutation
/// \return string representation of pi
std::string printPi(std::vector<unsigned short>& pi);

/// Simple depth-first-search implementation used to check whether a given subset of qubits is
/// connected on the given architecture
/// \param current index of current qubit
/// \param visited visited qubits
/// \param cm coupling map of architecture
void dfs(unsigned short current, std::set<unsigned short>& visited, const CouplingMap& rcm);

using filter_function = std::function<bool(const std::set<unsigned short>&)>;
std::vector<std::set<unsigned short>> subsets(const std::set<unsigned short>& input, int size, filter_function filter = nullptr);

void parse_line(const std::string& line, char separator, const std::set<char>& escape_chars,
                const std::set<char>& ignored_chars, std::vector<std::string>& result);
std::set<std::pair<unsigned short, unsigned short>>
getFullyConnectedMap(unsigned short nQubits);

std::string escapeChars(const std::string& s, const std::string& chars);
#endif //QMAP_UTILS_HPP
