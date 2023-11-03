//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "operations/Operation.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using Matrix      = std::vector<std::vector<double>>;
using Edge        = std::pair<std::uint16_t, std::uint16_t>;
using CouplingMap = std::set<Edge>;
using QubitSubset = std::set<std::uint16_t>;

struct Exchange {
  Exchange(const std::uint16_t f, const std::uint16_t s, const qc::OpType type)
      : first(f), second(s),
        middleAncilla(std::numeric_limits<decltype(middleAncilla)>::max()),
        op(type) {}
  Exchange(const std::uint16_t f, const std::uint16_t s,
           const std::uint16_t middleAnc, const qc::OpType type)
      : first(f), second(s), middleAncilla(middleAnc), op(type) {}
  std::uint16_t first;
  std::uint16_t second;
  std::uint16_t middleAncilla;
  qc::OpType    op;
};

class QMAPException : public std::runtime_error {
  std::string msg;

public:
  explicit QMAPException(std::string m)
      : std::runtime_error("QMAP Exception"), msg(std::move(m)) {}

  [[nodiscard]] const char* what() const noexcept override {
    return msg.c_str();
  }
};

class Dijkstra {
public:
  struct Node {
    /** true if the path contains any forward edge
        (i.e. an edge on which a CNOT with directionality q1->q2
        can be applied without the need of a CNOT reversal) */
    bool containsCorrectEdge = false;
    /** true if the node has already been expanded */
    bool visited = false;
    /** current qubit */
    std::optional<std::uint16_t> pos = std::nullopt;
    /** current cost of the path */
    double cost = -1.;
    /** current cost of the path with the last swap removed */
    double prevCost = -1.;
  };

  /**
   * @brief builds a distance table containing the minimal costs for moving
   * logical qubits from one physical qubit to/next to another (along the
   * cheapest path)
   *
   * @param n size of the distance table (i.e. number of qubits)
   * @param couplingMap coupling map specifying all edges in the architecture
   * @param distanceTable target table
   * @param edgeWeights matrix containing costs for swapping any two, connected
   * qubits (this might be uniform for all edges or different for each edge, as
   * e.g. in the case of fidelity-aware distances or distances on
   * mixed bi/unidirectional architectures)
   * @param reversalCost cost added if path consists only of back edges, i.e.
   * the 2-qubit gate to be mapped would need to be reversed in any case
   * (if reversal costs are handled elsewhere this should be set to 0)
   * @param removeLastEdge if set to true, the cost for the last swap on any
   * path is removed (i.e. distance is given for moving "next to" qubit)
   */
  static void buildTable(std::uint16_t n, const CouplingMap& couplingMap,
                         Matrix& distanceTable, const Matrix& edgeWeights,
                         double reversalCost, bool removeLastEdge);

protected:
  static void dijkstra(const CouplingMap& couplingMap, std::vector<Node>& nodes,
                       std::uint16_t start, const Matrix& edgeWeights,
                       double reversalCost);

  struct NodeComparator {
    bool operator()(const Node* x, const Node* y) {
      if (x->cost != y->cost) {
        return x->cost > y->cost;
      }
      return !x->containsCorrectEdge && y->containsCorrectEdge;
    }
  };
};

inline bool operator<(const Dijkstra::Node& x, const Dijkstra::Node& y) {
  if (x.cost != y.cost) {
    return x.cost < y.cost;
  }
  return x.containsCorrectEdge && !y.containsCorrectEdge;
}

/// Iterating routine through all combinations
/// \tparam Iterator iterator type
/// \param first iterator to beginning
/// \param k current iterator
/// \param last iterator to end
/// \return true if another combination was found
template <typename Iterator>
bool nextCombination(Iterator first, Iterator k, Iterator last) {
  /* Credits: Thomas Draper */
  if ((first == last) || (first == k) || (last == k)) {
    return false;
  }
  Iterator itr1 = first;
  Iterator itr2 = last;
  ++itr1;
  if (last == itr1) {
    return false;
  }
  itr1 = last;
  --itr1;
  itr1 = k;
  --itr2;
  while (first != itr1) {
    if (*--itr1 < *itr2) {
      Iterator j = k;
      while (!(*itr1 < *j)) {
        ++j;
      }
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
std::string printPi(std::vector<std::uint16_t>& pi);

/// Simple depth-first-search implementation used to check whether a given
/// subset of qubits is connected on the given architecture
///
/// \param current index of current qubit
/// \param visited visited qubits
/// \param cm coupling map of architecture
void dfs(std::uint16_t current, std::set<std::uint16_t>& visited,
         const CouplingMap& rcm);

using filter_function = std::function<bool(const QubitSubset&)>;
std::vector<QubitSubset> subsets(const QubitSubset& input, std::size_t size,
                                 const filter_function& filter = nullptr);

void        parseLine(const std::string& line, char separator,
                      const std::set<char>&     escapeChars,
                      const std::set<char>&     ignoredChars,
                      std::vector<std::string>& result);
CouplingMap getFullyConnectedMap(std::uint16_t nQubits);
