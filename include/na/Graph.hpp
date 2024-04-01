//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/mqt-qmap for more
// information.
//

#pragma once

#include "Definitions.hpp"
#include "DisjointSet.hpp"

#include <__algorithm/remove.h>
#include <__algorithm/remove_if.h>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <stack>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using Color = std::uint16_t;

namespace na {
template <class T> struct PairHash {
  size_t operator()(const std::pair<T, T>& x) const {
    return std::hash<T>{}(x.first) ^ std::hash<T>{}(x.second);
  }
};

template <class E> class Graph {
protected:
  // the adjecency matrix works with indices
  std::vector<std::vector<E>> adjacencyMatrix{};
  // the mapping of qubits to indices in the graph are stored in a map
  std::unordered_map<qc::Qubit, std::size_t> mapping;
  // the inverse mapping is used to get the qubit from the index
  std::unordered_map<std::size_t, qc::Qubit> invMapping;
  std::size_t                                nVertices = 0;

public:
  auto addVertex(qc::Qubit v) -> void {
    // check whether the vertex is already in the graph, if so do nothing
    if (mapping.find(v) == mapping.end()) {
      mapping[v]            = nVertices;
      invMapping[nVertices] = v;
      ++nVertices;
      for (auto& row : adjacencyMatrix) {
        row.emplace_back(nullptr);
      }
      adjacencyMatrix.emplace_back(1, nullptr);
    }
  }
  auto addEdge(qc::Qubit u, qc::Qubit v, E e) -> void {
    if (mapping.find(u) == mapping.end()) {
      addVertex(u);
    }
    if (mapping.find(v) == mapping.end()) {
      addVertex(v);
    }
    std::size_t const i = mapping.at(u);
    std::size_t const j = mapping.at(v);
    if (i < j) {
      adjacencyMatrix[i][j - i] = e;
    } else {
      adjacencyMatrix[j][i - j] = e;
    }
  }
  [[nodiscard]] auto getNVertices() const -> std::size_t { return nVertices; }
  [[nodiscard]] auto getNEdges() const -> std::size_t {
    return std::accumulate(
        adjacencyMatrix.cbegin(), adjacencyMatrix.cend(), 0UL,
        [](std::size_t acc, const auto& row) {
          return acc + static_cast<std::size_t>(std::count_if(
                           row.cbegin(), row.cend(),
                           [](const auto& edge) { return edge != nullptr; }));
        });
  }
  [[nodiscard]] auto getEdge(qc::Qubit v, qc::Qubit u) const -> E {
    std::size_t const i = mapping.at(v);
    std::size_t const j = mapping.at(u);
    if (i < j ? adjacencyMatrix[i][j - i] != nullptr
              : adjacencyMatrix[j][i - j] != nullptr) {
      return i < j ? adjacencyMatrix[i][j - i] : adjacencyMatrix[j][i - j];
    }
    std::stringstream ss;
    ss << "The edge (" << v << ", " << u << ") does not exist.";
    throw std::invalid_argument(ss.str());
  }
  [[nodiscard]] auto getDegree(qc::Qubit v) const -> std::size_t {
    if (mapping.find(v) == mapping.end()) {
      std::stringstream ss;
      ss << "The vertex " << v << " is not in the graph.";
      throw std::invalid_argument(ss.str());
    }
    std::size_t const i      = mapping.at(v);
    std::size_t       degree = 0;
    for (std::size_t j = 0; j < nVertices; ++j) {
      if ((i < j and adjacencyMatrix[i][j - i] != nullptr) or
          (j < i and adjacencyMatrix[j][i - j] != nullptr)) {
        ++degree;
      }
    }
    return degree;
  }
  [[nodiscard]] auto getVertices() const -> std::unordered_set<qc::Qubit> {
    return std::accumulate(mapping.cbegin(), mapping.cend(),
                           std::unordered_set<qc::Qubit>(),
                           [](auto& acc, const auto& v) {
                             acc.emplace(v.first);
                             return acc;
                           });
  }
  [[nodiscard]] auto isAdjacent(const qc::Qubit u, const qc::Qubit v) const
      -> bool {
    std::size_t const i = mapping.at(u);
    std::size_t const j = mapping.at(v);
    return (i < j and adjacencyMatrix[i][j - i] != nullptr) or
           (j < i and adjacencyMatrix[j][i - j] != nullptr);
  }
  [[nodiscard]] auto
  // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
  isAdjacentEdge(const std::pair<qc::Qubit, qc::Qubit>& e,
                 const std::pair<qc::Qubit, qc::Qubit>& f) const -> bool {
    return e.first == f.first or e.first == f.second or e.second == f.first or
           e.second == f.second;
  }
  [[nodiscard]] auto getMaxIndependentSet() const
      -> std::unordered_set<qc::Qubit> {
    std::unordered_set<qc::Qubit> result;
    std::vector<qc::Qubit>        queue = sortByDegreeDesc(getVertices());
    for (qc::Qubit v = queue.front(); !queue.empty(); v = queue.front()) {
      result.emplace(v);
      queue.erase(std::remove_if(queue.begin(), queue.end(),
                                 [&](const qc::Qubit& u) {
                                   return u == v or isAdjacent(u, v);
                                 }),
                  queue.end());
    }
    return result;
  }
  [[nodiscard]] auto coveredEdges(const std::unordered_set<qc::Qubit>& vs) const
      -> std::unordered_set<std::pair<qc::Qubit, qc::Qubit>,
                            PairHash<qc::Qubit>> {
    for (const auto& v : vs) {
      if (mapping.find(v) == mapping.end()) {
        throw std::invalid_argument("The set of qubits must be a subset of the "
                                    "domain of the mapping.");
      }
    }
    std::unordered_set<std::pair<qc::Qubit, qc::Qubit>, PairHash<qc::Qubit>>
        result;
    for (qc::Qubit const& v : vs) {
      std::size_t const i = mapping.at(v);
      for (std::size_t j = 0; j < nVertices; ++j) {
        // check whether indices i and j are adjacent
        if (i < j and adjacencyMatrix[i][j - i] != nullptr) {
          qc::Qubit const& u = invMapping.at(j);
          if (u < v) {
            result.emplace(u, v);
          } else {
            result.emplace(v, u);
          }
        } else if (j < i and adjacencyMatrix[j][i - j] != nullptr) {
          qc::Qubit const& u = invMapping.at(j);
          if (u < v) {
            result.emplace(u, v);
          } else {
            result.emplace(v, u);
          }
        }
      }
    }
    return result;
  }
  /**
   * @brief Get the Least Admissable Color for an edge.
   * @details For a coloring to be valid no two adjacent edges can have the same
   * color. Consequently, the least admissable color is one that is not used by
   * any adjacent edge.
   * Additionally, the color must be greater than the maximum color of any
   * adjacent edge that does not contain the vertex v to ensure that the
   * following constraint is satisfied:
   * Let two nodes u and v be in the minimum maximal independent set and the
   * node w and w' both adjacent to u and v. The edge (u, w) has a smaller
   * coloring than the edge (w, v) iff the edge (u, w') has a smaller coloring
   * than the edge (w', v), e.g.
   *
   *                (u)—– 0 —(w)
   *                  \        \
   *                   3        1
   *                    \        \
   *                    (w')— 4 —(v)
   *
   * A '2' instead of the '4' would not be allowed the other colors unchanged.
   *
   * @param coloring
   * @param e
   * @param v
   * @return auto
   */
  [[nodiscard]] auto getLeastAdmissableColor(
      const std::unordered_map<std::pair<qc::Qubit, qc::Qubit>, Color,
                               PairHash<qc::Qubit>>& coloring,
      const Color& maxColor, const std::pair<qc::Qubit, qc::Qubit>& e,
      const qc::Qubit& v) const -> Color {
    // compute the minimum admissable color as the maximum color +1 of adjacent
    // edges that do not contain the vertex v
    Color minAdmissableColor = 0;
    for (const auto& [f, k] : coloring) {
      if (isAdjacentEdge(e, f) and v != f.first and v != f.second) {
        minAdmissableColor =
            std::max(minAdmissableColor, static_cast<Color>(k + 1));
      }
    }
    std::set<Color> freeColors; // shall be sorted
    for (Color k = minAdmissableColor; k <= maxColor + 1; ++k) {
      freeColors.emplace(k);
    }
    for (const auto& [f, k] : coloring) {
      if (isAdjacentEdge(e, f)) {
        freeColors.erase(k);
      }
    }
    // The minAdmissableColor is now the minimum of the free colors
    minAdmissableColor = *freeColors.begin();
    return minAdmissableColor;
  }
  /**
   * @brief Colors all given edges starting with edges that are adjacent to the
   * first vertex in the queue.
   *
   * @param edges e.g. edges adjacent to a selected vertex
   * @param nodesQueue could for example be sorted by degree, highest first
   * @return
   */
  [[nodiscard]] auto
  colorEdges(const std::unordered_set<std::pair<qc::Qubit, qc::Qubit>,
                                      PairHash<qc::Qubit>>& edges,
             const std::vector<qc::Qubit>&                  nodesQueue) const
      -> std::unordered_map<std::pair<qc::Qubit, qc::Qubit>, Color,
                            PairHash<qc::Qubit>> {
    std::unordered_map<std::pair<qc::Qubit, qc::Qubit>, Color,
                       PairHash<qc::Qubit>>
          coloring{};
    Color maxColor = 0;
    // number of distinct colors of edges adjacent to an edge
    std::unordered_map<std::pair<qc::Qubit, qc::Qubit>, std::size_t,
                       PairHash<qc::Qubit>>
        nAdjColors;
    // the degree of the edge seen as a node
    std::unordered_map<std::pair<qc::Qubit, qc::Qubit>, std::size_t,
                       PairHash<qc::Qubit>>
        edgeDegree;
    for (const auto& e : edges) {
      nAdjColors[e] = 0;
      edgeDegree[e] = static_cast<std::size_t>(
          std::count_if(edges.cbegin(), edges.cend(),
                        [&](const std::pair<qc::Qubit, qc::Qubit>& f) {
                          return isAdjacentEdge(e, f);
                        }));
    }

    for (const auto& v : nodesQueue) {
      std::vector<std::pair<qc::Qubit, qc::Qubit>> adjacentEdges{};
      std::copy_if(edges.cbegin(), edges.cend(),
                   std::back_inserter(adjacentEdges),
                   [&](const std::pair<qc::Qubit, qc::Qubit>& e) {
                     return e.first == v or e.second == v;
                   });
      while (!adjacentEdges.empty()) {
        // calculate the maximum number of distinct colors of adjacent edges
        std::size_t maxNAdjColor = 0;
        for (const auto& e : adjacentEdges) {
          maxNAdjColor = std::max(maxNAdjColor, nAdjColors[e]);
        }
        // select the edges with the maximum number of distinct colors
        std::vector<std::pair<qc::Qubit, qc::Qubit>> selectedEdges{};
        std::copy_if(adjacentEdges.cbegin(), adjacentEdges.cend(),
                     std::back_inserter(selectedEdges),
                     [&](const std::pair<qc::Qubit, qc::Qubit>& e) {
                       return nAdjColors[e] == maxNAdjColor;
                     });
        // if there are multiple edges with the maximum number of distinct
        if (selectedEdges.size() > 1) {
          // calculate the maximum degree of the selected edges
          std::size_t maxDegree = 0;
          for (const auto& e : selectedEdges) {
            maxDegree = std::max(maxDegree, edgeDegree[e]);
          }
          // remove all edges with a degree less than the maximum degree
          selectedEdges.erase(
              std::remove_if(selectedEdges.begin(), selectedEdges.end(),
                             [&](const std::pair<qc::Qubit, qc::Qubit>& e) {
                               return edgeDegree[e] < maxDegree;
                             }),
              selectedEdges.end());
        }
        // pick any edge from the remaining selected edges
        auto& e = selectedEdges[0];
        // color the edge
        adjacentEdges.erase(
            std::remove(adjacentEdges.begin(), adjacentEdges.end(), e),
            adjacentEdges.end());
        coloring[e] = getLeastAdmissableColor(coloring, maxColor, e, v);
        maxColor    = std::max(maxColor, coloring[e]);
        // update the number of distinct colors of adjacent edges
        for (const auto& f : edges) {
          if (isAdjacentEdge(e, f)) {
            std::vector<bool> usedColors(maxColor + 1, false);
            for (const auto& g : edges) {
              if (isAdjacentEdge(f, g) and coloring.find(g) != coloring.end()) {
                usedColors[coloring[g]] = true;
              }
            }
            nAdjColors[f] = static_cast<std::size_t>(
                std::count(usedColors.cbegin(), usedColors.cend(),
                           [](const bool& b) { return b; }));
          }
        }
      }
    }
    return coloring;
  }
  [[nodiscard]] auto
  sortByDegreeDesc(const std::unordered_set<qc::Qubit>& unsorted) const
      -> std::vector<qc::Qubit> {
    std::vector<qc::Qubit> sorted(unsorted.begin(), unsorted.end());
    std::sort(sorted.begin(), sorted.end(),
              [&](const qc::Qubit& u, const qc::Qubit& v) {
                return getDegree(u) > getDegree(v);
              });
    return sorted;
  }
  /**
   * @brief Performs DFS to return the set of vertices in their topological
   * order with respective to the given neighboring function.
   * @details The topological order is a linear ordering of the vertices such
   * that for every directed edge u -> v, u comes before v in the ordering.
   * To check whether there is an edge from u to v, the function isEdge is
   * called. Note, that this function assumes directed edges. Otherwise there
   * cannot be a topological ordering. The function throws an exception if the
   * graph contains a cycle.
   * @param vertices
   * @param isEdge
   * @return std::vector<qc::Qubit>
   */
  [[nodiscard]] static auto
  topoOrder(const std::unordered_set<qc::Qubit>&             vertices,
            const std::function<bool(qc::Qubit, qc::Qubit)>& isEdge)
      -> std::vector<qc::Qubit> {
    auto n = vertices.size();
    // Calculate indegree of each node
    std::unordered_map<qc::Qubit, std::uint32_t> indegree{};
    for (const auto& v : vertices) {
      indegree[v] = 0;
      for (const auto& u : vertices) {
        if (isEdge(u, v)) {
          indegree[v]++;
        }
      }
    }
    std::stack<qc::Qubit>         stack{};
    std::vector<qc::Qubit>        result{};
    std::unordered_set<qc::Qubit> visited{};
    // Push nodes with 0 indegree onto the stack
    for (const auto& v : vertices) {
      if (indegree[v] == 0) {
        stack.push(v);
        visited.insert(v);
      }
    }
    // Perform DFS
    while (!stack.empty()) {
      const qc::Qubit u = stack.top();
      stack.pop();
      result.push_back(u);

      for (const auto& v : vertices) {
        if (isEdge(u, v)) {
          if (visited.find(v) == visited.end()) {
            indegree[v]--;
            if (indegree[v] == 0) {
              stack.push(v);
              visited.insert(v);
            }
          }
        }
      }
    }
    // Check for cycle
    if (result.size() != n) {
      throw std::invalid_argument("The graph contains a cycle.");
    }
    return result;
  }
  [[nodiscard]] auto computeSlackPositions(
      std::vector<qc::Qubit> const&                  moveable,
      std::vector<qc::Qubit> const&                  fixed,
      std::unordered_map<std::pair<qc::Qubit, qc::Qubit>, Color,
                         PairHash<qc::Qubit>> const& coloring) const
      -> std::vector<std::size_t> {
    const Color maxColor = std::accumulate(coloring.cbegin(), coloring.cend(),
                                           static_cast<Color>(0),
                                           [](Color acc, const auto& value) {
                                             return std::max(acc, value.second);
                                           });
    std::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t,
                       PairHash<std::size_t>>
        slack{};
    for (Color t = 0; t <= maxColor; ++t) {
      // slack positions for timestamp t
      std::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t,
                         PairHash<std::size_t>>
          tSlack{};
      // x-coordinates of movable vertices at timestamp t
      std::unordered_map<qc::Qubit, std::size_t> moveableXs{};
      for (const qc::Qubit v : moveable) {
        // get the neighbors of v
        const std::set<qc::Qubit> neighborsOfV = std::accumulate(
            coloring.cbegin(), coloring.cend(), std::set<qc::Qubit>(),
            [&](std::set<qc::Qubit> acc, auto value) {
              std::pair<qc::Qubit, qc::Qubit> e = value.first;
              if (value.second == t) {
                if (e.first == v) {
                  acc.emplace(e.second);
                } else if (e.second == v) {
                  acc.emplace(e.first);
                }
              }
              return acc;
            });
        if (!neighborsOfV.empty()) {
          assert(neighborsOfV.size() == 1);
          qc::Qubit const u = *neighborsOfV.begin();
          // get index of u in fixed which is the x-coordinate \wo slack
          const auto& it = std::find(fixed.cbegin(), fixed.cend(), u);
          assert(it != fixed.cend());
          // set x-pos of v to x-pos of u, which is the index of u in fixed
          moveableXs[v] =
              static_cast<std::size_t>(std::distance(fixed.cbegin(), it));
        }
      }
      // map the keys of moveableXs to their indices in moveable
      std::vector<std::size_t> moveableXsIds{};
      std::transform(moveableXs.cbegin(), moveableXs.cend(),
                     std::back_inserter(moveableXsIds), [&](const auto& value) {
                       const auto& it = std::find(moveable.cbegin(),
                                                  moveable.cend(), value.first);
                       return static_cast<std::size_t>(
                           std::distance(moveable.cbegin(), it));
                     });
      for (const qc::Qubit v : moveable) {
        if (moveableXs.find(v) == moveableXs.end()) {
          // index of v in moveable
          const auto& it = std::find(moveable.cbegin(), moveable.cend(), v);
          const auto  i =
              static_cast<std::size_t>(std::distance(moveable.cbegin(), it));
          // get the index of the left neighbor in the keys of moveableXs
          std::vector<std::size_t> leftNeighbors;
          std::copy_if(moveableXsIds.cbegin(), moveableXsIds.cend(),
                       std::back_inserter(leftNeighbors),
                       [&](const auto& j) { return j > i; });
          std::vector<std::size_t> rightNeighbors;
          std::copy_if(moveableXsIds.cbegin(), moveableXsIds.cend(),
                       std::back_inserter(rightNeighbors),
                       [&](const auto& j) { return j < i; });
          // if there is a left and right neighbor, we have to add a
          // slack position
          if (!leftNeighbors.empty() and !rightNeighbors.empty()) {
            // get min left and max right neighbor
            const auto& leftNeighbor =
                std::min_element(leftNeighbors.cbegin(), leftNeighbors.cend());
            const auto& rightNeighbor = std::max_element(
                rightNeighbors.cbegin(), rightNeighbors.cend());
            const auto& pair = std::pair<std::int64_t, std::int64_t>(
                moveableXs[moveable[*leftNeighbor]],
                moveableXs[moveable[*rightNeighbor]]);
            if (tSlack.find(pair) == tSlack.end()) {
              tSlack[pair] = 1;
            } else {
              tSlack[pair]++;
            }
          }
        }
      }
      // tSlack is completely computed until here
      for (const auto& elem : slack) {
        const auto& pair = elem.first;
        for (std::size_t i = 0; i < elem.second; ++i) {
          // check whether narrow slack positions are subsumed by wider slack
          // positions added from previous timestamps
          const auto& superPairs = std::accumulate(
              tSlack.cbegin(), tSlack.cend(),
              std::vector<std::pair<std::size_t, std::size_t>>(),
              [&](auto& acc, const auto& value) {
                if (value.first.first >= pair.first &&
                    value.first.second <= pair.second) {
                  acc.emplace_back(value.first);
                }
                return acc;
              });
          if (!superPairs.empty()) {
            const std::pair<std::size_t, std::size_t>& minSuperPair =
                *(std::min_element(superPairs.cbegin(), superPairs.cend(),
                                   [](const auto& a, const auto& b) {
                                     return a.second - a.first <
                                            b.second - b.first;
                                   }));
            tSlack[minSuperPair]--;
            if (tSlack[minSuperPair] == 0) {
              tSlack.erase(minSuperPair);
            }
          }
        }
      }
      for (const auto& [pair, s] : tSlack) {
        if (slack.find(pair) == slack.end()) {
          slack[pair] = s;
        } else {
          slack[pair] += s;
        }
      }
    }
    std::vector<std::size_t> slackPositions{};
    for (const auto& [pair, s] : slack) {
      for (std::size_t i = 0; i < s; ++i) {
        slackPositions.emplace_back(pair.first);
      }
    }
    std::sort(slackPositions.begin(), slackPositions.end());
    return slackPositions;
  }
  [[nodiscard]] auto
  groupByConnectedComponent(const std::vector<qc::Qubit>& sequence) const
      -> std::vector<qc::Qubit> {
    const auto&            vertices = getVertices();
    DisjointSet<qc::Qubit> ds(vertices.cbegin(), vertices.cend());
    for (const qc::Qubit& v : vertices) {
      for (const qc::Qubit& u : vertices) {
        if (isAdjacent(v, u)) {
          ds.unionSet(v, u);
        }
      }
    }
    std::vector<qc::Qubit> result{};
    for (const qc::Qubit& v : vertices) {
      if (ds.findSet(v) == v) {
        for (const qc::Qubit& u : sequence) {
          if (ds.findSet(u) == v) {
            result.emplace_back(u);
          }
        }
      }
    }
    return result;
  }
  /**
   * @brief Partitions the set of vertices into moveable and fixed vertices
   * with the aim to maximize the number of executable gates in one run
   * without reloading atoms.
   *
   * @return A triple of vectors containing the moveable vertices, the fixed
   * vertices and slack positions between the fixed vertices.
   */
  [[nodiscard]] auto computeSequence() const
      -> std::pair<std::vector<std::unordered_map<qc::Qubit, std::int64_t>>,
                   std::unordered_map<qc::Qubit, std::int64_t>> {
    auto        mis               = getMaxIndependentSet();
    const auto& sequenceUngrouped = sortByDegreeDesc(mis);
    const auto& sequence = groupByConnectedComponent(sequenceUngrouped);
    std::unordered_map<std::pair<qc::Qubit, qc::Qubit>, Color,
                       PairHash<qc::Qubit>> const coloring =
        colorEdges(coveredEdges(mis), sequence);
    // take the difference of all vertices and the mis
    std::unordered_set<qc::Qubit> const& difference = std::accumulate(
        mapping.cbegin(), mapping.cend(), std::unordered_set<qc::Qubit>(),
        [&](auto& acc, const auto& v) {
          if (mis.find(v.first) == mis.cend()) {
            acc.emplace(v.first);
          }
          return acc;
        });
    std::vector<
        qc::Qubit> const& fixed = topoOrder(difference, [&](const qc::Qubit& v,
                                                            const qc::Qubit&
                                                                u) {
      const std::set<qc::Qubit> neighborsOfV = std::accumulate(
          coloring.cbegin(), coloring.cend(), std::set<qc::Qubit>(),
          [&](std::set<qc::Qubit> acc, auto value) {
            std::pair<qc::Qubit, qc::Qubit> e = value.first;
            if (e.first == v) {
              acc.emplace(e.second);
            } else if (e.second == v) {
              acc.emplace(e.first);
            }
            return acc;
          });
      const std::set<qc::Qubit> neighborsOfU = std::accumulate(
          coloring.cbegin(), coloring.cend(), std::set<qc::Qubit>(),
          [&](std::set<qc::Qubit> acc, auto value) {
            std::pair<qc::Qubit, qc::Qubit> e = value.first;
            if (e.first == u) {
              acc.emplace(e.second);
            } else if (e.second == u) {
              acc.emplace(e.first);
            }
            return acc;
          });
      std::set<qc::Qubit> inter{};
      std::set_intersection(neighborsOfV.cbegin(), neighborsOfV.cend(),
                            neighborsOfU.cbegin(), neighborsOfU.cend(),
                            std::inserter(inter, inter.end()));
      if (!inter.empty()) {
        const qc::Qubit w = *inter.begin();
        const auto& vit = coloring.find(std::pair<qc::Qubit, qc::Qubit>(v, w));
        const Color vwColor =
            vit != coloring.end()
                ? vit->second
                : coloring.find(std::pair<qc::Qubit, qc::Qubit>(w, v))->second;
        const auto& uit = coloring.find(std::pair<qc::Qubit, qc::Qubit>(u, w));
        const Color uwColor =
            uit != coloring.end()
                ? uit->second
                : coloring.find(std::pair<qc::Qubit, qc::Qubit>(w, u))->second;
        if (vwColor + 1 == uwColor) {
          return true;
        }
      } else { // inter.empty() ==> no common neighbor
        const auto maxColor = std::accumulate(
            coloring.cbegin(), coloring.cend(), static_cast<Color>(0),
            [](Color acc, const auto& value) {
              return std::max(acc, value.second);
            });
        for (Color t = 0; t <= maxColor; ++t) {
          std::vector<std::size_t> enumerate(sequence.size());
          std::iota(enumerate.begin(), enumerate.end(), 0);
          const auto vMaxSeq = std::accumulate(
              enumerate.cbegin(), enumerate.cend(), sequence.size(),
              [&](std::size_t acc, const std::size_t i) {
                const auto value = sequence[i];
                if (neighborsOfV.find(value) != neighborsOfV.end() and
                    ((coloring.find(std::pair<qc::Qubit, qc::Qubit>(
                          v, value)) != coloring.end() and
                      coloring.find(std::pair<qc::Qubit, qc::Qubit>(v, value))
                              ->second == t) or
                     (coloring.find(std::pair<qc::Qubit, qc::Qubit>(
                          value, v)) != coloring.end() and
                      coloring.find(std::pair<qc::Qubit, qc::Qubit>(value, v))
                              ->second == t))) {
                  return std::min(acc, static_cast<std::size_t>(i));
                }
                return acc;
              });
          const auto uMaxSeq = std::accumulate(
              enumerate.cbegin(), enumerate.cend(), sequence.size(),
              [&](std::size_t acc, const std::size_t i) {
                const auto value = sequence[i];
                if (neighborsOfU.find(value) != neighborsOfU.end() and
                    ((coloring.find(std::pair<qc::Qubit, qc::Qubit>(
                          u, value)) != coloring.end() and
                      coloring.find(std::pair<qc::Qubit, qc::Qubit>(u, value))
                              ->second == t) or
                     (coloring.find(std::pair<qc::Qubit, qc::Qubit>(
                          value, u)) != coloring.end() and
                      coloring.find(std::pair<qc::Qubit, qc::Qubit>(value, u))
                              ->second == t))) {
                  return std::min(acc, static_cast<std::size_t>(i));
                }
                return acc;
              });
          if (vMaxSeq < sequence.size() and uMaxSeq < sequence.size() and
              vMaxSeq > uMaxSeq) {
            return true;
          }
        }
      }
      return false;
    });
    const auto& slack = computeSlackPositions(sequence, fixed, coloring);
    // compute relative x positions of fixed vertices
    std::unordered_map<qc::Qubit, std::int64_t> fixedPositions{};
    for (std::uint32_t x = 0, i = 0; x < fixed.size(); ++x) {
      fixedPositions[fixed[x]] = static_cast<std::int64_t>(x) + i;
      while (i < slack.size() and x == slack[i]) {
        ++i;
      }
    }
    // compute relative x positions of moveable vertices at every timestamp
    const Color maxColor = std::accumulate(coloring.cbegin(), coloring.cend(),
                                           static_cast<Color>(0),
                                           [](Color acc, const auto& value) {
                                             return std::max(acc, value.second);
                                           });
    std::vector<std::unordered_map<qc::Qubit, std::int64_t>> moveablePositions(
        maxColor + 1);
    for (Color t = 0; t <= maxColor; ++t) {
      for (const qc::Qubit v : sequence) {
        // get the neighbors of v at timestamp t
        const std::set<qc::Qubit>& neighborsOfV = std::accumulate(
            coloring.cbegin(), coloring.cend(), std::set<qc::Qubit>(),
            [&](std::set<qc::Qubit> acc, auto value) {
              if (value.second == t) {
                std::pair<qc::Qubit, qc::Qubit> e = value.first;
                if (e.first == v) {
                  acc.emplace(e.second);
                } else if (e.second == v) {
                  acc.emplace(e.first);
                }
              }
              return acc;
            });
        if (!neighborsOfV.empty()) {
          assert(neighborsOfV.size() == 1);
          qc::Qubit const u = *neighborsOfV.begin();
          // get x-coordinate of u which is a fixed vertex
          const auto& it = fixedPositions.find(u);
          assert(it != fixed.cend());
          moveablePositions[t][v] = it->second;
        }
      }
      for (std::size_t i = 0; i < sequence.size(); ++i) {
        const qc::Qubit v = sequence[i];
        if (moveablePositions[t].find(v) == moveablePositions[t].end()) {
          if (i > 0) {
            // right neighbor already has a position
            const auto rightNeighbor = i - 1;
            const auto rightNeighborX =
                moveablePositions[t][sequence[rightNeighbor]];
            const std::int64_t minX = std::min(rightNeighborX - 1, -1LL);
            assert(rightNeighborX - minX >= 0);
            std::vector<std::int64_t> xrange(
                static_cast<std::size_t>(rightNeighborX - minX));
            std::iota(xrange.begin(), xrange.end(), minX);
            std::vector<std::int64_t> freeX;
            std::copy_if(xrange.cbegin(), xrange.cend(),
                         std::back_inserter(freeX), [&](const std::int64_t x) {
                           assert(x >= 0);
                           return std::all_of(fixedPositions.cbegin(),
                                              fixedPositions.cend(),
                                              [&](const auto& value) {
                                                return x != value.second;
                                              });
                         });
            const auto& maxIt = std::max_element(freeX.cbegin(), freeX.cend());
            moveablePositions[t][v] = *maxIt;
          } else {
            // it is the rightmost atom (w/o position)
            const auto leftNeighbor = std::accumulate(
                moveablePositions[t].cbegin(), moveablePositions[t].cend(),
                std::make_pair(0UL, 0LL),
                [](const std::pair<const qc::Qubit, std::int64_t>& acc,
                   const auto&                                     value) {
                  if (value.second > acc.second) {
                    return value;
                  }
                  return acc;
                });
            // k is the number of spots between the left positioned neighbor and
            // the leftmost atom (including itself)
            const auto         k = static_cast<std::size_t>(std::distance(
                sequence.cbegin(), std::find(sequence.cbegin(), sequence.cend(),
                                                     leftNeighbor.first)));
            const std::int64_t maxX =
                std::accumulate(fixedPositions.cbegin(), fixedPositions.cend(),
                                0LL, [](const auto& acc, const auto& value) {
                                  return std::max(acc, value.second);
                                });
            std::vector<std::int64_t> xrange(
                static_cast<std::size_t>(maxX - leftNeighbor.second));
            std::iota(xrange.begin(), xrange.end(), leftNeighbor.second + 1);
            std::vector<std::int64_t> freeX;
            std::copy_if(xrange.cbegin(), xrange.cend(),
                         std::back_inserter(freeX), [&](const std::int64_t x) {
                           return std::all_of(fixedPositions.cbegin(),
                                              fixedPositions.cend(),
                                              [&](const auto& value) {
                                                return x != value.second;
                                              });
                         });
            if (k <= freeX.size()) {
              moveablePositions[t][v] = freeX[k - 1];
            } else {
              moveablePositions[t][v] = maxX + static_cast<std::int64_t>(k) -
                                        static_cast<std::int64_t>(freeX.size());
            }
          }
        }
      }
    }
    return std::make_pair(moveablePositions, fixedPositions);
  }
  // Outputs a string representation of the graph in the DOT format
  [[nodiscard]] auto toString() const -> std::string {
    std::stringstream ss;
    ss << "graph {\n";
    for (const auto& [v, i] : mapping) {
      ss << "  " << i << " [label=\"" << v << "\"];\n";
    }
    for (std::size_t i = 0; i < nVertices; ++i) {
      for (std::size_t j = i + 1; j < nVertices; ++j) {
        if (adjacencyMatrix[i][j - i] != nullptr) {
          ss << "  " << i << " -- " << j << ";\n";
        }
      }
    }
    ss << "}\n";
    return ss.str();
  }
  friend auto operator<<(std::ostream& os, const Graph& g) -> std::ostream& {
    os << g.toString(); // Using toString() method
    return os;
  }
};
} // namespace na
