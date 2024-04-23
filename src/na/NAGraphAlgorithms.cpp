//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/mqt-qmap for more
// information.
//

#include "NAGraphAlgorithms.hpp"

#include "Definitions.hpp"
#include "DirectedAcyclicGraph.hpp"
#include "DisjointSet.hpp"
#include "Layer.hpp"
#include "UndirectedGraph.hpp"

#include <__algorithm/remove_if.h>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace na {

auto NAGraphAlgorithms::getMaxIndependentSet(const InteractionGraph& g)
    -> std::unordered_set<qc::Qubit> {
  std::unordered_set<qc::Qubit> result;
  const auto&                   vertices = g.getVertices();
  std::vector                   queue(vertices.cbegin(), vertices.cend());
  // sort the vertices by degree in descending order
  std::sort(queue.begin(), queue.end(), [&](const auto& u, const auto& v) {
    return g.getDegree(u) > g.getDegree(v);
  });
  for (qc::Qubit v = queue.front(); !queue.empty(); v = queue.front()) {
    result.emplace(v);
    queue.erase(std::remove_if(queue.begin(), queue.end(),
                               [&](const qc::Qubit& u) {
                                 return u == v or g.isAdjacent(u, v);
                               }),
                queue.end());
  }
  return result;
}

auto NAGraphAlgorithms::coveredEdges(const InteractionGraph&              g,
                                     const std::unordered_set<qc::Qubit>& vs)
    -> std::unordered_set<Edge, qc::PairHash<qc::Qubit>> {
  std::unordered_set<Edge, qc::PairHash<qc::Qubit>> result;
  for (qc::Qubit const& v : vs) {
    try {
      result.merge(g.getAdjacentEdges(v));
    } catch (const std::invalid_argument&) {
      std::stringstream ss;
      ss << "Vertices must be a subset of all verticies in the graph, in "
            "particular "
         << v << " is not in the graph.";
      throw std::invalid_argument(ss.str());
    }
  }
  return result;
}

auto NAGraphAlgorithms::getLeastAdmissableColor(
    const InteractionGraph&                                         g,
    const std::unordered_map<Edge, Color, qc::PairHash<qc::Qubit>>& coloring,
    const Color& maxColor, const Edge& e, const qc::Qubit& v,
    const std::vector<qc::Qubit>&              sequence,
    const qc::DirectedAcyclicGraph<qc::Qubit>& partialOrder,
    std::unordered_map<qc::Qubit, std::unordered_map<Color, std::size_t>> ranks)
    -> Color {
  // compute the minimum admissable color as the maximum color +1 of adjacent
  // edges that do not contain the vertex v
  Color minAdmissableColor = 0;
  for (const auto& [f, k] : coloring) {
    if (g.isAdjacentEdge(e, f) and v != f.first and v != f.second) {
      minAdmissableColor =
          std::max(minAdmissableColor, static_cast<Color>(k + 1));
    }
  }
  std::set<Color> freeColors; // shall be sorted
  for (Color k = minAdmissableColor; k <= maxColor + 1; ++k) {
    freeColors.emplace(k);
  }
  for (const auto& [f, k] : coloring) {
    if (g.isAdjacentEdge(e, f)) {
      freeColors.erase(k);
    }
  }
  // The minAdmissableColor is now the minimum of the free colors
  // that does not generate a cycle in the graph induced by the partial order
  // of SLM traps.

  // u is the SLM trap that is adjacent to v in the edge e
  const qc::Qubit u = v == e.first ? e.second : e.first;
  // get the rank of u, i.e., the index of v in the sequence of AOD traps
  // (not contained in ranks since the edge e is not colored yet)
  const auto rankOfU = static_cast<std::size_t>(std::distance(
      sequence.cbegin(), std::find(sequence.cbegin(), sequence.cend(), v)));
  for (const auto leastAdmissableColor : freeColors) {
    bool isAdmissable = true;
    for (const auto& [f, k] : coloring) {
      if (f.first == v or f.second == v) {
        const qc::Qubit w = f.first == v ? f.second : f.first;
        if (k > leastAdmissableColor) {
          if (partialOrder.isReachable(w, u)) {
            isAdmissable = false;
            break;
          }
        } else if (k < leastAdmissableColor) {
          if (partialOrder.isReachable(u, w)) {
            throw std::logic_error("Coloring cannot be completed to a valid "
                                   "one (cycle is unavoidable).");
          }
        }
      } else if (k == leastAdmissableColor) {
        // get the SLM atom from the edge f
        if (const qc::Qubit w = std::find(sequence.cbegin(), sequence.cend(),
                                          f.first) == sequence.end()
                                    ? f.first
                                    : f.second;
            rankOfU > ranks[w][k]) {
          if (partialOrder.isReachable(w, u)) {
            isAdmissable = false;
            break;
          }
        } else if (rankOfU < ranks[w][k]) {
          if (partialOrder.isReachable(u, w)) {
            isAdmissable = false;
            break;
          }
        } else {
          throw std::logic_error("Ranks are not consistent.");
        }
      }
    }
    if (isAdmissable) {
      return leastAdmissableColor;
    }
  }
  throw std::logic_error("No admissable color found (should never occur).");
}

auto NAGraphAlgorithms::colorEdges(
    const InteractionGraph&                                  g,
    const std::unordered_set<Edge, qc::PairHash<qc::Qubit>>& edges,
    const std::vector<qc::Qubit>&                            nodesQueue)
    -> std::pair<std::unordered_map<Edge, Color, qc::PairHash<qc::Qubit>>,
                 qc::DirectedAcyclicGraph<qc::Qubit>> {
  std::unordered_map<Edge, Color, qc::PairHash<qc::Qubit>> coloring{};
  Color                                                    maxColor = 0;
  // number of distinct colors of edges adjacent to an edge
  std::unordered_map<Edge, std::size_t, qc::PairHash<qc::Qubit>> nAdjColors;
  // the degree of the edge seen as a node
  std::unordered_map<Edge, std::size_t, qc::PairHash<qc::Qubit>> edgeDegree;
  for (const auto& e : edges) {
    nAdjColors[e] = 0;
    edgeDegree[e] = static_cast<std::size_t>(
        std::count_if(edges.cbegin(), edges.cend(),
                      [&](const Edge& f) { return g.isAdjacentEdge(e, f); }));
  }
  // represent the partial order on the SLM atoms as a directed acyclic graph
  qc::DirectedAcyclicGraph<qc::Qubit> partialOrder;
  const std::unordered_set aodTraps(nodesQueue.cbegin(), nodesQueue.cend());
  for (const auto& v : g.getVertices()) {
    if (aodTraps.find(v) == aodTraps.end()) {
      // only add SLM traps, AOD traps are the one in the queue
      partialOrder.addVertex(v);
    }
  }
  // for every SLM trap and color, calculate its rank, i.e., the index of the
  // interaction partner in the sequence of AOD traps
  std::unordered_map<qc::Qubit, std::unordered_map<Color, std::size_t>> ranks;

  for (const auto& v : nodesQueue) {
    std::vector<Edge> adjacentEdges{};
    std::copy_if(edges.cbegin(), edges.cend(),
                 std::back_inserter(adjacentEdges),
                 [&](const Edge& e) { return e.first == v or e.second == v; });
    std::sort(adjacentEdges.begin(), adjacentEdges.end(),
              [&](const Edge& a, const Edge& b) {
                const auto u = a.first == v ? a.second : a.first;
                const auto w = b.first == v ? b.second : b.first;
                return partialOrder.isReachable(u, w) or
                       (!partialOrder.isReachable(w, u) and
                        (nAdjColors[a] > nAdjColors[b] or
                         (nAdjColors[a] == nAdjColors[b] and
                          edgeDegree[a] > edgeDegree[b])));
              });
    for (const auto& e : adjacentEdges) {
      // color the edge
      coloring[e] = getLeastAdmissableColor(g, coloring, maxColor, e, v,
                                            nodesQueue, partialOrder, ranks);
      // update partial order
      const qc::Qubit u     = e.first == v ? e.second : e.first;
      ranks[u][coloring[e]] = static_cast<std::size_t>(
          std::distance(nodesQueue.cbegin(),
                        std::find(nodesQueue.cbegin(), nodesQueue.cend(), v)));
      for (const auto& [f, k] : coloring) {
        if (f.first == v or f.second == v) {
          const qc::Qubit w = f.first == v ? f.second : f.first;
          if (k < coloring[e]) {
            partialOrder.addEdge(w, u);
          } else if (k > coloring[e]) {
            partialOrder.addEdge(u, w);
          }
        } else if (k == coloring[e]) {
          const qc::Qubit w = std::find(nodesQueue.cbegin(), nodesQueue.cend(),
                                        f.first) == nodesQueue.cend()
                                  ? f.first
                                  : f.second;
          if (ranks[u][k] < ranks[w][k]) {
            partialOrder.addEdge(w, u);
          } else if (ranks[u][k] > ranks[w][k]) {
            partialOrder.addEdge(u, w);
          } else {
            throw std::logic_error("Coloring is not valid.");
          }
        }
      }

      maxColor = std::max(maxColor, coloring[e]);
      // update the number of distinct colors of adjacent edges
      for (const auto& f : edges) {
        if (g.isAdjacentEdge(e, f)) {
          std::vector<bool> usedColors(maxColor + 1, false);
          for (const auto& h : edges) {
            if (g.isAdjacentEdge(f, h) and coloring.find(h) != coloring.end()) {
              usedColors[coloring[h]] = true;
            }
          }
          nAdjColors[f] = static_cast<std::size_t>(
              std::count(usedColors.cbegin(), usedColors.cend(),
                         [](const bool& b) { return b; }));
        }
      }
    }
  }
  return {coloring, partialOrder};
}

auto NAGraphAlgorithms::computeRestingPositions(
    const std::vector<qc::Qubit>& moveable, const std::vector<qc::Qubit>& fixed,
    const std::unordered_map<Edge, Color, qc::PairHash<qc::Qubit>>& coloring)
    -> std::vector<std::size_t> {
  const Color maxColor = std::accumulate(
      coloring.cbegin(), coloring.cend(), static_cast<Color>(0),
      [](Color acc, const auto& value) { return std::max(acc, value.second); });
  std::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t,
                     qc::PairHash<std::size_t>>
      resting{};
  for (Color t = 0; t <= maxColor; ++t) {
    // resting positions for timestamp t
    std::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t,
                       qc::PairHash<std::size_t>>
        tResting{};
    // x-coordinates of movable vertices at timestamp t
    std::unordered_map<qc::Qubit, std::size_t> moveableXs{};
    for (const qc::Qubit v : moveable) {
      // get the neighbors of v
      const std::set<qc::Qubit> neighborsOfV = std::accumulate(
          coloring.cbegin(), coloring.cend(), std::set<qc::Qubit>(),
          [&](std::set<qc::Qubit> acc, auto value) {
            Edge e = value.first;
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
        // get index of u in fixed which is the x-coordinate \wo resting
        const auto& it = std::find(fixed.cbegin(), fixed.cend(), u);
        assert(it != fixed.cend());
        // set x-pos of v to x-pos of u, which is the index of u in fixed
        moveableXs[v] =
            static_cast<std::size_t>(std::distance(fixed.cbegin(), it));
      }
    }
    // map the keys of moveableXs to their indices in moveable
    std::vector<std::size_t> moveableXsIds{};
    std::transform(
        moveableXs.cbegin(), moveableXs.cend(),
        std::back_inserter(moveableXsIds), [&](const auto& value) {
          const auto& it =
              std::find(moveable.cbegin(), moveable.cend(), value.first);
          return static_cast<std::size_t>(std::distance(moveable.cbegin(), it));
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
        // resting position
        if (!leftNeighbors.empty() and !rightNeighbors.empty()) {
          // get min left and max right neighbor
          const auto& leftNeighbor =
              std::min_element(leftNeighbors.cbegin(), leftNeighbors.cend());
          const auto& rightNeighbor =
              std::max_element(rightNeighbors.cbegin(), rightNeighbors.cend());
          const auto& pair = std::pair<std::int64_t, std::int64_t>(
              moveableXs[moveable[*leftNeighbor]],
              moveableXs[moveable[*rightNeighbor]]);
          if (tResting.find(pair) == tResting.end()) {
            tResting[pair] = 1;
          } else {
            tResting[pair]++;
          }
        }
      }
    }
    // tResting is completely computed until here
    // merge tResting with resting
    std::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t,
                       qc::PairHash<std::size_t>>
        newResting;
    for (const auto& elem : resting) {
      const auto& pair = elem.first;
      for (std::size_t i = 0; i < elem.second; ++i) {
        // check whether existing resting positions overlap with newly required
        // ones
        const auto& overlappingPairs =
            std::accumulate(tResting.cbegin(), tResting.cend(),
                            std::vector<std::pair<std::size_t, std::size_t>>(),
                            [&](auto& acc, const auto& value) {
                              if (value.first.first < pair.second &&
                                  pair.first < value.first.second) {
                                acc.emplace_back(value.first);
                              }
                              return acc;
                            });
        if (!overlappingPairs.empty()) {
          // select the smallest one to not loose a lot flexibility
          const std::pair<std::size_t, std::size_t>& minSuperPair =
              *(std::min_element(
                  overlappingPairs.cbegin(), overlappingPairs.cend(),
                  [](const auto& a, const auto& b) {
                    return a.second - a.first < b.second - b.first;
                  }));
          --tResting[minSuperPair];
          if (tResting[minSuperPair] == 0) {
            tResting.erase(minSuperPair);
          }
          const std::pair<std::size_t, std::size_t> newPair{
              std::max(pair.first, minSuperPair.first),
              std::min(pair.second, minSuperPair.second)};
          if (newResting.find(newPair) == newResting.end()) {
            newResting[newPair] = 1;
          } else {
            ++newResting[newPair];
          }
        } else {
          if (newResting.find(pair) == newResting.end()) {
            newResting[pair] = 1;
          } else {
            ++newResting[pair];
          }
        }
      }
    }
    for (const auto& [pair, s] : tResting) {
      if (newResting.find(pair) == newResting.end()) {
        newResting[pair] = s;
      } else {
        newResting[pair] += s;
      }
    }
    resting.clear();
    resting = newResting;
  }
  std::vector<std::size_t> restingPositions{};
  for (const auto& [pair, s] : resting) {
    for (std::size_t i = 0; i < s; ++i) {
      restingPositions.emplace_back(pair.first);
    }
  }
  std::sort(restingPositions.begin(), restingPositions.end());
  return restingPositions;
}

auto NAGraphAlgorithms::groupByConnectedComponent(
    const InteractionGraph& g, const std::vector<qc::Qubit>& sequence)
    -> std::vector<qc::Qubit> {
  const auto&                vertices = g.getVertices();
  qc::DisjointSet<qc::Qubit> ds(vertices.cbegin(), vertices.cend());
  for (const qc::Qubit& v : vertices) {
    for (const qc::Qubit& u : vertices) {
      if (g.isAdjacent(v, u)) {
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

auto NAGraphAlgorithms::computeSequence(const InteractionGraph& g)
    -> std::pair<std::vector<std::unordered_map<qc::Qubit, std::int64_t>>,
                 std::unordered_map<qc::Qubit, std::int64_t>> {
  const auto& mis = getMaxIndependentSet(g);
  std::vector sequenceUngrouped(mis.cbegin(), mis.cend());
  // sort the vertices by degree in descending order
  std::sort(sequenceUngrouped.begin(), sequenceUngrouped.end(),
            [&](const auto& u, const auto& v) {
              return g.getDegree(u) > g.getDegree(v);
            });
  const auto& sequence = groupByConnectedComponent(g, sequenceUngrouped);
  const auto& [coloring, partialOrder] =
      colorEdges(g, coveredEdges(g, mis), sequence);
  const auto& fixed   = partialOrder.orderTopologically();
  const auto& resting = computeRestingPositions(sequence, fixed, coloring);
  // compute relative x positions of fixed vertices
  std::unordered_map<qc::Qubit, std::int64_t> fixedPositions{};
  for (std::uint32_t x = 0, i = 0; x < fixed.size(); ++x) {
    fixedPositions[fixed[x]] = static_cast<std::int64_t>(x) + i;
    while (i < resting.size() and x == resting[i]) {
      ++i;
    }
  }
  // compute relative x positions of moveable vertices at every timestamp
  const Color maxColor = std::accumulate(
      coloring.cbegin(), coloring.cend(), static_cast<Color>(0),
      [](Color acc, const auto& value) { return std::max(acc, value.second); });
  std::vector<std::unordered_map<qc::Qubit, std::int64_t>> moveablePositions(
      maxColor + 1);
  for (Color t = 0; t <= maxColor; ++t) {
    for (const qc::Qubit v : sequence) {
      // get the neighbors of v at timestamp t
      const std::set<qc::Qubit>& neighborsOfV = std::accumulate(
          coloring.cbegin(), coloring.cend(), std::set<qc::Qubit>(),
          [&](std::set<qc::Qubit> acc, auto value) {
            if (value.second == t) {
              Edge e = value.first;
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
        assert(it != fixedPositions.end());
        moveablePositions[t][v] = it->second;
      }
    }
    for (std::size_t i = 0; i < sequence.size(); ++i) {
      if (const qc::Qubit v = sequence[i];
          moveablePositions[t].find(v) == moveablePositions[t].end()) {
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
} // namespace na
