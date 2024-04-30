//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/mqt-qmap for more
// information.
//

#pragma once

#include "Definitions.hpp"
#include "datastructures/DirectedAcyclicGraph.hpp"
#include "datastructures/Layer.hpp"
#include "datastructures/UndirectedGraph.hpp"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace na {

using Color            = std::uint16_t;
using Edge             = std::pair<qc::Qubit, qc::Qubit>;
using InteractionGraph = qc::Layer::InteractionGraph;

class NAGraphAlgorithms {
public:
  [[nodiscard]] static auto getMaxIndependentSet(const InteractionGraph& g)
      -> std::unordered_set<qc::Qubit>;

  [[nodiscard]] static auto
  coveredEdges(const InteractionGraph&              g,
               const std::unordered_set<qc::Qubit>& vs)
      -> std::unordered_set<Edge, qc::PairHash<qc::Qubit, qc::Qubit>>;

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
  [[nodiscard]] static auto getLeastAdmissableColor(
      const InteractionGraph& g,
      const std::unordered_map<Edge, Color, qc::PairHash<qc::Qubit, qc::Qubit>>&
                   coloring,
      const Color& maxColor, const Edge& e, const qc::Qubit& v,
      const std::vector<qc::Qubit>&              sequence,
      const qc::DirectedAcyclicGraph<qc::Qubit>& partialOrder,
      std::unordered_map<qc::Qubit, std::unordered_map<Color, std::size_t>>
          ranks) -> Color;

  /**
   * @brief Colors all given edges starting with edges that are adjacent to the
   * first vertex in the queue.
   *
   * @param edges e.g. edges adjacent to a selected vertex
   * @param nodesQueue could for example be sorted by degree, highest first
   * @return
   */
  [[nodiscard]] static auto colorEdges(
      const InteractionGraph&                                             g,
      const std::unordered_set<Edge, qc::PairHash<qc::Qubit, qc::Qubit>>& edges,
      const std::vector<qc::Qubit>& nodesQueue)
      -> std::pair<
          std::unordered_map<Edge, Color, qc::PairHash<qc::Qubit, qc::Qubit>>,
          qc::DirectedAcyclicGraph<qc::Qubit>>;

  [[nodiscard]] static auto computeRestingPositions(
      const std::vector<qc::Qubit>& moveable,
      const std::vector<qc::Qubit>& fixed,
      const std::unordered_map<Edge, Color, qc::PairHash<qc::Qubit, qc::Qubit>>&
          coloring) -> std::vector<std::size_t>;

  [[nodiscard]] static auto
  groupByConnectedComponent(const InteractionGraph&       g,
                            const std::vector<qc::Qubit>& sequence)
      -> std::vector<qc::Qubit>;

  /**
   * @brief Partitions the set of vertices into moveable and fixed vertices
   * with the aim to maximize the number of executable gates in one run
   * without reloading atoms.
   *
   * @return A triple of vectors containing the moveable vertices, the fixed
   * vertices and resting positions between the fixed vertices.
   */
  [[nodiscard]] static auto computeSequence(const InteractionGraph& g)
      -> std::pair<std::vector<std::unordered_map<qc::Qubit, std::int64_t>>,
                   std::unordered_map<qc::Qubit, std::int64_t>>;
};
} // namespace na
