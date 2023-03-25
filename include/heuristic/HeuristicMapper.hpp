//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "Mapper.hpp"
#include "heuristic/UniquePriorityQueue.hpp"

#pragma once

/**
 * number of two-qubit gates acting on pairs of logical qubits in some layer
 * where the keys correspond to logical qubit pairs ({q1, q2}, with q1<=q2)
 * and the values to the number of gates acting on a pair in each direction
 * (the first number with control=q1, target=q2 and the second the reverse).
 *
 * e.g., with multiplicity {{0,1},{2,3}} there are 2 gates with logical
 * qubit 0 as control and qubit 1 as target, and 3 gates with 1 as control
 * and 0 as target.
 */
using TwoQubitMultiplicity =
    std::map<Edge, std::pair<std::uint16_t, std::uint16_t>>;

class HeuristicMapper : public Mapper {
public:
  using Mapper::Mapper; // import constructors from parent class

  /**
   * @brief map the circuit passed at initialization to the architecture
   *
   * @param config the settings for this mapping run (controls e.g. layering
   * methods, pre- and post-optimizations, etc.)
   */
  void map(const Configuration& configuration) override;

  /**
   * @brief struct representing one node in the A* search containing info about
   * swaps, mappings and costs
   */
  struct Node {
    /** current fixed cost (for non-fidelity-aware mapping cost of all swaps
     * already added) */
    double costFixed = 0;
    /** heuristic cost expected for future swaps needed in current circuit layer
     */
    double costHeur = 0.;
    /** heuristic cost expected for future swaps needed in later circuit layers
     * (further layers contribute less) */
    double lookaheadPenalty = 0.;
    /**
     * containing the logical qubit currently mapped to each physical qubit.
     * `qubits[physical_qubit] = logical_qubit`
     *
     * The inverse of `locations`
     */
    std::array<std::int16_t, MAX_DEVICE_QUBITS> qubits{};
    /**
     * containing the logical qubit currently mapped to each physical qubit.
     * `locations[logical_qubit] = physical_qubit`
     *
     * The inverse of `qubits`
     */
    std::array<std::int16_t, MAX_DEVICE_QUBITS> locations{};
    /** true if all qubit pairs are mapped next to each other on the
     * architecture */
    bool done = true;
    /** swaps used to get from mapping after last layer to the current mapping;
     * each search node begins a new entry in the outer vector */
    std::vector<std::vector<Exchange>> swaps = {};
    /** number of swaps used to get from mapping after last layer to the current
     * mapping */
    std::size_t nswaps = 0;

    Node() = default;
    Node(const std::array<std::int16_t, MAX_DEVICE_QUBITS>& q,
         const std::array<std::int16_t, MAX_DEVICE_QUBITS>& loc,
         const std::vector<std::vector<Exchange>>&          sw            = {},
         const double                                       initCostFixed = 0)
        : costFixed(initCostFixed) {
      std::copy(q.begin(), q.end(), qubits.begin());
      std::copy(loc.begin(), loc.end(), locations.begin());
      std::copy(sw.begin(), sw.end(), std::back_inserter(swaps));
    }

    /**
     * @brief returns costFixed + costHeur + lookaheadPenalty
     */
    [[nodiscard]] double getTotalCost() const {
      return costFixed + costHeur + lookaheadPenalty;
    }

    /**
     * @brief returns costFixed + lookaheadPenalty
     */
    [[nodiscard]] double getTotalFixedCost() const {
      return costFixed + lookaheadPenalty;
    }

    /**
     * @brief applies an in-place swap of 2 qubits in `qubits` and `locations`
     * of the node
     */
    void applySWAP(const Edge& swap, Architecture& arch);

    /**
     * @brief applies an in-place teleportation of 2 qubits in `qubits` and
     * `locations` of the node
     */
    void applyTeleportation(const Edge& swap, Architecture& arch);

    /**
     * @brief recalculates the fixed cost of the node from current mapping and
     * swaps
     *
     * @param arch the architecture for calculating distances between physical
     * qubits and supplying qubit information such as fidelity
     */
    void recalculateFixedCost(const Architecture& arch);

    /**
     * @brief calculates the heuristic cost of the current mapping in the node
     * for some given layer and writes it to `Node::costHeur` additional
     * `Node::done` is set to true if all qubits shared by a gate in the layer
     * are mapped next to each other
     *
     * @param arch the architecture for calculating distances between physical
     * qubits and supplying qubit information such as fidelity
     * @param twoQubitGateMultiplicity number of two qubit gates acting on pairs
     * of logical qubits in the current layer
     * @param admissibleHeuristic controls if the heuristic should be calculated
     * such that it is admissible (i.e. A*-search should yield the optimal
     * solution using this heuristic)
     */
    void
    updateHeuristicCost(const Architecture&         arch,
                        const TwoQubitMultiplicity& twoQubitGateMultiplicity,
                        bool                        admissibleHeuristic);

    std::ostream& print(std::ostream& out) const {
      out << "{\n";
      out << "\t\"done\": " << done << ",\n";
      out << "\t\"cost\": {\n";
      out << "\t\t\"fixed\": " << costFixed << ",\n";
      out << "\t\t\"heuristic\": " << costHeur << ",\n";
      out << "\t\t\"lookahead_penalty\": " << lookaheadPenalty << "\n";
      out << "\t},\n";
      out << "\t\"nswaps\": " << nswaps << "\n}\n";
      return out;
    }
  };

protected:
  UniquePriorityQueue<Node> nodes{};

  /**
   * @brief creates an initial mapping of logical qubits to physical qubits with
   * different methods depending on `Mapper::results.config.initialLayout`
   */
  virtual void createInitialMapping();

  /**
   * @brief statically creates an initial mapping of logical qubits to physical
   * qubits by considering qubits that share a gate in the first layer and
   * mapping those to any free connected qubit pair in the architecture. The
   * remaining qubits are then just mapped by order of index.
   */
  virtual void staticInitialMapping();

  /**
   * @brief returns distance of the given logical qubit pair according to the
   * current mapping
   */
  double distanceOnArchitectureOfLogicalQubits(std::uint16_t control,
                                               std::uint16_t target) {
    return architecture.distance(
        static_cast<std::uint16_t>(locations.at(control)),
        static_cast<std::uint16_t>(locations.at(target)));
  }

  /**
   * @brief returns distance of the given physical qubit pair on the
   * architecture
   */
  double distanceOnArchitectureOfPhysicalQubits(std::uint16_t control,
                                                std::uint16_t target) {
    return architecture.distance(control, target);
  }

  /**
   * @brief map the logical qubit `target` to a free physical qubit, that is
   * nearest to the physical qubit `source` is mapped to
   *
   * @param source an already mapped logical qubit, which should be mapped near
   * to `target`
   * @param target an unmapped logical qubit
   */
  virtual void mapToMinDistance(std::uint16_t source, std::uint16_t target);

  /**
   * @brief maps any yet unmapped qubits, which are acted on in a given layer,
   * to a physical qubit.
   *
   * @param twoQubitGateMultiplicity number of two qubit gates acting on pairs
   * of logical qubits in the current layer
   */
  virtual void
  mapUnmappedGates(const TwoQubitMultiplicity& twoQubitGateMultiplicity);

  /**
   * @brief search for an optimal mapping/set of swaps using A*-search and the
   * heuristic specified in `HeuristicMapper::Node::updateHeuristicCost`
   *
   * uses `HeuristicMapper::nodes` as a priority queue for the A*-search,
   * assumed to be empty (or at least containing only nodes compliant with the
   * current layer in their fields `costHeur` and `done`)
   *
   * @param layer index of the current circuit layer
   */
  virtual Node aStarMap(std::size_t layer);

  /**
   * @brief expand the given node by calling `expand_node_add_one_swap` for all
   * possible swaps, which creates new search nodes and adds them to
   * `HeuristicMapper::nodes`
   *
   * @param consideredQubits set of all qubits that are acted on by a
   * 2-qubit-gate in the respective layer
   * @param node current search node
   * @param layer index of current circuit layer
   * @param twoQubitGateMultiplicity number of two qubit gates acting on pairs
   * of logical qubits in the current layer
   */
  void expandNode(const std::unordered_set<std::uint16_t>& consideredQubits,
                  Node& node, std::size_t layer,
                  const TwoQubitMultiplicity& twoQubitGateMultiplicity);

  /**
   * @brief creates a new node with a swap on the given edge and adds it to
   * `HeuristicMapper::nodes`
   *
   * @param swap edge on which to perform a swap
   * @param node current search node
   * @param layer index of current circuit layer
   * @param twoQubitGateMultiplicity number of two qubit gates acting on pairs
   * of logical qubits in the current layer
   */
  void
  expandNodeAddOneSwap(const Edge& swap, Node& node, std::size_t layer,
                       const TwoQubitMultiplicity& twoQubitGateMultiplicity);

  /**
   * @brief calculates the heuristic cost for the following layers and saves it
   * in the node as `lookaheadPenalty`
   *
   * @param layer index of current circuit layer
   * @param node search node for which to calculate lookahead penalty
   */
  void lookahead(std::size_t layer, Node& node);

  double heuristicAddition(const double currentCost, const double newCost) {
    if (results.config.admissibleHeuristic) {
      return std::max(currentCost, newCost);
    }
    return currentCost + newCost;
  }
};

inline bool operator<(const HeuristicMapper::Node& x,
                      const HeuristicMapper::Node& y) {
  auto itx = x.qubits.begin(); // NOLINT (readability-qualified-auto)
  auto ity = y.qubits.begin(); // NOLINT (readability-qualified-auto)
  while (itx != x.qubits.end() && ity != y.qubits.end()) {
    if (*itx != *ity) {
      return *itx < *ity;
    }
    ++itx;
    ++ity;
  }
  return false;
}

inline bool operator>(const HeuristicMapper::Node& x,
                      const HeuristicMapper::Node& y) {
  const auto xcost = x.getTotalCost();
  const auto ycost = y.getTotalCost();
  if (std::abs(xcost - ycost) > 1e-6) {
    return xcost > ycost;
  }

  if (x.done) {
    return false;
  }
  if (y.done) {
    return true;
  }

  const auto xheur = x.costHeur + x.lookaheadPenalty;
  const auto yheur = y.costHeur + y.lookaheadPenalty;
  if (std::abs(xheur - yheur) > 1e-6) {
    return xheur > yheur;
  }
  return x < y;
}
