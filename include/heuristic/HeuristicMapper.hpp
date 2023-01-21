//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "Mapper.hpp"
#include "heuristic/UniquePriorityQueue.hpp"

#pragma once

using EdgeMultiplicity =
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
    /** cost of all swaps in the node */
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
         const double                                       initCostFixed = 0) {
      std::copy(q.begin(), q.end(), qubits.begin());
      std::copy(loc.begin(), loc.end(), locations.begin());
      std::copy(sw.begin(), sw.end(), std::back_inserter(swaps));
      costFixed = initCostFixed;
    }

    /**
     * @brief applies an in-place swap of 2 qubits in `qubits` and `locations`
     * of the node
     */
    void applySWAP(const Edge& swap, Architecture& arch,
                   const std::vector<std::uint16_t>& singleQubitGateMultiplicity,
                   bool                             considerFidelity) {
      const auto q1 = qubits.at(swap.first);
      const auto q2 = qubits.at(swap.second);

      qubits.at(swap.first)  = q2;
      qubits.at(swap.second) = q1;

      if (q1 != -1) {
        locations.at(static_cast<std::size_t>(q1)) =
            static_cast<std::int16_t>(swap.second);
      }
      if (q2 != -1) {
        locations.at(static_cast<std::size_t>(q2)) =
            static_cast<std::int16_t>(swap.first);
      }

      if (arch.getCouplingMap().find(swap) != arch.getCouplingMap().end() ||
          arch.getCouplingMap().find(Edge{swap.second, swap.first}) !=
              arch.getCouplingMap().end()) {
        swaps.back().emplace_back(swap.first, swap.second, qc::SWAP);
      } else {
        throw QMAPException("Something wrong in applySWAP.");
      }

      if (considerFidelity) {
        // accounting for fidelity difference of single qubit gates (two qubit
        // gates are handled in the heuristic)
        costFixed +=
            ((singleQubitGateMultiplicity.at(static_cast<std::size_t>(q2)) -
               singleQubitGateMultiplicity.at(static_cast<std::size_t>(q1))) *
                 arch.getSingleQubitFidelityCost(swap.first) +
             (singleQubitGateMultiplicity.at(static_cast<std::size_t>(q1)) -
              singleQubitGateMultiplicity.at(static_cast<std::size_t>(q2))) *
                 arch.getSingleQubitFidelityCost(swap.second));
        // adding cost of the swap gate itself
        costFixed += arch.getSwapFidelityCost(swap.first, swap.second);
      } else {
        if (arch.bidirectional()) {
          costFixed += COST_BIDIRECTIONAL_SWAP;
        } else {
          costFixed += COST_UNIDIRECTIONAL_SWAP;
        }
      }
    }

    /**
     * @brief applies an in-place teleportation of 2 qubits in `qubits` and
     * `locations` of the node
     */
    void applyTeleportation(const Edge& swap, Architecture& arch,
                            bool considerFidelity) {
      const auto q1 = qubits.at(swap.first);
      const auto q2 = qubits.at(swap.second);

      qubits.at(swap.first)  = q2;
      qubits.at(swap.second) = q1;

      if (q1 != -1) {
        locations.at(static_cast<std::size_t>(q1)) =
            static_cast<std::int16_t>(swap.second);
      }
      if (q2 != -1) {
        locations.at(static_cast<std::size_t>(q2)) =
            static_cast<std::int16_t>(swap.first);
      }

      std::uint16_t middleAnc = std::numeric_limits<decltype(middleAnc)>::max();
      for (const auto& qpair : arch.getTeleportationQubits()) {
        if (swap.first == qpair.first || swap.second == qpair.first) {
          middleAnc = static_cast<std::uint16_t>(qpair.second);
        } else if (swap.first == qpair.second || swap.second == qpair.second) {
          middleAnc = static_cast<std::uint16_t>(qpair.first);
        }
      }

      if (middleAnc == std::numeric_limits<decltype(middleAnc)>::max()) {
        throw QMAPException("Teleportation between seemingly wrong qubits: " +
                            std::to_string(swap.first) + " <--> " +
                            std::to_string(swap.second));
      }

      std::uint16_t source = std::numeric_limits<decltype(source)>::max();
      std::uint16_t target = std::numeric_limits<decltype(target)>::max();
      if (arch.getCouplingMap().find({swap.first, middleAnc}) !=
              arch.getCouplingMap().end() ||
          arch.getCouplingMap().find({middleAnc, swap.first}) !=
              arch.getCouplingMap().end()) {
        source = swap.first;
        target = swap.second;
      } else {
        source = swap.second;
        target = swap.first;
      }

      if (source == middleAnc || target == middleAnc) {
        std::clog << "FAIL: TELE " << source << " -(" << middleAnc << ")-> "
                  << target << "\n";
        throw QMAPException("Overlap between source/target and middle "
                            "ancillary in teleportation.");
      }

      swaps.back().emplace_back(source, target, middleAnc, qc::Teleportation);

      if (considerFidelity) {
        throw QMAPException(
            "Teleportation currently not supported for noise-aware mapping");
        // TODO: check this cost implementation and add teleportation support in
        // noise-aware heuristic accounting for fidelity difference of single
        // qubit gates (two qubit gates are handled in the heuristic)
        /*costFixed += (-(singleQubitGateMultiplicity.at(q2) -
        singleQubitGateMultiplicity.at(q1)) *
                         std::log2(singleQubitFidelities.at(swap.first))
                      -(singleQubitGateMultiplicity.at(q1) -
        singleQubitGateMultiplicity.at(q2)) *
                        std::log2(singleQubitFidelities.at(swap.second))
                      );
        // adding cost oft the swap gate
        costFixed += (-    std::log2(twoQubitFidelities.at(source).at(middleAnc))
                      -    std::log2(twoQubitFidelities.at(middleAnc).at(target))
                      -    std::log2(singleQubitFidelities.at(source))
                      -    std::log2(singleQubitFidelities.at(middleAnc))
                      -2 * std::log2(singleQubitFidelities.at(target))
                      -    std::log2(1.0 -
        arch.getProperties().readoutErrorRate.get(source))
                      -    std::log2(1.0 -
        arch.getProperties().readoutErrorRate.get(middleAnc)));*/
      } else {
        costFixed += COST_TELEPORTATION;
      }
    }

    /**
     * @brief recalculates the fixed cost of the node from current mapping and
     * swaps
     *
     * @param arch the architecture for calculating distances between physical
     * qubits and supplying qubit information such as fidelity
     * @param singleQubitGateMultiplicity vector containing the number of gates
     * acting on each logical qubit
     * @param considerFidelity controls if qubit fidelity should be taken into
     * account
     */
    void recalculateFixedCost(
        const Architecture&              arch,
        const std::vector<std::uint16_t>& singleQubitGateMultiplicity,
        bool                             considerFidelity) {

      costFixed = 0;
      if (considerFidelity) {
        // adding costs of single qubit gates
        for (std::uint16_t i = 0U; i < arch.getNqubits(); ++i) {
          if (singleQubitGateMultiplicity.at(i) == 0) {
            continue;
          }
          costFixed += singleQubitGateMultiplicity.at(i) *
                       arch.getSingleQubitFidelityCost(locations.at(i));
        }
        // adding cost of the swap gates
        for (auto& swapNode : swaps) {
          for (auto& swap : swapNode) {
            if (swap.op == qc::SWAP) {
              costFixed += arch.getSwapFidelityCost(swap.first, swap.second);
            } else if (swap.op == qc::Teleportation) {
              throw QMAPException("Teleportation currently not supported for "
                                  "noise-aware mapping");
            }
          }
        }
      } else {
        for (auto& swapNode : swaps) {
          for (auto& swap : swapNode) {
            if (swap.op == qc::SWAP) {
              if (arch.bidirectional()) {
                costFixed += COST_BIDIRECTIONAL_SWAP;
              } else {
                costFixed += COST_UNIDIRECTIONAL_SWAP;
              }
            } else if (swap.op == qc::Teleportation) {
              costFixed += COST_TELEPORTATION;
            }
          }
        }
      }
    }

    /**
     * @brief calculates the heuristic cost of the current mapping in the node
     * for some given layer and writes it to `Node::costHeur` additional
     * `Node::done` is set to true if all qubits shared by a gate in the layer
     * are mapped next to each other
     *
     * @param arch the architecture for calculating distances between physical
     * qubits and supplying qubit information such as fidelity
     * @param singleQubitGateMultiplicity vector containing the number of gates
     * acting on each logical qubit
     * @param twoQubitGateMultiplicity number of two qubit gates acting on each
     * logical qubit edge in the current layer where the first number in the
     * value pair corresponds to the number of edges having their gates given as
     * (control, target) in the key, and the second with all gates in reverse to
     * that
     * @param admissibleHeuristic controls if the heuristic should be calculated
     * such that it is admissible (i.e. A*-search should yield the optimal
     * solution using this heuristic)
     * @param considerFidelity controls if the heuristic should consider
     * fidelity data of the architecture
     */
    void updateHeuristicCost(
        const Architecture&              arch,
        const std::vector<std::uint16_t>& singleQubitGateMultiplicity,
        const EdgeMultiplicity&           twoQubitGateMultiplicity,
        bool admissibleHeuristic, bool considerFidelity) {

      costHeur = 0.;
      done     = true;
      // single qubit gate savings potential by moving it to another physical qubit with higher fidelity
      double savingsPotential = 0.;
      if (considerFidelity) {
        for (std::uint16_t log_qbit = 0U; log_qbit < arch.getNqubits();
             ++log_qbit) {
          if (singleQubitGateMultiplicity.at(log_qbit) == 0) {
            continue;
          }
          double qbitSavings     = 0;
          double currFidelity = arch.getSingleQubitFidelityCost(locations.at(log_qbit));
          for (std::uint16_t phys_qbit = 0U; phys_qbit < arch.getNqubits();
               ++phys_qbit) {
            if (arch.getSingleQubitFidelityCost(phys_qbit) >= currFidelity) {
              continue;
            }
            double curSavings =
                singleQubitGateMultiplicity.at(log_qbit) *
                    (currFidelity -
                     arch.getSingleQubitFidelityCost(phys_qbit)) -
                arch.fidelityDistance(locations.at(log_qbit), phys_qbit);
            qbitSavings = std::max(qbitSavings, curSavings);
          }
          savingsPotential += qbitSavings;
        }
      }

      // iterating over all virtual qubit pairs, that share a gate on the
      // current layer
      for (const auto& edgeMultiplicity : twoQubitGateMultiplicity) {
        const auto& q1                   = edgeMultiplicity.first.first;
        const auto& q2                   = edgeMultiplicity.first.second;
        const auto& straightMultiplicity = edgeMultiplicity.second.first;
        const auto& reverseMultiplicity  = edgeMultiplicity.second.second;
        const auto& totalMultiplicity =
            straightMultiplicity + reverseMultiplicity;

        // only if all qubit pairs are mapped next to each other the mapping
        // is complete
        if (done &&
            arch.getCouplingMap().find(
                {static_cast<std::uint16_t>(locations.at(q1)),
                 static_cast<std::uint16_t>(locations.at(q2))}) ==
                arch.getCouplingMap().end() &&
            arch.getCouplingMap().find(
                {static_cast<std::uint16_t>(locations.at(q2)),
                 static_cast<std::uint16_t>(locations.at(q1))}) ==
                arch.getCouplingMap().end()) {
          done = false;
        }

        if (considerFidelity) {
          // find the optimal edge, to which to remap the given virtual qubit
          // pair and take the cost of moving it there via swaps plus the
          // fidelity cost  of executing all their shared gates on that edge
          // as the qubit pairs cost
          double swapCost = std::numeric_limits<double>::max();
          for (const auto& edge : arch.getCouplingMap()) {
            swapCost = std::min(
                swapCost,
                straightMultiplicity * arch.getTwoQubitFidelityCost(edge.first, edge.second) +
                reverseMultiplicity * arch.getTwoQubitFidelityCost(edge.second, edge.first) +
                arch.fidelityDistance(locations.at(q1), edge.first) +
                arch.fidelityDistance(locations.at(q2),edge.second));
            swapCost = std::min(
                swapCost,
                straightMultiplicity * arch.getTwoQubitFidelityCost(edge.second, edge.first) +
                reverseMultiplicity * arch.getTwoQubitFidelityCost(edge.first, edge.second) +
                arch.fidelityDistance(locations.at(q2), edge.first) +
                arch.fidelityDistance(locations.at(q1),edge.second));
          }

          if (admissibleHeuristic) {
            costHeur = std::max(costHeur, swapCost);
          } else {
            costHeur += swapCost;
          }
        } else {
          double swapCost1 =
              arch.distance(static_cast<std::uint16_t>(locations.at(q1)),
                            static_cast<std::uint16_t>(locations.at(q2)));
          double swapCost2 =
              arch.distance(static_cast<std::uint16_t>(locations.at(q2)),
                            static_cast<std::uint16_t>(locations.at(q1)));

          if (admissibleHeuristic) {
            if (straightMultiplicity > 0) {
              costHeur = std::max(costHeur, swapCost1);
            }
            if (reverseMultiplicity > 0) {
              costHeur = std::max(costHeur, swapCost2);
            }
          } else {
            costHeur += swapCost1 * straightMultiplicity +
                        swapCost2 * reverseMultiplicity;
          }
        }
      }
      costHeur -= savingsPotential;
    }

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
   * @brief gathers all qubits that are acted on by a 2-qubit-gate in the given
   * layer in `consideredQubits`, and maps any of them that are not yet mapped
   * to a physical qubit.
   *
   * All gates are mapped in order of their index in the layer. The qubits are
   * mapped to any 2 qubits with minimal distance on the architecture.
   *
   * @param layer index of the circuit layer to consider
   * @param consideredQubits vector in which to gather all relevant qubits of
   * this layer
   */
  virtual void mapUnmappedGates(std::size_t                 layer,
                                std::vector<std::uint16_t>& consideredQubits);

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
   */
  void expandNode(const std::vector<std::uint16_t>& consideredQubits,
                  Node& node, std::size_t layer,
                  const std::vector<std::uint16_t>& singleQubitGateMultiplicity,
                  const EdgeMultiplicity&           twoQubitGateMultiplicity);

  /**
   * @brief creates a new node with a swap on the given edge and adds it to
   * `HeuristicMapper::nodes`
   *
   * @param swap edge on which to perform a swap
   * @param node current search node
   * @param layer index of current circuit layer
   */
  void expandNodeAddOneSwap(
      const Edge& swap, Node& node, std::size_t layer,
      const std::vector<std::uint16_t>& singleQubitGateMultiplicity,
      const EdgeMultiplicity&           twoQubitGateMultiplicity);

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
  const auto xcost = static_cast<double>(x.costFixed) + x.lookaheadPenalty;
  const auto ycost = static_cast<double>(y.costFixed) + y.lookaheadPenalty;
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
