/*
 * This file is part of the MQT QMAP library which is released under the MIT license.
 * See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
 */

#include "Mapper.hpp"
#include "heuristic/unique_priority_queue.hpp"

#ifndef QMAP_HEURISTICMAPPER_HPP
    #define QMAP_HEURISTICMAPPER_HPP

class HeuristicMapper: public Mapper {
public:
    using Mapper::Mapper; // import constructors from parent class

    /**
     * @brief map the circuit passed at initialization to the architecture
     *
     * @param config the settings for this mapping run (controls e.g. layering methods, pre- and post-optimizations, etc.)
     */
    void map(const Configuration& config) override;

    /**
     * @brief struct representing one node in the A* search containing info about swaps, mappings and costs
     */
    struct Node {
        /** cost of all swaps in the node */
        unsigned long costFixed = 0;
        /** heuristic cost expected for future swaps needed in current circuit layer */
        double costHeur = 0.;
        /** heuristic cost expected for future swaps needed in later circuit layers (further layers contribute less) */
        double lookaheadPenalty = 0.;
        double costTotal        = 0.;
        /**
         * containing the logical qubit currently mapped to each physical qubit.
         * `qubits[physical_qubit] = logical_qubit`
         *
         * The inverse of `locations`
         */
        std::array<short, MAX_DEVICE_QUBITS> qubits{};
        /**
         * containing the logical qubit currently mapped to each physical qubit.
         * `locations[logical_qubit] = physical_qubit`
         *
         * The inverse of `qubits`
         */
        std::array<short, MAX_DEVICE_QUBITS> locations{};
        /** true if all qubit pairs are mapped next to each other on the architecture */
        bool done = true;
        /** swaps used to get from mapping after last layer to the current mapping;
         * each search node begins a new entry in the outer vector */
        std::vector<std::vector<Exchange>> swaps = {};
        /** number of swaps used to get from mapping after last layer to the current mapping */
        unsigned long nswaps = 0;

        Node() = default;
        Node(const std::array<short, MAX_DEVICE_QUBITS>& q, const std::array<short, MAX_DEVICE_QUBITS>& loc, const std::vector<std::vector<Exchange>>& sw = {}) {
            std::copy(q.begin(), q.end(), qubits.begin());
            std::copy(loc.begin(), loc.end(), locations.begin());
            std::copy(sw.begin(), sw.end(), std::back_inserter(swaps));
        }

        /**
         * @brief applies an in-place swap of 2 qubits in `qubits` and `locations` of the node
         */
        void applySWAP(const Edge& swap, Architecture& arch) {
            short q1 = qubits.at(swap.first);
            short q2 = qubits.at(swap.second);

            qubits.at(swap.first)  = q2;
            qubits.at(swap.second) = q1;

            if (q1 != -1) {
                locations.at(q1) = static_cast<short>(swap.second);
            }
            if (q2 != -1) {
                locations.at(q2) = static_cast<short>(swap.first);
            }

            if (arch.getCouplingMap().find(swap) != arch.getCouplingMap().end() || arch.getCouplingMap().find(Edge{swap.second, swap.first}) != arch.getCouplingMap().end()) {
                swaps.back().emplace_back(swap.first, swap.second, qc::SWAP);
            } else {
                throw QMAPException("Something wrong in applySWAP.");
            }
        }

        /**
         * @brief applies an in-place teleportation of 2 qubits in `qubits` and `locations` of the node
         */
        void applyTeleportation(const Edge& swap, Architecture& arch) {
            short q1 = qubits.at(swap.first);
            short q2 = qubits.at(swap.second);

            qubits.at(swap.first)  = q2;
            qubits.at(swap.second) = q1;

            if (q1 != -1) {
                locations.at(q1) = static_cast<short>(swap.second);
            }
            if (q2 != -1) {
                locations.at(q2) = static_cast<short>(swap.first);
            }

            unsigned short middle_anc = std::numeric_limits<decltype(middle_anc)>::max();
            for (const auto& qpair: arch.getTeleportationQubits()) {
                if (swap.first == qpair.first || swap.second == qpair.first) {
                    middle_anc = qpair.second;
                } else if (swap.first == qpair.second || swap.second == qpair.second) {
                    middle_anc = qpair.first;
                }
            }

            if (middle_anc == std::numeric_limits<decltype(middle_anc)>::max()) {
                throw QMAPException("Teleportation between seemingly wrong qubits: " + std::to_string(swap.first) + " <--> " + std::to_string(swap.second));
            }

            unsigned short source = std::numeric_limits<decltype(source)>::max();
            unsigned short target = std::numeric_limits<decltype(target)>::max();
            if (arch.getCouplingMap().find({swap.first, middle_anc}) != arch.getCouplingMap().end() || arch.getCouplingMap().find({middle_anc, swap.first}) != arch.getCouplingMap().end()) {
                source = swap.first;
                target = swap.second;
            } else {
                source = swap.second;
                target = swap.first;
            }

            if (source == middle_anc || target == middle_anc) {
                std::clog << "FAIL: TELE " << source << " -(" << middle_anc << ")-> " << target << "\n";
                throw QMAPException("Overlap between source/target and middle ancillary in teleportation.");
            }

            swaps.back().emplace_back(source, target, middle_anc, qc::Teleportation);
        }

        void updateHeuristicCost(const Architecture& arch, const Gate& gate, bool admissibleHeuristic) {
            auto cost = arch.distance(locations.at(gate.control), locations.at(gate.target));
            if (admissibleHeuristic) {
                costHeur = std::max(costHeur, cost);
            } else {
                costHeur += cost;
            }
        }

        /**
         * @brief checks if the qubits of the given gate are mapped next to each other. If not, set `done` to false
         */
        void checkUnfinished(const Architecture& arch, const Gate& gate) {
            if (arch.distance(locations.at(gate.control), locations.at(gate.target)) > COST_DIRECTION_REVERSE) {
                done = false;
            }
        }

        std::ostream& print(std::ostream& out) const {
            out << "{\n";
            out << "\t\"done\": " << done << ",\n";
            out << "\t\"cost\": {\n";
            out << "\t\t\"fixed\": " << costFixed << ",\n";
            out << "\t\t\"heuristic\": " << costHeur << ",\n";
            out << "\t\t\"total\": " << costTotal << ",\n";
            out << "\t\t\"lookahead_penalty\": " << lookaheadPenalty << "\n";
            out << "\t},\n";
            out << "\t\"nswaps\": " << nswaps << "\n}\n";
            return out;
        }
    };

protected:
    unique_priority_queue<Node> nodes{};

    /**
     * @brief creates an initial mapping of logical qubits to physical qubits with different methods depending on
     * `Mapper::results.config.initialLayout`
     */
    virtual void createInitialMapping();

    /**
     * @brief statically creates an initial mapping of logical qubits to physical qubits by considering qubits that
     * share a gate in the first layer and mapping those to any free connected qubit pair in the architecture.
     * The remaining qubits are then just mapped by order of index.
     */
    virtual void staticInitialMapping();

    /**
     * @brief returns distance of the given logical qubit pair according to the current mapping
     */
    double distanceOnArchitectureOfLogicalQubits(unsigned short control, unsigned short target) {
        return architecture.distance(locations.at(control), locations.at(target));
    }

    /**
     * @brief returns distance of the given physical qubit pair on the architecture
     */
    double distanceOnArchitectureOfPhysicalQubits(unsigned short control, unsigned short target) {
        return architecture.distance(control, target);
    }

    /**
     * @brief map the logical qubit `target` to a free physical qubit, that is nearest to the physical qubit `source` is mapped to
     *
     * @param source an already mapped logical qubit, which should be mapped near to `target`
     * @param target an unmapped logical qubit
     */
    virtual void mapToMinDistance(unsigned short source, unsigned short target);

    /**
     * @brief gathers all qubits that are acted on by a 2-qubit-gate in the given layer in `consideredQubits`,
     * and maps any of them that are not yet mapped to a physical qubit.
     *
     * All gates are mapped in order of their index in the layer. The qubits are mapped to any 2 qubits with minimal distance on the architecture.
     *
     * Additionally sets the fields `costHeur` (maximum distance between any 2 qubits which share a gate) and `done` (all qubit considered pairs
     * are mapped next to each other) in the current search node.
     *
     * @param layer index of the circuit layer to consider
     * @param node current AStar search node
     * @param consideredQubits vector in which to gather all relevant qubits of this layer
     */
    virtual void mapUnmappedGates(long layer, Node& node, std::vector<unsigned short>& consideredQubits);

    /**
     * @brief search for an optimal mapping/set of swaps using A*-search and the heuristic specified in `HeuristicMapper::Node::updateHeuristicCost`
     *
     * uses `HeuristicMapper::nodes` as a priority queue for the A*-search, assumed to be empty (or at least containing only nodes
     * compliant with the current layer in their fields `costHeur` and `done`)
     *
     * @param layer index of the current circuit layer
     */
    virtual Node AstarMap(long layer);

    /**
     * @brief expand the given node by calling `expand_node_add_one_swap` for all possible swaps, which creates new search nodes and adds them to `HeuristicMapper::nodes`
     *
     * @param consideredQubits set of all qubits that are acted on by a 2-qubit-gate in the respective layer
     * @param node current search node
     * @param layer index of current circuit layer
     */
    void expandNode(const std::vector<unsigned short>& consideredQubits, Node& node, long layer);

    /**
     * @brief creates a new node with a swap on the given edge and adds it to `HeuristicMapper::nodes`
     *
     * @param swap edge on which to perform a swap
     * @param node current search node
     * @param layer index of current circuit layer
     */
    void expand_node_add_one_swap(const Edge& swap, Node& node, long layer);

    /**
     * @brief calculates the heuristic cost for the following layers and saves it in the node as `lookaheadPenalty`
     *
     * @param layer index of current circuit layer
     * @param node search node for which to calculate lookahead penalty
     */
    void lookahead(long layer, Node& node);

    // TODO: also use in `HeuristicMapper::mapUnmappedGates` and `HeuristicMapper::Node::updateHeuristicCost`
    double heuristicCost(double currentCost, double newCost) {
        if (results.config.admissibleHeuristic) {
            return std::max(currentCost, newCost);
        } else {
            return currentCost + newCost;
        }
    }
};

inline bool operator<(const HeuristicMapper::Node& x, const HeuristicMapper::Node& y) {
    auto itx = x.qubits.begin();
    auto ity = y.qubits.begin();
    while (itx != x.qubits.end() && ity != y.qubits.end()) {
        if (*itx != *ity) {
            return *itx < *ity;
        }
        ++itx;
        ++ity;
    }
    return false;
}

inline bool operator>(const HeuristicMapper::Node& x, const HeuristicMapper::Node& y) {
    const auto xcost = x.costTotal + static_cast<double>(x.costFixed) + x.lookaheadPenalty;
    const auto ycost = y.costTotal + static_cast<double>(y.costFixed) + y.lookaheadPenalty;
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
    } else {
        return x < y;
    }
}

#endif //QMAP_HEURISTICMAPPER_HPP
