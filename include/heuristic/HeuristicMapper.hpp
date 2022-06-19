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

    void map(const Configuration& config) override;

    struct Node {
        unsigned long                        costFixed        = 0;
        double                               costHeur         = 0.;
        double                               lookaheadPenalty = 0.;
        double                               costTotal        = 0.;
        std::array<short, MAX_DEVICE_QUBITS> qubits{};    // get qubit at specific location
        std::array<short, MAX_DEVICE_QUBITS> locations{}; // get location of specific qubit
        bool                                 done   = true;
        std::vector<std::vector<Exchange>>   swaps  = {};
        unsigned long                        nswaps = 0;

        Node() = default;
        Node(const std::array<short, MAX_DEVICE_QUBITS>& q, const std::array<short, MAX_DEVICE_QUBITS>& loc, const std::vector<std::vector<Exchange>>& sw = {}) {
            std::copy(q.begin(), q.end(), qubits.begin());
            std::copy(loc.begin(), loc.end(), locations.begin());
            std::copy(sw.begin(), sw.end(), std::back_inserter(swaps));
        }

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

    virtual void createInitialMapping();

    double distanceOnArchitectureOfLogicalQubits(unsigned short control, unsigned short target) {
        return architecture.distance(locations.at(control), locations.at(target));
    }

    double distanceOnArchitectureOfPhysicalQubits(unsigned short control, unsigned short target) {
        return architecture.distance(control, target);
    }

    virtual void mapToMinDistance(unsigned short source, unsigned short target);

    virtual void mapUnmappedGates(long layer, Node& node, std::vector<unsigned short>& consideredQubits);

    virtual Node AstarMap(long layer);

    void expandNode(const std::vector<unsigned short>& consideredQubits, Node& node, long layer);
    void expand_node_add_one_swap(const Edge& swap, Node& node, long layer);

    void lookahead(long layer, Node& node);

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
    auto xcost = x.costTotal + static_cast<double>(x.costFixed) + x.lookaheadPenalty;
    auto ycost = y.costTotal + static_cast<double>(y.costFixed) + y.lookaheadPenalty;
    if (xcost != ycost) {
        return xcost > ycost;
    }

    if (x.done) {
        return false;
    }
    if (y.done) {
        return true;
    }

    auto xheur = x.costHeur + x.lookaheadPenalty;
    auto yheur = y.costHeur + y.lookaheadPenalty;
    if (xheur != yheur) {
        return xheur > yheur;
    } else {
        return x < y;
    }
}

#endif //QMAP_HEURISTICMAPPER_HPP
