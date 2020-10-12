/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#include "Mapper.hpp"
#include "heuristic/unique_priority_queue.hpp"

#ifndef QMAP_HEURISTICMAPPER_HPP
#define QMAP_HEURISTICMAPPER_HPP

class HeuristicMapper: public Mapper {
public:
	using Mapper::Mapper; // import constructors from parent class

	void map(const MappingSettings& ms) override;

	struct Node {
		unsigned long costFixed = 0;
		double costHeur = 0.;
		double lookaheadPenalty = 0.;
		double costTotal = 0.;
		std::array<short, MAX_DEVICE_QUBITS> qubits{}; // get qubit at specific location
		std::array<short, MAX_DEVICE_QUBITS> locations{}; // get location of specific qubit
		bool done = true;
		std::vector<std::vector<Edge>> swaps = {};
		unsigned long nswaps = 0;

		Node() = default;
		Node(const std::array<short, MAX_DEVICE_QUBITS>& q, const std::array<short, MAX_DEVICE_QUBITS>& loc, const std::vector<std::vector<Edge>>& sw = {}) {
			std::copy(q.begin(), q.end(), qubits.begin());
			std::copy(loc.begin(), loc.end(), locations.begin());
			std::copy(sw.begin(), sw.end(), std::back_inserter(swaps));
		}

		void applySWAP(const Edge& swap) {
			short q1 = qubits.at(swap.first);
			short q2 = qubits.at(swap.second);

			qubits.at(swap.first) = q2;
			qubits.at(swap.second) = q1;

			if (q1 != -1) {
				locations.at(q1) = swap.second;
			}
			if (q2 != -1) {
				locations.at(q2) = swap.first;
			}
			swaps.back().emplace_back(swap);
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
	unique_priority_queue<Node> nodes {};

	void initResults() override;

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
		if (settings.admissibleHeuristic) {
			return std::max(currentCost, newCost);
		} else {
			return currentCost + newCost;
		}
	}
};

inline bool operator<(const HeuristicMapper::Node& x, const HeuristicMapper::Node& y) {
	auto itx = x.qubits.begin();
	auto ity = y.qubits.begin();
	while(itx != x.qubits.end() && ity != y.qubits.end()) {
		if (*itx != *ity) {
			return *itx < *ity;
		}
		++itx;
		++ity;
	}
	return false;
}

inline bool operator>(const HeuristicMapper::Node& x, const HeuristicMapper::Node& y) {
	auto xcost = x.costTotal + x.costFixed + x.lookaheadPenalty;
	auto ycost = y.costTotal + y.costFixed + y.lookaheadPenalty;
	if (xcost != ycost) {
		return xcost > ycost;
	}

	if(x.done){
		return false;
	}
	if(y.done) {
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
