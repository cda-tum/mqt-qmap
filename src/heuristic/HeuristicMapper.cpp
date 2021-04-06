/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#include "heuristic/HeuristicMapper.hpp"

#include <chrono>

void HeuristicMapper::map(const MappingSettings& ms) {
	settings = ms;
	if (settings.layeringStrategy == LayeringStrategy::OddGates || settings.layeringStrategy == LayeringStrategy::QubitTriangle) {
		std::cerr << "Layering strategy " << toString(settings.layeringStrategy) << " not suitable for heuristic mapper!" << std::endl;
		return;
	}
	const auto start = std::chrono::steady_clock::now();
	qc.stripIdleQubits(true, false);
	initResults();

	createLayers();
	if (ms.verbose) {
        std::clog << "Teleportation qubits: " << ms.teleportation_qubits << "\n";
		printLayering(std::clog);
	}

	createInitialMapping();
	if (ms.verbose) {
		printLocations(std::clog );
		printQubits(std::clog );
	}

	unsigned long gateidx = 0;
	std::vector<unsigned long> gatesToAdjust{};
	for (size_t i = 0; i < layers.size(); ++i) {
		Node result = AstarMap(i);

		qubits = result.qubits;
		locations = result.locations;

		if (ms.verbose) {
			printLocations(std::clog );
			printQubits(std::clog );
		}

		// initial layer needs no swaps
		if(i != 0) {
			for (const auto& swaps: result.swaps) {
				for (const auto& swap: swaps) {
				    if (swap.op == qc::SWAP) {
                        qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(), std::vector<qc::Control>{ }, swap.first, swap.second, qc::SWAP);
                        if (architecture.bidirectional()) {
                            results.output_gates += GATES_OF_BIDIRECTIONAL_SWAP;
                        } else {
                            results.output_gates += GATES_OF_UNIDIRECTIONAL_SWAP;
                        }
                        results.output_swaps++;
				    } else if (swap.op == qc::Teleportation){
                        qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(), std::vector<qc::Control>{}, std::vector{swap.first, swap.second, swap.middle_ancilla}, qc::Teleportation);
                        results.output_gates += COST_TELEPORTATION;
                        results.output_teleportations++;
				    }
					std::swap(qcMapped.outputPermutation.at(swap.first), qcMapped.outputPermutation.at(swap.second));

					gateidx++;
				}
			}
		}

		// add gates of the layer to circuit
		for(const auto& gate: layers.at(i)) {
			auto op = dynamic_cast<qc::StandardOperation*>(gate.op);
			if (!op) {
				throw QMAPException("Cast to StandardOperation not possible during mapping. Check that circuit contains only StandardOperations");
			}

			if (gate.singleQubit()) {
				if (locations.at(gate.target) == DEFAULT_POSITION) {
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(),
					                                             gate.target,
					                                             op->getType(),
					                                             op->getParameter().at(0),
					                                             op->getParameter().at(1),
					                                             op->getParameter().at(2));
					gatesToAdjust.push_back(gateidx);
					results.output_gates++;
					results.output_singlequbitgates++;
					gateidx++;
				} else {
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(),
                                                                 locations.at(gate.target),
                                                                 op->getType(),
                                                                 op->getParameter().at(0),
                                                                 op->getParameter().at(1),
                                                                 op->getParameter().at(2));
					results.output_gates++;
					results.output_singlequbitgates++;
					gateidx++;
				}
			} else {
				Edge cnot = {locations.at(gate.control), locations.at(gate.target)};
				if (architecture.getCouplingMap().find(cnot) == architecture.getCouplingMap().end()) {
					Edge reverse = {cnot.second, cnot.first};
					if (architecture.getCouplingMap().find(reverse) == architecture.getCouplingMap().end()) {
						throw QMAPException("Invalid CNOT: " + std::to_string(reverse.first) + "-" + std::to_string(reverse.second));
					}
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(), reverse.first, qc::H);
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(), reverse.second, qc::H);
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(), qc::Control(reverse.first), reverse.second, qc::X);
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(), reverse.second, qc::H);
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(), reverse.first, qc::H);

					results.output_direction_reverse++;
					results.output_cnots++;
					results.output_gates += 5;
					gateidx += 5;
				} else {
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(), qc::Control(cnot.first), cnot.second, qc::X);
					results.output_cnots++;
					results.output_gates++;
					gateidx++;
				}
			}
		}
	}

	// fix single qubit gates
	if (!gatesToAdjust.empty()) {
		gateidx--; // index of last operation
		for(auto it = qcMapped.rbegin(); it != qcMapped.rend(); ++it, --gateidx) {
			auto op = dynamic_cast<qc::StandardOperation*>(it->get());
			if (!op) {
				throw QMAPException("Cast to StandardOperation not possible during mapping. Check that circuit contains only StandardOperations");
			}
			if (op->getType() == qc::SWAP) {
				short q0 = qubits.at(op->getTargets().at(0));
				short q1 = qubits.at(op->getTargets().at(1));
				qubits.at(op->getTargets().at(0)) = q1;
				qubits.at(op->getTargets().at(1)) = q0;

				if (q0 != DEFAULT_POSITION) {
					locations.at(q0) = op->getTargets().at(1);
				}
				if (q1 != DEFAULT_POSITION) {
					locations.at(q1) = op->getTargets().at(0);
				}
			}
			if (gatesToAdjust.back() == gateidx) {
				gatesToAdjust.pop_back();
				auto target = op->getTargets().at(0);
				auto targetLocation = locations.at(target);

				if (targetLocation == -1) {
					// qubit only occurs in single qubit gates, can be mapped to an arbitrary free qubit
					unsigned short loc = 0;
					while(qubits.at(loc) != DEFAULT_POSITION) {
						++loc;
					}
					locations.at(target) = loc;
					op->setTargets({loc});
					qcMapped.initialLayout.at(target) = loc;
					qc.outputPermutation.at(target) = loc;
				} else {
					op->setTargets({ static_cast<unsigned short>(targetLocation)});
				}
			}
		}
	}

	for (unsigned short i=0; i<qubits.size(); ++i) {
		if (qubits[i] == -1) {
			qcMapped.outputPermutation.erase(i);
		}
	}

	const auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> diff = end - start;
	results.time = diff.count();
	results.timeout = false;
}

void HeuristicMapper::initResults() {
	Mapper::initResults();
	results.method = Method::Heuristic;

	results.seed = settings.teleportation_seed;
	results.input_teleportation_qubits = settings.teleportation_qubits;
	results.output_teleportation_qubits = settings.teleportation_qubits;
	results.output_teleportation_fake = settings.teleportation_fake;
}

void HeuristicMapper::createInitialMapping() {
	if (layers.empty())
		return;

	if (settings.teleportation_qubits > 0) {
        std::mt19937_64 mt;
        if (settings.teleportation_seed == 0) {
            std::array<std::mt19937_64::result_type, std::mt19937_64::state_size> random_data{};
            std::random_device rd;
            std::generate(std::begin(random_data), std::end(random_data), [&rd]() { return rd(); });
            std::seed_seq seeds(std::begin(random_data), std::end(random_data));
            mt.seed(seeds);
        } else {
            mt.seed(settings.teleportation_seed);
        }

        std::uniform_int_distribution<> dis(0, architecture.getNqubits() - 1);

        for (int i = 0; i < settings.teleportation_qubits; i += 2) {
            Edge e;
            do {
                auto it = std::begin(architecture.getCouplingMap());
                std::advance(it, dis(mt));
                e = *it;
            } while (qubits[e.first] != -1 || qubits[e.second] != -1);
            locations[qc.getNqubits() + i] = e.first;
            locations[qc.getNqubits() + i + 1] = e.second;
            qubits[e.first] = qc.getNqubits() + i;
            qubits[e.second] = qc.getNqubits() + i + 1;
        }

        if (settings.teleportation_fake) {
            settings.teleportation_qubits = 0;
            results.output_teleportation_qubits = settings.teleportation_qubits;
        }
    }

	switch (settings.initialLayoutStrategy) {
		case InitialLayoutStrategy::Identity:
			for (unsigned short i = 0; i < architecture.getNqubits(); ++i) {
				if (qc.initialLayout.count(i)) {
					locations.at(i) = i;
					qubits.at(i) = i;
				}
			}
			break;
		case InitialLayoutStrategy::Static:
			for (const auto& gate: layers.at(0)) {
				if (gate.singleQubit())
					continue;

				for (const auto& edge: architecture.getCouplingMap()) {
					if (qubits.at(edge.first) == DEFAULT_POSITION && qubits.at(edge.second) == DEFAULT_POSITION) {
						qubits.at(edge.first) = gate.control;
						qubits.at(edge.second) = gate.target;
						locations.at(gate.control) = edge.first;
						locations.at(gate.target) = edge.second;
						qcMapped.initialLayout.at(edge.first) = gate.control;
						qcMapped.initialLayout.at(edge.second) = gate.target;
						qcMapped.outputPermutation.at(edge.first) = gate.control;
						qcMapped.outputPermutation.at(edge.second) = gate.target;
						break;
					}
				}
			}

			// assign remaining logical qubits
			for (unsigned short i = 0; i < architecture.getNqubits(); ++i) {
				if (qc.initialLayout.count(i)) {
					if (locations.at(i) == DEFAULT_POSITION) {
						for (unsigned short j = 0; j < architecture.getNqubits(); ++j) {
							if (qubits.at(j) == DEFAULT_POSITION) {
								locations.at(i) = j;
								qubits.at(j) = i;
								qcMapped.initialLayout.at(j) = i;
								qcMapped.outputPermutation.at(j) = i;
								break;
							}
						}
					}
				}
			}

			break;
		case InitialLayoutStrategy::Dynamic:
		case InitialLayoutStrategy::None:
			// nothing to be done here
			break;

		// TODO: Design strategy that maps most used qubit to most connected qubits on architecture
	}
}

void HeuristicMapper::mapUnmappedGates(long layer, HeuristicMapper::Node& node, std::vector<unsigned short>& consideredQubits) {
	for (const auto& gate: layers.at(layer)) {
		if (gate.singleQubit())
			continue;

		consideredQubits.push_back(gate.control);
		consideredQubits.push_back(gate.target);

		auto controlLocation = locations.at(gate.control);
		auto targetLocation = locations.at(gate.target);

		if (controlLocation == DEFAULT_POSITION && targetLocation == DEFAULT_POSITION) {
			std::set<Edge> possibleEdges{};
			for (const auto& edge: architecture.getCouplingMap()) {
				if (qubits.at(edge.first) == DEFAULT_POSITION && qubits.at(edge.second) == DEFAULT_POSITION) {
					possibleEdges.emplace(edge);
				}
			}

			if (possibleEdges.empty()) {
                double bestScore = std::numeric_limits<int>::max();
                std::pair<unsigned short, unsigned short> bestMapping;
                for(int i=0; i < architecture.getNqubits(); i++) {
                    for (int j = i + 1; j < architecture.getNqubits(); j++) {
                        if (qubits.at(i) == DEFAULT_POSITION && qubits.at(j) == DEFAULT_POSITION) {
                            double dist = architecture.distance(i, j);
                            if (dist < bestScore) {
                                bestScore = dist;
                                bestMapping = std::make_pair(i, j);
                            }
                        }
                    }
                }
                auto& chosenEdge = bestMapping;
                locations.at(gate.control) = chosenEdge.first;
                locations.at(gate.target) = chosenEdge.second;
                qubits.at(chosenEdge.first) = gate.control;
                qubits.at(chosenEdge.second) = gate.target;
                qc::QuantumComputation::findAndSWAP(gate.control, chosenEdge.first, qcMapped.initialLayout);
                qc::QuantumComputation::findAndSWAP(gate.target, chosenEdge.second, qcMapped.initialLayout);
                qc::QuantumComputation::findAndSWAP(gate.control, chosenEdge.first, qcMapped.outputPermutation);
                qc::QuantumComputation::findAndSWAP(gate.target, chosenEdge.second, qcMapped.outputPermutation);
			} else {
                // TODO: Consider fidelity here if available. The best available edge should be chosen
                auto& chosenEdge = *possibleEdges.begin();
                locations.at(gate.control) = chosenEdge.first;
                locations.at(gate.target) = chosenEdge.second;
                qubits.at(chosenEdge.first) = gate.control;
                qubits.at(chosenEdge.second) = gate.target;
                qc::QuantumComputation::findAndSWAP(gate.control, chosenEdge.first, qcMapped.initialLayout);
                qc::QuantumComputation::findAndSWAP(gate.target, chosenEdge.second, qcMapped.initialLayout);
                qc::QuantumComputation::findAndSWAP(gate.control, chosenEdge.first, qcMapped.outputPermutation);
                qc::QuantumComputation::findAndSWAP(gate.target, chosenEdge.second, qcMapped.outputPermutation);
			}
		} else if (controlLocation == DEFAULT_POSITION) {
			mapToMinDistance(gate.target, gate.control);
		} else if (targetLocation == DEFAULT_POSITION) {
			mapToMinDistance(gate.control, gate.target);
		}
		node.costHeur = std::max(node.costHeur, distanceOnArchitectureOfLogicalQubits(gate.control, gate.target));
	}

	node.done = node.costHeur <= COST_DIRECTION_REVERSE;
}

void HeuristicMapper::mapToMinDistance(unsigned short source, unsigned short target) {
	auto min = std::numeric_limits<double>::max();
	unsigned short pos = DEFAULT_POSITION;
	for (int i = 0; i < architecture.getNqubits(); ++i) {
		if (qubits.at(i) == DEFAULT_POSITION) {
			// TODO: Consider fidelity here if available
			auto distance = distanceOnArchitectureOfPhysicalQubits(locations.at(source), i);
			if (distance < min) {
				min = distance;
				pos = i;
			}
		}
	}
	qubits.at(pos) = target;
	locations.at(target) = pos;
	qc::QuantumComputation::findAndSWAP(target, pos, qcMapped.initialLayout);
	qc::QuantumComputation::findAndSWAP(target, pos, qcMapped.outputPermutation);
}

HeuristicMapper::Node HeuristicMapper::AstarMap(long layer) {
	std::vector<unsigned short> consideredQubits{};
	Node node{};

	mapUnmappedGates(layer, node, consideredQubits);

	node.locations = locations;
	node.qubits = qubits;

	nodes.push(node);

	while (!nodes.top().done) {
		Node current = nodes.top();
		nodes.pop();
		expandNode(consideredQubits, current, layer);
	}

	Node result = nodes.top();
	nodes.pop();

	while (!nodes.empty()) {
		nodes.pop();
	}

	return result;

}

void HeuristicMapper::expandNode(const std::vector<unsigned short>& consideredQubits, Node& node, long layer) {
	std::vector<std::vector<bool>> used_swaps;
	used_swaps.reserve(architecture.getNqubits());
	for (int p=0; p<architecture.getNqubits(); ++p) {
		used_swaps.emplace_back(architecture.getNqubits());
	}

    std::set<Edge> perms = architecture.getCouplingMap();
    architecture.getCurrentTeleportations().clear();
    architecture.teleportation_qubits.clear();
    for(int i = 0; i < settings.teleportation_qubits; i+=2) {
        architecture.teleportation_qubits.emplace_back(node.locations[qc.getNqubits() + i], node.locations[qc.getNqubits() + i + 1]);
        Edge e;
        for(auto const& g : architecture.getCouplingMap()) {
            if(g.first == node.locations[qc.getNqubits() + i] && g.second != node.locations[qc.getNqubits() + i + 1]) {
                e.first = g.second;
                e.second = node.locations[qc.getNqubits() + i + 1];
                architecture.getCurrentTeleportations().insert(e);
                perms.insert(e);
            }
            if(g.second == node.locations[qc.getNqubits() + i] && g.first != node.locations[qc.getNqubits() + i + 1]) {
                e.first = g.first;
                e.second = node.locations[qc.getNqubits() + i + 1];
                architecture.getCurrentTeleportations().insert(e);
                perms.insert(e);
            }
            if(g.first == node.locations[qc.getNqubits() + i + 1] && g.second != node.locations[qc.getNqubits() + i]) {
                e.first = g.second;
                e.second = node.locations[qc.getNqubits() + i];
                architecture.getCurrentTeleportations().insert(e);
                perms.insert(e);
            }
            if(g.second == node.locations[qc.getNqubits() + i + 1] && g.first != node.locations[qc.getNqubits() + i]) {
                e.first = g.first;
                e.second = node.locations[qc.getNqubits() + i];
                architecture.getCurrentTeleportations().insert(e);
                perms.insert(e);
            }
        }
    }

	for (const auto& q: consideredQubits) {
		for (const auto& edge: perms) {
			if (edge.first == node.locations.at(q) || edge.second == node.locations.at(q)) {
				auto q1 = node.qubits.at(edge.first);
				auto q2 = node.qubits.at(edge.second);
				if  (q2 == -1 || q1 == -1) {
					expand_node_add_one_swap(edge, node, layer);
				} else if (!used_swaps.at(q1).at(q2)){
					used_swaps.at(q1).at(q2) = true;
					used_swaps.at(q2).at(q1) = true;
					expand_node_add_one_swap(edge, node, layer);
				}
			}
		}
	}
}

void HeuristicMapper::expand_node_add_one_swap(const Edge &swap, Node& node, long layer) {
	auto& currentLayer = layers.at(layer);

	Node new_node = Node(node.qubits, node.locations, node.swaps);
	new_node.nswaps++;

	new_node.swaps.emplace_back();
	if(architecture.getCouplingMap().find(swap) != architecture.getCouplingMap().end() ||
	    architecture.getCouplingMap().find(Edge {swap.second, swap.first}) != architecture.getCouplingMap().end()) {

	    if (architecture.bidirectional()) {
            new_node.costFixed = node.costFixed + COST_BIDIRECTIONAL_SWAP;
        } else {
            new_node.costFixed = node.costFixed + COST_UNIDIRECTIONAL_SWAP;
        }

        new_node.applySWAP(swap, architecture);
    } else {
        new_node.costFixed = node.costFixed + COST_TELEPORTATION;
	    new_node.applyTeleportation(swap, architecture);
	}
    new_node.costTotal = new_node.costFixed;
    new_node.done = true;

	for (const auto& gate: currentLayer) {
		if (gate.singleQubit())
			continue;

		new_node.updateHeuristicCost(architecture, gate, settings.admissibleHeuristic);
		new_node.checkUnfinished(architecture, gate);
	}

	// calculate heuristics for the cost of the following layers
	if (settings.lookahead) {
		lookahead(getNextLayer(layer), new_node);
	}

	nodes.push(new_node);
}

void HeuristicMapper::lookahead(long layer, HeuristicMapper::Node& node) {
	auto nextLayer = layer;
	double factor = settings.firstLookaheadFactor;

	for (int i = 0; i < settings.nrLookaheads; ++i) {
		if (nextLayer == -1) {
			break;
		}

		double penalty = 0.;
		for (const auto& gate: layers.at(nextLayer)) {
			if (gate.singleQubit())
				continue;

			auto loc1 = node.locations.at(gate.control);
			auto loc2 = node.locations.at(gate.target);
			if (loc1 == DEFAULT_POSITION && loc2 == DEFAULT_POSITION) {
				// no penalty
			} else if (loc1 == DEFAULT_POSITION) {
				auto min = std::numeric_limits<double>::max();
				for (int j = 0; j < architecture.getNqubits(); ++j) {
					if (node.qubits.at(j) == DEFAULT_POSITION) {
						// TODO: Consider fidelity here if available
						min = std::min(min, distanceOnArchitectureOfPhysicalQubits(j, node.locations.at(gate.target)));
					}
				}
				penalty = heuristicCost(penalty, min);
			} else if (loc2 == DEFAULT_POSITION) {
				auto min = std::numeric_limits<double>::max();
				for (int j = 0; j < architecture.getNqubits(); ++j) {
					if (node.qubits.at(j) == DEFAULT_POSITION) {
						// TODO: Consider fidelity here if available
						min = std::min(min, distanceOnArchitectureOfPhysicalQubits(node.locations.at(gate.control), j));
					}
				}
				penalty = heuristicCost(penalty, min);
			} else {
				auto cost = architecture.distance(node.locations.at(gate.control), node.locations.at(gate.target));
				penalty = heuristicCost(penalty, cost);
			}
		}

		node.lookaheadPenalty += factor * penalty;
		factor *= settings.lookaheadFactor;
		nextLayer = getNextLayer(nextLayer); // TODO: consider single qubits here for better fidelity lookahead
	}

}
