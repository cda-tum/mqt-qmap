#include "mapper.hpp"

#include <string.h>

/**
 * maps all qubits to locations
 */
void initial_mapping(circuit_properties& properties) {
	int* qubits    = properties.qubits;
	int* locations = properties.locations;
#if VERIFICATION
	for(unsigned int i = 0; i < nqubits; i++) {
		locations[i] = i;
		qubits[i]    = i;
	}
#else
	for (std::vector<QASMparser::gate>::iterator it = layers[0].begin(); it != layers[0].end(); it++) {
		QASMparser::gate g = *it;
		if (g.control != -1) {
			for(std::set<edge>::iterator it = arch.graph.begin(); it != arch.graph.end(); it++) {
				if(qubits[it->v1] == -1 && qubits[it->v2] == -1) {
					qubits[it->v1]       = g.control;
					qubits[it->v2]       = g.target;
					locations[g.control] = it->v1;
					locations[g.target]  = it->v2;
					break;
				}
			}
		}
	}
	for(unsigned int i = 0; i < nqubits; i++) {
		if(locations[i] == -1) {
			int j = 0;
			while(qubits[j] != -1){
				j++;
			}
			locations[i] = j;
			qubits[j]    = i;
		}
	}
#endif
}

/**
 * maps the qubit to the physical location that is not yet mapped and has minimum distance
 */
void map_to_min_distance(int* map, int* loc, const int source, const int target) {
	int min     = 1000;
	int min_pos = -1;
	for(int i = 0; i < arch.positions; i++) {
		if(map[i] == -1 && arch.dist[loc[source]][i] < min) {
			min     = arch.dist[loc[source]][i];
			min_pos = i;
		}
	}
	map[min_pos] = target;
	loc[target]  = min_pos;
}

/**
 * maps the unmapped gates
 */
void map_unmapped_gates(const int layer, circuit_properties& p, node& n, std::vector<int>& considered_qubits) {
	int* map = p.qubits;
	int* loc = p.locations;

	//Find a mapping for all logical qubits in the CNOTs of the layer that are not yet mapped
	for (std::vector<QASMparser::gate>::iterator it = layers[layer].begin(); it != layers[layer].end(); it++) {
		QASMparser::gate g = *it;
		if (g.control != -1) {
			considered_qubits.push_back(g.control);
			considered_qubits.push_back(g.target);

			if(loc[g.control] == -1 && loc[g.target] == -1) {
                std::set<edge> possible_edges;
				for(std::set<edge>::iterator it = arch.graph.begin(); it != arch.graph.end(); it++) {
					if(map[it->v1] == -1 && map[it->v2] == -1) {
						possible_edges.insert(*it);
					}
				}
				if(!possible_edges.empty()) {
					edge e = *possible_edges.begin();
					loc[g.control] = e.v1;
					map[e.v1] = g.control;
					loc[g.target] = e.v2;
					map[e.v2] = g.target;
				} else {
                    std::cerr << "no edge available!";
                    exit(1);
				}
			} else if(loc[g.control] == -1) {
				map_to_min_distance(map, loc, g.target, g.control);
			} else if(loc[g.target] == -1) {
				map_to_min_distance(map, loc, g.control, g.target);
			}
			n.cost_heur = std::max(n.cost_heur, arch.dist[loc[g.control]][loc[g.target]]);
		}
	}
}

/**
 * fix the position of the single qubit gates
 */
void fix_positions_of_single_qubit_gates(int* locations, int* qubits, std::vector<QASMparser::gate>& all_gates) {
	for(std::vector<QASMparser::gate>::reverse_iterator it = all_gates.rbegin(); it != all_gates.rend(); it++) {
		if(strcmp(it->type, "SWP") == 0) {
			int tmp_qubit1 = qubits[it->control];
			int tmp_qubit2 = qubits[it->target];

			qubits[it->control] = tmp_qubit2;
			qubits[it->target]  = tmp_qubit1;

			if(tmp_qubit1 != -1) {
				locations[tmp_qubit1] = it->target;
			}
			if(tmp_qubit2 != -1) {
				locations[tmp_qubit2] = it->control;
			}
		}
		if(it->target < 0) {
			int target = -(it->target +1);
			it->target = locations[target];
			if(locations[target] == -1) {
				//This qubit occurs only in single qubit gates -> it can be mapped to an arbirary physical qubit
				int loc = 0;
				while(qubits[loc] != -1) {
					loc++;
				}
				locations[target] = loc;
			}
		}
	}
}

/**
 * generate the circuit based on the gates
 */
void generate_circuit(std::vector<std::vector<QASMparser::gate>>& mapped_circuit, const std::vector<QASMparser::gate>& all_gates) {
	int *last_layer = new int[arch.positions];
	for(int i = 0; i < arch.positions; i++) {
		last_layer[i] = -1;
	}

	//build resulting circuit
	for(std::vector<QASMparser::gate>::const_iterator it = all_gates.begin(); it != all_gates.end(); it++) {
		if(strcmp(it->type, "SWP") == 0) {
			continue;
		}
		if(it->control == -1) {
			//single qubit gate
			QASMparser::gate g = *it;
			unsigned int layer = last_layer[g.target] + 1;

			if (mapped_circuit.size() <= layer) {
				mapped_circuit.push_back(std::vector<QASMparser::gate>());
			}
			mapped_circuit[layer].push_back(g);
			last_layer[g.target] = layer;
		} else {
			QASMparser::gate g = *it;
			unsigned int layer = std::max(last_layer[g.control], last_layer[g.target]) + 1;
			if (mapped_circuit.size() <= layer) {
				mapped_circuit.push_back(std::vector<QASMparser::gate>());
			}
			mapped_circuit[layer].push_back(g);

			last_layer[g.target]  = layer;
			last_layer[g.control] = layer;
		}
	}
	
	delete[] last_layer;
}

void map_to_inital_permutation(std::vector<QASMparser::gate>& all_gates, circuit_properties& properties) { // add swaps so that each logical qubit is mapped to the pysical qubit with the same index
	int locations[nqubits];
	int qubits[arch.positions];
	
	memcpy(locations, properties.locations, sizeof(int) * nqubits);
	memcpy(qubits,    properties.qubits,    sizeof(int) * arch.positions);
	
	std::cout << std::endl;
	for(int i = 0; i < (int)nqubits; i++) {
		if(locations[i] != i) {
			QASMparser::gate cnot;
			QASMparser::gate h1;
			QASMparser::gate h2;

			cnot.control = locations[qubits[i]];
			cnot.target  = locations[i];

			strcpy(cnot.type, "CX");
			strcpy(h1.type,   "U3(pi/2,0,pi)");
			strcpy(h2.type,   "U3(pi/2,0,pi)");

			h1.control = h2.control = -1;
			h1.target  = i;
			h2.target  = locations[i];

			all_gates.push_back(cnot);
			all_gates.push_back(h1);
			all_gates.push_back(h2);
			all_gates.push_back(cnot);
			all_gates.push_back(h1);
			all_gates.push_back(h2);
			all_gates.push_back(cnot);	

			int tmp_qubit1 = qubits[cnot.control];
			int tmp_qubit2 = qubits[cnot.target];

			qubits[cnot.target]  = tmp_qubit1;
			qubits[cnot.control] = tmp_qubit2;

			locations[tmp_qubit1] = cnot.target;
			locations[tmp_qubit2] = cnot.control;
		}
	}
	for(int i = 0; i < (int)nqubits; i++) {
		//std::cout << "  q" << i            << " mapped to x" << qubits[i]            << std::endl;
		//std::cout << "  q" << i            << " mapped to Q" << locations[i]         << std::endl;
		//std::cout << "  Q" << locations[i] << " mapped to q" << qubits[locations[i]] << std::endl;		
		assert(i            == locations[i]);
		assert(locations[i] == qubits[locations[i]]);
	}
	/*
	for(unsigned int i = 0; i < nqubits; i++) {
		locations[i] = i;
		qubits[i]    = i;
	}
	for(int i = 0; i < 5; i++) {
		//std::cout << "  q" << i            << " mapped to x" << qubits[i]            << std::endl;
		std::cout << "  q" << i            << " mapped to Q" << locations[i]         << std::endl;
		std::cout << "  Q" << locations[i] << " mapped to q" << qubits[locations[i]] << std::endl;		
		//assert(i            == locations[i]);
		//assert(locations[i] == qubits[locations[i]]);
	}
	for(std::vector<QASMparser::gate>::const_iterator it = all_gates.begin(); it != all_gates.end(); it++) {
		if(strcmp(it->type, "SWP") == 0) {
			std::cout << "  q" << qubits[it->control] << "[" << it->control << "] <-> q" <<  qubits[it->target] << "[" << it->target << "]" << std::endl;
			int tmp_qubit1 = qubits[it->control];
			int tmp_qubit2 = qubits[it->target];

			qubits[it->control] = tmp_qubit2;
			qubits[it->target]  = tmp_qubit1;

			if(tmp_qubit1 != -1) {
				locations[tmp_qubit1] = it->target;
			}
			if(tmp_qubit2 != -1) {
				locations[tmp_qubit2] = it->control;
			}
		
		}
	}
	for(int i = 0; i < 5; i++) {
		//std::cout << "  q" << i            << " mapped to x" << qubits[i]            << std::endl;
		std::cout << "  q" << i            << " mapped to Q" << properties.locations[i]         << std::endl;
		std::cout << "  Q" << properties.locations[i] << " mapped to q" <<properties.qubits[properties.locations[i]] << std::endl;		
		//assert(i            == locations[i]);
		//assert(locations[i] == qubits[locations[i]]);
	}

	
	for(int i = 0; i < 5; i++) {
		//std::cout << "  q" << i            << " mapped to x" << qubits[i]            << std::endl;
		std::cout << "  q" << i            << " mapped to Q" << locations[i]         << std::endl;
		std::cout << "  Q" << locations[i] << " mapped to q" << qubits[locations[i]] << std::endl;		
		//assert(i            == locations[i]);
		//assert(locations[i] == qubits[locations[i]]);
	}
	for(std::vector<QASMparser::gate>::reverse_iterator it = all_gates.rbegin(); it != all_gates.rend(); it++) {
		if(strcmp(it->type, "SWP") == 0) {
			std::cout << "  q" << qubits[it->control] << "[" << it->control << "] <-> q" <<  qubits[it->target] << "[" << it->target << "]" << std::endl;
			int tmp_qubit1 = qubits[it->control];
			int tmp_qubit2 = qubits[it->target];

			qubits[it->control] = tmp_qubit2;
			qubits[it->target]  = tmp_qubit1;

			if(tmp_qubit1 != -1) {
				locations[tmp_qubit1] = it->target;
			}
			if(tmp_qubit2 != -1) {
				locations[tmp_qubit2] = it->control;
			}
		
		}
	}
	
	*/
}