#include "mapper.hpp"

/**
 * initializes the circuit properties
 * 	- locations
 *  - qubits
 * if special optimization is used -> also 
 *  - depths
 *  - workload
 *  - fidelities
 */
circuit_properties create_circuit_properties() {
    circuit_properties p;
    p.locations  = new int[nqubits];
	p.qubits     = new int[arch.positions];
#if SPECIAL_OPT
	p.depths     = new int[arch.positions]();
	p.workload   = new int[arch.positions](); 
	p.fidelities = new double[arch.positions];
#endif

	//Initially, no physical qubit is occupied
	for (int i = 0; i < arch.positions; i++) {
		p.qubits[i] = -1;
#if SPECIAL_OPT
		p.fidelities[i] = arch.initial_fidelities[i];
#endif
	}

	//Initially, no logical qubit is mapped to a physical one
	for(unsigned i = 0; i < nqubits; i++) {
		p.locations[i] = -1;
	}

    return p;
}

void adapt_circuit_properties(circuit_properties& p, const node& n) {
	delete_circuit_properties(p);
	p.locations  = n.locations;
	p.qubits     = n.qubits;
#if SPECIAL_OPT
	p.depths     = n.depths;
	p.workload   = n.workload;
	p.fidelities = n.fidelities;
#endif
}

/**
 * adapts the properties of the current qubits by considering all gates of the specified layer
 */
void update_properties(circuit_properties& p, const int layer) {
#if SPECIAL_OPT
	for(std::vector<QASMparser::gate>::iterator it = layers[layer].begin(); it != layers[layer].end(); it++) {
	    QASMparser::gate g = *it;
		int pt = p.locations[g.target];
		if(g.control != -1) {
			int pc        = p.locations[g.control];
			int max_depth = std::max(p.depths[pc], p.depths[pt]) + DEPTH_GATE;
            p.depths[pc]  = max_depth;
			p.depths[pt]  = max_depth;
            p.workload[pt]   += WORKLOAD_CNOT;
            p.workload[pc]   += WORKLOAD_CNOT;
            p.fidelities[pt] *= arch.fidelity_dist[pc][pt];
            p.fidelities[pc] *= arch.fidelity_dist[pc][pt];

			edge e(pc, pt);
			if (arch.graph.find(e) == arch.graph.end()) {
				p.depths[pt]     += DEPTH_GATE    << 1; 
				p.depths[pc]     += DEPTH_GATE    << 1;
				p.workload[pt]   += WORKLOAD_GATE << 1;
				p.workload[pc]   += WORKLOAD_GATE << 1; 
				p.fidelities[pt] *= arch.singlequbit_fidelities[pt];
				p.fidelities[pc] *= arch.singlequbit_fidelities[pc];
			}
#if USE_INITIAL_MAPPING
		} else {
#else
		} else if(pt >= 0) {
#endif
			p.depths[pt]     += DEPTH_GATE;
        	p.workload[pt]   += WORKLOAD_GATE;
        	p.fidelities[pt] *= arch.singlequbit_fidelities[pt];
		}
	}
#endif
}

/**
 * clean up the circuit properties
 */
void delete_circuit_properties(circuit_properties& p) {
    delete[] p.locations;
	delete[] p.qubits;	
#if SPECIAL_OPT
	delete[] p.depths;	
	delete[] p.workload;
	delete[] p.fidelities;
#endif
}
