#include "mapper.hpp"

#include <math.h>

/**
 * returns the maximal depth based on depths
 */
int get_maximal_depth(int const * const depths) {	
	// calculate max 
    int max_depth = 0;
	for(int i = 0; i < arch.positions; i++) {
		max_depth = std::max(max_depth, depths[i]);
	}
    return max_depth;
}

/**
 * calculates the total workload cost based on the different workload of gates of the qubits
 * currently the standard deviation of all workloads is used
 */
long long workload_cost(int const * const workload) {
	int       avg              = 0;
	int       nqubits_not_null = 0;
	long long variance         = 0;
	for(int i = 0; i < arch.positions; i++) {  // calcualte average
		if(workload[i] != 0) {
			avg              += workload[i];
			nqubits_not_null ++;
		}
	}
	avg /= nqubits_not_null;
	for(int i = 0; i < arch.positions; i++) {
		if(workload[i] != 0) {
			long long diff = workload[i] - avg;
			//variance += abs(diff); 
			variance += diff * diff;
		}
	}
	
	return sqrt(variance/nqubits_not_null);
}

/**
 * calculates the total fidelity cost based on the different fidelities of gates of the qubits
 */
double fidelity_cost(double const * const fidelities) {
	/*
	double min_fidelity = 1;
	for(int i = 0; i < arch.positions; i++) {
		min_fidelity = std::min(min_fidelity, fidelities[i]);
	}
    return 1 / min_fidelity;

	*/

	/*
	double  avg             = 0;
	double  nqubits_not_one = 0;
	double  variance        = 0;
	for(int i = 0; i < arch.positions; i++) {  // calcualte average
		if(fidelities[i] != 1) {
			avg += fidelities[i];
			nqubits_not_one ++;
		}
	}
	avg /= nqubits_not_one;
	for(int i = 0; i < arch.positions; i++) {
		if(fidelities[i] != 1) {
			double diff = fidelities[i] - avg; 
			variance += diff * diff;
		}
	}
	
	return sqrt(variance/(nqubits_not_one + 0.000000001));
	*/
	
	double  nqubits_not_one = 0;
	double  variance        = 0;
	
	for(int i = 0; i < arch.positions; i++) {
		if(fidelities[i] != 1) {
			double diff = 1 - fidelities[i]; 
			variance   += diff * diff;
			nqubits_not_one ++;
		}
	}
	
	return sqrt(variance/(nqubits_not_one + 0.000000001));
}

/**
 * calculates the heuristic cost for a certain dijkstra node based on the path length
 */
double calculate_heuristic_cost(dijkstra_node const * const node) {
	int path_length = node->cost - 1;
	if(node->contains_correct_edge) {
#if SPECIAL_OPT
		return path_length;
#else 
		return path_length * COST_SWAP;
#endif
	}
#if SPECIAL_OPT
	return path_length + INVERSE; 
#else
	return path_length * COST_SWAP + 4;
#endif
}

/**
 * calculates the total cost of a node based on the cost
 * if special opt is enabled also the workload and the depth is considered according 
 * to their factors
 */
double get_total_cost(const node& n) {
#if SPECIAL_OPT
	return (fidelity_cost(n.fidelities)                       * FIDELITY_NORM)    +
		   (workload_cost(n.workload)                         * WORKLOAD_NORM)    + 
		   ((get_maximal_depth(n.depths) - current_depth) /((double) DEPTH_SWAP) * DEPTH_PERCENTAGE) + 
		   (n.cost_fixed/((double)  COST_SWAP)                * COST_PERCENTAGE);
#else
    return n.cost_fixed;
#endif
}

/**
 * calculates the total lookahead cost of a node based on cost of the special opt 
 */
double get_total_lookahead_cost(int const * const depths, int const * const workload, double const * const fidelities) {
	return ((get_maximal_depth(depths) - current_depth) / ((double) DEPTH_SWAP) * DEPTH_PERCENTAGE) + 
	  	   (workload_cost(workload)                           * WORKLOAD_NORM)    + 
		   (fidelity_cost(fidelities)                         * FIDELITY_NORM);
}

/**
 * combines the heuristic values of the old value and the new value
 */
double heuristic_function(const double old_heur, const double new_heur) {
#if HEURISTIC_ADMISSIBLE
	return std::max(old_heur, new_heur);
#else
	return old_heur + new_heur;
#endif
}

/**
 * returns the heuristic cost for a certain node and considering the current cost
 */
double get_heuristic_cost(const double cost_heur, const node& n, const QASMparser::gate& g) {
	return heuristic_function(cost_heur, arch.dist[n.locations[g.control]][n.locations[g.target]]);
}