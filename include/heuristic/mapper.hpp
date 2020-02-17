#ifndef MAPPER_H_
#define MAPPER_H_


#include <vector>


/**
 * Own Includes
 */
#include "parser.hpp"
#include "unique_priority_queue.hpp"


/**
 * General Defines
 */
#define SUCCESS 0
#define ERROR   1

#define ARCH_LINEAR_N      0
#define ARCH_IBM_QX5       1
#define ARCH_IBM_MELBOURNE 2

#define UNUSED_FUNCTION __attribute__ ((unused))

#define ARCH ARCH_IBM_MELBOURNE

/**
 * Control Defines
 */
#define VERIFICATION         0 // maps all logical qubits to physical qubits with the same index and adds a swap layer at the the end
							   // this swap layer ignores constraints

#define LOOK_AHEAD           1 // enables the lookahead; is additionally controlled by the constants below
#define USE_INITIAL_MAPPING  1 // enables initial mapping, it is automatically enabled when using special_opt
#define HEURISTIC_ADMISSIBLE 1 // enables the admissible heuristic approach
#define ONE_SWAP_PER_EXPAND  1 // decides whether whole permutations or only one swap should be considered for a expansion step
#define SPECIAL_OPT          1 // enables special optimizations like depth and workload; is additionally controlled by the constants below

#if SPECIAL_OPT || VERIFCATION // force initial mapping -> because of keeping track of unmapped gates
#undef  USE_INITIAL_MAPPING
#define USE_INITIAL_MAPPING  1
#endif

#ifndef ARCH
#define ARCH ARCH_LINEAR_N     // assume default architecture
#endif

/*
 * Constants
 */
// cost
const int COST_GATE     = 1;
const int COST_SWAP     = 7 * COST_GATE;

// fidelity
const int FIDELITY_GATE = 1;
const int FIDELITY_CNOT = 5;
const int FIDELITY_SWAP = 2 * FIDELITY_GATE + 3 * FIDELITY_CNOT;

// depth
const int DEPTH_GATE    = 1;
const int DEPTH_SWAP    = 5 * DEPTH_GATE;

// workload
const int WORKLOAD_GATE = 1;
const int WORKLOAD_CNOT = 5;
const int WORKLOAD_SWAP = 2 * WORKLOAD_GATE + 3 * WORKLOAD_CNOT;

// special opt. factors
const double COST_PERCENTAGE  = 1;
const double DEPTH_PERCENTAGE = 1 - COST_PERCENTAGE;
const double WORKLOAD_FACTOR  = 0;
const double WORKLOAD_NORM    = WORKLOAD_FACTOR / 1000;
const double FIDELITY_FACTOR  = 0; // good values: 10, 100
const double FIDELITY_NORM    = FIDELITY_FACTOR / 1;
const double INVERSE          = DEPTH_PERCENTAGE * (((double)2) * DEPTH_GATE / DEPTH_SWAP) + COST_PERCENTAGE  * 0.57; // additional cost if no edge is in the correct direction

// lookahead
const int    N_LOOK_AHEADS             = 15;
const double FIRST_LOOK_AHEAD_FACTOR   = 0.75;
const double GENERAL_LOOK_AHEAD_FACTOR = 0.5;

const bool   SPECIAL_OPT_VALUES_SET    = (DEPTH_PERCENTAGE != 0 || WORKLOAD_NORM != 0 || FIDELITY_NORM != 0);

/** 
 * Global variables - standard
 */
extern unsigned long ngates;
extern unsigned long current_depth;
extern unsigned int  nqubits;

/**
 * Types 
 */
struct edge {
	int    v1;
	int    v2;
	double fidelity;
 
	edge() = default;
	edge(int v1_, int v2_, double fidelity_ = 1) : v1(v1_), v2(v2_), fidelity(fidelity_) {}
};

inline bool operator<(const edge& lhs, const edge& rhs) {
	if (lhs.v1 != rhs.v1) {
		return lhs.v1 < rhs.v1;
	}
	return lhs.v2 < rhs.v2;
}

#if ONE_SWAP_PER_EXPAND
typedef edge                   SWAP_TYPE;
#else
typedef std::vector<edge>      SWAP_TYPE;
#endif

typedef std::vector<SWAP_TYPE> SWAP_LIST_TYPE;   


struct architecture {
	int            positions;
	double**       dist;
	double**       fidelity_dist;
	double*        initial_fidelities;
	double*        singlequbit_fidelities;
	std::set<edge> graph;
};
extern architecture arch;


struct node {
	int    cost_fixed;
	double cost_heur;
	double lookahead_penalty;
	double total_cost;
	int*   qubits;    // get qubit of location -> -1 indicates that there is "no" qubit at a certain location
	int*   locations; // get location of qubits -> -1 indicates that a qubit does not have a location -> shall only occur for i > nqubits
#if SPECIAL_OPT
	int*    depths;
	int*    workload;
	double* fidelities;
#endif
	int    nswaps;
	int    done;
	SWAP_LIST_TYPE swaps;
};

struct node_func_less {
	// true iff x < y
	bool operator()(const node& x, const node& y) const {
		for(int i=0; i < arch.positions; i++) {
			if (x.qubits[i] != y.qubits[i]) {
				return x.qubits[i] < y.qubits[i];
			}
		}
		return false;
	}
};

struct node_cost_greater {
	// true iff x > y
	bool operator()(const node& x, const node& y) const {
		if ((x.total_cost + x.cost_heur + x.lookahead_penalty) != (y.total_cost + y.cost_heur + y.lookahead_penalty)) {
			return (x.total_cost + x.cost_heur + x.lookahead_penalty) > (y.total_cost + y.cost_heur + y.lookahead_penalty);
		}

		if(x.done == 1) {
			return false;
		}
		if(y.done == 1) {
			return true;
		}

		if (x.cost_heur + x.lookahead_penalty != y.cost_heur + y.lookahead_penalty) {
			return x.cost_heur + x.lookahead_penalty > y.cost_heur + y.lookahead_penalty;
		} else {
			return node_func_less{}(x, y);
		}

	}
};

struct cleanup_node {
	void operator()(const node& n) {
		delete[] n.qubits;
		delete[] n.locations;
#if SPECIAL_OPT
		delete[] n.depths;
		delete[] n.workload;
		delete[] n.fidelities;
#endif
	}
};

// circuit properties
struct circuit_properties {
	int* locations;
	int* qubits;
#if SPECIAL_OPT
	int* depths;
	int* workload;
	double* fidelities;
#endif
};

// dijkstra
struct dijkstra_node {
	bool   contains_correct_edge;
	bool   visited;
	int    pos;
    double cost;
};

/** 
 * Global variables - special types
 */
extern std::vector<std::vector<QASMparser::gate>>                                   layers;
extern unique_priority_queue<node, cleanup_node, node_cost_greater, node_func_less> nodes;


/**
 * Functions
 */
// architecture_handling
bool create_architecture_properties(const std::string input);
void delete_architecture_properties();

// cost
int       get_maximal_depth(int const * const depths);
long long workload_cost(int const * const workload);
double    fidelity_cost(double const * const fidelities);
double    calculate_heuristic_cost(dijkstra_node const * const node);
double    get_total_lookahead_cost(int const * const depths, int const * const workload, double const * const fidelities);
double    heuristic_function(const double old_heur, const double new_heur);
double    get_total_cost(const node& n);
double    heuristic_function(const double old_heur, const double new_heur);
double    get_heuristic_cost(const double cost_heur, const node& n, const QASMparser::gate& g);

// node_handling
node create_node();
node create_node(const node& base, edge const * const new_swaps, const int nswaps = 1);
void update_node(node& n, const circuit_properties& p);
void check_if_not_done(node& n, const int value);
void delete_node(const node& n);

// layer_handling
std::vector<std::vector<QASMparser::gate>> init_layers(const std::vector<QASMparser::gate> &gates);
unsigned int                     get_next_layer(const unsigned int layer);
unsigned int                     calculate_max_layer_width();

// circuit_property_handling
circuit_properties create_circuit_properties();
void               delete_circuit_properties(circuit_properties& p);
void               adapt_circuit_properties(circuit_properties& p, const node& n);
void 			   update_properties(circuit_properties& p, const int layer);

// mapping.cpp
void mapping(const std::vector<QASMparser::gate>& gates, std::vector<std::vector<QASMparser::gate>>& mapped_circuit, 
			std::vector<QASMparser::gate>& all_gates, int &total_swaps, circuit_properties& properties);

// util.cpp
void initial_mapping(circuit_properties& properties);
void map_to_min_distance(int* map, int* loc, const int source, const int target);
void map_unmapped_gates(const int layer, circuit_properties& p, node& n, std::vector<int>& considered_qubits);
void fix_positions_of_single_qubit_gates(int* locations, int* qubits, std::vector<QASMparser::gate>& all_gates);
void generate_circuit(std::vector<std::vector<QASMparser::gate>>& mapped_circuit, const std::vector<QASMparser::gate>& all_gates);
void map_to_inital_permutation(std::vector<QASMparser::gate>& all_gates, circuit_properties& properties); // add swaps so that each logical qubit is mapped to the pysical qubit with the same index


#endif /* MAPPER_H_ */
