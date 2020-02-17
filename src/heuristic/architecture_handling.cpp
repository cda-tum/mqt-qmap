#include "mapper.hpp"

/**
 * static function declarations 
 */
static void init_arch_arrays(int positions);
template<class DijkstraCmp>
static void set_dijkstra_node(dijkstra_node* nodes, std::priority_queue<dijkstra_node*, std::vector<dijkstra_node*>, DijkstraCmp>& queue,
					          const edge e, const int parent, const int pos, const bool contains_correct_edge);
template<class DijkstraCmp>
static void dijkstra(const std::set<edge>& graph, dijkstra_node* nodes, const int start);
template<class DijkstraCmp>
static void build_dijkstra_table(const std::set<edge>& graph, double*** distance_table);

static bool build_graph_from_file(const std::string input);
static void build_graph_linear(int nqubits)           UNUSED_FUNCTION;
static void build_graph_QX5()                         UNUSED_FUNCTION;
static void build_graph_melbourne()                   UNUSED_FUNCTION;

/**
 * File specific types
 */
class dijkstra_distance_cmp {
	public:
		static bool compare_parameters(double x_cost, double y_cost, bool x_contains_correct_edge, bool y_contains_correct_edge) {
			if(x_cost != y_cost) {
				return x_cost < y_cost;
			}

			return x_contains_correct_edge && !y_contains_correct_edge;
		}

		static void set_initial_cost(dijkstra_node* x) {
			x->cost = 0;
		}
 
		static double get_cost(dijkstra_node const * const x) {
			return calculate_heuristic_cost(x);
		}

		static double step(dijkstra_node const * const x, edge e) {
			return x->cost + 1;
		}

		bool operator()(dijkstra_node const * const x, dijkstra_node const * const y) const {
			return !compare_parameters(x->cost, y->cost, x->contains_correct_edge, y->contains_correct_edge);
		}
};

/*
class dijkstra_fidelity_cmp {
	public:
		static bool compare_parameters(double x_cost, double y_cost, bool x_contains_correct_edge, bool y_contains_correct_edge) {
			if(x_cost != y_cost) {
				return x_cost > y_cost;
			}

			return x_contains_correct_edge && !y_contains_correct_edge;
		}

		static void set_initial_cost(dijkstra_node* x) {
			x->cost = 1;
		}

		static double get_cost(dijkstra_node const * const x) {
			return x->cost;
		}

		static double step(dijkstra_node const * const x, edge e) {
			return x->cost * e.fidelity;
		}

		bool operator()(dijkstra_node const * const x, dijkstra_node const * const y) const {
			return !compare_parameters(x->cost, y->cost, x->contains_correct_edge, y->contains_correct_edge);
		}
};
*/

static void init_arch_arrays(int positions) {	
	arch.positions              = positions;
	arch.initial_fidelities     = new double[arch.positions];
	arch.singlequbit_fidelities = new double[arch.positions];

	for (int i = 0; i < arch.positions; i++) {
		arch.initial_fidelities[i]     = 1;
		arch.singlequbit_fidelities[i] = 1;
	}
}

/**
 * generates a graph from the file with the filename input, 
 * if the name is empty a graph according to ARCH is generated
 */
bool create_architecture_properties(const std::string input) {
    if(input.empty()) {
#if ARCH == ARCH_LINEAR_N
		build_graph_linear(nqubits);
#elif ARCH == ARCH_IBM_QX5
		build_graph_QX5();
#elif ARCH == ARCH_IBM_MELBOURNE
		build_graph_melbourne();
#else
    	static_assert(false, "No architecture specified!");
#endif
    } else {
        if(!build_graph_from_file(input)) {
			return false;
		}
    }


	build_dijkstra_table<dijkstra_distance_cmp>(arch.graph, &arch.dist);
	//build_dijkstra_table<dijkstra_fidelity_cmp>(graph, &arch.fidelity_dist);

	// setup fidelity_table
	arch.fidelity_dist = new double*[arch.positions];

	for (int i = 0; i < arch.positions; i++) {
		arch.fidelity_dist[i] = new double[arch.positions]();
	}
	for (std::set<edge>::iterator it = arch.graph.begin(); it != arch.graph.end(); it++) {
		arch.fidelity_dist[it->v1][it->v2] = it->fidelity;
	}
	for (int i = 0; i < arch.positions; i++) {
		for (int j = 0; j < arch.positions; j++) {
			if(arch.fidelity_dist[i][j] == 0 && arch.fidelity_dist[j][i] != 0) {
				arch.fidelity_dist[i][j] = arch.fidelity_dist[j][i]; 
			}
		}
	}
	std::cout << "Finished building of the dist table" << std::endl;

	return true;
}

/**
 * Deletes the allocated architecture properties
 */
void delete_architecture_properties() {
	delete[] arch.dist;
	delete[] arch.fidelity_dist;
	delete[] arch.initial_fidelities;
}


/**
 * sets the properties of the dijkstra node for position pos
 */
template<class DijkstraCmp>
static void set_dijkstra_node(dijkstra_node* nodes, std::priority_queue<dijkstra_node*, std::vector<dijkstra_node*>, DijkstraCmp>& queue,
					          const edge e, const int parent, const int pos, const bool contains_correct_edge) {
	if(nodes[pos].visited) {
		return;
	}
	double new_cost = DijkstraCmp::step(nodes + parent, e);
	if(nodes[pos].cost < 0 || DijkstraCmp::compare_parameters(new_cost, nodes[pos].cost, contains_correct_edge, nodes[pos].contains_correct_edge)) {
		nodes[pos].contains_correct_edge = contains_correct_edge;
		nodes[pos].cost                  = new_cost;

		queue.push(nodes + pos);
	} 
}

/**
 * the dijkstra algorithm calculates the distance from one node to all others
 */
template<class DijkstraCmp>
static void dijkstra(const std::set<edge>& graph, dijkstra_node* nodes, const int start) {
	std::priority_queue<dijkstra_node*, std::vector<dijkstra_node*>, DijkstraCmp> queue;
	queue.push(nodes + start);
	while(!queue.empty()) {
		dijkstra_node* current = queue.top();
		current->visited = true;
		queue.pop();
		int cur = current->pos;
		for (std::set<edge>::iterator it = graph.begin(); it != graph.end(); it++) {
			edge e = *it;
			if (cur == e.v1) { 
				set_dijkstra_node<DijkstraCmp>(nodes, queue, e, e.v1, e.v2, true);	
			} else if (cur == e.v2) {
				set_dijkstra_node<DijkstraCmp>(nodes, queue, e, e.v2, e.v1, current->contains_correct_edge);
			}
		}
	}
}

template<class DijkstraCmp>
static void build_dijkstra_table(const std::set<edge>& graph, double*** distance_table) {
	*distance_table = new double*[arch.positions];

	for (int i = 0; i < arch.positions; i++) {
		(*distance_table)[i] = new double[arch.positions];
	}

	for (int i = 0; i < arch.positions; i++) {
		dijkstra_node* nodes = new dijkstra_node[arch.positions];
		for (int j = 0; j < arch.positions; j++) {
			nodes[j].contains_correct_edge = false;
			nodes[j].visited               = false;
			nodes[j].pos                   = j;
			nodes[j].cost                  = -1;
		}
		DijkstraCmp::set_initial_cost(nodes + i);
		

		dijkstra<DijkstraCmp>(graph, nodes, i);

		for (int j = 0; j < arch.positions; j++) {
			if (i != j) {
				(*distance_table)[i][j] = DijkstraCmp::get_cost(nodes + j);
			} else {
				(*distance_table)[i][j] = 0;
			}
		}

		delete[] nodes;
	}
}

/**
 * builds the graph by parsing the configuration file
 */
static bool build_graph_from_file(const std::string input) {
	std::ifstream infile(input, std::ifstream::in);
	if(infile.fail()) {
		std::cerr << "Failed to open file '" << input << "'!" << std::endl;
		return false;
	} 

	std::string line;
	arch.graph.clear();
	

	if(std::getline(infile, line)) {
		int positions;
		sscanf(line.c_str(), "Positions: %d", &positions);
		init_arch_arrays(positions);
		std::cout << "Positions: " << arch.positions << std::endl;
	} else {
		std::cout << "First Line has to be: Positions: [0-9]*" << std::endl;
		return false;
	}
	if(std::getline(infile, line)) {
		int e1, e2, nread;
		double e3, e4;
		edge e;

		if(line.compare("QUBITS") == 0) {
			while(getline(infile, line)) {
				nread = sscanf(line.c_str(), "q%d: %lf,%lf", &e1, &e3, &e4);
				if(nread < 1) {
					break;
				}

				if(e1 >= arch.positions) {
					std::cerr << "File error: qubit out of range " << e1 << " >= " << arch.positions << std::endl;
				}					
				if(nread >= 2) {
					arch.singlequbit_fidelities[e1] = e4;
				}
				if(nread >= 1) {
					arch.initial_fidelities[e1] = e3;
				}	
			}
		}
		do {
			nread = sscanf(line.c_str(), "[%d,%d,%lf]", &e1, &e2, &e3);
			if (nread >= 3) {
				e.fidelity = e3;
				//std::cout << "Edge: " << "Fidelity " << e3 << std::endl;
			}
			if (nread >= 2) {
				//std::cout << "Edge: " << e1 << "  " << e2 << "  detected." << std::endl;
				e.v1 = e1;
				e.v2 = e2;
				arch.graph.insert(e);
			}
		} while (getline(infile, line));
	}
	infile.close();
	std::cout << "Finished reading of the coupling graph" << std::endl;
	return true;
}


/**
 * builds a arch.graph representing the coupling map of a linear architecture
 */
static void build_graph_linear(int nqubits) {
	arch.graph.clear();
    
	init_arch_arrays(nqubits);
	for(int i = 0; i < nqubits - 1; i++) {
		arch.graph.emplace(i,   i+1);
        arch.graph.emplace(i+1, i);
	}
}

/**
 * builds a graph representing the coupling map of IBM QX5
 */
static void build_graph_QX5() {
	arch.graph.clear();
	init_arch_arrays(16);
	arch.graph.emplace( 1,  0);
	arch.graph.emplace( 1,  2);
	arch.graph.emplace( 2,  3);
	arch.graph.emplace( 3, 14);
	arch.graph.emplace( 3,  4);
	arch.graph.emplace( 5,  4);
	arch.graph.emplace( 6,  5);
	arch.graph.emplace( 6, 11);
	arch.graph.emplace( 6,  7);
	arch.graph.emplace( 7, 10);
	arch.graph.emplace( 8,  7);
	arch.graph.emplace( 9,  8);
	arch.graph.emplace( 9, 10);
	arch.graph.emplace(11, 10);
	arch.graph.emplace(12,  5);
	arch.graph.emplace(12, 11);
	arch.graph.emplace(12, 13);
	arch.graph.emplace(13,  4);
	arch.graph.emplace(13, 14);
	arch.graph.emplace(15,  0);
	arch.graph.emplace(15, 14);
	arch.graph.emplace(15,  2);
}

/**
 * builds a graph representing the coupling map of IBM QX5
 */
static void build_graph_melbourne() {
	arch.graph.clear();
	init_arch_arrays(16);
	arch.graph.emplace( 0,  1, 1 - 0.05008); arch.graph.emplace( 1,  0, 1 - 0.05008);
	arch.graph.emplace( 1,  2, 1 - 0.02242); arch.graph.emplace( 2,  1, 1 - 0.02242);
	arch.graph.emplace( 2,  3, 1 - 0.03372); arch.graph.emplace( 3,  2, 1 - 0.03372);
	arch.graph.emplace( 3,  4, 1 - 0.02215); arch.graph.emplace( 4,  3, 1 - 0.02215);
	arch.graph.emplace( 4,  5, 1 - 0.03099); arch.graph.emplace( 5,  4, 1 - 0.03099);
	arch.graph.emplace( 5,  6, 1 - 0.04057); arch.graph.emplace( 6,  5, 1 - 0.04057);

	arch.graph.emplace( 0, 14, 1 - 0.05392); arch.graph.emplace(14,  0, 1 - 0.05392);
	arch.graph.emplace( 1, 13, 1 - 0.06918); arch.graph.emplace(13,  1, 1 - 0.06918);
	arch.graph.emplace( 2, 12, 1 - 0.05196); arch.graph.emplace(12,  2, 1 - 0.05196);
	arch.graph.emplace( 3, 11, 1 - 0.02423); arch.graph.emplace(11,  3, 1 - 0.02423);
	arch.graph.emplace( 4, 10, 1 - 0.03644); arch.graph.emplace(10,  4, 1 - 0.03644);
	arch.graph.emplace( 5,  9, 1 - 0.06800); arch.graph.emplace( 9,  5, 1 - 0.06800);
	//arch.graph.emplace( 6,  8, 1 - 1      ); arch.graph.emplace( 8,  6, 1 - 1      );
	
	//artificial
	arch.graph.emplace( 6,  8, 1 - 0.06800); arch.graph.emplace( 8,  6, 1 - 0.06800);


	arch.graph.emplace(14, 13, 1 - 0.07199); arch.graph.emplace(13, 14, 1 - 0.07199);
	arch.graph.emplace(13, 12, 1 - 0.03570); arch.graph.emplace(12, 13, 1 - 0.03570);
	arch.graph.emplace(12, 11, 1 - 0.02706); arch.graph.emplace(11, 12, 1 - 0.02706);
	arch.graph.emplace(11, 10, 1 - 0.02300); arch.graph.emplace(10, 11, 1 - 0.02300);
	arch.graph.emplace(10,  9, 1 - 0.03827); arch.graph.emplace( 9, 10, 1 - 0.03827);
	//arch.graph.emplace( 9,  8, 1 - 1      ); arch.graph.emplace( 8,  9, 1 - 1      );
	//arch.graph.emplace( 8,  7, 1 - 1      ); arch.graph.emplace( 7,  8, 1 - 1      );
	
	//artificial
	arch.graph.emplace( 9,  8, 1 - 0.03827); arch.graph.emplace( 8,  9, 1 - 0.03827);
	arch.graph.emplace( 8,  7, 1 - 0.03570); arch.graph.emplace( 7,  8, 1 - 0.03570);
	arch.graph.emplace(15, 14, 1 - 0.03570); arch.graph.emplace(14,  15, 1 - 0.03570);

	arch.singlequbit_fidelities[ 0] = 1 - 0.001861414;
	arch.singlequbit_fidelities[ 1] = 1 - 0.000842034;
	arch.singlequbit_fidelities[ 2] = 1 - 0.002173244;
	arch.singlequbit_fidelities[ 3] = 1 - 0.000505368;
	arch.singlequbit_fidelities[ 4] = 1 - 0.001262586;
	arch.singlequbit_fidelities[ 5] = 1 - 0.002429494;
	arch.singlequbit_fidelities[ 6] = 1 - 0.000876746;
	arch.singlequbit_fidelities[ 7] = 1 - 0.001687877;
	arch.singlequbit_fidelities[ 8] = 1 - 0.000344825;
	arch.singlequbit_fidelities[ 9] = 1 - 0.00250304;
	arch.singlequbit_fidelities[10] = 1 - 0.00105976;
	arch.singlequbit_fidelities[11] = 1 - 0.000608612;
	arch.singlequbit_fidelities[12] = 1 - 0.005052505;
	arch.singlequbit_fidelities[13] = 1 - 0.001711347;
	arch.singlequbit_fidelities[14] = 1 - 0.000669206;
	arch.singlequbit_fidelities[15] = 1 - 0.000669206; //artificial
}

/*
static void build_graph_QX5() {
	arch.graph.clear();
	init_arch_arrays(16);
	edge e (1, 0);
	arch.graph.insert(e);
	e.v1 = 1;
	e.v2 = 2;
	arch.graph.insert(e);
	e.v1 = 2;
	e.v2 = 3;
	arch.graph.insert(e);
	e.v1 = 3;
	e.v2 = 14;
	arch.graph.insert(e);
	e.v1 = 3;
	e.v2 = 4;
	arch.graph.insert(e);
	e.v1 = 5;
	e.v2 = 4;
	arch.graph.insert(e);
	e.v1 = 6;
	e.v2 = 5;
	arch.graph.insert(e);
	e.v1 = 6;
	e.v2 = 11;
	arch.graph.insert(e);
	e.v1 = 6;
	e.v2 = 7;
	arch.graph.insert(e);
	e.v1 = 7;
	e.v2 = 10;
	arch.graph.insert(e);
	e.v1 = 8;
	e.v2 = 7;
	arch.graph.insert(e);
	e.v1 = 9;
	e.v2 = 8;
	arch.graph.insert(e);
	e.v1 = 9;
	e.v2 = 10;
	arch.graph.insert(e);
	e.v1 = 11;
	e.v2 = 10;
	arch.graph.insert(e);
	e.v1 = 12;
	e.v2 = 5;
	arch.graph.insert(e);
	e.v1 = 12;
	e.v2 = 11;
	arch.graph.insert(e);
	e.v1 = 12;
	e.v2 = 13;
	arch.graph.insert(e);
	e.v1 = 13;
	e.v2 = 4;
	arch.graph.insert(e);
	e.v1 = 13;
	e.v2 = 14;
	arch.graph.insert(e);
	e.v1 = 15;
	e.v2 = 0;
	arch.graph.insert(e);
	e.v1 = 15;
	e.v2 = 14;
	arch.graph.insert(e);
	e.v1 = 15;
	e.v2 = 2;
	arch.graph.insert(e);
}
*/