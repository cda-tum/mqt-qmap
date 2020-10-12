/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#include <utils.hpp>

void Dijkstra::build_table(unsigned short n, const std::set<Edge>& couplingMap, Matrix& distanceTable, const std::function<double(const Node&)>& cost) {
	distanceTable.clear();
	distanceTable.resize(n, std::vector<double>(n, -1.));

	for (unsigned short i=0; i<n; ++i) {
		std::vector<Dijkstra::Node> nodes(n);
		for (unsigned short j=0; j<n; ++j) {
			nodes.at(j).contains_correct_edge = false;
			nodes.at(j).visited = false;
			nodes.at(j).pos = j;
			nodes.at(j).cost = -1.;
		}

		nodes.at(i).cost = 0.;

		dijkstra(couplingMap, nodes, i);

		if (VERBOSE) {
			for (const auto& node: nodes) {
				std::cout << node.cost << " ";
			}
			std::cout << std::endl;
		}

		for (int j=0; j<n; ++j) {
			if (i == j) {
				distanceTable.at(i).at(j) = 0;
			} else {
				distanceTable.at(i).at(j) = cost(nodes.at(j));
			}
		}
	}
}

void Dijkstra::dijkstra(const CouplingMap& couplingMap, std::vector<Node>& nodes, unsigned short start) {
	std::priority_queue<Node*> queue{};
	queue.push(&nodes.at(start));
	while (!queue.empty()) {
		auto current = queue.top();
		current->visited = true;
		queue.pop();
		auto pos = current->pos;

		for (const auto& edge: couplingMap) {
			short to = -1;
			bool correctEdge = false;
			if (pos == edge.first) {
				to = edge.second;
				correctEdge = true;
			} else if (pos == edge.second) {
				to = edge.first;
			}
			if (to != -1) {
				if (nodes.at(to).visited)
					continue;

				Node new_node;
				new_node.cost = current->cost + 1.0;
				new_node.pos = to;
				new_node.contains_correct_edge = correctEdge;
				if (nodes.at(to).cost < 0 || new_node < nodes.at(to)) {
					nodes.at(to) = new_node;
					queue.push(&nodes.at(to));
				}
			}
		}
	}
}



/// Create a string representation of a given permutation
/// \param pi permutation
/// \return string representation of pi
std::string printPi(std::vector<unsigned short>& pi){
	if (std::is_sorted(pi.begin(),pi.end())) {
		return "( )";
	}

	std::stringstream perm{};
	perm << '(';
	for (unsigned long i = 0; i < pi.size()-1; i++) {
		perm << pi[i] << ',';
	}
	perm << pi[pi.size()-1] << ')';

	return perm.str();
}

/// Simple depth-first-search implementation used to check whether a given subset of qubits is
/// connected on the given architecture
/// \param current index of current qubit
/// \param visited visited qubits
/// \param cm coupling map of architecture
void dfs(unsigned short current, std::set<unsigned short>& visited, CouplingMap& rcm) {
	for (auto edge: rcm) {
		if (edge.first == current) {
			if(!visited.count(edge.second)) {
				visited.insert(edge.second);
				dfs(edge.second, visited, rcm);
			}
		} else if (edge.second == current) {
			if (!visited.count(edge.first)) {
				visited.insert(edge.first);
				dfs(edge.first, visited, rcm);
			}
		}
	}
}

/// Helper function returning correct 1D array index for 3D array
/// \param k first index
/// \param i second index
/// \param j third index
/// \return index in 1D array
unsigned long idx(unsigned int k, unsigned short i, unsigned short j, const std::set<unsigned short>& iValues, const std::set<unsigned short>& jValues) {
	unsigned short counti = 0;
	for (unsigned short iVal : iValues) {
		if (iVal == i) break;
		counti++;
	}
	unsigned short countj = 0;
	for (unsigned short jVal : jValues) {
		if (jVal == j) break;
		countj++;
	}

	return k*jValues.size()*iValues.size() + counti*jValues.size() + countj;
}

unsigned long idx(unsigned int k, unsigned short i, unsigned short j, const std::set<unsigned short>& iValues, unsigned short nj) {
	unsigned short counti = 0;
	for (unsigned short iVal : iValues) {
		if (iVal == i) break;
		counti++;
	}

	return k*nj*iValues.size() + counti*nj + j;
}
