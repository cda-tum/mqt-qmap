/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#ifndef QMAP_UTILS_HPP
#define QMAP_UTILS_HPP

#include <stdexcept>
#include <algorithm>
#include <set>
#include <vector>
#include <string>
#include <queue>
#include <functional>
#include <iostream>
#include <sstream>
#include <unordered_set>
#include <unordered_map>

using Matrix = std::vector<std::vector<double>>;
using Edge = std::pair<unsigned short, unsigned short>;
using CouplingMap = std::set<Edge>;

class QMAPException : public std::runtime_error {
	std::string msg;
public:
	explicit QMAPException(std::string  msg) : std::runtime_error("QMAP Exception"), msg(std::move(msg)) { }

	const char *what() const noexcept override {
		return msg.c_str();
	}
};

/// Computes n! recursively
/// \param n interger to compute factorial of
/// \return n!
static inline
unsigned long factorial(unsigned long n) {
	if (n == 1)
		return 1;
	else
		return n * factorial(n - 1);
}

class Dijkstra {
	static constexpr bool VERBOSE = false;

public:
	struct Node {
		bool   contains_correct_edge = false;
		bool   visited = false;
		int    pos = -1;
		double cost = -1.;
	};

	static void build_table(unsigned short n, const std::set<Edge>& graph, Matrix& distanceTable, const std::function<double(const Node&)>& cost);

protected:

	static void dijkstra(const CouplingMap& couplingMap, std::vector<Node>& nodes, unsigned short start);

};

inline bool operator<(const Dijkstra::Node& x, const Dijkstra::Node& y) {
	if (x.cost != y.cost) {
		return x.cost < y.cost;
	}
	return x.contains_correct_edge && ! y.contains_correct_edge;
}

/// Iterating routine through all combinations
/// \tparam Iterator iterator type
/// \param first iterator to beginning
/// \param k current iterator
/// \param last iterator to end
/// \return true if another combination was found
template<typename Iterator> bool next_combination(Iterator first, Iterator k, Iterator last) {
	/* Credits: Thomas Draper */
	if ((first == last) || (first == k) || (last == k))
		return false;
	Iterator itr1 = first;
	Iterator itr2 = last;
	++itr1;
	if (last == itr1)
		return false;
	itr1 = last;
	--itr1;
	itr1 = k;
	--itr2;
	while (first != itr1) {
		if (*--itr1 < *itr2) {
			Iterator j = k;
			while (!(*itr1 < *j)) ++j;
			std::iter_swap(itr1, j);
			++itr1;
			++j;
			itr2 = k;
			std::rotate(itr1, j, last);
			while (last != j) {
				++j;
				++itr2;
			}
			std::rotate(k, itr2, last);
			return true;
		}
	}
	std::rotate(first, k, last);
	return false;
}

/// Create a string representation of a given permutation
/// \param pi permutation
/// \return string representation of pi
std::string printPi(std::vector<unsigned short>& pi);

/// Simple depth-first-search implementation used to check whether a given subset of qubits is
/// connected on the given architecture
/// \param current index of current qubit
/// \param visited visited qubits
/// \param cm coupling map of architecture
void dfs(unsigned short current, std::set<unsigned short>& visited, CouplingMap& rcm);

/// Helper function returning correct 1D array index for 3D array
/// \param k first index
/// \param i second index
/// \param j third index
/// \return index in 1D array
unsigned long idx(unsigned int k, unsigned short i, unsigned short j, const std::set<unsigned short>& iValues, const std::set<unsigned short>& jValues);
unsigned long idx(unsigned int k, unsigned short i, unsigned short j, const std::set<unsigned short>& iValues, unsigned short nj);

#endif //QMAP_UTILS_HPP
