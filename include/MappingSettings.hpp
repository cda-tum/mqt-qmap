/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#ifndef QMAP_MAPPINGSETTINGS_HPP
#define QMAP_MAPPINGSETTINGS_HPP

#include "nlohmann/json.hpp"

/// Identity: q_i -> Q_i
/// Static: first layer is mapped q_c -> Q_c and q_t -> Q_t
/// Dynamic: Layout is generated on demand upon encountering a specific gate
enum class InitialLayoutStrategy {
	Identity, Static, Dynamic, None
};
static std::string toString(const InitialLayoutStrategy strategy) {
	switch (strategy) {
		case InitialLayoutStrategy::Identity:
			return "identity";
		case InitialLayoutStrategy::Static:
			return "static";
		case InitialLayoutStrategy::Dynamic:
			return "dynamic";
		case InitialLayoutStrategy::None:
			return "none";
	}
	return " ";
}
// map InitialLayoutStrategy values to JSON as strings
NLOHMANN_JSON_SERIALIZE_ENUM( InitialLayoutStrategy, {
	{InitialLayoutStrategy::None, "none"},
	{InitialLayoutStrategy::Identity, "identity"},
	{InitialLayoutStrategy::Static, "static"},
	{InitialLayoutStrategy::Dynamic, "dynamic"},
})

enum class LayeringStrategy {
	IndividualGates, DisjointQubits, OddGates, QubitTriangle, None
};
static std::string toString(const LayeringStrategy strategy) {
	switch (strategy) {
		case LayeringStrategy::IndividualGates:
			return "individual_gates";
		case LayeringStrategy::DisjointQubits:
			return "disjoint_qubits";
		case LayeringStrategy::OddGates:
			return "odd_gates";
		case LayeringStrategy::QubitTriangle:
			return "qubit_triangle";
		case LayeringStrategy::None:
			return "none";
	}
	return " ";
}
// map LayeringStrategy values to JSON as strings
NLOHMANN_JSON_SERIALIZE_ENUM( LayeringStrategy, {
	{LayeringStrategy::None, "none"},
	{LayeringStrategy::IndividualGates, "individual_gates"},
	{LayeringStrategy::DisjointQubits, "disjoint_qubits"},
	{LayeringStrategy::OddGates, "odd_gates"},
	{LayeringStrategy::QubitTriangle, "qubit_triangle"},
})

struct MappingSettings {
	MappingSettings() = default;

	unsigned int timeout = 3600000; // 60min timeout
	void setTimeout(unsigned int sec) { timeout = sec;}

	LayeringStrategy layeringStrategy = LayeringStrategy::None;

	/// Settings for heuristic approach
	InitialLayoutStrategy initialLayoutStrategy = InitialLayoutStrategy::None;

	bool admissibleHeuristic = true;
	bool verbose = false;

	bool lookahead = true;
	int nrLookaheads = 15;
	double firstLookaheadFactor = 0.75;
	double lookaheadFactor = 0.5;

};

#endif //QMAP_MAPPINGSETTINGS_HPP
