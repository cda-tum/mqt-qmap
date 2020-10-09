//
// Created by Lukas Burgholzer on 29.02.20.
//

#ifndef QMAP_MAPPINGSETTINGS_HPP
#define QMAP_MAPPINGSETTINGS_HPP

/// Identity: q_i -> Q_i
/// Static: first layer is mapped q_c -> Q_c and q_t -> Q_t
/// Dynamic: Layout is generated on demand upon encountering a specific gate
enum InitialLayoutStrategy {
	Identity, Static, Dynamic
};

enum LayeringStrategy {
	IndividualGates, DisjointQubits, OddGates, QubitTriangle
};

struct MappingSettings {
	MappingSettings() = default;

	unsigned int timeout = 3600000; // 60min timeout
	void setTimeout(unsigned int sec) { timeout = sec;}

	LayeringStrategy layeringStrategy = IndividualGates;

	/// Settings for heuristic approach
	InitialLayoutStrategy initialLayoutStrategy = Dynamic;

	bool admissibleHeuristic = true;
	bool verbose = false;

	bool lookahead = true;
	int nrLookaheads = 15;
	double firstLookaheadFactor = 0.75;
	double lookaheadFactor = 0.5;

};

#endif //QMAP_MAPPINGSETTINGS_HPP
