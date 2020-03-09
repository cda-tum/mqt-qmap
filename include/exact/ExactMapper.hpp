/*
Minimal Mapping of Quantum Circuits to IBM QX Architectures by JKU Linz, Austria

Developer: Robert Wille, Lukas Burgholzer, Alwin Zulehner

For more information, please visit http://iic.jku.at/eda/research/ibm_qx_mapping

If you have any questions feel free to contact us using
robert.wille@jku.at, lukas.burgholzer@jku.at or alwin.zulehner@jku.at

If you use the compiler for your research, we would be thankful if you referred to it
by citing the following publication:

@inproceedings{wille2019mapping,
    title={Mapping Quantum Circuits to {IBM QX} Architectures Using the Minimal Number of {SWAP} and {H} Operations},
    author={Wille, Robert and Burgholzer, Lukas and Zulehner, Alwin},
    booktitle={Design Automation Conference},
    year={2019}
}
*/

#ifndef EXACT_MAPPER_hpp
#define EXACT_MAPPER_hpp

#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <functional>
#include <chrono>
#include <set>
#include <unordered_set>

#include <z3++.h>
#include "Mapper.hpp"

using namespace z3;

/// Main structure representing the circuit and mapping functionality
class ExactMapper : public Mapper{
	using Mapper::Mapper;

	static constexpr bool VERBOSE = true;

protected:
	// inputs
	std::vector<unsigned long> reducedLayerIndices{};
	std::vector<std::vector<std::pair<unsigned short, unsigned short>>> mappingSwaps{};
	void coreMappingRoutine(const std::unordered_set<unsigned short>& qubitChoice, const CouplingMap& rcm, MappingResults& choiceResults, std::vector<std::vector<std::pair<unsigned short, unsigned short>>>& swaps);
	void initResults() override;

public:
	void map(const MappingSettings& settings) override;
};

#endif /* EXACT_MAPPER_hpp */
