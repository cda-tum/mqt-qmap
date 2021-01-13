/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
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
using matrix = std::vector<expr_vector>;

/// Main structure representing the circuit and mapping functionality
class ExactMapper : public Mapper{
	using Mapper::Mapper;

protected:
	// inputs
	std::vector<unsigned long> reducedLayerIndices{};
	std::vector<std::vector<std::pair<unsigned short, unsigned short>>> mappingSwaps{};
	void coreMappingRoutine(const std::set<unsigned short>& qubitChoice, const CouplingMap& rcm, MappingResults& choiceResults, std::vector<std::vector<std::pair<unsigned short, unsigned short>>>& swaps);
	void initResults() override;

public:
	void map(const MappingSettings& settings) override;
};

#endif /* EXACT_MAPPER_hpp */
