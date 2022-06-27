/*
 * This file is part of the MQT QMAP library which is released under the MIT license.
 * See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
 */

#ifndef EXACT_MAPPER_hpp
#define EXACT_MAPPER_hpp

#include "Encodings.hpp"
#include "LogicBlock/LogicBlock.hpp"
#include "Mapper.hpp"

#include <algorithm>
#include <bitset>
#include <chrono>
#include <cmath>
#include <functional>
#include <set>
#include <unordered_set>
#include <z3++.h>

using namespace z3;
using matrix      = std::vector<expr_vector>;
using Swaps       = std::vector<std::pair<unsigned short, unsigned short>>;
using QubitChoice = std::set<unsigned short>;

/// Main structure representing the circuit and mapping functionality
class ExactMapper: public Mapper {
    using Mapper::Mapper;

protected:
    // inputs
    std::vector<unsigned long> reducedLayerIndices{};
    std::vector<Swaps>         mappingSwaps{};
    void                       coreMappingRoutine(const QubitChoice&  qubitChoice,
                                                  const CouplingMap&  rcm,
                                                  MappingResults&     choiceResults,
                                                  std::vector<Swaps>& swaps,
                                                  long unsigned int   limit,
                                                  unsigned int        timeout);

public:
    void map(const Configuration& settings) override;
};

#endif /* EXACT_MAPPER_hpp */
