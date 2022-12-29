//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "Mapper.hpp"

#include <algorithm>
#include <bitset>
#include <chrono>
#include <cmath>
#include <functional>
#include <set>
#include <unordered_set>

using Swap        = std::pair<std::uint16_t, std::uint16_t>;
using Swaps       = std::vector<Swap>;
using QubitChoice = std::set<std::uint16_t>;

/// Main structure representing the circuit and mapping functionality
class ExactMapper : public Mapper {
  using Mapper::Mapper;

protected:
  // inputs
  std::vector<std::size_t> reducedLayerIndices{};
  std::vector<Swaps>       mappingSwaps{};
  void                     coreMappingRoutine(const QubitChoice& qubitChoice,
                                              const CouplingMap& rcm, MappingResults& choiceResults,
                                              std::vector<Swaps>& swaps, std::size_t limit,
                                              unsigned int timeout);

public:
  void map(const Configuration& settings) override;
};
