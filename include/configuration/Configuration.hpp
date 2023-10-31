//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "CommanderGrouping.hpp"
#include "Encoding.hpp"
#include "InitialLayout.hpp"
#include "Layering.hpp"
#include "Method.hpp"
#include "SwapReduction.hpp"
#include "nlohmann/json.hpp"

#include <set>

struct Configuration {
  Configuration() = default;

  // which method to use
  Method method = Method::Heuristic;

  bool preMappingOptimizations  = true;
  bool postMappingOptimizations = true;

  bool addMeasurementsToMappedCircuit = true;
  bool swapOnFirstLayer               = false;

  bool        verbose = false;
  bool        debug   = false;
  std::string dataLoggingPath;

  // map to particular subgraph of architecture (in exact mapper)
  std::set<std::uint16_t> subgraph{};

  // how to cluster the gates into layers
  Layering layering = Layering::None;

  // initial layout to use for heuristic approach
  InitialLayout initialLayout = InitialLayout::None;

  // lookahead scheme settings
  bool        lookahead            = true;
  bool        admissibleHeuristic  = true;
  bool        considerFidelity     = false;
  std::size_t nrLookaheads         = 15;
  double      firstLookaheadFactor = 0.75;
  double      lookaheadFactor      = 0.5;

  // teleportation settings
  bool          useTeleportation    = false;
  std::size_t   teleportationQubits = 0;
  std::uint64_t teleportationSeed   = 0;
  bool          teleportationFake   = false;

  // timeout merely affects exact mapper
  std::size_t timeout = 3600000; // 60min timeout

  // encoding of at most and exactly one constraints in exact mapper
  Encoding          encoding          = Encoding::Commander;
  CommanderGrouping commanderGrouping = CommanderGrouping::Fixed3;

  // use qubit subsets in exact mapper
  bool useSubsets = true;

  // include WCNF file in results of exact mapper
  bool includeWCNF = false;

  // limit the number of considered swaps
  bool          enableSwapLimits = true;
  SwapReduction swapReduction    = SwapReduction::CouplingLimit;
  std::size_t   swapLimit        = 0;
  bool          useBDD           = false;

  [[nodiscard]] nlohmann::json json() const;
  [[nodiscard]] std::string    toString() const { return json().dump(2); }

  [[nodiscard]] bool dataLoggingEnabled() const {
    return !dataLoggingPath.empty();
  }

  void               setTimeout(const std::size_t sec) { timeout = sec; }
  [[nodiscard]] bool swapLimitsEnabled() const {
    return (swapReduction != SwapReduction::None) && enableSwapLimits;
  }
};
