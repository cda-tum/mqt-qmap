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

struct Configuration {
  Configuration() = default;

  // which method to use
  Method method = Method::Heuristic;

  bool preMappingOptimizations  = true;
  bool postMappingOptimizations = true;

  bool addMeasurementsToMappedCircuit = true;

  bool verbose = false;

  // map to particular subgraph of architecture (in exact mapper)
  std::set<std::uint16_t> subgraph{};

  // how to cluster the gates into layers
  Layering layering = Layering::None;

  // initial layout to use for heuristic approach
  InitialLayout initialLayout = InitialLayout::None;

  // lookahead scheme settings
  bool        lookahead            = true;
  bool        admissibleHeuristic  = true;
  bool        considerFidelity     = true;
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
  Encoding          encoding          = Encoding::Naive;
  CommanderGrouping commanderGrouping = CommanderGrouping::Halves;

  // use qubit subsets in exact mapper
  bool useSubsets = true;

  // include WCNF file in results of exact mapper
  bool includeWCNF = false;

  // limit the number of considered swaps
  bool          enableSwapLimits = true;
  SwapReduction swapReduction    = SwapReduction::CouplingLimit;
  std::size_t   swapLimit        = 0;
  bool          useBDD           = false;

  [[nodiscard]] nlohmann::json json() const {
    nlohmann::json config{};
    config["method"]            = ::toString(method);
    config["layering_strategy"] = ::toString(layering);
    if (!subgraph.empty()) {
      config["subgraph"] = subgraph;
    }
    config["pre_mapping_optimizations"]  = preMappingOptimizations;
    config["post_mapping_optimizations"] = postMappingOptimizations;
    config["add_measurements_to_mapped_circuit"] =
        addMeasurementsToMappedCircuit;
    config["verbose"] = verbose;

    if (method == Method::Heuristic) {
      auto& heuristic             = config["settings"];
      heuristic["initial_layout"] = ::toString(initialLayout);
      if (lookahead) {
        auto& lookaheadSettings                   = heuristic["lookahead"];
        lookaheadSettings["admissible_heuristic"] = admissibleHeuristic;
        lookaheadSettings["consider_fidelity"]    = considerFidelity;
        lookaheadSettings["lookaheads"]           = nrLookaheads;
        lookaheadSettings["first_factor"]         = firstLookaheadFactor;
        lookaheadSettings["factor"]               = lookaheadFactor;
      }
      if (useTeleportation) {
        auto& teleportation     = heuristic["teleportation"];
        teleportation["qubits"] = teleportationQubits;
        teleportation["seed"]   = teleportationSeed;
        teleportation["fake"]   = teleportationFake;
      }
    }

    if (method == Method::Exact) {
      auto& exact       = config["settings"];
      exact["timeout"]  = timeout;
      exact["encoding"] = ::toString(encoding);
      if (encoding == Encoding::Commander || encoding == Encoding::Bimander) {
        exact["commander_grouping"] = ::toString(commanderGrouping);
      }
      exact["include_WCNF"] = includeWCNF;
      exact["use_subsets"]  = useSubsets;
      if (enableSwapLimits) {
        auto& limits             = exact["limits"];
        limits["swap_reduction"] = ::toString(swapReduction);
        if (swapLimit > 0) {
          limits["swap_limit"] = swapLimit;
        }
        limits["use_bdd"] = useBDD;
      }
    }

    return config;
  }
  [[nodiscard]] std::string toString() const { return json().dump(2); }

  void setTimeout(const std::size_t sec) { timeout = sec; }
};
