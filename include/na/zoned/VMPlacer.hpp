/*
 * Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
 * Copyright (c) 2025 Munich Quantum Software Company GmbH
 * All rights reserved.
 *
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License
 */

#pragma once

#include "ir/Definitions.hpp"
#include "na/zoned/Architecture.hpp"
#include "na/zoned/Types.hpp"

#include <cstddef>
#include <functional>
#include <nlohmann/json_fwd.hpp>
#include <optional>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na::zoned {
/// class to find a qubit layout based on vertex matching of a weighted
/// bipartite graph
class VMPlacer {
  friend class VMPlacerTest_MinimumWeightFullBipartiteMatching1_Test;
  friend class VMPlacerTest_MinimumWeightFullBipartiteMatching2_Test;
  friend class VMPlacerTest_MinimumWeightFullBipartiteMatchingExceptions_Test;
  friend class VMPlacerTest_MinimumWeightFullBipartiteMatchingEmpty_Test;

  std::reference_wrapper<const Architecture> architecture_;
  /// If true, during the initial placement the atoms are placed starting in the
  /// last row instead of the first row in the first SLM
  bool reverseInitialPlacement_ = false;

  /// this flag indicates whether the  placement should use a window when
  /// selecting potential free sites
  bool useWindow_ = true;
  size_t windowSize_ = 10;

  /// this flag indicates whether the placement between gates is dynamic, i.e.,
  /// if this flag is false, the initial placement is used after all gates
  bool dynamicPlacement_ = true;

  // todo: Why is that?
  constexpr static double costAtomTransfer_ = 0.9999;

public:
  /// Create a VMPlacer based on the given architecture and configuration
  VMPlacer(const Architecture& architecture, const nlohmann::json& config);
  /// generate qubit placement based on minimum weight matching
  [[nodiscard]] auto
  place(size_t nQubits,
        const std::vector<TwoQubitGateLayer>& twoQubitGateLayers,
        const std::vector<std::unordered_set<qc::Qubit>>& reuseQubits)
      -> std::vector<Placement>;

private:
  /// generate qubit initial layout
  auto makeInitialPlacement(size_t nQubits) const -> Placement;

  /// @note implemented following pseudocode in
  /// https://www2.eecs.berkeley.edu/Pubs/TechRpts/1978/ERL-m-78-67.pdf
  [[nodiscard]] static auto minimumWeightFullBipartiteMatching(
      const std::vector<std::vector<std::optional<double>>>& costMatrix)
      -> std::vector<size_t>;

  [[nodiscard]] auto
  computeMovementCostBetweenPlacements(const Placement& placementBefore,
                                       const Placement& placementAfter) const
      -> double;

  [[nodiscard]] auto
  computeLayersMovementCost(const Placement& placementBefore,
                            const Placement& placementBetween,
                            const Placement& placementAfter) const -> double;

  [[nodiscard]] auto filterMapping(
      const Placement& previousGatePlacement,
      const std::pair<Placement, Placement>& placementsWithoutReuse,
      const std::pair<Placement, Placement>& placementsWithReuse) const
      -> std::pair<Placement, Placement>;
  /// generate gate mapping based on minimum weight matching for the first
  /// layer of gates
  [[nodiscard]] auto
  placeGatesInEntanglementZone(const Placement& previousQubitPlacement,
                               const std::unordered_set<qc::Qubit>& reuseQubits,
                               const TwoQubitGateLayer& twoQubitGates,
                               const TwoQubitGateLayer& nextTwoQubitGates,
                               bool reuse) const -> Placement;

  /// Generate qubit mapping based on minimum weight matching.
  auto placeAtomsInStorageZone(const Placement& initialPlacement,
                               const Placement& previousGatePlacement,
                               const std::unordered_set<qc::Qubit>& reuseQubits,
                               const TwoQubitGateLayer& nextTwoQubitGates,
                               bool reuse) const -> Placement;
};
} // namespace na::zoned
