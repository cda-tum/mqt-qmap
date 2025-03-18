#pragma once

#include "Definitions.hpp"
#include "na/azac/Architecture.hpp"

#include <array>
#include <cstddef>
#include <functional>
#include <nlohmann/json_fwd.hpp>
#include <optional>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na::azac {
/// class to find a qubit layout based on vertex matching of a weighted
/// bipartite graph
class VMPlacer {
  friend class VMPlacerTest_MinimumWeightFullBipartiteMatching1_Test;
  friend class VMPlacerTest_MinimumWeightFullBipartiteMatching2_Test;
  friend class VMPlacerTest_MinimumWeightFullBipartiteMatchingExceptions_Test;

  std::reference_wrapper<const Architecture> architecture_;
  /// If true, during the initial placement the atoms are placed starting in the
  /// last row instead of the first row in the first SLM
  bool reverseInitialPlacement_ = false;

  /// this flag indicates whether the  placement should use a window when
  /// selecting potential free sites
  bool useWindow_ = true;
  size_t windowSize_ = 0;

  /// this flag indicates whether the placement between gates is dynamic, i.e.,
  /// if this flag is false, the initial placement is used after all gates
  bool dynamicPlacement_ = true;

  // todo: Why is that?
  constexpr static double costAtomTransfer_ = 0.9999;

public:
  VMPlacer(const Architecture& architecture, const nlohmann::json& config);
  [[nodiscard]] auto
  place(const size_t nQubits,
        const std::vector<std::vector<std::array<qc::Qubit, 2>>>&
            twoQubitGateLayers,
        const std::vector<std::unordered_set<qc::Qubit>>& reuseQubits)
      -> std::vector<std::vector<
          std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>>;
  /// generate qubit initial layout
  auto makeInitialPlacement(size_t nQubits) const -> std::vector<
      std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>;

private:
  /// @note implemented following pseudocode in
  /// https://www2.eecs.berkeley.edu/Pubs/TechRpts/1978/ERL-m-78-67.pdf
  [[nodiscard]] static auto minimumWeightFullBipartiteMatching(
      const std::vector<std::vector<std::optional<double>>>& costMatrix)
      -> std::vector<size_t>;

  [[nodiscard]] auto computeMovementCostBetweenPlacements(
      const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                   size_t>>& placementBefore,
      const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                   size_t>>& placementAfter) const -> double;

  [[nodiscard]] auto computeLayersMovementCost(
      const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                   size_t>>& placementBefore,
      const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                   size_t>>& placementBetween,
      const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                   size_t>>& placementAfter) const -> double;

  [[nodiscard]] auto filterMapping(
      const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                   size_t>>& previousGatePlacement,
      const std::pair<std::vector<std::tuple<std::reference_wrapper<const SLM>,
                                             size_t, size_t>>,
                      std::vector<std::tuple<std::reference_wrapper<const SLM>,
                                             size_t, size_t>>>&
          placementsWithoutReuse,
      const std::pair<std::vector<std::tuple<std::reference_wrapper<const SLM>,
                                             size_t, size_t>>,
                      std::vector<std::tuple<std::reference_wrapper<const SLM>,
                                             size_t, size_t>>>&
          placementsWithReuse) const
      -> std::pair<std::vector<std::tuple<std::reference_wrapper<const SLM>,
                                          size_t, size_t>>,
                   std::vector<std::tuple<std::reference_wrapper<const SLM>,
                                          size_t, size_t>>>;
  /// generate gate mapping based on minimum weight matching for the first
  /// layer of gates
  [[nodiscard]] auto placeGatesInEntanglementZone(
      const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                   size_t>>& previousQubitPlacement,
      const std::unordered_set<qc::Qubit>& reuseQubits,
      const std::vector<std::array<qc::Qubit, 2>>& twoQubitGates,
      const std::vector<std::array<qc::Qubit, 2>>& nextTwoQubitGates,
      bool reuse) const
      -> std::vector<
          std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>;

  /// Generate qubit mapping based on minimum weight matching.
  auto placeQubitsInStorageZone(
      const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                   size_t>>& initialPlacement,
      const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                   size_t>>& previousGatePlacement,
      const std::unordered_set<qc::Qubit>& reuseQubits,
      const std::vector<std::array<qc::Qubit, 2>>& nextTwoQubitGates,
      bool reuse) const
      -> std::vector<
          std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>;
};
} // namespace na::azac
