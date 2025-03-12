#pragma once

#include "Definitions.hpp"
#include "na/azac/Architecture.hpp"

#include <cstddef>
#include <functional>
#include <nlohmann/json_fwd.hpp>
#include <optional>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na {
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
        const std::vector<std::vector<std::pair<qc::Qubit, qc::Qubit>>>&
            twoQubitGateLayers,
        const std::vector<std::unordered_set<qc::Qubit>>& reuseQubits)
      -> std::vector<std::vector<
          std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>> {
    std::vector<std::vector<
        std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>>
        placement;
    placement.reserve((2 * twoQubitGateLayers.size()) + 1);
    placement.emplace_back(makeInitialPlacement(nQubits));
    // early return if no two-qubit gates are present
    if (twoQubitGateLayers.empty()) {
      return placement;
    }
    placement.emplace_back(placeGatesInEntanglementZone(
        placement.front(), std::unordered_set<qc::Qubit>{},
        twoQubitGateLayers.front(),
        twoQubitGateLayers.size() > 1
            ? twoQubitGateLayers[1]
            : std::vector<std::pair<qc::Qubit, qc::Qubit>>{},
        false));
    for (size_t layer = 0; layer < twoQubitGateLayers.size(); ++layer) {
      // first compute the next qubit and gate placement without reusing atoms
      std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>
          qubitPlacementWithoutReuse;
      if (dynamicPlacement_) {
        qubitPlacementWithoutReuse = placeQubitsInStorageZone(
            placement.front(), placement.back(), reuseQubits[layer],
            twoQubitGateLayers.size() > layer + 1
                ? twoQubitGateLayers[layer + 1]
                : std::vector<std::pair<qc::Qubit, qc::Qubit>>{},
            false);
      } else {
        // keep the initial mapping for static placement
        qubitPlacementWithoutReuse = placement.front();
      }
      if (layer + 1 < twoQubitGateLayers.size()) {
        const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                     size_t>>& gatePlacementWithoutReuse =
            placeGatesInEntanglementZone(
                qubitPlacementWithoutReuse, reuseQubits[layer],
                twoQubitGateLayers[layer + 1],
                twoQubitGateLayers.size() > layer + 2
                    ? twoQubitGateLayers[layer + 2]
                    : std::vector<std::pair<qc::Qubit, qc::Qubit>>{},
                false);
        // then compute the next qubit and gate placement with reusing atoms
        if (!reuseQubits[layer].empty()) {
          std::vector<
              std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>
              qubitPlacementWithReuse;
          if (dynamicPlacement_) {
            qubitPlacementWithReuse = placeQubitsInStorageZone(
                placement.front(), placement.back(), reuseQubits[layer],
                twoQubitGateLayers.size() > layer + 1
                    ? twoQubitGateLayers[layer + 1]
                    : std::vector<std::pair<qc::Qubit, qc::Qubit>>{},
                true);
          } else {
            // keep the initial mapping for static placement
            qubitPlacementWithReuse = placement.front();
            for (const auto q : reuseQubits[layer]) {
              qubitPlacementWithReuse[q] = placement.back()[q];
            }
          }
          const std::vector<
              std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>&
              gatePlacementWithReuse = placeGatesInEntanglementZone(
                  qubitPlacementWithoutReuse, reuseQubits[layer],
                  twoQubitGateLayers[layer + 1],
                  twoQubitGateLayers.size() > layer + 2
                      ? twoQubitGateLayers[layer + 2]
                      : std::vector<std::pair<qc::Qubit, qc::Qubit>>{},
                  true);
          // keep the mapping with shorter distance
          const auto& [gatePlacement, qubitPlacement] = filterMapping(
              placement.back(),
              std::pair{qubitPlacementWithoutReuse, gatePlacementWithoutReuse},
              std::pair{qubitPlacementWithReuse, gatePlacementWithReuse});
          placement.emplace_back(gatePlacement);
          placement.emplace_back(qubitPlacement);
        } else {
          placement.emplace_back(qubitPlacementWithoutReuse);
          placement.emplace_back(gatePlacementWithoutReuse);
        }
      } else {
        placement.emplace_back(qubitPlacementWithoutReuse);
      }
    }
    return placement;
  }
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
      const std::vector<std::pair<qc::Qubit, qc::Qubit>>& twoQubitGates,
      const std::vector<std::pair<qc::Qubit, qc::Qubit>>& nextTwoQubitGates,
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
      const std::vector<std::pair<qc::Qubit, qc::Qubit>>& nextTwoQubitGates,
      bool reuse) const
      -> std::vector<
          std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>;
};
} // namespace na
