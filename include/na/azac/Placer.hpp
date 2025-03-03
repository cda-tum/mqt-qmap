#pragma once

#include "na/azac/Architecture.hpp"
#include "na/azac/Utils.hpp"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace na {

/// class to find a qubit layout
template <typename T> class Placer {
protected:
  /// generate qubit initial layout
  auto placeQubitInitial() -> void {
    const auto t_p = std::chrono::system_clock::now();
    if (static_cast<T*>(this)->getGivenInitialMapping()) {
      static_cast<T*>(this)->getQubitMapping().emplace_back(
          *static_cast<T*>(this)->getGivenInitialMapping());
    } else {
      if (static_cast<T*>(this)->isTrivialPlacement()) {
        placeTrivial();
      } else {
        throw std::invalid_argument(
            "Initial placement via simulated annealing is not implemented");
      }
    }
    static_cast<T*>(this)->getRuntimeAnalysis().initialPlacement =
        std::chrono::system_clock::now() - t_p;
  }

  /// generate qubit initial layout
  auto placeQubitIntermediate() -> void {
    const auto t_p = std::chrono::system_clock::now();
    VertexMatchingPlacer intermediatePlacer(
        static_cast<T*>(this)->getArchitecture());
    intermediatePlacer.run(static_cast<T*>(this)->getQubitMapping().front(),
                           static_cast<T*>(this)->getGateScheduling(),
                           static_cast<T*>(this)->isDynamicPlacement(),
                           static_cast<T*>(this)->getReuseQubits());
    static_cast<T*>(this)->getQubitMapping() = intermediatePlacer.getMapping();

    static_cast<T*>(this)->getRuntimeAnalysis().intermediatePlacement =
        std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::system_clock::now() - t_p);
  }

private:
  auto placeTrivial() -> void {
    auto slmIt = static_cast<T*>(this)->getArchitecture().storageZones.cbegin();
    const SLM* slm = slmIt->get();
    // decide whether to begin with row 0 or row n
    const double dis1 = static_cast<T*>(this)
                            ->getArchitecture()
                            .nearestEntanglementSiteDistance(slm, 0, 0);
    const double dis2 =
        static_cast<T*>(this)
            ->getArchitecture()
            .nearestEntanglementSiteDistance(slm, slm->nRows - 1, 0);
    std::size_t c = 0;
    std::int64_t r =
        dis1 < dis2 ? 0 : static_cast<std::int64_t>(slm->nRows) - 1;
    const std::int64_t step = dis1 < dis2 ? 1 : -1;
    std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>
        listPossiblePosition{};
    for (std::size_t i = 0; i < static_cast<T*>(this)->getNQubits(); ++i) {
      listPossiblePosition.emplace_back(slm, r, c);
      ++c;
      if (c % slm->nCols == 0) {
        r += step;
        c = 0;
        if (r == static_cast<std::int64_t>(slm->nRows)) {
          ++slmIt;
          slm = slmIt->get();
          if (step > 0) {
            r = static_cast<std::int64_t>(slm->nRows) - 1;
          } else {
            r = 0;
          }
        }
      }
    }
    static_cast<T*>(this)->getQubitMapping().emplace_back(listPossiblePosition);
  }

  /// class to find a qubit placement via vertex matching
  class VertexMatchingPlacer {
  private:
    const Architecture& architecture;
    std::vector<std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>
        mapping{};
    bool l2 = false;
    double costAtomTransfer = 0.9999;
    std::size_t nQubit = 0;
    std::vector<std::unordered_set<std::size_t>> listReuseQubits{};

  public:
    VertexMatchingPlacer(const Architecture& architecture,
                         const bool l2 = false)
        : architecture(architecture), l2(l2) {}

    [[nodiscard]] auto getMapping() const -> const std::vector<
        std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>& {
      return mapping;
    }

    auto
    run(const std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>&
            initialMapping,
        const std::vector<std::vector<const std::pair<qc::Qubit, qc::Qubit>*>>&
            listGate,
        const bool dynamicPlacement,
        const std::vector<std::unordered_set<std::size_t>>& reuseQubits)
        -> void {
      listReuseQubits = reuseQubits;
      nQubit = initialMapping.size();
      mapping.clear();
      mapping.emplace_back(initialMapping);
      std::cout << "[INFO] AZAC: Minimum-weight-full-matching-based "
                   "intermediate placement: Start\n";
      placeGate(initialMapping, listGate, false);
      for (std::size_t layer = 0; layer < listGate.size(); ++layer) {
        if (dynamicPlacement) {
          placeQubit(listGate, layer, false);
        } else {
          mapping.emplace_back(
              mapping[0]); // keep the initial mapping for static placement
        }
        if (layer + 1 < listGate.size()) {
          placeGate(mapping[mapping.size() - 2], mapping.back(), listGate,
                    layer + 1, false);
        }
        if (!reuseQubits[layer].empty()) {
          if (dynamicPlacement) {
            placeQubit(listGate, layer, true);
          } else {
            mapping.emplace_back(
                mapping[0]); // keep the initial mapping for static placement
            for (const auto q : reuseQubits[layer]) {
              // todo: check if this is correct
              mapping.back()[q] = mapping[mapping.size() - 4][q];
            }
          }
          if (layer + 1 < listGate.size()) {
            placeGate(mapping[mapping.size() - 4], mapping.back(), listGate,
                      layer + 1, true);
            // keep the mapping with shorter distance
            filterMapping(layer);
          }
        }
      }
      std::cout << "[INFO] AZAC: Minimum-weight-full-matching-based "
                   "intermediate placement: Finish\n";
    }

    auto filterMapping(const std::size_t layer) -> void {
      const auto& lastGateMapping = mapping[mapping.size() - 5];
      auto& qubitMapping = mapping[mapping.size() - 4];
      auto& gateMapping = mapping[mapping.size() - 3];
      double costNoReuse = 0;
      std::unordered_map<
          std::tuple<const SLM*, std::size_t, const SLM*, std::size_t>, double>
          movementParallelMovement1;
      std::unordered_map<
          std::tuple<const SLM*, std::size_t, const SLM*, std::size_t>, double>
          movementParallelMovement2;

      for (std::size_t q = 0; q < lastGateMapping.size(); ++q) {
        if (lastGateMapping[q] != gateMapping[q]) {
          const auto* slm1 = std::get<0>(lastGateMapping[q]);
          if (slm1->isEntanglement()) {
            slm1 = slm1->entanglementZone->front().get();
          }
          const auto* slm2 = std::get<0>(qubitMapping[q]);
          if (slm2->isEntanglement()) {
            slm2 = slm2->entanglementZone->front().get();
          }
          const std::tuple key{slm1, std::get<1>(lastGateMapping[q]), slm2,
                               std::get<1>(qubitMapping[q])};
          const double dis = architecture.distance(
              std::get<0>(lastGateMapping[q]), std::get<1>(lastGateMapping[q]),
              std::get<2>(lastGateMapping[q]), std::get<0>(qubitMapping[q]),
              std::get<1>(qubitMapping[q]), std::get<2>(qubitMapping[q]));
          if (const auto it = movementParallelMovement1.find(key);
              it != movementParallelMovement1.end()) {
            movementParallelMovement1[key] = std::max(it->second, dis);
          } else {
            movementParallelMovement1[key] = dis;
          }
        }
        if (qubitMapping[q] != gateMapping[q]) {
          const auto* slm1 = std::get<0>(gateMapping[q]);
          if (slm1->isEntanglement()) {
            slm1 = slm1->entanglementZone->front().get();
          }
          const auto* slm2 = std::get<0>(qubitMapping[q]);
          if (slm2->isEntanglement()) {
            slm2 = slm2->entanglementZone->front().get();
          }
          const std::tuple key{slm2, std::get<1>(qubitMapping[q]), slm1,
                               std::get<1>(gateMapping[q])};
          const double dis = architecture.distance(
              std::get<0>(qubitMapping[q]), std::get<1>(qubitMapping[q]),
              std::get<2>(qubitMapping[q]), std::get<0>(gateMapping[q]),
              std::get<1>(gateMapping[q]), std::get<2>(gateMapping[q]));
          if (const auto it = movementParallelMovement2.find(key);
              it != movementParallelMovement2.end()) {
            movementParallelMovement2[key] = std::max(it->second, dis);
          } else {
            movementParallelMovement2[key] = dis;
          }
        }
      }
      for (const auto& [_, value] : movementParallelMovement1) {
        costNoReuse += std::sqrt(value);
      }
      for (const auto& [_, value] : movementParallelMovement2) {
        costNoReuse += std::sqrt(value);
      }
      // now calculate cost with reuse
      gateMapping = mapping[mapping.size() - 1];
      qubitMapping = mapping[mapping.size() - 2];
      double costReuse = 0;
      movementParallelMovement1.clear();
      movementParallelMovement2.clear();
      for (std::size_t q = 0; q < lastGateMapping.size(); ++q) {
        if (lastGateMapping[q] != gateMapping[q]) {
          const auto* slm1 = std::get<0>(lastGateMapping[q]);
          if (slm1->isEntanglement()) {
            slm1 = slm1->entanglementZone->front().get();
          }
          const auto* slm2 = std::get<0>(qubitMapping[q]);
          if (slm2->isEntanglement()) {
            slm2 = slm2->entanglementZone->front().get();
          }
          const std::tuple key{slm1, std::get<1>(lastGateMapping[q]), slm2,
                               std::get<1>(qubitMapping[q])};
          double dis = architecture.distance(
              std::get<0>(lastGateMapping[q]), std::get<1>(lastGateMapping[q]),
              std::get<2>(lastGateMapping[q]), std::get<0>(qubitMapping[q]),
              std::get<1>(qubitMapping[q]), std::get<2>(qubitMapping[q]));
          if (const auto it = movementParallelMovement1.find(key);
              it != movementParallelMovement1.end()) {
            movementParallelMovement1[key] =
                std::max(movementParallelMovement1[key], dis);
          } else {
            movementParallelMovement1[key] = dis;
          }
        }
        if (qubitMapping[q] != gateMapping[q]) {
          const auto* slm1 = std::get<0>(gateMapping[q]);
          if (slm1->isEntanglement()) {
            slm1 = slm1->entanglementZone->front().get();
          }
          const auto* slm2 = std::get<0>(qubitMapping[q]);
          if (slm2->isEntanglement()) {
            slm2 = slm2->entanglementZone->front().get();
          }
          const std::tuple key{slm2, std::get<1>(qubitMapping[q]), slm1,
                               std::get<1>(gateMapping[q])};
          const double dis = architecture.distance(
              std::get<0>(qubitMapping[q]), std::get<1>(qubitMapping[q]),
              std::get<2>(qubitMapping[q]), std::get<0>(gateMapping[q]),
              std::get<1>(gateMapping[q]), std::get<2>(gateMapping[q]));
          if (const auto it = movementParallelMovement2.find(key);
              it != movementParallelMovement2.end()) {
            movementParallelMovement2[key] =
                std::max(movementParallelMovement2[key], dis);
          } else {
            movementParallelMovement2[key] = dis;
          }
        }
      }
      for (const auto& [_, value] : movementParallelMovement1) {
        costReuse += std::sqrt(value);
      }
      for (const auto& [_, value] : movementParallelMovement2) {
        costReuse += std::sqrt(value);
      }
      if (costAtomTransfer * pow((1 - costNoReuse / 1.5e6), nQubit) >
          pow((1 - costReuse / 1.5e6), nQubit)) {
        listReuseQubits[layer] = {};
        mapping.pop_back();
        mapping.pop_back();
      } else {
        mapping.erase(mapping.end() - 3);
        mapping.erase(mapping.end() - 3);
      }
    }
    /// generate gate mapping based on minimum weight matching for the first
    /// layer of gates
    auto placeGate(
        const std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>&
            qubitMapping,
        const std::vector<std::vector<const std::pair<qc::Qubit, qc::Qubit>*>>&
            listTwoGateLayer,
        const bool testReuse) -> void {
      const auto& listGate = listTwoGateLayer.front();
      std::unordered_map<std::size_t, std::size_t> dictReuseQubitNeighbor;
      if (listTwoGateLayer.size() > 1 and testReuse) {
        for (const auto q : listReuseQubits.front()) {
          for (const auto gate : listTwoGateLayer[1]) {
            if (q == gate->first) {
              dictReuseQubitNeighbor[q] = gate->second;
              break;
            }
            if (q == gate->second) {
              dictReuseQubitNeighbor[q] = gate->first;
              break;
            }
          }
        }
      }
      std::unordered_map<std::tuple<const SLM*, std::size_t, std::size_t>,
                         std::size_t>
          siteRydbergToIdx;
      std::vector<std::tuple<const SLM*, std::size_t, std::size_t>> listRydberg;
      // WARNING: The role of listColCoo and listRowCoo is swapped compared to
      // the original implementation because our matching algorithm only
      // supports a single direction, i.e., fewer rows than columns.
      std::vector<std::size_t> listRowCoo;
      std::vector<std::size_t> listColCoo;
      std::vector<double> listData;
      const auto expandFactor =
          static_cast<std::size_t>(std::ceil(std::sqrt(listGate.size() / 2)));
      for (std::size_t i = 0; i < listGate.size(); ++i) {
        const auto& [q1, q2] = *listGate[i];
        std::unordered_set<std::tuple<const SLM*, std::size_t, std::size_t>>
            setNearbySite;
        const auto* slm = std::get<0>(qubitMapping[q1]);
        std::unordered_set<std::tuple<const SLM*, std::size_t, std::size_t>>
            nearestSites{};
        nearestSites.emplace(architecture.nearestEntanglementSite(
            std::get<0>(qubitMapping[q1]), std::get<1>(qubitMapping[q1]),
            std::get<2>(qubitMapping[q1]), std::get<0>(qubitMapping[q2]),
            std::get<1>(qubitMapping[q2]), std::get<2>(qubitMapping[q2])));
        nearestSites.emplace(architecture.nearestEntanglementSite(
            std::get<0>(qubitMapping[q1]), 0, std::get<2>(qubitMapping[q1]),
            std::get<0>(qubitMapping[q2]), 0, std::get<2>(qubitMapping[q2])));
        nearestSites.emplace(architecture.nearestEntanglementSite(
            std::get<0>(qubitMapping[q1]), slm->nRows - 1,
            std::get<2>(qubitMapping[q1]), std::get<0>(qubitMapping[q2]),
            slm->nRows - 1, std::get<2>(qubitMapping[q2])));
        for (const auto& nearestSite : nearestSites) {
          setNearbySite.emplace(nearestSite);
          const auto& [slm_idx, slm_r, slm_c] = nearestSite;
          auto low_r = slm_r > expandFactor ? slm_r - expandFactor : 0;
          auto high_r = std::min(slm->nRows, slm_r + expandFactor + 1);
          auto low_c = slm_c > expandFactor ? slm_c - expandFactor : 0;
          auto high_c = std::min(slm->nCols, slm_c + expandFactor + 1);
          if (high_c - low_c < 2 * expandFactor) {
            const auto heightGap = static_cast<std::size_t>(
                std::ceil(listGate.size() / (high_c - low_c)) - expandFactor);
            low_r = low_r > (heightGap / 2) ? low_r - (heightGap / 2) : 0;
            high_r = std::min(slm->nRows, low_r + heightGap + expandFactor);
          }
          if (high_r - low_r < 2 * expandFactor) {
            const auto widthGap = static_cast<std::size_t>(
                std::ceil(listGate.size() / (high_r - low_r)) - expandFactor);
            low_c = low_c > (widthGap / 2) ? low_c - (widthGap / 2) : 0;
            high_c = std::min(slm->nCols, low_c + widthGap + expandFactor);
          }
          for (std::size_t r = low_r; r < high_r; ++r) {
            for (std::size_t c = low_c; c < high_c; ++c) {
              setNearbySite.emplace(slm_idx, r, c);
            }
          }
        }
        for (const auto& site : setNearbySite) {
          if (siteRydbergToIdx.find(site) == siteRydbergToIdx.end()) {
            siteRydbergToIdx[site] = listRydberg.size();
            listRydberg.emplace_back(site);
          }
          const auto idxRydberg = siteRydbergToIdx[site];
          double dis1 = architecture.distance(
              std::get<0>(qubitMapping[q1]), std::get<1>(qubitMapping[q1]),
              std::get<2>(qubitMapping[q1]), std::get<0>(site),
              std::get<1>(site), std::get<2>(site));
          double dis2 = architecture.distance(
              std::get<0>(qubitMapping[q2]), std::get<1>(qubitMapping[q2]),
              std::get<2>(qubitMapping[q2]), std::get<0>(site),
              std::get<1>(site), std::get<2>(site));
          double dis3 = 0;
          std::optional<std::size_t> q3;
          if (dictReuseQubitNeighbor.find(q1) != dictReuseQubitNeighbor.end()) {
            q3 = dictReuseQubitNeighbor[q1];
          } else if (dictReuseQubitNeighbor.find(q2) !=
                     dictReuseQubitNeighbor.end()) {
            q3 = dictReuseQubitNeighbor[q2];
          }
          if (q3) {
            dis3 = architecture.distance(
                std::get<0>(qubitMapping[*q3]), std::get<1>(qubitMapping[*q3]),
                std::get<2>(qubitMapping[*q3]), std::get<0>(site),
                std::get<1>(site), std::get<2>(site));
          }
          listColCoo.emplace_back(idxRydberg);
          listRowCoo.emplace_back(i);
          if (std::get<1>(qubitMapping[q1]) == std::get<1>(qubitMapping[q2]) &&
              std::get<0>(qubitMapping[q1]) == std::get<0>(qubitMapping[q2])) {
            listData.emplace_back(std::sqrt(std::max(dis1, dis2)) +
                                  std::sqrt(dis3));
          } else {
            listData.emplace_back(std::sqrt(dis1) + std::sqrt(dis2) +
                                  std::sqrt(dis3));
          }
        }
      }

      if (listRydberg.size() < listGate.size()) {
        std::stringstream ss;
        ss << "[Error] AZAC: Minimum-weight-full-matching-based intermediate "
              "placement: No enough sites for gates ("
           << listRydberg.size() << " vs " << listGate.size() << ").";
        throw std::invalid_argument(ss.str());
      }

      std::vector matrix(
          listGate.size(),
          std::vector<std::optional<double>>(listRydberg.size(), std::nullopt));
      for (std::size_t i = 0; i < listRowCoo.size(); ++i) {
        matrix[listRowCoo[i]][listColCoo[i]] = listData[i];
      }
      const auto& matching = minimumWeightFullBipartiteMatching(matrix);
      double cost = 0;
      for (std::size_t i = 0; i < matching.size(); ++i) {
        cost = cost + matrix[i][matching[i]].value();
      }
      auto tmpMapping = qubitMapping;
      for (std::size_t idxGate = 0; idxGate < matching.size(); ++idxGate) {
        const auto q0 = listGate[idxGate]->first;
        const auto q1 = listGate[idxGate]->second;
        const auto& site = listRydberg[matching[idxGate]];
        if (std::get<2>(qubitMapping[q0]) < std::get<2>(qubitMapping[q1])) {
          tmpMapping[q0] = site;
          tmpMapping[q1] = {std::get<0>(site)->entanglementZone->back().get(),
                            std::get<1>(site), std::get<2>(site)};
        } else {
          tmpMapping[q0] = {std::get<0>(site)->entanglementZone->back().get(),
                            std::get<1>(site), std::get<2>(site)};
          tmpMapping[q1] = site;
        }
      }
      mapping.emplace_back(tmpMapping);
    }
    /// generate gate mapping based on minimum weight matching for all layers
    /// of gates except the first one
    auto placeGate(
        const std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>&
            gateMapping,
        const std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>&
            qubitMapping,
        const std::vector<std::vector<const std::pair<qc::Qubit, qc::Qubit>*>>&
            listTwoGateLayer,
        const std::size_t layer, const bool testReuse) -> void {
      const auto& listGate = listTwoGateLayer[layer];
      std::unordered_map<std::size_t, std::size_t> dictReuseQubitNeighbor;
      if (listTwoGateLayer.size() > layer + 1 and testReuse) {
        for (const auto q : listReuseQubits[layer]) {
          for (const auto gate : listTwoGateLayer[layer + 1]) {
            if (q == gate->first) {
              dictReuseQubitNeighbor[q] = gate->second;
              break;
            }
            if (q == gate->second) {
              dictReuseQubitNeighbor[q] = gate->first;
              break;
            }
          }
        }
      }
      std::unordered_map<std::tuple<const SLM*, std::size_t, std::size_t>,
                         std::size_t>
          siteRydbergToIdx;
      std::vector<std::tuple<const SLM*, std::size_t, std::size_t>> listRydberg;
      // WARNING: The role of listColCoo and listRowCoo is swapped compared to
      // the original implementation because our matching algorithm only
      // supports a single direction, i.e., fewer rows than columns.
      std::vector<std::size_t> listRowCoo;
      std::vector<std::size_t> listColCoo;
      std::vector<double> listData;
      const auto expandFactor =
          static_cast<std::size_t>(std::ceil(std::sqrt(listGate.size() / 2)));
      for (std::size_t i = 0; i < listGate.size(); ++i) {
        const auto& [q1, q2] = *listGate[i];
        std::unordered_set<std::tuple<const SLM*, std::size_t, std::size_t>>
            setNearbySite;

        if (testReuse && (listReuseQubits[layer - 1].find(q1) !=
                          listReuseQubits[layer - 1].end())) {
          const auto location = gateMapping[q1];
          const auto slmIdx =
              std::get<0>(location)->entanglementZone->front().get();
          setNearbySite.emplace(slmIdx, std::get<1>(location),
                                std::get<2>(location));
        } else if (testReuse && (listReuseQubits[layer - 1].find(q2) !=
                                 listReuseQubits[layer - 1].end())) {
          const auto location = gateMapping[q2];
          const auto slmIdx =
              std::get<0>(location)->entanglementZone->front().get();
          setNearbySite.emplace(slmIdx, std::get<1>(location),
                                std::get<2>(location));
        } else {
          const auto* slm = std::get<0>(qubitMapping[q1]);
          std::unordered_set<std::tuple<const SLM*, std::size_t, std::size_t>>
              nearestSites{};
          nearestSites.emplace(architecture.nearestEntanglementSite(
              std::get<0>(qubitMapping[q1]), std::get<1>(qubitMapping[q1]),
              std::get<2>(qubitMapping[q1]), std::get<0>(qubitMapping[q2]),
              std::get<1>(qubitMapping[q2]), std::get<2>(qubitMapping[q2])));
          nearestSites.emplace(architecture.nearestEntanglementSite(
              std::get<0>(qubitMapping[q1]), 0, std::get<2>(qubitMapping[q1]),
              std::get<0>(qubitMapping[q2]), 0, std::get<2>(qubitMapping[q2])));
          nearestSites.emplace(architecture.nearestEntanglementSite(
              std::get<0>(qubitMapping[q1]), slm->nRows - 1,
              std::get<2>(qubitMapping[q1]), std::get<0>(qubitMapping[q2]),
              slm->nRows - 1, std::get<2>(qubitMapping[q2])));
          for (const auto& nearestSite : nearestSites) {
            setNearbySite.emplace(nearestSite);
            const auto& [slm_idx, slm_r, slm_c] = nearestSite;
            auto low_r = slm_r > expandFactor ? slm_r - expandFactor : 0;
            auto high_r = std::min(slm->nRows, slm_r + expandFactor + 1);
            auto low_c = slm_c > expandFactor ? slm_c - expandFactor : 0;
            auto high_c = std::min(slm->nCols, slm_c + expandFactor + 1);
            if (high_c - low_c < 2 * expandFactor) {
              const auto heightGap = static_cast<std::size_t>(
                  std::ceil(listGate.size() / (high_c - low_c)) - expandFactor);
              low_r = low_r > (heightGap / 2) ? low_r - (heightGap / 2) : 0;
              high_r = std::min(slm->nRows, low_r + heightGap + expandFactor);
            }
            if (high_r - low_r < 2 * expandFactor) {
              const auto widthGap = static_cast<std::size_t>(
                  std::ceil(listGate.size() / (high_r - low_r)) - expandFactor);
              low_c = low_c > (widthGap / 2) ? low_c - (widthGap / 2) : 0;
              high_c = std::min(slm->nCols, low_c + widthGap + expandFactor);
            }
            for (std::size_t r = low_r; r < high_r; ++r) {
              for (std::size_t c = low_c; c < high_c; ++c) {
                setNearbySite.emplace(slm_idx, r, c);
              }
            }
          }
        }
        for (const auto& site : setNearbySite) {
          if (siteRydbergToIdx.find(site) == siteRydbergToIdx.end()) {
            siteRydbergToIdx[site] = listRydberg.size();
            listRydberg.emplace_back(site);
          }
          const auto idxRydberg = siteRydbergToIdx[site];
          double dis1 = architecture.distance(
              std::get<0>(qubitMapping[q1]), std::get<1>(qubitMapping[q1]),
              std::get<2>(qubitMapping[q1]), std::get<0>(site),
              std::get<1>(site), std::get<2>(site));
          double dis2 = architecture.distance(
              std::get<0>(qubitMapping[q2]), std::get<1>(qubitMapping[q2]),
              std::get<2>(qubitMapping[q2]), std::get<0>(site),
              std::get<1>(site), std::get<2>(site));
          double dis3 = 0;
          std::optional<std::size_t> q3;
          if (dictReuseQubitNeighbor.find(q1) != dictReuseQubitNeighbor.end()) {
            q3 = dictReuseQubitNeighbor[q1];
          } else if (dictReuseQubitNeighbor.find(q2) !=
                     dictReuseQubitNeighbor.end()) {
            q3 = dictReuseQubitNeighbor[q2];
          }
          if (q3) {
            dis3 = architecture.distance(
                std::get<0>(qubitMapping[*q3]), std::get<1>(qubitMapping[*q3]),
                std::get<2>(qubitMapping[*q3]), std::get<0>(site),
                std::get<1>(site), std::get<2>(site));
          }
          listColCoo.emplace_back(idxRydberg);
          listRowCoo.emplace_back(i);
          // if row and slm of both qubits coinicde, then the distance is the
          // maximum of the two distances
          if (std::get<1>(qubitMapping[q1]) == std::get<1>(qubitMapping[q2]) &&
              std::get<0>(qubitMapping[q1]) == std::get<0>(qubitMapping[q2])) {
            listData.emplace_back(std::sqrt(std::max(dis1, dis2)) +
                                  std::sqrt(dis3));
          }
          // otherwise, the distance is the sum
          else {
            listData.emplace_back(std::sqrt(dis1) + std::sqrt(dis2) +
                                  std::sqrt(dis3));
          }
        }
      }

      if (listRydberg.size() < listGate.size()) {
        std::stringstream ss;
        ss << "[Error] AZAC: Minimum-weight-full-matching-based intermediate "
              "placement: No enough sites for gates ("
           << listRydberg.size() << " vs " << listGate.size() << ").";
        throw std::invalid_argument(ss.str());
      }

      std::vector matrix(
          listGate.size(),
          std::vector<std::optional<double>>(listRydberg.size(), std::nullopt));
      for (std::size_t i = 0; i < listRowCoo.size(); ++i) {
        matrix[listRowCoo[i]][listColCoo[i]] = listData[i];
      }
      const auto& matching = minimumWeightFullBipartiteMatching(matrix);
      double cost = 0;
      for (std::size_t i = 0; i < matching.size(); ++i) {
        cost = cost + matrix[i][matching[i]].value();
      }
      auto tmpMapping = qubitMapping;
      for (std::size_t idxGate = 0; idxGate < matching.size(); ++idxGate) {
        const auto q0 = listGate[idxGate]->first;
        const auto q1 = listGate[idxGate]->second;
        const auto& site = listRydberg[matching[idxGate]];
        if (testReuse && (listReuseQubits[layer - 1].find(q0) !=
                          listReuseQubits[layer - 1].end())) {
          tmpMapping[q0] = gateMapping[q0];
          if (site == gateMapping[q0]) {
            tmpMapping[q1] = {std::get<0>(site)->entanglementZone->back().get(),
                              std::get<1>(site), std::get<2>(site)};
          } else {
            tmpMapping[q1] = site;
          }
        } else if (testReuse && (listReuseQubits[layer - 1].find(q1) !=
                                 listReuseQubits[layer - 1].end())) {
          tmpMapping[q1] = gateMapping[q1];
          if (site == gateMapping[q1]) {
            tmpMapping[q0] = {std::get<0>(site)->entanglementZone->back().get(),
                              std::get<1>(site), std::get<2>(site)};
          } else {
            tmpMapping[q0] = site;
          }
        } else {
          if (std::get<2>(qubitMapping[q0]) < std::get<2>(qubitMapping[q1])) {
            tmpMapping[q0] = site;
            tmpMapping[q1] = {std::get<0>(site)->entanglementZone->back().get(),
                              std::get<1>(site), std::get<2>(site)};
          } else {
            tmpMapping[q0] = {std::get<0>(site)->entanglementZone->back().get(),
                              std::get<1>(site), std::get<2>(site)};
            tmpMapping[q1] = site;
          }
        }
      }
      mapping.emplace_back(tmpMapping);
    }

    /// Generate qubit mapping based on minimum weight matching.
    /// @param listGate List of 2-qubit gates as pairs of qubits
    /// - @code listGate[layer]@endcode: gates to be executed in the current
    ///   Rydberg stage
    /// - @code listGate[i], i > layer@endcode: is the list of yet unexecuted
    ///   gates
    /// @param layer the current Rydberg stage
    /// @param testReuse whether to test the reuse of qubits or not
    auto placeQubit(
        const std::vector<std::vector<const std::pair<qc::Qubit, qc::Qubit>*>>&
            listGate,
        const std::size_t layer, const bool testReuse) -> void {
      // the very initial placement of qubits
      const auto& qubitMapping = mapping.front();
      // the placement of qubits after the last gate
      const auto& lastGateMapping =
          mapping[mapping.size() - (testReuse ? 3 : 1)];
      // for each storage SLM array, construct a matrix indicating occupancy of
      // sites
      std::unordered_map<const SLM*, std::vector<std::vector<bool>>>
          isEmptyStorageSite{};
      // for each SLM array, initialize the site as empty
      for (const std::unique_ptr<SLM>& slmId : architecture.storageZones) {
        isEmptyStorageSite.emplace(
            slmId.get(),
            std::vector(slmId->nRows, std::vector(slmId->nCols, true)));
      }
      // qubits that need to be placed, in particular, this does not include
      // qubits that are reused NOTE: the indices stored in qubitToPlace are
      // the indices of the lastGateMapping and do not refer to the actual
      // index of a qubit
      std::vector<qc::Qubit> qubitToPlace{};
      // go through the placement of qubits after the last gate
      for (std::size_t q = 0; q < lastGateMapping.size(); ++q) {
        const auto& lastMapping = lastGateMapping[q];
        const auto arrayId = std::get<0>(lastMapping);
        if (isEmptyStorageSite.find(arrayId) != isEmptyStorageSite.end()) {
          // the mapped qubit is in the storage zone, set the site as occupied
          isEmptyStorageSite[arrayId][std::get<1>(lastMapping)]
                            [std::get<2>(lastMapping)] = false;
        } else if (!testReuse || listReuseQubits[layer].find(q) ==
                                     listReuseQubits[layer].end()) {
          // the mapped qubit is in the entangling zone and must be placed (if
          // not reused)
          qubitToPlace.emplace_back(q);
        }
      }
      // occupied sites in the initial mapping that are currently unoccupied;
      // those are used as candidate sites for qubit placement as it is
      // guaranteed that they exist
      std::unordered_set<std::tuple<const SLM*, std::size_t, std::size_t>>
          commonSite{};
      for (const std::tuple<const SLM*, std::size_t, std::size_t>& m :
           mapping[0]) {
        const auto* arrayId = std::get<0>(m);
        if (isEmptyStorageSite.find(arrayId) != isEmptyStorageSite.end()) {
          if (isEmptyStorageSite[arrayId][std::get<1>(m)][std::get<2>(m)]) {
            // add site to common sites if it is currently unoccupied
            commonSite.emplace(arrayId, std::get<1>(m), std::get<2>(m));
          }
        } else {
          // site is in an entangling zone and not yet added to
          // isEmptyStorageSite
          isEmptyStorageSite[arrayId] =
              std::vector(arrayId->nRows, std::vector(arrayId->nCols, true));
          commonSite.emplace(arrayId, std::get<1>(m), std::get<2>(m));
        }
      }
      // dictionary to store the qubit interactions with other qubits
      std::unordered_map<qc::Qubit, std::vector<qc::Qubit>>
          dictQubitInteraction{};
      for (const auto q : qubitToPlace) {
        dictQubitInteraction.emplace(q, std::vector<qc::Qubit>{});
      }
      if (listGate.size() > layer + 1) {
        for (const auto* gate : listGate[layer + 1]) {
          if (dictQubitInteraction.find(gate->first) !=
                  dictQubitInteraction.end() &&
              ((!testReuse) || (listReuseQubits[layer].find(gate->second) ==
                                listReuseQubits[layer].end()))) {
            dictQubitInteraction[gate->first].emplace_back(gate->second);
          }
          if (dictQubitInteraction.find(gate->second) !=
                  dictQubitInteraction.end() &&
              ((!testReuse) || (listReuseQubits[layer].find(gate->first) ==
                                listReuseQubits[layer].end()))) {
            dictQubitInteraction[gate->second].emplace_back(gate->first);
          }
        }
      }
      const std::size_t expandFactor = 1;

      std::unordered_map<std::tuple<const SLM*, std::size_t, std::size_t>,
                         std::size_t>
          siteStorageToIdx{};
      std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>
          listStorage{};
      std::vector<std::size_t> listColCoo{};
      std::vector<std::size_t> listRowCoo{};
      std::vector<double> listData{};

      for (std::size_t i = 0; i < qubitToPlace.size(); ++i) {
        const auto q = qubitToPlace[i];
        std::unordered_map<const SLM*, std::tuple<std::size_t, std::size_t,
                                                  std::size_t, std::size_t>>
            dictBoudingbox{};
        const auto* slm = std::get<0>(qubitMapping[q]);
        auto lowerRow = std::get<1>(qubitMapping[q]);
        auto upperRow = lowerRow;
        auto leftCol = std::get<2>(qubitMapping[q]);
        auto rightCol = leftCol;
        const std::pair<size_t, size_t>& exactLocQ = std::apply(
            [&](auto&&... args) {
              return architecture.exactSlmLocation(
                  std::forward<decltype(args)>(args)...);
            },
            qubitMapping[q]);
        const std::pair<size_t, size_t>& exactLocGate = std::apply(
            [&](auto&&... args) {
              return architecture.exactSlmLocation(
                  std::forward<decltype(args)>(args)...);
            },
            mapping[0][q]);
        if (exactLocGate.second < exactLocQ.second) {
          lowerRow = 0;
        } else {
          upperRow = slm->nRows;
        }

        dictBoudingbox.emplace(
            std::get<0>(qubitMapping[q]),
            std::tuple{lowerRow, upperRow, leftCol, rightCol});
        for (const qc::Qubit neighbor_q : dictQubitInteraction[q]) {
          const SLM* tmpSlmIdx = std::get<0>(lastGateMapping[neighbor_q]);
          const std::tuple<const SLM*, size_t, size_t>& neighborQLocation =
              tmpSlmIdx->isEntanglement()
                  ? std::apply(
                        [&](auto&&... args)
                            -> std::tuple<const SLM*, size_t, size_t> {
                          return architecture.nearestStorageSite(
                              std::forward<decltype(args)>(args)...);
                        },
                        lastGateMapping[neighbor_q])
                  : lastGateMapping[neighbor_q];
          if (const auto& it =
                  dictBoudingbox.find(std::get<0>(neighborQLocation));
              it != dictBoudingbox.end()) {
            dictBoudingbox.insert_or_assign(
                std::get<0>(neighborQLocation),
                std::tuple{std::min(std::get<1>(neighborQLocation),
                                    std::get<0>(it->second)),
                           std::max(std::get<1>(neighborQLocation),
                                    std::get<1>(it->second)),
                           std::min(std::get<2>(neighborQLocation),
                                    std::get<2>(it->second)),
                           std::max(std::get<2>(neighborQLocation),
                                    std::get<3>(it->second))});
          } else {
            const auto* slmId = std::get<0>(neighborQLocation);
            lowerRow = std::get<1>(neighborQLocation);
            upperRow = std::get<1>(neighborQLocation);
            leftCol = std::get<2>(neighborQLocation);
            rightCol = std::get<2>(neighborQLocation);
            const std::pair<std::size_t, std::size_t>& exactLocNeighborQ =
                std::apply(
                    [&](auto&&... args) {
                      return architecture.exactSlmLocation(
                          std::forward<decltype(args)>(args)...);
                    },
                    neighborQLocation);
            if (exactLocGate.second < exactLocNeighborQ.second) {
              lowerRow = 0;
            } else {
              upperRow = slmId->nRows;
            }
            dictBoudingbox.emplace(
                std::get<0>(neighborQLocation),
                std::tuple{lowerRow, upperRow, leftCol, rightCol});
          }
        }
        const auto& gateLocation = lastGateMapping[q];
        const std::tuple<const SLM*, std::size_t, std::size_t>
            nearestStorageSite = std::apply(
                [&](auto&&... args) {
                  return architecture.nearestStorageSite(
                      std::forward<decltype(args)>(args)...);
                },
                gateLocation);
        // todo: what is ratio?
        const std::size_t ratio = 3;
        if (const auto it =
                dictBoudingbox.find(std::get<0>(nearestStorageSite));
            it != dictBoudingbox.end()) {
          dictBoudingbox.insert_or_assign(
              std::get<0>(nearestStorageSite),
              std::tuple{std::min(std::get<1>(nearestStorageSite) - ratio,
                                  std::get<0>(it->second)),
                         std::max(std::get<1>(nearestStorageSite) + ratio,
                                  std::get<1>(it->second)),
                         std::min(std::get<2>(nearestStorageSite) - ratio,
                                  std::get<2>(it->second)),
                         std::max(std::get<2>(nearestStorageSite) + ratio,
                                  std::get<3>(it->second))});
        } else {
          dictBoudingbox.emplace(
              std::get<0>(nearestStorageSite),
              std::tuple{std::get<1>(nearestStorageSite) - ratio,
                         std::get<1>(nearestStorageSite) + ratio,
                         std::get<2>(nearestStorageSite) - ratio,
                         std::get<2>(nearestStorageSite) + ratio});
        }
        auto setNearbySite = commonSite;
        if (isEmptyStorageSite[std::get<0>(qubitMapping[q])][std::get<1>(
                qubitMapping[q])][std::get<2>(qubitMapping[q])]) {
          setNearbySite.emplace(qubitMapping[q]);
        }

        for (const auto& [slmId, _] : dictBoudingbox) {
          auto& boundingBox = dictBoudingbox[slmId];
          boundingBox = {std::get<0>(boundingBox) > expandFactor
                             ? static_cast<std::size_t>(
                                   std::get<0>(boundingBox) - expandFactor)
                             : 0,
                         std::min(std::get<1>(boundingBox) + expandFactor + 1,
                                  slmId->nRows),
                         std::get<2>(boundingBox) > expandFactor
                             ? static_cast<std::size_t>(
                                   std::get<2>(boundingBox) - expandFactor)
                             : 0,
                         std::min(std::get<3>(boundingBox) + expandFactor + 1,
                                  slmId->nCols)};
          for (auto r = std::get<0>(boundingBox); r < std::get<1>(boundingBox);
               ++r) {
            for (auto c = std::get<2>(boundingBox);
                 c < std::get<3>(boundingBox); ++c) {
              if (isEmptyStorageSite.find(slmId) == isEmptyStorageSite.end() ||
                  isEmptyStorageSite[slmId][r][c]) {
                setNearbySite.emplace(slmId, r, c);
              }
            }
          }
        }

        for (const auto& site : setNearbySite) {
          if (siteStorageToIdx.find(site) == siteStorageToIdx.end()) {
            siteStorageToIdx.emplace(site, listStorage.size());
            listStorage.emplace_back(site);
          }
          const auto idxStorage = siteStorageToIdx[site];
          const double dis = architecture.distance(
              std::get<0>(gateLocation), std::get<1>(gateLocation),
              std::get<2>(gateLocation), std::get<0>(site), std::get<1>(site),
              std::get<2>(site));
          double lookaheadCost = 0;
          for (const auto& neighborQ : dictQubitInteraction[q]) {
            const auto& siteNeighborQ = lastGateMapping[neighborQ];
            if (!std::get<0>(siteNeighborQ)->entanglementZone) {
              lookaheadCost += architecture.nearestEntanglementSiteDistance(
                  std::get<0>(site), std::get<1>(site), std::get<2>(site),
                  std::get<0>(siteNeighborQ), std::get<1>(siteNeighborQ),
                  std::get<2>(siteNeighborQ));
            } else {
              const std::pair<std::size_t, std::size_t>& exactLocNeighborQ =
                  std::apply(
                      [&](auto&&... args) {
                        return architecture.exactSlmLocation(
                            std::forward<decltype(args)>(args)...);
                      },
                      lastGateMapping[neighborQ]);
              const auto dx =
                  static_cast<std::int64_t>(exactLocNeighborQ.first) -
                  static_cast<std::int64_t>(exactLocQ.first);
              const auto dy =
                  static_cast<std::int64_t>(exactLocNeighborQ.second) -
                  static_cast<std::int64_t>(exactLocQ.second);
              lookaheadCost += std::sqrt(std::sqrt((dx * dx) + (dy * dy)));
            }
          }
          const double cost = std::sqrt(dis) + (0.1 * lookaheadCost);
          listColCoo.emplace_back(idxStorage);
          listRowCoo.emplace_back(i);
          listData.emplace_back(cost);
        }
      }

      // to construct from three arrays:
      // - listData the entries of the matrix, in any order
      // - listRowCoo the row indices of the matrix entries
      // - listColCoo the column indices of the matrix entries
      // Where A[listRowCoo[k], listColCoo[k]] = listData[k].
      std::vector costMatrix(
          qubitToPlace.size(),
          std::vector<std::optional<double>>(listStorage.size(), std::nullopt));
      for (std::size_t k = 0; k < listData.size(); ++k) {
        costMatrix[listRowCoo[k]][listColCoo[k]] = listData[k];
      }
      const auto& matching = minimumWeightFullBipartiteMatching(costMatrix);
      auto tmpMapping = lastGateMapping;
      for (std::size_t j = 0; j < matching.size(); ++j) {
        tmpMapping[qubitToPlace[j]] = listStorage[matching[j]];
      }
      mapping.emplace_back(tmpMapping);
    }
  };
};

/// class to find a qubit layout
template <typename T> class AStarPlacer {
protected:
  /// generate qubit initial layout
  auto placeQubitInitial() -> void {
    const auto t_p = std::chrono::system_clock::now();
    if (static_cast<T*>(this)->getGivenInitialMapping()) {
      static_cast<T*>(this)->getQubitMapping().emplace_back(
          *static_cast<T*>(this)->getGivenInitialMapping());
    } else {
      if (static_cast<T*>(this)->isTrivialPlacement()) {
        placeTrivial();
      } else {
        throw std::invalid_argument(
            "Initial placement via simulated annealing is not implemented");
      }
    }
    static_cast<T*>(this)->getRuntimeAnalysis().initialPlacement =
        std::chrono::system_clock::now() - t_p;
  }

  /// generate qubit initial layout
  auto placeQubitIntermediate() -> void {
    const auto t_p = std::chrono::system_clock::now();
    /* Architecture =        */ static_cast<T*>(this)->getArchitecture();
    /* Initial Mapping =     */ static_cast<T*>(this)
        ->getQubitMapping()
        .front();
    /* Gate Scheduling =     */ static_cast<T*>(this)->getGateScheduling();
    /* Dynamic Placement =   */ static_cast<T*>(this)->isDynamicPlacement();
    /* Reuse Qubits =        */ static_cast<T*>(this)->getReuseQubits();
    /* Final Qubit Mapping = */ static_cast<T*>(this)->getQubitMapping();

    static_cast<T*>(this)->getRuntimeAnalysis().intermediatePlacement =
        std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::system_clock::now() - t_p);
  }

private:
  auto placeTrivial() -> void {
    auto slmIt = static_cast<T*>(this)->getArchitecture().storageZones.cbegin();
    const SLM* slm = slmIt->get();
    // decide whether to begin with row 0 or row n
    const double dis1 = static_cast<T*>(this)
                            ->getArchitecture()
                            .nearestEntanglementSiteDistance(slm, 0, 0);
    const double dis2 =
        static_cast<T*>(this)
            ->getArchitecture()
            .nearestEntanglementSiteDistance(slm, slm->nRows - 1, 0);
    std::size_t c = 0;
    std::int64_t r =
        dis1 < dis2 ? 0 : static_cast<std::int64_t>(slm->nRows) - 1;
    const std::int64_t step = dis1 < dis2 ? 1 : -1;
    std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>
        listPossiblePosition{};
    for (std::size_t i = 0; i < static_cast<T*>(this)->getNQubits(); ++i) {
      listPossiblePosition.emplace_back(slm, r, c);
      ++c;
      if (c % slm->nCols == 0) {
        r += step;
        c = 0;
        if (r == static_cast<std::int64_t>(slm->nRows)) {
          ++slmIt;
          slm = slmIt->get();
          if (step > 0) {
            r = static_cast<std::int64_t>(slm->nRows) - 1;
          } else {
            r = 0;
          }
        }
      }
    }
    static_cast<T*>(this)->getQubitMapping().emplace_back(listPossiblePosition);
  }

  /// A node representing one stage in the process of placing all atoms
  /// that must be moved for the next stage starting from the last mapping
  /// until a new mapping is found satisfying all constraints of the next
  /// stage
  struct Node {
    /// the maximum distance an already placed atom must travel to its
    /// target location
    double maxDistanceOfPlacedAtom = 0.0;
    /// a binary search tree representing the horizontal groups
    /// @see getNeighbors for more details
    std::vector<std::map<size_t, size_t>> hGroups;
    /// the maximum distance of placed atoms in every group to their
    /// target location
    std::vector<double> maxDistancesOfPlacedAtomsPerHGroup;
    /// @see hGroups
    std::vector<std::map<size_t, size_t>> vGroups;
    /// @see maxDistancesOfPlacedAtomsPerHGroup
    std::vector<double> maxDistancesOfPlacedAtomsPerVGroup;
  };

  std::vector<std::unique_ptr<Node>> nodes;

  /// @brief Returns the cost of a node, i.e., the total cost to reach that node
  /// from the start node.
  /// @details Different groups cannot be rearranged concurrently in one step.
  /// Hence, we add up the time it takes to perform the rearrangements of one
  /// group in one step and sum it up over all groups.
  /// This will not resemble the exact time to rearrange all costs because at
  /// this point it is not clear how the horizontal and vertical groups can be
  /// combined.
  auto getCost(Node current) -> double {
    double cost = 0.0;
    for (const auto d : current->maxDistancesOfPlacedAtomsPerHGroup) {
      cost += std::sqrt(d);
    }
    for (const auto d : current->maxDistancesOfPlacedAtomsPerVGroup) {
      cost += std::sqrt(d);
    }
    return cost;
  }

  /// @brief Return the estimated cost still required to reach a goal node.
  /// @details To yield an optimal results, the heuristic must be admissible,
  /// i.e., never overestimating the cost.
  /// The heuristic returns the estimated costs that are still added to the
  /// current actual cost to reach a goal node.
  /// Hence, the heuristic must always be less or equal to the additional cost
  /// needed to reach a goal.
  /// In the best case, all atoms that are not placed yet are compatible with
  /// an existing group and can just be added to that group.
  /// Hence, the sum in the cost function does not get an additional summand
  /// just the existing summands may increase.
  /// In the case of minimal increase in the overall cost, only one summand
  /// increases its value.
  /// This increase is bounded from below by the maximal distance of an atom to
  /// its nearest potential target site minus the maximum distance already
  /// placed atoms must travel to their determined target site.
  double getHeuristic(Node current) {
    double maxDistanceOfUnplacedAtom = 0.0;
    for (const auto& atom : /* unplaced atoms */) {
      maxDistanceOfUnplacedAtom = std::max(maxDistanceOfUnplacedAtom,
                                           /* distance to nearest site */)
    }
    return maxDistanceOfUnplacedAtom - current->maxDistanceOfPlacedAtom;
  }

  /// @brief Return pointers to all neighbors of the given node.
  /// @details When calling this function, the neighbors are allocated
  /// permanently such that (1) the returned pointers remain valid when the
  /// execution returned from this function and (2) not all nodes in the tree
  /// have to be created before they are needed.
  /// Hence, nodes are only created on demand in this function.
  /// Consequently, this function must only be called once per node.
  /// Otherwise, neighbors for the same node are created twice.
  /// @par
  /// When creating a new node, the horizontal and vertical groups are checked
  /// whether the new corresponding placement is compatible with any of the
  /// existing groups.
  /// If yes, the new placement is added to the respective group and otherwise,
  /// a new group is formed with the new placement.
  auto getNeighbors(const Node* node) -> std::vector<const Node*> {
    size_t atomToBePlacedNext = node->partialPlacement.size();
    std::vector<const Node*> neighbors;
    for (const auto site : /* all possible target sites for the atom */) {
      // assume nodes is of type std::vector<std::unique_ptr<Node>>
      // make a copy of node, the parent of neighbor
      Node* neighbor = nodes.emplace_back(std::make_unique<Node>(*node)).get();
      neighbor->maxDistanceOfPlacedAtom =
            std::max(node->maxDistanceOfPlacedAtom, /* distance for
                current atom from its current site to `site` */);
      // check whether the current placement is compatible with any
      // existing group
      const size_t key = /* the atom's row it is currently in */;
      const size_t value = /* the atom's row it should be placed in */;
      size_t i = 0;
      for (auto& hGroup : neighbor->hGroups) {
        auto it = hGroup.lower_bound(key);
        if (it != hGroup.end()) {
          // an assignment for this key already exists in this group
          const auto& [upperKey, upperValue] = *it;
          if (upperKey == key) {
            if (upperValue == value) {
              // new placement is compatible with this group
              break;
            }
          } else { // if (upperKey > key)
            if (it != hGroup.begin()) {
              // it can be safely decremented
              --it;
              const auto& [_, lowerValue] = *it;
              if (lowerValue < value && value < upperValue) {
                // new placement is compatible with this group
                break;
              }
            } else { // if (it == hGroup.begin())
              if (value < upperValue) {
                // new placement is compatible with this group
                break;
              }
            }
          }
        } else { // if (it == hGroup.end())
          // it can be safely decremented because group must contain
          // at least one element
          --it;
          const auto& [_, lowerValue] = *it;
          if (lowerValue < value) {
            // new placement is compatible with this group
            break;
          }
        }
        ++i;
      }
      if (i == neighbor->hGroups.size()) {
        // no compatible group could be found and a new group is created
        neighbor->hGroups.emplace_back();
        neighbor->maxDistancesOfPlacedAtomsPerHGroup.emplace_back(0.0);
      }
      neighbor->hGroups[i].emplace(key, value);
      neighbor->maxDistancesOfPlacedAtomsPerHGroup[i] =
                std::max(neighbor->maxDistancesOfPlacedAtomsPerHGroup[i],
                         /* distance for current atom from its current site
                            to `site` */);
      // [ do the same for the vertical group... ]
      neighbors.emplace_back(neighbor);
    }
  }
};
} // namespace na
