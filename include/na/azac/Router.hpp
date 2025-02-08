#pragma once

#include "na/azac/CompilerBase.hpp"

#include <algorithm>

namespace na {

template <typename T> class Router {
private:
  /// constant, the distance of AOD row and col to some trap. We use 1µm here.
  static constexpr std::size_t parkingDist = 1;

  std::priority_queue<std::pair<double, std::size_t>,
                      std::vector<std::pair<double, std::size_t>>,
                      std::greater<std::pair<double, std::size_t>>>
      aodEndTime;
  std::vector<std::size_t> aodDependency;
  std::vector<std::size_t> rydbergDependency;
  std::vector<std::size_t> qubitDependency;
  std::unordered_map<std::tuple<const SLM*, std::size_t, std::size_t>,
                     std::size_t>
      siteDependency;

protected:
  /// generate rearrangement layers between two Rydberg layers
  auto routeQubit() -> void {
    for (std::size_t i = 0;
         i < static_cast<T*>(this)->getArchitecture().aods.size(); ++i) {
      aodEndTime.emplace(0, i);
    }
    aodDependency = std::vector<std::size_t>(
        static_cast<T*>(this)->getArchitecture().aods.size(), 0);
    rydbergDependency = std::vector<std::size_t>(
        static_cast<T*>(this)->getArchitecture().entanglementZones.size(), 0);
    std::chrono::microseconds timeMis{};
    qubitDependency =
        std::vector<std::size_t>(static_cast<T*>(this)->getNQubits(), 0);
    writeInitialInstruction();

    for (std::size_t layer = 0;
         layer < static_cast<T*>(this)->getGateScheduling().size(); ++layer) {
      // extract sets of movement that can be performed simultaneously
      const auto t_s = std::chrono::system_clock::now();
      routeQubitMis(layer);
      timeMis += (std::chrono::system_clock::now() - t_s);
      std::cout << "[INFO] AZAC: Solve for Rydberg stage " << (layer + 1) << "/"
                << static_cast<T*>(this)->getGateScheduling().size()
                << ". mis time="
                << std::chrono::duration_cast<std::chrono::microseconds>(
                       timeMis)
                       .count()
                << "\n";
    }
    static_cast<T*>(this)->getRuntimeAnalysis().routing = timeMis;
  }

private:
  /// process layers of movement from storage zone to Rydberg and back to
  /// storage zone
  auto routeQubitMis(const std::size_t layer) -> void {
    const std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>&
        initialMapping = static_cast<T*>(this)->getQubitMapping()[2 * layer];
    const std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>&
        gateMapping = static_cast<T*>(this)->getQubitMapping()[2 * layer + 1];
    const std::optional<
        std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>&
        finalMapping =
            layer + 2 < static_cast<T*>(this)->getQubitMapping().size()
                ? std::make_optional(
                      static_cast<T*>(this)->getQubitMapping()[2 * layer + 2])
                : std::nullopt;

    // sort remainGraph based on qubit distance if using maximal is
    std::vector<std::size_t> remainGraph; // consist qubits to be moved
    for (const std::pair<qc::Qubit, qc::Qubit>* gate :
         static_cast<T*>(this)->getGateScheduling()[layer]) {
      if (initialMapping[gate->first] != gateMapping[gate->first]) {
        remainGraph.emplace_back(gate->first);
      }
      if (initialMapping[gate->second] != gateMapping[gate->second]) {
        remainGraph.emplace_back(gate->second);
      }
    }

    if (static_cast<T*>(this)->getRoutingStrategy() !=
        CompilerBase::RoutingStrategy::MaximalIs) {
      std::sort(remainGraph.begin(), remainGraph.end(),
                [&](const std::size_t a, const std::size_t b) {
                  const auto [aInitX, aInitY] = std::apply(
                      [&](auto&&... args) {
                        return static_cast<T*>(this)
                            ->getArchitecture()
                            .exactSlmLocation(
                                std::forward<decltype(args)>(args)...);
                      },
                      initialMapping[a]);
                  const auto [bInitX, bInitY] = std::apply(
                      [&](auto&&... args) {
                        return static_cast<T*>(this)
                            ->getArchitecture()
                            .exactSlmLocation(
                                std::forward<decltype(args)>(args)...);
                      },
                      initialMapping[b]);
                  const auto [aGateX, aGateY] = std::apply(
                      [&](auto&&... args) {
                        return static_cast<T*>(this)
                            ->getArchitecture()
                            .exactSlmLocation(
                                std::forward<decltype(args)>(args)...);
                      },
                      gateMapping[a]);
                  const auto [bGateX, bGateY] = std::apply(
                      [&](auto&&... args) {
                        return static_cast<T*>(this)
                            ->getArchitecture()
                            .exactSlmLocation(
                                std::forward<decltype(args)>(args)...);
                      },
                      gateMapping[b]);
                  const auto aDx =
                      static_cast<double>(aGateX) - static_cast<double>(aInitX);
                  const auto aDy =
                      static_cast<double>(aGateY) - static_cast<double>(aInitY);
                  const auto bDx =
                      static_cast<double>(bGateX) - static_cast<double>(bInitX);
                  const auto bDy =
                      static_cast<double>(bGateY) - static_cast<double>(bInitY);
                  return aDx * aDx + aDy * aDy > bDx * bDx + bDy * bDy;
                });
    }
    const auto idLayerStart =
        static_cast<T*>(this)->getResult().instructions.size();
    std::size_t batch = 0;
    while (!remainGraph.empty()) {
      // graph construction
      const auto& tuples =
          graphConstruction(remainGraph, initialMapping, gateMapping);
      // collect violation
      const auto& violations = collectViolation(tuples);
      // solve MIS
      const auto& movedQubits = maximalIsSolve(tuples.size(), violations);

      std::unordered_set<std::size_t>
          setAod{}; // use to record aods per movement layer
      setAod.reserve(movedQubits.size());
      for (const auto i : movedQubits) {
        setAod.emplace(remainGraph[i]);
      }
      processMovementLayer(setAod, initialMapping, gateMapping);
      std::vector<std::size_t> tmp;
      for (const auto q : remainGraph) {
        if (setAod.find(q) == setAod.end()) {
          tmp.emplace_back(q);
        }
      }
      remainGraph = tmp;
      batch = batch + 1;
    }

    // append a layer for gate execution
    processGateLayer(layer, gateMapping);
    // move qubit back to the final location
    if (finalMapping) {
      if (static_cast<T*>(this)->isDynamicPlacement() ||
          static_cast<T*>(this)->isReuse()) {
        remainGraph.clear(); // consist qubits to be moved
        for (const std::pair<qc::Qubit, qc::Qubit>* gate :
             static_cast<T*>(this)->getGateScheduling()[layer]) {
          if ((*finalMapping)[gate->first] != gateMapping[gate->first]) {
            remainGraph.emplace_back(gate->first);
          }
          if ((*finalMapping)[gate->second] != gateMapping[gate->second]) {
            remainGraph.emplace_back(gate->second);
          }
        }

        if (static_cast<T*>(this)->getRoutingStrategy() !=
            CompilerBase::RoutingStrategy::MaximalIs) {
          std::sort(remainGraph.begin(), remainGraph.end(),
                    [&](const std::size_t a, const std::size_t b) {
                      const auto [aFinalX, aFinalY] = std::apply(
                          [&](auto&&... args) {
                            return static_cast<T*>(this)
                                ->getArchitecture()
                                .exactSlmLocation(
                                    std::forward<decltype(args)>(args)...);
                          },
                          (*finalMapping)[a]);
                      const auto [bFinalX, bFinalY] = std::apply(
                          [&](auto&&... args) {
                            return static_cast<T*>(this)
                                ->getArchitecture()
                                .exactSlmLocation(
                                    std::forward<decltype(args)>(args)...);
                          },
                          (*finalMapping)[b]);
                      const auto [aGateX, aGateY] = std::apply(
                          [&](auto&&... args) {
                            return static_cast<T*>(this)
                                ->getArchitecture()
                                .exactSlmLocation(
                                    std::forward<decltype(args)>(args)...);
                          },
                          gateMapping[a]);
                      const auto [bGateX, bGateY] = std::apply(
                          [&](auto&&... args) {
                            return static_cast<T*>(this)
                                ->getArchitecture()
                                .exactSlmLocation(
                                    std::forward<decltype(args)>(args)...);
                          },
                          gateMapping[b]);
                      const auto aDx = static_cast<double>(aGateX) -
                                       static_cast<double>(aFinalX);
                      const auto aDy = static_cast<double>(aGateY) -
                                       static_cast<double>(aFinalY);
                      const auto bDx = static_cast<double>(bGateX) -
                                       static_cast<double>(bFinalX);
                      const auto bDy = static_cast<double>(bGateY) -
                                       static_cast<double>(bFinalY);
                      return aDx * aDx + aDy * aDy > bDx * bDx + bDy * bDy;
                    });
        }
        while (!remainGraph.empty()) {
          // graph construction
          const auto& vectors =
              graphConstruction(remainGraph, *finalMapping, gateMapping);
          // collect violation
          const auto& violations = collectViolation(vectors);

          const auto& movedQubits = maximalIsSolve(vectors.size(), violations);
          // todo: add layer
          std::unordered_set<std::size_t>
              setAod; // use to record aods per movement layer
          setAod.reserve(movedQubits.size());
          for (const auto i : movedQubits) {
            setAod.emplace(remainGraph[i]);
          }

          processMovementLayer(setAod, gateMapping, *finalMapping);

          std::vector<std::size_t> tmp;
          for (const auto q : remainGraph) {
            if (setAod.find(q) == setAod.end()) {
              tmp.emplace_back(q);
            }
          }
          remainGraph = tmp;
          ++batch;
        }
      } else {
        // construct reverse layer
        constructReverseLayer(idLayerStart, gateMapping, *finalMapping);
      }
      aodAssignment(idLayerStart);
    }
  }

  auto graphConstruction(
      const std::vector<std::size_t>& remainGraph,
      const std::vector<std::tuple<const SLM*, size_t, size_t>>& initialMapping,
      const std::vector<std::tuple<const SLM*, size_t, size_t>>& finalMapping)
      -> std::vector<
          std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>> {
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>>
        vectors{};
    const std::size_t vectorLength =
        static_cast<T*>(this)->isUseWindow()
            ? std::min(static_cast<T*>(this)->getWindowSize(),
                       remainGraph.size())
            : remainGraph.size();
    for (std::size_t i = 0; i < vectorLength; ++i) {
      const auto q = remainGraph[i];
      const auto& [q_x, q_y] = std::apply(
          [&](auto&&... args) {
            return static_cast<T*>(this)->getArchitecture().exactSlmLocation(
                std::forward<decltype(args)>(args)...);
          },
          initialMapping[q]);
      const auto& [site_x, site_y] = std::apply(
          [&](auto&&... args) {
            return static_cast<T*>(this)->getArchitecture().exactSlmLocation(
                std::forward<decltype(args)>(args)...);
          },
          finalMapping[q]);
      vectors.emplace_back(q_x, site_x, q_y, site_y);
    }
    return vectors;
  }

  auto collectViolation(
      const std::vector<std::tuple<std::size_t, std::size_t, std::size_t,
                                   std::size_t>>& vectors)
      -> std::vector<std::pair<std::size_t, std::size_t>> {
    std::vector<std::pair<std::size_t, std::size_t>> violations{};
    for (std::size_t i = 0; i < vectors.size(); ++i) {
      for (std::size_t j = i + 1; j < vectors.size(); ++j) {
        if (!compatible2D(vectors[i], vectors[j])) {
          violations.emplace_back(i, j);
        }
      }
    }
    return violations;
  }

  /// solve maximal independent set
  auto
  maximalIsSolve(const std::size_t n,
                 const std::vector<std::pair<std::size_t, std::size_t>>& edges)
      -> std::vector<std::size_t> {
    // assume the vertices are sorted based on qubit distance
    std::vector isNodeConflict(n, false);
    std::vector nodeNeighbors(n, std::vector<std::size_t>{});
    for (const auto& edge : edges) {
      nodeNeighbors[edge.first].emplace_back(edge.second);
      nodeNeighbors[edge.second].emplace_back(edge.first);
    }
    std::vector<std::size_t> result{};
    for (std::size_t i = 0; i < isNodeConflict.size(); ++i) {
      if (!isNodeConflict[i]) {
        result.emplace_back(i);
        for (const auto j : nodeNeighbors[i]) {
          isNodeConflict[j] = true;
        }
      }
    }
    return result;
  }

  /// check if move a and b can be performed simultaneously
  auto
  compatible2D(std::tuple<std::size_t, std::size_t, std::size_t, std::size_t> a,
               std::tuple<std::size_t, std::size_t, std::size_t, std::size_t> b)
      -> bool {
    if (std::get<0>(a) == std::get<0>(b) && std::get<1>(a) != std::get<1>(b)) {
      return false;
    }
    if (std::get<1>(a) == std::get<1>(b) && std::get<0>(a) != std::get<0>(b)) {
      return false;
    }
    if (std::get<0>(a) < std::get<0>(b) && std::get<1>(a) >= std::get<1>(b)) {
      return false;
    }
    if (std::get<0>(a) > std::get<0>(b) && std::get<1>(a) <= std::get<1>(b)) {
      return false;
    }
    if (std::get<2>(a) == std::get<2>(b) && std::get<3>(a) != std::get<3>(b)) {
      return false;
    }
    if (std::get<3>(a) == std::get<3>(b) && std::get<2>(a) != std::get<2>(b)) {
      return false;
    }
    if (std::get<2>(a) < std::get<2>(b) && std::get<3>(a) >= std::get<3>(b)) {
      return false;
    }
    if (std::get<2>(a) > std::get<2>(b) && std::get<3>(a) <= std::get<3>(b)) {
      return false;
    }
    return true;
  }

  auto writeInitialInstruction() -> void {
    static_cast<T*>(this)->getResult().instructions.clear();
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>>
        initialLocs{};
    initialLocs.reserve(static_cast<T*>(this)->getNQubits());
    for (std::size_t i = 0; i < static_cast<T*>(this)->getNQubits(); ++i) {
      initialLocs.emplace_back(
          i,
          std::get<0>(static_cast<T*>(this)->getQubitMapping().front()[i])
              ->id,
          std::get<1>(static_cast<T*>(this)->getQubitMapping().front()[i]),
          std::get<2>(static_cast<T*>(this)->getQubitMapping().front()[i]));
    }
    static_cast<T*>(this)->getResult().instructions.emplace_back(
        nlohmann::json{{"type", "init"},
                       {"id", 0},
                       {"begin_time", 0.},
                       {"end_time", 0.},
                       {"init_locs", initialLocs}});

    // process single-qubit gates
    std::unordered_set<std::size_t> setQubitDependency{};
    const std::size_t instIdx =
        static_cast<T*>(this)->getResult().instructions.size();
    const std::vector<qc::StandardOperation>& list1qgate =
        static_cast<T*>(this)->getDictG1QParent().at(nullptr);
    std::vector<nlohmann::json> resultGate{};
    for (const auto& gateInfo : list1qgate) {
      // collect qubit dependency
      const auto qubit = gateInfo.getTargets().front();
      setQubitDependency.emplace(qubitDependency[qubit]);
      qubitDependency[qubit] = instIdx;
      resultGate.emplace_back(
          nlohmann::json{{"name", gateInfo.getName()}, {"q", qubit}});
    }
    nlohmann::json dependency;
    dependency["qubit"] = setQubitDependency;
    if (!resultGate.empty()) {
      write1QGateInstruction(
          instIdx, resultGate, dependency,
          static_cast<T*>(this)->getQubitMapping().front());
      static_cast<T*>(this)->getResult().instructions.back()["begin_time"] = 0.;
      static_cast<T*>(this)->getResult().instructions.back()["end_time"] =
          static_cast<T*>(this)->getArchitecture().time1QGate *
          static_cast<double>(resultGate.size()); // due to sequential execution
    }
  }

  /// generate layers for row-by-row based atom transfer
  auto processMovementLayer(
      const std::unordered_set<std::size_t>& setAodQubit,
      const std::vector<std::tuple<const SLM*, size_t, size_t>>& initialMapping,
      const std::vector<std::tuple<const SLM*, size_t, size_t>>& finalMapping)
      -> void {
    // seperate qubits in listAodqubit into multiple lists where qubits in one
    // list can pick up simultaneously we use row-based pick up
    std::unordered_map<std::size_t, std::vector<std::size_t>>
        pickupDict; // key: row, value: a list of qubit in the same row
    for (const auto q : setAodQubit) {
      const auto& [x, y] = std::apply(
          [&](auto&&... args) {
            return static_cast<T*>(this)->getArchitecture().exactSlmLocation(
                std::forward<decltype(args)>(args)...);
          },
          initialMapping[q]);
      pickupDict.try_emplace(y, std::vector<std::size_t>{});
      pickupDict[y].emplace_back(q);
    }
    std::vector<std::vector<std::size_t>> listAodqubits{};
    std::vector<std::vector<std::tuple<size_t, const SLM*, size_t, size_t>>>
        listendLocation{};
    std::vector<std::vector<std::tuple<size_t, const SLM*, size_t, size_t>>>
        listBeginlocation{};
    nlohmann::json dependency = {{"qubit", std::vector<std::size_t>()},
                                 {"site", std::vector<std::size_t>()}};
    // process aod dependency
    const std::size_t instIdx =
        static_cast<T*>(this)->getResult().instructions.size();

    std::unordered_set<std::size_t> setQubitdependency;
    std::unordered_set<std::size_t> setSiteDependency;
    for (const auto& [dictKey, dictValue] : pickupDict) {
      // collect set of aod qubits to pick up
      listAodqubits.emplace_back(dictValue);
      std::vector<std::tuple<std::size_t, const SLM*, std::size_t, std::size_t>>
          rowBeginLocation{};
      std::vector<std::tuple<std::size_t, const SLM*, std::size_t, std::size_t>>
          rowendLocation{};
      for (const auto q : dictValue) {
        // collect qubit begin location
        rowBeginLocation.emplace_back(q, std::get<0>(initialMapping[q]),
                                      std::get<1>(initialMapping[q]),
                                      std::get<2>(initialMapping[q]));
        // collect qubit end location
        rowendLocation.emplace_back(q, std::get<0>(finalMapping[q]),
                                    std::get<1>(finalMapping[q]),
                                    std::get<2>(finalMapping[q]));
        const auto& siteKey = finalMapping[q];
        if (siteDependency.find(siteKey) != siteDependency.end()) {
          setSiteDependency.emplace(siteDependency[siteKey]);
        }
        siteDependency[initialMapping[q]] = instIdx;

        // collect qubit dependency
        setQubitdependency.emplace(qubitDependency[q]);
        qubitDependency[q] = instIdx;
      }
      listBeginlocation.emplace_back(std::move(rowBeginLocation));
      listendLocation.emplace_back(std::move(rowendLocation));
    }
    dependency["qubit"] = setQubitdependency;
    dependency["site"] = setSiteDependency;
    writeRearrangementInstruction(instIdx, listAodqubits, listBeginlocation,
                                  listendLocation, dependency);
  }

  auto writeRearrangementInstruction(
      const std::size_t instIdx,
      const std::vector<std::vector<size_t>>& aodQubits,
      const std::vector<std::vector<
          std::tuple<std::size_t, const SLM*, std::size_t, std::size_t>>>&
          beginLocation,
      const std::vector<std::vector<
          std::tuple<std::size_t, const SLM*, std::size_t, std::size_t>>>&
          endLocation,
      nlohmann::json dependency) -> void {
    std::vector<std::vector<
        std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>>>
        beginLocationId{};
    beginLocationId.reserve(beginLocation.size());
    for (const auto& row : beginLocation) {
      std::vector<
          std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>>
          rowTuples{};
      rowTuples.reserve(row.size());
      for (const auto& qubit : row) {
        rowTuples.emplace_back(std::get<0>(qubit), std::get<1>(qubit)->id,
                               std::get<2>(qubit), std::get<3>(qubit));
      }
      beginLocationId.emplace_back(std::move(rowTuples));
    }
    std::vector<std::vector<
        std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>>>
        endLocationId{};
    endLocationId.reserve(endLocation.size());
    for (const auto& row : endLocation) {
      std::vector<
          std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>>
          rowTuples{};
      rowTuples.reserve(row.size());
      for (const auto& qubit : row) {
        rowTuples.emplace_back(std::get<0>(qubit), std::get<1>(qubit)->id,
                               std::get<2>(qubit), std::get<3>(qubit));
      }
      endLocationId.emplace_back(std::move(rowTuples));
    }
    nlohmann::json inst = {{"type", "rearrangeJob"},
                           {"id", instIdx},
                           {"aod_id", -1},
                           {"aod_qubits", aodQubits},
                           {"begin_locs", beginLocationId},
                           {"end_locs", endLocationId},
                           {"dependency", dependency}};
    inst["insts"] = expandArrangement(inst, beginLocation, endLocation);
    static_cast<T*>(this)->getResult().instructions.emplace_back(inst);
  }

  auto flattenRearrangmentInstruction() -> void {
    for (nlohmann::json& inst :
         static_cast<T*>(this)->getResult().instructions) {
      if (inst["type"] == "rearrangeJob") {
        nlohmann::json flattened{};
        for (const auto& row : inst["aod_qubits"]) {
          flattened.insert(flattened.end(), row.cbegin(), row.cend());
        }
        inst["aod_qubits"] = flattened;
        flattened.clear();
        for (const auto& row : inst["begin_locs"]) {
          flattened.insert(flattened.end(), row.cbegin(), row.cend());
        }
        inst["begin_locs"] = flattened;
        flattened.clear();
        for (const auto& row : inst["end_locs"]) {
          flattened.insert(flattened.end(), row.cbegin(), row.cend());
        }
        inst["end_locs"] = flattened;
      }
    }
  }

  /// generate a layer for gate execution
  auto processGateLayer(
      const std::size_t layer,
      const std::vector<std::tuple<const SLM*, size_t, size_t>>& gateMapping)
      -> void {
    const std::vector<std::size_t>& listGateidx =
        static_cast<T*>(this)->getGateSchedulingIdx()[layer];
    const std::vector<const std::pair<qc::Qubit, qc::Qubit>*>& listGate =
        static_cast<T*>(this)->getGateScheduling()[layer];
    const std::vector<const qc::StandardOperation*>& list1qgate =
        static_cast<T*>(this)->getGate1QScheduling()[layer];
    std::unordered_map<std::size_t, std::vector<std::size_t>> dictGateZone{};
    for (std::size_t i = 0; i < listGate.size(); ++i) {
      const auto* slmIdx = std::get<0>(gateMapping[listGate[i]->first]);
      const auto zoneIdx = *slmIdx->entanglementId;
      dictGateZone.try_emplace(zoneIdx, std::vector<std::size_t>{});
      dictGateZone[zoneIdx].emplace_back(i);
    }
    for (const auto& [rydbergIdx, gateIdxs] : dictGateZone) {
      std::vector<nlohmann::json> resultGate{};
      for (const auto i : gateIdxs) {
        resultGate.emplace_back(nlohmann::json{{"id", listGateidx[i]},
                                               {"q0", listGate[i]->first},
                                               {"q1", listGate[i]->second}});
      }
      std::unordered_set<std::size_t> setQubitdependency{};
      const std::size_t instIdx =
          static_cast<T*>(this)->getResult().instructions.size();
      for (const auto gateIdx : gateIdxs) {
        const auto gate = listGate[gateIdx];
        // collect qubit dependency
        setQubitdependency.emplace(qubitDependency[gate->first]);
        qubitDependency[gate->first] = instIdx;
        setQubitdependency.emplace(qubitDependency[gate->second]);
        qubitDependency[gate->second] = instIdx;
      }
      nlohmann::json dependency = {{"qubit", std::vector<std::size_t>{}},
                                   {"rydberg", rydbergDependency[rydbergIdx]}};
      rydbergDependency[rydbergIdx] = instIdx;
      dependency["qubit"] = setQubitdependency;
      writeGateInstruction(instIdx, rydbergIdx, resultGate, dependency);
    }

    // process single-qubit gates
    const std::size_t instIdx =
        static_cast<T*>(this)->getResult().instructions.size();
    std::vector<nlohmann::json> resultGate;
    std::unordered_set<std::size_t> setQubitdependency;
    for (const auto* gateInfo : list1qgate) {
      // collect qubit dependency
      const auto qubit = gateInfo->getTargets().front();
      setQubitdependency.emplace(qubitDependency[qubit]);
      qubitDependency[qubit] = instIdx;
      resultGate.emplace_back(
          nlohmann::json{{"name", gateInfo->getName()}, {"q", qubit}});
    }
    nlohmann::json dependency = {"qubit", setQubitdependency};
    if (!resultGate.empty()) {
      write1QGateInstruction(instIdx, resultGate, dependency, gateMapping);
    }
  }

  auto writeGateInstruction(const std::size_t instIdx,
                            const std::size_t rydbergIdx,
                            const std::vector<nlohmann::json>& resultGate,
                            const nlohmann::json& dependency) -> void {
    static_cast<T*>(this)->getResult().instructions.emplace_back(
        nlohmann::json{{"type", "rydberg"},
                       {"id", instIdx},
                       {"zone_id", rydbergIdx},
                       {"gates", resultGate},
                       {"dependency", dependency}});
  }

  auto write1QGateInstruction(
      std::size_t instIdx, const std::vector<nlohmann::json>& resultGate,
      const nlohmann::json& dependency,
      const std::vector<std::tuple<const SLM*, size_t, size_t>>& gateMapping)
      -> void {
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>>
        locs{};
    for (const auto& gate : resultGate) {
      locs.emplace_back(gate["q"], std::get<0>(gateMapping[gate["q"]])->id,
                        std::get<1>(gateMapping[gate["q"]]),
                        std::get<2>(gateMapping[gate["q"]]));
    }

    static_cast<T*>(this)->getResult().instructions.emplace_back(
        nlohmann::json{{"type", "1qGate"},
                       {"unitary", "u3"},
                       {"id", instIdx},
                       {"locs", locs},
                       {"gates", resultGate},
                       {"dependency", dependency}});
  }

  /// construct reverse movement layer by processing the forward movement
  auto constructReverseLayer(
      const std::size_t idLayerStart,
      const std::vector<std::tuple<const SLM*, size_t, size_t>>& initialMapping,
      const std::vector<std::tuple<const SLM*, size_t, size_t>>& finalMapping)
      -> void {
    const std::size_t idLayerend =
        static_cast<T*>(this)->getResult().instructions.size();
    for (std::size_t layer = idLayerStart; layer < idLayerend; ++layer) {
      if (static_cast<T*>(this)->getResult().instructions[layer]["type"] ==
          "rydberg") {
        break;
      } else {
        // process a rearrangement layer
        const std::size_t instIdx =
            static_cast<T*>(this)->getResult().instructions.size();
        nlohmann::json dependency = {{"qubit", std::vector<std::size_t>{}},
                                     {"site", std::vector<std::size_t>{}}};
        // process aod dependency
        std::unordered_set<std::size_t> setQubitdependency{};
        std::unordered_set<std::size_t> setSitedependency{};
        nlohmann::json listAodQubits = static_cast<T*>(this)
                                           ->getResult()
                                           .instructions[layer]["aod_qubits"];
        std::vector<std::vector<
            std::tuple<std::size_t, const SLM*, std::size_t, std::size_t>>>
            listEndLocation{};
        std::vector<std::vector<
            std::tuple<std::size_t, const SLM*, std::size_t, std::size_t>>>
            listBeginLocation{};
        for (const auto& subListQubits : listAodQubits) {
          std::vector<
              std::tuple<std::size_t, const SLM*, std::size_t, std::size_t>>
              rowBeginLocation{};
          std::vector<
              std::tuple<std::size_t, const SLM*, std::size_t, std::size_t>>
              rowendLocation{};
          for (const auto& q : subListQubits) {
            rowBeginLocation.emplace_back(q, std::get<0>(initialMapping[q]),
                                          std::get<1>(initialMapping[q]),
                                          std::get<2>(initialMapping[q]));
            rowendLocation.emplace_back(q, std::get<0>(finalMapping[q]),
                                        std::get<1>(finalMapping[q]),
                                        std::get<2>(finalMapping[q]));
            // process site dependency
            std::tuple siteKey{std::get<0>(finalMapping[q]),
                               std::get<1>(finalMapping[q]),
                               std::get<2>(finalMapping[q])};
            if (const auto& it = siteDependency.find(siteKey);
                it != siteDependency.end()) {
              setSitedependency.emplace(it->second);
            }
            siteKey = {std::get<0>(initialMapping[q]),
                       std::get<1>(initialMapping[q]),
                       std::get<2>(initialMapping[q])};
            siteDependency[siteKey] = instIdx;
            // collect qubit dependency
            setQubitdependency.emplace(qubitDependency[q]);
            qubitDependency[q] = instIdx;
          }

          listBeginLocation.emplace_back(rowBeginLocation);
          listEndLocation.emplace_back(rowendLocation);
        }
        dependency["qubit"] = setQubitdependency;
        dependency["site"] = setSitedependency;
        writeRearrangementInstruction(instIdx, listAodQubits, listBeginLocation,
                                      listEndLocation, dependency);
      }
    }
  }

  /// processs the aod assignment between two Rydberg stages
  auto aodAssignment(const std::size_t idLayerStart) -> void {
    std::vector listInstructionDuration(
        2, std::vector<std::pair<double, std::size_t>>{});
    const std::size_t idLayerEnd =
        static_cast<T*>(this)->getResult().instructions.size();
    std::size_t durationIdx = 0;
    std::vector<std::size_t> listGateLayerIdx{};
    for (std::size_t idx = idLayerStart; idx < idLayerEnd; ++idx) {
      if (static_cast<T*>(this)->getResult().instructions[idx]["type"] !=
          "rearrangeJob") {
        durationIdx = 1;
        listGateLayerIdx.emplace_back(idx);
        continue;
      }
      const double duration =
          getDuration(static_cast<T*>(this)->getResult().instructions[idx]);
      listInstructionDuration[durationIdx].emplace_back(duration, idx);
    }
    std::sort(listInstructionDuration[0].begin(),
              listInstructionDuration[0].end(),
              std::greater<std::pair<double, std::size_t>>());
    std::sort(listInstructionDuration[1].begin(),
              listInstructionDuration[1].end(),
              std::greater<std::pair<double, std::size_t>>());
    // assign instruction according to the duration in descending order
    for (const std::size_t i : {0UL, 1UL}) {
      for (const auto& item : listInstructionDuration[i]) {
        const double duration = item.first;
        nlohmann::json& inst =
            static_cast<T*>(this)->getResult().instructions[item.second];
        auto [beginTime, aodId] = aodEndTime.top();
        aodEndTime.pop();
        beginTime =
            std::max(beginTime, getBeginTime(item.second, inst["dependency"]));
        const double endTime = beginTime + duration;
        inst["dependency"]["aod"] = aodDependency[aodId];
        aodDependency[aodId] = item.second;
        inst["begin_time"] = beginTime;
        inst["end_time"] = endTime;
        inst["aod_id"] = aodId;
        aodEndTime.emplace(endTime, aodId);
        for (auto& detailInst : inst["insts"]) {
          detailInst["begin_time"] = detailInst["begin_time"].get<double>() + beginTime;
          detailInst["end_time"] = detailInst["end_time"].get<double>() + beginTime;
        }
        if (static_cast<T*>(this)->getResult().runtime < endTime) {
          static_cast<T*>(this)->getResult().runtime = endTime;
        }
      }
      if (i == 0) {
        for (const auto gateLayeridx : listGateLayerIdx) {
          // laser scheduling
          nlohmann::json& inst =
              static_cast<T*>(this)->getResult().instructions[gateLayeridx];
          const auto beginTime = getBeginTime(gateLayeridx, inst["dependency"]);
          const double endTime =
              inst["type"] == "rydberg"
                  ? beginTime +
                        static_cast<T*>(this)->getArchitecture().timeRydberg
                  : beginTime +
                        (static_cast<T*>(this)->getArchitecture().time1QGate *
                         inst["gates"].size()); // for sequential gate execution
          if (static_cast<T*>(this)->getResult().runtime < endTime) {
            static_cast<T*>(this)->getResult().runtime = endTime;
          }
          inst["begin_time"] = beginTime;
          inst["end_time"] = endTime;
        }
      }
    }
  }

  auto getBeginTime(const std::size_t curInstIdx, const nlohmann::json& dependency)
      -> double {
    double beginTime = 0;
    for (const auto& dependencyType : dependency.items()) {
      if (dependencyType.value().is_number_integer()) {
        const auto instIdx = dependencyType.value().get<std::size_t>();
        if (beginTime < static_cast<T*>(this)
                            ->getResult()
                            .instructions[instIdx]["end_time"]
                            .template get<double>()) {
          beginTime = static_cast<T*>(this)
                          ->getResult()
                          .instructions[instIdx]["end_time"];
        }
      } else {
        if (dependencyType.key() == "site") {
          for (const std::size_t instIdx : dependencyType.value()) {
            if (static_cast<T*>(this)
                    ->getResult()
                    .instructions[instIdx]["type"] == "rearrangeJob") {
              // find the time that the instruction finish atom transfer
              double atomTransferFinishTime = 0;
              for (const nlohmann::json& detailInst :
                   static_cast<T*>(this)->getResult().instructions[instIdx]
                                                                   ["insts"]) {
                std::string instType = detailInst["type"];
                instType = instType.substr(0, instType.find(":"));
                if (instType == "activate") {
                  atomTransferFinishTime =
                      std::max(detailInst["end_time"].get<double>(),
                               atomTransferFinishTime);
                }
              }
              // find the time until dropping of the qubits
              double atomTransferBeginTime = 0;
              for (const auto& detailInst :
                   static_cast<T*>(this)->getResult().instructions[curInstIdx]
                                                                   ["insts"]) {
                std::string instType = detailInst["type"];
                instType = instType.substr(0, instType.find(":"));
                if (instType == "deactivate") {
                  atomTransferBeginTime =
                      std::max(detailInst["begin_time"].template get<double>(),
                               atomTransferBeginTime);
                  break;
                }
              }
              const double tmpBegintime =
                  atomTransferFinishTime - atomTransferBeginTime;
              if (beginTime < tmpBegintime) {
                beginTime = tmpBegintime;
              }
            } else {
              if (beginTime < static_cast<T*>(this)
                                  ->getResult()
                                  .instructions[instIdx]["end_time"]
                                  .template get<double>()) {
                beginTime = static_cast<T*>(this)
                                ->getResult()
                                .instructions[instIdx]["end_time"];
              }
            }
          }
        } else {
          for (const std::size_t instIdx : dependencyType.value()) {
            if (beginTime < static_cast<T*>(this)
                                ->getResult()
                                .instructions[instIdx]["end_time"].template get<double>()) {
              beginTime = static_cast<T*>(this)
                              ->getResult()
                              .instructions[instIdx]["end_time"];
            }
          }
        }
      }
    }
    return beginTime;
  }

  auto getDuration(nlohmann::json& inst) -> double {
    auto& listDetailinst = inst["insts"];
    double duration = 0;

    for (auto& detailInst : listDetailinst) {
      std::string instType = detailInst["type"];
      instType = instType.substr(0, instType.find(":"));
      detailInst["begin_time"] = duration;
      if (instType == "activate" || instType == "deactivate") {
        duration +=
            static_cast<T*>(this)->getArchitecture().timeAtomTransfer;
        detailInst["end_time"] = duration;
      } else if (instType == "move") {
        double moveDuration = 0;
        for (std::size_t r = 0; r < detailInst["row_y_begin"].size(); ++r) {
          const auto& rowBegin = detailInst["row_y_begin"][r];
          const auto& rowEnd = detailInst["row_y_end"][r];
          for (std::size_t c = 0; c < detailInst["col_x_begin"].size(); ++c) {
            const auto& colBegin = detailInst["col_x_begin"][c];
            const auto& colEnd = detailInst["col_x_end"][c];
            const double tmp =
                static_cast<T*>(this)->getArchitecture().movementDuration(
                    colBegin, rowBegin, colEnd, rowEnd);
            if (moveDuration < tmp) {
              moveDuration = tmp;
            }
          }
        }
        detailInst["end_time"] = moveDuration + duration;
        duration += moveDuration;
      } else {
        throw std::invalid_argument(
            "[ERROR] Invalid instruction type in `get_duration`, must be one "
            "of 'activate', 'deactivate', 'move'");
      }
    }

    return duration;
  }

  auto expandArrangement(
      const nlohmann::json& inst,
      const std::vector<
          std::vector<std::tuple<size_t, const SLM*, size_t, size_t>>>&
          beginLocation,
      const std::vector<
          std::vector<std::tuple<size_t, const SLM*, size_t, size_t>>>&
          endLocation) -> nlohmann::json {
    nlohmann::json details{};

    // ---------------------- find out number of cols ----------------------
    std::vector<std::size_t> allColX{};   // all the x coord of qubits
    std::vector<nlohmann::json> coords{}; // coords of qubits
    // these coords are going to be updated as we construct the detail insts

    for (const auto& locs : beginLocation) {
      std::vector<nlohmann::json> coordsRow{};
      for (const auto& [q, slm, r, c] : locs) {
        const auto& [x, y] =
            static_cast<T*>(this)->getArchitecture().exactSlmLocation(slm, r, c);
        coordsRow.emplace_back(nlohmann::json{{"id", q}, {"x", x}, {"y", y}});
        allColX.emplace_back(x);
      }
      coords.emplace_back(std::move(coordsRow));
    }
    const auto initCoords = coords;

    std::sort(allColX.begin(), allColX.end());

    // assign AOD column ids based on all x coords needed
    std::unordered_map<std::size_t, std::size_t> colXToId{};
    colXToId.reserve(allColX.size());
    for (std::size_t i = 0; i < allColX.size(); ++i) {
      colXToId.emplace(allColX[i], i);
    }
    // ---------------------------------------------------------------------

    // -------------------- activation and parking -------------------------
    std::unordered_set<std::size_t> allColIdxSoFar{}; // which col has been activated
    for (std::size_t rowId = 0; rowId < beginLocation.size(); ++rowId) {
      // for each row
      const auto& locs = beginLocation[rowId];
      const std::size_t row_y =
          static_cast<T*>(this)
              ->getArchitecture()
              .exactSlmLocation(std::get<1>(locs.front()),
                                  std::get<2>(locs.front()),
                                  std::get<3>(locs.front()))
              .second;
      const std::pair rowLoc{std::get<1>(locs.front())->id,
                             std::get<2>(locs.front())};
      // before activation, adjust column position. This is necessary
      // whenever cols are parked (the `parking` movement below).
      nlohmann::json shiftBack = {
          {"type", "move"},
          {"move_type", "before"},
          {"row_id", std::vector<std::size_t>{}},
          {"row_y_begin", std::vector<std::size_t>{}},
          {"row_y_end", std::vector<std::size_t>{}},
          {"row_loc_begin", std::vector<std::size_t>{}},
          {"row_loc_end", std::vector<std::size_t>{}},
          {"col_id", std::vector<std::size_t>{}},
          {"col_x_begin", std::vector<std::size_t>{}},
          {"col_x_end", std::vector<std::size_t>{}},
          {"col_loc_begin",
           std::vector<std::pair<std::int64_t, std::int64_t>>{}},
          {"col_loc_end", std::vector<std::pair<std::int64_t, std::int64_t>>{}},
          {"begin_coord", coords},
          {"end_coord", std::vector<std::size_t>{}}};

      // activate one row and some columns
      nlohmann::json activate = {
          {"type", "activate"},
          {"row_id", std::vector{rowId}},
          {"row_y", std::vector{row_y}},
          {"row_loc", std::vector{rowLoc}},
          {"col_id", std::vector<std::size_t>{}},
          {"col_x", std::vector<std::size_t>{}},
          {"col_loc", std::vector<std::pair<std::size_t, std::size_t>>{}},
      };

      for (std::size_t j = 0; j < locs.size(); ++j) {
        const auto& [q, slm, r, c] = locs[j];
        const std::size_t colX = static_cast<T*>(this)
                                      ->getArchitecture()
                                      .exactSlmLocation(slm, r, c)
                                      .first;
        const std::pair colLoc{slm->id, c};
        const auto colId = colXToId[colX];
        if (allColIdxSoFar.find(colId) == allColIdxSoFar.end()) {
          // the col hasn't been activated, so there's no shift back
          // and we need to activate it at `col_x`.`
          allColIdxSoFar.emplace(colId);
          activate["col_id"].emplace_back(colId);
          activate["col_x"].emplace_back(colX);
          activate["col_loc"].emplace_back(colLoc);
        } else {
          // the col has been activated, thus parked previously and we
          // need the shift back, but we do not activate again.
          shiftBack["col_id"].emplace_back(colId);
          shiftBack["col_x_begin"].emplace_back(colX + parkingDist);
          shiftBack["col_x_end"].emplace_back(colX);
          shiftBack["col_loc_begin"].emplace_back(std::pair{-1, -1});
          shiftBack["col_loc_end"].emplace_back(colLoc);
          // since there's a shift, update the coords of the qubit
          coords[rowId][j]["x"] = colX;
        }
      }

      shiftBack["end_coord"] = coords;

      if (!shiftBack["col_id"].empty()) {
        details.emplace_back(shiftBack);
      }
      details.emplace_back(activate);

      if (rowId < inst["begin_locs"].size() - 1) {
        // parking movement after the activation
        // parking is required if we have activated some col, and there is
        // some qubit we don't want to pick up at the intersection of this
        // col and some future row to activate. We just always park here.
        // the last parking is not needed since there's a big move after it.
        nlohmann::json parking = {
            {"type", "move"},
            {"move_type", "after"},
            {"row_id", std::vector{rowId}},
            {"row_y_begin", std::vector{row_y}},
            {"row_y_end", std::vector{row_y + parkingDist}},
            {"row_loc_begin", std::vector{rowLoc}},
            {"row_loc_end", std::vector{std::pair{-1, -1}}},
            {"col_id", std::vector<std::size_t>{}},
            {"col_x_begin", std::vector<std::size_t>{}},
            {"col_x_end", std::vector<std::size_t>{}},
            {"col_loc_begin", std::vector<std::size_t>{}},
            {"col_loc_end", std::vector<std::size_t>{}},
            {"begin_coord", coords},
            {"end_coord", std::vector<std::size_t>{}}};
        for (std::size_t j = 0; j < locs.size(); ++j) {
          const auto& [q, slm, r, c] = locs[j];
          const std::size_t col_x = static_cast<T*>(this)
                                        ->getArchitecture()
                                        .exactSlmLocation(slm, r, c)
                                        .first;
          const std::pair colLoc{slm, c};
          const auto colId = colXToId[col_x];
          // all columns used in this row are parked after the activation
          parking["col_id"].emplace_back(colId);
          parking["col_x_begin"].emplace_back(col_x);
          parking["col_x_end"].emplace_back(col_x + parkingDist);
          parking["col_loc_begin"].emplace_back(
              std::pair{colLoc.first->id, colLoc.second});
          parking["col_loc_end"].emplace_back(std::pair{-1, -1});
          coords[rowId][j]["x"] = parking["col_x_end"].back();
          coords[rowId][j]["y"] = parking["row_y_end"].front();
        }
        parking["end_coord"] = coords;
        details.emplace_back(parking);
      }
    }
    // ---------------------------------------------------------------------

    // ------------------------- big move ----------------------------------
    nlohmann::json bigMove = {
        {"type", "move:big"},
        {"move_type", "big"},
        {"row_id", std::vector<std::size_t>{}},
        {"row_y_begin", std::vector<std::size_t>{}},
        {"row_y_end", std::vector<std::size_t>{}},
        {"row_loc_begin", std::vector<std::size_t>{}},
        {"row_loc_end", std::vector<std::size_t>{}},
        {"col_id", std::vector<std::size_t>{}},
        {"col_x_begin", std::vector<std::size_t>{}},
        {"col_x_end", std::vector<std::size_t>{}},
        {"col_loc_begin", std::vector<std::size_t>{}},
        {"col_loc_end", std::vector<std::size_t>{}},
        {"begin_coord", coords},
        {"end_coord", std::vector<std::size_t>{}},
    };

    for (std::size_t rowId = 0; rowId < beginLocation.size(); ++rowId) {
      const auto& beginLocs = beginLocation[rowId];
      const auto& endLocs = endLocation[rowId];
      bigMove["row_id"].emplace_back(rowId);
      bigMove["row_y_begin"].emplace_back(coords[rowId][0]["y"]);
      if (initCoords[rowId][0]["y"] == coords[rowId][0]["y"]) {
        // AOD row is align with SLM row
        bigMove["row_loc_begin"].emplace_back(std::pair{
            std::get<1>(beginLocs[0])->id, std::get<2>(beginLocs[0])});
      } else {
        bigMove["row_loc_begin"].emplace_back(std::pair{-1, -1});
      }
      bigMove["row_y_end"].emplace_back(
          static_cast<T*>(this)
              ->getArchitecture()
              .exactSlmLocation(std::get<1>(endLocs.front()),
                                  std::get<2>(endLocs.front()),
                                  std::get<3>(endLocs.front()))
              .second);
      bigMove["row_loc_end"].emplace_back(std::pair{
          std::get<1>(endLocs.front())->id, std::get<2>(endLocs.front())});

      for (std::size_t j = 0; j < beginLocs.size(); ++j) {
        const auto& beginLoc = beginLocs[j];
        const auto& endLoc = endLocs[j];
        const auto col_x = static_cast<T*>(this)
                               ->getArchitecture()
                               .exactSlmLocation(std::get<1>(beginLoc),
                                                   std::get<2>(beginLoc),
                                                   std::get<3>(beginLoc))
                               .first;
        const std::size_t colId = colXToId[col_x];

        if (std::find(bigMove["col_id"].cbegin(), bigMove["col_id"].cend(),
                      colId) == bigMove["col_id"].cend()) {
          // the movement of this rol has not been recorded before
          bigMove["col_id"].emplace_back(colId);
          bigMove["col_x_begin"].emplace_back(coords[rowId][j]["x"]);
          if (initCoords[rowId][j]["x"] == coords[rowId][j]["x"]) {
            // AOD col is aligned with SLM col
            bigMove["col_loc_begin"].emplace_back(
                std::pair{std::get<1>(beginLoc)->id, std::get<3>(beginLoc)});
          } else {
            bigMove["col_loc_begin"].emplace_back(std::pair{-1, -1});
          }
          bigMove["col_x_end"].emplace_back(
              static_cast<T*>(this)
                  ->getArchitecture()
                  .exactSlmLocation(std::get<1>(endLoc), std::get<2>(endLoc),
                                      std::get<3>(endLoc))
                  .first);
          bigMove["col_loc_end"].emplace_back(
              std::pair{std::get<1>(endLoc)->id, std::get<3>(endLoc)});
        }
        // whether or not the movement of this col has been considered
        // before, we need to update the coords of the qubit.
        coords[rowId][j]["x"] =
            static_cast<T*>(this)
                ->getArchitecture()
                .exactSlmLocation(std::get<1>(endLoc), std::get<2>(endLoc),
                                  std::get<3>(endLoc))
                .first;
        coords[rowId][j]["y"] =
            static_cast<T*>(this)
                ->getArchitecture()
                .exactSlmLocation(std::get<1>(endLocs.front()),
                                  std::get<2>(endLocs.front()),
                                  std::get<3>(endLocs.front()))
                .second;
      }
    }
    bigMove["end_coord"] = coords;
    details.emplace_back(bigMove);
    // ---------------------------------------------------------------------

    // --------------------------- deactivation ----------------------------
    details.emplace_back(
        nlohmann::json{{"type", "deactivate"},
                       {"row_id", std::vector<std::size_t>{}},
                       {"col_id", std::vector<std::size_t>{}}});
    for (std::size_t i = 0; i < beginLocation.size(); ++i) {
      details.back()["row_id"].emplace_back(i);
    }
    for (std::size_t i = 0; i < allColX.size(); ++i) {
      details.back()["col_id"].emplace_back(i);
    }
    // ---------------------------------------------------------------------

    for (std::size_t instCounter = 0; instCounter < details.size();
         ++instCounter) {
      auto& detailInst = details[instCounter];
      detailInst["id"] = instCounter;
    }
    return details;
  }
};

} // namespace na
