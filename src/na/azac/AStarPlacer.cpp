#include "na/azac/AStarPlacer.hpp"

#include "Definitions.hpp"
#include "na/azac/Architecture.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <nlohmann/json_fwd.hpp>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na::azac {
template <class Node>
auto AStarPlacer::aStarTreeSearch(
    const Node& start,
    const std::function<std::vector<std::reference_wrapper<const Node>>(
        const Node&)>& getNeighbors,
    const std::function<bool(const Node&)>& isGoal,
    const std::function<double(const Node&)>& getCost,
    const std::function<double(const Node&)>& getHeuristic)
    -> std::vector<std::reference_wrapper<const Node>> {
  //===--------------------------------------------------------------------===//
  // Setup open set structure
  //===--------------------------------------------------------------------===//
  // struct for items in the open set
  struct Item {
    double priority;  //< sum of cost and heuristic
    const Node* node; //< pointer to the node
    // pointer to the parent item to reconstruct the path in the end
    Item* parent;

    Item(const double priority, const Node& node, Item* parent)
        : priority(priority), node(&node), parent(parent) {
      assert(!std::isnan(priority));
    }
  };
  // compare function for the open set
  struct ItemCompare {
    bool operator()(const Item* a, const Item* b) const {
      // this way, the item with the lowest priority is on top of the heap
      return a->priority > b->priority;
    }
  };
  // vector of items to store all items and keep them alive also after they
  // are popped from the open set. they are required alive to reconstruct the
  // path in the end.
  std::vector<std::unique_ptr<Item>> items;
  // open list of nodes to be evaluated as a minimum heap based on the
  // priority. whenever an item is placed in the queue it is created in the
  // vector `items` before and only a reference is placed in the queue
  std::priority_queue<Item*, std::vector<Item*>, ItemCompare> openSet;
  openSet.emplace(items
                      .emplace_back(std::make_unique<Item>(getHeuristic(start),
                                                           start, nullptr))
                      .get());
  //===--------------------------------------------------------------------===//
  // Perform A* search
  //===--------------------------------------------------------------------===//
  while (!openSet.empty()) {
    Item* itm = openSet.top();
    openSet.pop();
    // if a goal is reached, that is the shortest path to a goal under the
    // assumption that the heuristic is admissible
    if (isGoal(*itm->node)) {
      // reconstruct the path from the goal to the start and then reverse it
      std::vector<std::reference_wrapper<const Node>> path;
      for (; itm != nullptr; itm = itm->parent) {
        path.emplace_back(*itm->node);
      }
      std::reverse(path.begin(), path.end());
      return path;
    }
    // expand the current node by adding all neighbors to the open set
    const auto& neighbors = getNeighbors(*itm->node);
    if (!neighbors.empty()) {
      for (const auto& neighbor : neighbors) {
        // getCost returns the total cost to reach the current node
        const auto cost = getCost(neighbor);
        const auto heuristic = getHeuristic(neighbor);
        openSet.emplace(items
                            .emplace_back(std::make_unique<Item>(
                                cost + heuristic, neighbor, itm))
                            .get());
      }
    }
  }
  throw std::runtime_error("No path from start to any goal found.");
}

auto AStarPlacer::discretizePlacementOfAtoms(
    const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                 size_t>>& placement,
    const std::vector<qc::Qubit>& atoms) const
    -> std::pair<
        std::unordered_map<std::pair<std::reference_wrapper<const SLM>, size_t>,
                           uint8_t, std::hash<std::pair<const SLM&, size_t>>,
                           std::equal_to<std::pair<const SLM&, size_t>>>,
        std::unordered_map<std::pair<std::reference_wrapper<const SLM>, size_t>,
                           uint8_t, std::hash<std::pair<const SLM&, size_t>>,
                           std::equal_to<std::pair<const SLM&, size_t>>>> {
  std::map<size_t, std::unordered_set<
                       std::pair<std::reference_wrapper<const SLM>, size_t>,
                       std::hash<std::pair<const SLM&, size_t>>,
                       std::equal_to<std::pair<const SLM&, size_t>>>>
      rows;
  std::map<size_t, std::unordered_set<
                       std::pair<std::reference_wrapper<const SLM>, size_t>,
                       std::hash<std::pair<const SLM&, size_t>>,
                       std::equal_to<std::pair<const SLM&, size_t>>>>
      columns;
  for (const auto atom : atoms) {
    const auto& [slm, r, c] = placement[atom];
    const auto& [x, y] = architecture_.get().exactSlmLocation(slm, r, c);
    rows.try_emplace(y).first->second.emplace(std::tie(slm, r));
    columns.try_emplace(x).first->second.emplace(std::tie(slm, c));
  }
  std::unordered_map<std::pair<std::reference_wrapper<const SLM>, size_t>,
                     uint8_t, std::hash<std::pair<const SLM&, size_t>>,
                     std::equal_to<std::pair<const SLM&, size_t>>>
      rowIndices;
  uint8_t rowIndex = 0;
  for (const auto& [_, sites] : rows) {
    for (const auto& site : sites) {
      rowIndices.emplace(site, rowIndex);
    }
    ++rowIndex;
  }
  std::unordered_map<std::pair<std::reference_wrapper<const SLM>, size_t>,
                     uint8_t, std::hash<std::pair<const SLM&, size_t>>,
                     std::equal_to<std::pair<const SLM&, size_t>>>
      columnIndices;
  uint8_t columnIndex = 0;
  for (const auto& [_, sites] : columns) {
    for (const auto& site : sites) {
      columnIndices.emplace(site, columnIndex);
    }
    ++columnIndex;
  }
  return std::pair{rowIndices, columnIndices};
}

auto AStarPlacer::discretizeNonOccupiedStorageSites(
    const std::unordered_set<
        std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>,
        std::hash<std::tuple<const SLM&, size_t, size_t>>,
        std::equal_to<std::tuple<const SLM&, size_t, size_t>>>& occupiedSites)
    const -> std::pair<
        std::unordered_map<std::pair<std::reference_wrapper<const SLM>, size_t>,
                           uint8_t, std::hash<std::pair<const SLM&, size_t>>,
                           std::equal_to<std::pair<const SLM&, size_t>>>,
        std::unordered_map<std::pair<std::reference_wrapper<const SLM>, size_t>,
                           uint8_t, std::hash<std::pair<const SLM&, size_t>>,
                           std::equal_to<std::pair<const SLM&, size_t>>>> {
  std::map<size_t, std::pair<std::reference_wrapper<const SLM>, size_t>> rows;
  std::map<size_t, std::pair<std::reference_wrapper<const SLM>, size_t>>
      columns;
  for (const auto& slm : architecture_.get().storageZones) {
    // find rows with free sites
    for (size_t r = 0; r < slm->nRows; ++r) {
      for (size_t c = 0; c < slm->nCols; ++c) {
        if (occupiedSites.find(std::tie(*slm, r, c)) == occupiedSites.end()) {
          // free site in row r found at column c
          rows.emplace(slm->location.second + (slm->siteSeparation.second * r),
                       std::tie(*slm, r));
          break;
        }
      }
    }
    // find columns with free sites
    for (size_t c = 0; c < slm->nCols; ++c) {
      for (size_t r = 0; r < slm->nRows; ++r) {
        if (occupiedSites.find(std::tie(*slm, r, c)) == occupiedSites.end()) {
          // free site in column c found at row r
          columns.emplace(slm->location.first + (slm->siteSeparation.first * c),
                          std::tie(*slm, c));
          break;
        }
      }
    }
  }
  std::unordered_map<std::pair<std::reference_wrapper<const SLM>, size_t>,
                     uint8_t, std::hash<std::pair<const SLM&, size_t>>,
                     std::equal_to<std::pair<const SLM&, size_t>>>
      rowIndices;
  uint8_t rowIndex = 0;
  for (const auto& [_, site] : rows) {
    rowIndices.emplace(site, rowIndex++);
  }
  std::unordered_map<std::pair<std::reference_wrapper<const SLM>, size_t>,
                     uint8_t, std::hash<std::pair<const SLM&, size_t>>,
                     std::equal_to<std::pair<const SLM&, size_t>>>
      columnIndices;
  uint8_t columnIndex = 0;
  for (const auto& [_, site] : columns) {
    columnIndices.emplace(site, columnIndex++);
  }
  return std::pair{rowIndices, columnIndices};
}

auto AStarPlacer::discretizeNonOccupiedEntanglementSites(
    const std::unordered_set<
        std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>,
        std::hash<std::tuple<const SLM&, size_t, size_t>>,
        std::equal_to<std::tuple<const SLM&, size_t, size_t>>>& occupiedSites)
    const -> std::pair<
        std::unordered_map<std::pair<std::reference_wrapper<const SLM>, size_t>,
                           uint8_t, std::hash<std::pair<const SLM&, size_t>>,
                           std::equal_to<std::pair<const SLM&, size_t>>>,
        std::unordered_map<std::pair<std::reference_wrapper<const SLM>, size_t>,
                           uint8_t, std::hash<std::pair<const SLM&, size_t>>,
                           std::equal_to<std::pair<const SLM&, size_t>>>> {
  std::map<size_t, std::unordered_set<
                       std::pair<std::reference_wrapper<const SLM>, size_t>,
                       std::hash<std::pair<const SLM&, size_t>>,
                       std::equal_to<std::pair<const SLM&, size_t>>>>
      rows;
  std::map<size_t, std::unordered_set<
                       std::pair<std::reference_wrapper<const SLM>, size_t>,
                       std::hash<std::pair<const SLM&, size_t>>,
                       std::equal_to<std::pair<const SLM&, size_t>>>>
      columns;
  for (const auto& zone : architecture_.get().entanglementZones) {
    for (const auto& slm : *zone) {
      // find rows with free sites
      for (size_t r = 0; r < slm.nRows; ++r) {
        for (size_t c = 0; c < slm.nCols; ++c) {
          if (occupiedSites.find(std::tie(slm, r, c)) == occupiedSites.end()) {
            // free site in row r found at column c
            rows.try_emplace(slm.location.second +
                             (slm.siteSeparation.second * r))
                .first->second.emplace(std::tie(slm, r));
            break;
          }
        }
      }
      // find columns with free sites
      for (size_t c = 0; c < slm.nCols; ++c) {
        for (size_t r = 0; r < slm.nRows; ++r) {
          if (occupiedSites.find(std::tie(slm, r, c)) == occupiedSites.end()) {
            // free site in column c found at row r
            columns
                .try_emplace(slm.location.first +
                             (slm.siteSeparation.first * c))
                .first->second.emplace(std::tie(slm, c));
            break;
          }
        }
      }
    }
  }
  std::unordered_map<std::pair<std::reference_wrapper<const SLM>, size_t>,
                     uint8_t, std::hash<std::pair<const SLM&, size_t>>,
                     std::equal_to<std::pair<const SLM&, size_t>>>
      rowIndices;
  uint8_t rowIndex = 0;
  for (const auto& [_, sites] : rows) {
    for (const auto& site : sites) {
      rowIndices.emplace(site, rowIndex);
    }
    ++rowIndex;
  }
  std::unordered_map<std::pair<std::reference_wrapper<const SLM>, size_t>,
                     uint8_t, std::hash<std::pair<const SLM&, size_t>>,
                     std::equal_to<std::pair<const SLM&, size_t>>>
      columnIndices;
  uint8_t columnIndex = 0;
  for (const auto& [_, sites] : columns) {
    for (const auto& site : sites) {
      columnIndices.emplace(site, columnIndex);
    }
    ++columnIndex;
  }
  return std::pair{rowIndices, columnIndices};
}

auto AStarPlacer::makeInitialPlacement(const size_t nQubits) const
    -> std::vector<
        std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>> {
  auto slmIt = architecture_.get().storageZones.cbegin();
  std::size_t c = 0;
  std::int64_t r = reverseInitialPlacement_
                       ? static_cast<std::int64_t>((*slmIt)->nRows) - 1
                       : 0;
  const std::int64_t step = reverseInitialPlacement_ ? -1 : 1;
  std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>
      initialPlacement;
  initialPlacement.reserve(nQubits);
  for (qc::Qubit qubit = 0; qubit < nQubits; ++qubit) {
    initialPlacement.emplace_back(**slmIt, r, c++);
    if (c == (*slmIt)->nCols) {
      // the end of the row reached, go to the next row
      r += step;
      c = 0;
      if (r == static_cast<std::int64_t>((*slmIt)->nRows)) {
        // the end of the slm reached, go to the next slm
        ++slmIt;
        if (step > 0) {
          r = static_cast<std::int64_t>((*slmIt)->nRows) - 1;
        } else {
          r = 0;
        }
      }
    }
  }
  return initialPlacement;
}

auto AStarPlacer::makeIntermediatePlacement(
    const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                 size_t>>& previousPlacement,
    const std::unordered_set<qc::Qubit>& previousReuseQubits,
    const std::unordered_set<qc::Qubit>& reuseQubits,
    const std::vector<std::array<qc::Qubit, 2>>& twoQubitGates)
    -> std::pair<std::vector<std::tuple<std::reference_wrapper<const SLM>,
                                        size_t, size_t>>,
                 std::vector<std::tuple<std::reference_wrapper<const SLM>,
                                        size_t, size_t>>> {
  const auto& gatePlacement = placeGatesInEntanglementZone(
      previousPlacement, previousReuseQubits, twoQubitGates);
  return {gatePlacement,
          placeQubitsInStorageZone(gatePlacement, reuseQubits, twoQubitGates)};
}

auto AStarPlacer::addGateOption(
    const std::unordered_map<
        std::pair<std::reference_wrapper<const SLM>, size_t>, uint8_t,
        std::hash<std::pair<const SLM&, size_t>>,
        std::equal_to<std::pair<const SLM&, size_t>>>& discreteTargetRows,
    const std::unordered_map<
        std::pair<std::reference_wrapper<const SLM>, size_t>, uint8_t,
        std::hash<std::pair<const SLM&, size_t>>,
        std::equal_to<std::pair<const SLM&, size_t>>>& discreteTargetColumns,
    const SLM& leftSlm, const size_t leftRow, const size_t leftCol,
    const SLM& rightSlm, const size_t rightRow, const size_t rightCol,
    const SLM& nearestSlm, const size_t r, const size_t c, GateJob& job) const
    -> void {
  //                  other
  //         ┌─┐       ┌─┐ <-- Entanglement sites
  //         └┬┘       └┬┘
  //          │╲dis2   ╱│
  //     dis1 │  ╲   ╱  │
  //          │    ╳    │
  //          │  ╱   ╲  │ dis4
  //          │╱dis3   ╲│
  //         ┌┴┐       ┌┴┐ <-- Storage sites
  //         └─┘       └─┘
  //          ^         ^
  //        atom1     atom2
  const auto& [otherSlm, otherRow, otherCol] =
      architecture_.get().otherEntanglementSite(nearestSlm, r, c);
  const auto dis1 = static_cast<float>(architecture_.get().distance(
      leftSlm, leftRow, leftCol, nearestSlm, r, c));
  const auto dis2 = static_cast<float>(architecture_.get().distance(
      rightSlm, rightRow, rightCol, nearestSlm, r, c));
  const auto dis3 = static_cast<float>(architecture_.get().distance(
      leftSlm, leftRow, leftCol, otherSlm, otherRow, otherCol));
  const auto dis4 = static_cast<float>(architecture_.get().distance(
      rightSlm, rightRow, rightCol, otherSlm, otherRow, otherCol));
  if (dis1 + dis4 <= dis2 + dis3) {
    job.options.emplace_back(GateJob::Option{
        std::array{
            std::array{discreteTargetRows.at(std::tie(nearestSlm, r)),
                       discreteTargetColumns.at(std::tie(nearestSlm, c))},
            std::array{discreteTargetRows.at(std::tie(otherSlm, otherRow)),
                       discreteTargetColumns.at(std::tie(otherSlm, otherCol))}},
        std::array{dis1, dis4}});
  } else {
    job.options.emplace_back(GateJob::Option{
        std::array{
            std::array{discreteTargetRows.at(std::tie(otherSlm, otherRow)),
                       discreteTargetColumns.at(std::tie(otherSlm, otherCol))},
            std::array{discreteTargetRows.at(std::tie(nearestSlm, r)),
                       discreteTargetColumns.at(std::tie(nearestSlm, c))}},
        std::array{dis2, dis3}});
  }
}

auto AStarPlacer::placeGatesInEntanglementZone(
    const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                 size_t>>& previousPlacement,
    const std::unordered_set<qc::Qubit>& reuseQubits,
    const std::vector<std::array<qc::Qubit, 2>>& twoQubitGates)
    -> std::vector<
        std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>> {
  // Duplicate the previous placement as a starting point for the current
  std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>
      currentPlacement = previousPlacement;
  //===------------------------------------------------------------------===//
  // Find gates and atoms that must be placed
  //===------------------------------------------------------------------===//
  std::set<std::pair<double, std::array<qc::Qubit, 2>>, std::greater<>>
      gatesToPlace;
  std::vector<qc::Qubit> atomsToPlace;
  for (const auto& gate : twoQubitGates) {
    const auto& [first, second] = gate;
    if (const auto firstQubitReuse =
            reuseQubits.find(first) != reuseQubits.end();
        !firstQubitReuse && reuseQubits.find(second) == reuseQubits.end()) {
      const auto& [storageSlm1, storageRow1, storageCol1] =
          previousPlacement[first];
      const auto& [storageSlm2, storageRow2, storageCol2] =
          previousPlacement[second];
      const auto& [nearestSlm, nearestRow, nearestCol] =
          architecture_.get().nearestEntanglementSite(storageSlm1, storageRow1,
                                                      storageCol1, storageSlm2,
                                                      storageRow2, storageCol2);
      const auto& [otherSlm, otherRow, otherCol] =
          architecture_.get().otherEntanglementSite(nearestSlm, nearestRow,
                                                    nearestCol);
      // Calculate the various distances to the entanglement sites
      //
      // Example:
      //       nearest    other
      //         ┌─┐       ┌─┐ <-- Entanglement sites
      //         └┬┘       └┬┘
      //          │╲dis2   ╱│
      //     dis1 │  ╲   ╱  │
      //          │    ╳    │
      //          │  ╱   ╲  │ dis4
      //          │╱dis3   ╲│
      //         ┌┴┐       ┌┴┐ <-- Storage sites
      //         └─┘       └─┘
      //          ^         ^
      //   gate->first   gate->second
      const auto dis1 =
          architecture_.get().distance(storageSlm1, storageRow1, storageCol1,
                                       nearestSlm, nearestRow, nearestCol);
      const auto dis2 =
          architecture_.get().distance(storageSlm2, storageRow2, storageCol2,
                                       nearestSlm, nearestRow, nearestCol);
      const auto dis3 = architecture_.get().distance(
          storageSlm1, storageRow1, storageCol1, otherSlm, otherRow, otherCol);
      const auto dis4 = architecture_.get().distance(
          storageSlm2, storageRow2, storageCol2, otherSlm, otherRow, otherCol);
      if (dis1 + dis4 <= dis2 + dis3) {
        // if the situation is as depicted in the example above
        gatesToPlace.emplace(std::max(dis1, dis4), gate);
      } else {
        // otherwise, either the entanglement sites or storage sites are
        // flipped
        gatesToPlace.emplace(std::max(dis2, dis3), gate);
      }
      atomsToPlace.emplace_back(first);
      atomsToPlace.emplace_back(second);
    } else {
      // one of the qubits is reused, so no need to place a gate
      if (firstQubitReuse) {
        const auto& [slm, r, c] = previousPlacement[first];
        currentPlacement[second] =
            architecture_.get().otherEntanglementSite(slm, r, c);
      } else {
        // second qubit is reused
        const auto& [slm, r, c] = previousPlacement[second];
        currentPlacement[first] =
            architecture_.get().otherEntanglementSite(slm, r, c);
      }
    }
  }
  if (gatesToPlace.empty()) {
    return currentPlacement;
  }
  //===------------------------------------------------------------------===//
  // Discretize the previous placement of the atoms to be placed
  //===------------------------------------------------------------------===//
  const auto& [discreteRows, discreteColumns] =
      discretizePlacementOfAtoms(previousPlacement, atomsToPlace);
  //===------------------------------------------------------------------===//
  // Extract occupied entanglement sites from the previous placement
  //===------------------------------------------------------------------===//
  // This set will only contain the first SLM in a pair of entanglement
  // SLMs and represents the occupied pair of sites
  std::unordered_set<
      std::tuple<std::reference_wrapper<const SLM>, std::size_t, std::size_t>,
      std::hash<std::tuple<const SLM&, size_t, size_t>>,
      std::equal_to<std::tuple<const SLM&, size_t, size_t>>>
      occupiedEntanglementSites{};
  for (const auto qubit : reuseQubits) {
    const auto& [slm, r, c] = previousPlacement[qubit];
    assert(slm.get().isEntanglement());
    occupiedEntanglementSites.emplace(slm, r, c);
    // also add the interaction partner's site to the set
    const auto& [otherSlm, otherRow, otherCol] =
        architecture_.get().otherEntanglementSite(slm, r, c);
    occupiedEntanglementSites.emplace(otherSlm, otherRow, otherCol);
  }
  //===------------------------------------------------------------------===//
  // Discretize the free sites for the atoms to be placed
  //===------------------------------------------------------------------===//
  const auto& [discreteTargetRows, discreteTargetColumns] =
      discretizeNonOccupiedEntanglementSites(occupiedEntanglementSites);
  //===------------------------------------------------------------------===//
  // Initialize the gate jobs
  //===------------------------------------------------------------------===//
  /// The number of either atom or gate jobs that must be performed
  const size_t nJobs = gatesToPlace.size();
  /// a list of all gates that must be placed in the entanglement zone before a
  /// rydberg layer
  std::vector<GateJob> gateJobs;
  gateJobs.reserve(nJobs);
  for (const auto& [_, gate] : gatesToPlace) {
    const auto& [leftAtom, rightAtom] = gate;
    const auto& [leftSlm, leftRow, leftCol] = previousPlacement[leftAtom];
    const auto& [rightSlm, rightRow, rightCol] = previousPlacement[rightAtom];
    const auto& [nearestSlm, nearestRow, nearestCol] =
        architecture_.get().nearestEntanglementSite(
            leftSlm, leftRow, leftCol, rightSlm, rightRow, rightCol);
    auto& job = gateJobs.emplace_back();
    job.qubits = gate;
    job.currentSites = std::array{
        std::array{discreteRows.at(std::tie(leftSlm, leftRow)),
                   discreteColumns.at(std::tie(leftSlm, leftCol))},
        std::array{discreteRows.at(std::tie(rightSlm, rightRow)),
                   discreteColumns.at(std::tie(rightSlm, rightCol))}};
    size_t rLow = 0;
    size_t rHigh = nearestSlm.get().nRows;
    size_t cLow = 0;
    size_t cHigh = nearestSlm.get().nCols;
    if (useWindow_) {
      rLow = nearestRow > windowMinHeight_ / 2
                 ? nearestRow - (windowMinHeight_ / 2)
                 : 0;
      rHigh = std::min(nearestRow + (windowMinHeight_ / 2) + 1,
                       nearestSlm.get().nRows);
      cLow = nearestCol > windowMinWidth_ / 2
                 ? nearestCol - (windowMinWidth_ / 2)
                 : 0;
      cHigh = std::min(nearestCol + (windowMinWidth_ / 2) + 1,
                       nearestSlm.get().nCols);
    }
    for (size_t r = rLow; r < rHigh; ++r) {
      for (size_t c = cLow; c < cHigh; ++c) {
        if (occupiedEntanglementSites.find(std::tie(nearestSlm, r, c)) ==
            occupiedEntanglementSites.end()) {
          addGateOption(discreteTargetRows, discreteTargetColumns, leftSlm,
                        leftRow, leftCol, rightSlm, rightRow, rightCol,
                        nearestSlm, r, c, job);
        }
      }
    }
    size_t expansion = 0;
    while (useWindow_ && static_cast<double>(job.options.size()) <
                             windowShare_ * static_cast<double>(nJobs)) {
      // window does not contain enough options, so expand it
      ++expansion;
      size_t windowWidth = 0;
      size_t windowHeight = 0;
      if (windowRatio_ < 1.0) {
        // landscpe ==> expand width and adjust height
        windowWidth = windowMinWidth_ + expansion;
        windowHeight = static_cast<size_t>(
            std::round(windowRatio_ * static_cast<double>(windowWidth)));
      } else {
        // portrait ==> expand height and adjust width
        windowHeight = windowMinHeight_ + expansion;
        windowWidth = static_cast<size_t>(
            std::round(static_cast<double>(windowHeight) / windowRatio_));
      }
      auto rLowNew =
          nearestRow > windowHeight / 2 ? nearestRow - (windowHeight / 2) : 0;
      auto rHighNew =
          std::min(nearestRow + (windowHeight / 2) + 1, nearestSlm.get().nRows);
      auto cLowNew =
          nearestCol > windowWidth / 2 ? nearestCol - (windowWidth / 2) : 0;
      auto cHighNew =
          std::min(nearestCol + (windowWidth / 2) + 1, nearestSlm.get().nCols);
      if (rLowNew < rLow) {
        assert(rLow - rLowNew == 1);
        for (size_t c = cLowNew; c < cHighNew; ++c) {
          if (occupiedEntanglementSites.find(std::tie(
                  nearestSlm, rLowNew, c)) == occupiedEntanglementSites.end()) {
            addGateOption(discreteTargetRows, discreteTargetColumns, leftSlm,
                          leftRow, leftCol, rightSlm, rightRow, rightCol,
                          nearestSlm, rLowNew, c, job);
          }
        }
      }
      if (rHighNew > rHigh) {
        assert(rHighNew - rHigh == 1);
        for (size_t c = cLowNew; c < cHighNew; ++c) {
          // NOTE: we have to use rHighNew - 1 here, which is equal to rHigh
          if (occupiedEntanglementSites.find(std::tie(nearestSlm, rHigh, c)) ==
              occupiedEntanglementSites.end()) {
            addGateOption(discreteTargetRows, discreteTargetColumns, leftSlm,
                          leftRow, leftCol, rightSlm, rightRow, rightCol,
                          nearestSlm, rHigh, c, job);
          }
        }
      }
      if (cLowNew < cLow) {
        assert(cLow - cLowNew == 1);
        for (size_t r = rLow; r < rHigh; ++r) {
          if (occupiedEntanglementSites.find(std::tie(
                  nearestSlm, r, cLowNew)) == occupiedEntanglementSites.end()) {
            addGateOption(discreteTargetRows, discreteTargetColumns, leftSlm,
                          leftRow, leftCol, rightSlm, rightRow, rightCol,
                          nearestSlm, r, cLowNew, job);
          }
        }
      }
      if (cHighNew > cHigh) {
        assert(cHighNew - cHigh == 1);
        for (size_t r = rLow; r < rHigh; ++r) {
          // NOTE: we have to use cHighNew - 1 here, which is equal to cHigh
          if (occupiedEntanglementSites.find(std::tie(nearestSlm, r, cHigh)) ==
              occupiedEntanglementSites.end()) {
            addGateOption(discreteTargetRows, discreteTargetColumns, leftSlm,
                          leftRow, leftCol, rightSlm, rightRow, rightCol,
                          nearestSlm, r, cHigh, job);
          }
        }
      }
      rLow = rLowNew;
      rHigh = rHighNew;
      cLow = cLowNew;
      cHigh = cHighNew;
    }
    std::sort(
        job.options.begin(), job.options.end(),
        [](const GateJob::Option& lhs, const GateJob::Option& rhs) -> bool {
          return lhs.distance < rhs.distance;
        });
  }
  //===------------------------------------------------------------------===//
  // Get the extent of discrete source and target
  //===------------------------------------------------------------------===//
  const uint8_t maxDiscreteSourceRow =
      std::max_element(discreteRows.begin(), discreteRows.end(),
                       [](const auto& lhs, const auto& rhs) {
                         return lhs.second < rhs.second;
                       })
          ->second;
  const uint8_t maxDiscreteSourceColumn =
      std::max_element(discreteColumns.begin(), discreteColumns.end(),
                       [](const auto& lhs, const auto& rhs) {
                         return lhs.second < rhs.second;
                       })
          ->second;
  const uint8_t maxDiscreteTargetRow =
      std::max_element(discreteTargetRows.begin(), discreteTargetRows.end(),
                       [](const auto& lhs, const auto& rhs) {
                         return lhs.second < rhs.second;
                       })
          ->second;
  const uint8_t maxDiscreteTargetColumn =
      std::max_element(discreteTargetColumns.begin(),
                       discreteTargetColumns.end(),
                       [](const auto& lhs, const auto& rhs) {
                         return lhs.second < rhs.second;
                       })
          ->second;
  const std::array<float, 2> scaleFactors{
      std::min(1.F, static_cast<float>(1 + maxDiscreteTargetRow) /
                        static_cast<float>(1 + maxDiscreteSourceRow)),
      std::min(1.F, static_cast<float>(1 + maxDiscreteTargetColumn) /
                        static_cast<float>(1 + maxDiscreteSourceColumn))};
  //===------------------------------------------------------------------===//
  // Run the A* algorithm
  //===------------------------------------------------------------------===//
  /// A list of all nodes that have been created so far.
  /// This list is dynamically extended when new nodes are created.
  /// This happens when a node is expanded by calling getNeighbors.
  std::deque<std::unique_ptr<GateNode>> nodes;
  // make the root node
  nodes.emplace_back(std::make_unique<GateNode>());
  const auto deepeningFactor = deepeningFactor_;
  const auto& path = aStarTreeSearch<GateNode>(
      *nodes.front(),
      [&nodes, &gateJobs](const auto& node) {
        return getGatePlacementNeighbors(nodes, gateJobs, std::move(node));
      },
      [nJobs](const auto& node) { return isGoal(2 * nJobs, std::move(node)); },
      [](const auto& node) { return getCost(std::move(node)); },
      [&gateJobs, deepeningFactor, &scaleFactors](const auto& node) {
        return getGatePlacementHeuristic(gateJobs, deepeningFactor,
                                         scaleFactors, std::move(node));
      });
  //===------------------------------------------------------------------===//
  // Extract the final mapping
  //===------------------------------------------------------------------===//
  std::unordered_map<
      uint8_t,
      std::unordered_map<uint8_t, std::tuple<std::reference_wrapper<const SLM>,
                                             size_t, size_t>>>
      targetSites;
  for (const auto& [row, r] : discreteTargetRows) {
    const SLM& slm = row.first.get();
    auto& targetSitesForThisRow = targetSites.try_emplace(r).first->second;
    for (const auto& [column, c] : discreteTargetColumns) {
      if (&slm == &column.first.get()) {
        targetSitesForThisRow.emplace(c,
                                      std::tie(slm, row.second, column.second));
      }
    }
  }
  assert(!targetSites.empty());
  assert(path.size() == nJobs + 1);
  for (size_t i = 0; i < nJobs; ++i) {
    const auto& job = gateJobs[i];
    const auto& option = *path[i + 1].get().option;
    for (size_t j = 0; j < 2; ++j) {
      const auto atom = job.qubits[j];
      const auto& [row, col] = option.sites[j];
      currentPlacement[atom] = targetSites.at(row).at(col);
    }
  }
  return currentPlacement;
}

auto AStarPlacer::placeQubitsInStorageZone(
    const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                 size_t>>& previousPlacement,
    const std::unordered_set<qc::Qubit>& reuseQubits,
    const std::vector<std::array<qc::Qubit, 2>>& twoQubitGates)
    -> std::vector<
        std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>> {
  // Duplicate the previous placement as a starting point for the current
  std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>
      currentPlacement = previousPlacement;
  //===------------------------------------------------------------------===//
  // Find atoms that must be placed
  //===------------------------------------------------------------------===//
  std::set<std::pair<double, qc::Qubit>, std::greater<>> atomsToPlaceSet;
  for (const auto& gate : twoQubitGates) {
    for (const auto qubit : gate) {
      if (reuseQubits.find(qubit) == reuseQubits.end()) {
        const auto& [slm, r, c] = previousPlacement[qubit];
        const auto& [nearestSlm, nearestRow, nearestCol] =
            architecture_.get().nearestStorageSite(slm, r, c);
        const auto distance = architecture_.get().distance(
            slm, r, c, nearestSlm, nearestRow, nearestCol);
        atomsToPlaceSet.emplace(distance, qubit);
      }
    }
  }
  if (atomsToPlaceSet.empty()) {
    return currentPlacement;
  }
  //===------------------------------------------------------------------===//
  // Discretize the previous placement of the atoms to be placed
  //===------------------------------------------------------------------===//
  std::vector<qc::Qubit> atomsToPlace;
  atomsToPlace.reserve(atomsToPlace.size());
  for (const auto& [_, atom] : atomsToPlaceSet) {
    atomsToPlace.emplace_back(atom);
  }
  // Place the atoms with the longest distance first, but the place atoms in
  // increasing order of distance to the first atom
  const auto& frontPlacement = previousPlacement[atomsToPlace.front()];
  std::set<std::pair<double, qc::Qubit>> atomsWithoutFirstAtom;
  for (auto atomIt = atomsToPlace.cbegin() + 1; atomIt != atomsToPlace.cend();
       ++atomIt) {
    const auto& [frontSlm, frontRow, frontCol] = frontPlacement;
    const auto& [slm, r, c] = previousPlacement[*atomIt];
    const auto distance =
        architecture_.get().distance(slm, r, c, frontSlm, frontRow, frontCol);
    atomsWithoutFirstAtom.emplace(distance, *atomIt);
  }
  std::transform(atomsWithoutFirstAtom.cbegin(), atomsWithoutFirstAtom.cend(),
                 atomsToPlace.begin() + 1,
                 [](const auto& pair) { return pair.second; });
  // Discretize the previous placement of the atoms to be placed that are
  // ordered now
  const auto& [discreteRows, discreteColumns] =
      discretizePlacementOfAtoms(previousPlacement, atomsToPlace);
  //===------------------------------------------------------------------===//
  // Extract occupied storage sites from the previous placement
  //===------------------------------------------------------------------===//
  std::unordered_set<
      std::tuple<std::reference_wrapper<const SLM>, std::size_t, std::size_t>,
      std::hash<std::tuple<const SLM&, size_t, size_t>>,
      std::equal_to<std::tuple<const SLM&, size_t, size_t>>>
      occupiedStorageSites{};
  for (const auto& [slm, r, c] : previousPlacement) {
    if (slm.get().isStorage()) {
      occupiedStorageSites.emplace(slm, r, c);
    }
  }
  //===------------------------------------------------------------------===//
  // Discretize the free sites for the atoms to be placed
  //===------------------------------------------------------------------===//
  const auto& [discreteTargetRows, discreteTargetColumns] =
      discretizeNonOccupiedStorageSites(occupiedStorageSites);
  //===------------------------------------------------------------------===//
  // Initialize atom jobs
  //===------------------------------------------------------------------===//
  /// The number of either atom or gate jobs that must be performed
  size_t nJobs = atomsToPlace.size();
  /// a list of all atoms that must be placed in the storage zone after a
  /// rydberg layer
  std::vector<AtomJob> atomJobs;
  atomJobs.reserve(nJobs);
  uint8_t minDiscreteColumnOfNearestSite = std::numeric_limits<uint8_t>::max();
  uint8_t maxDiscreteColumnOfNearestSite = 0;
  for (const auto atom : atomsToPlace) {
    const auto& [previousSlm, previousRow, previousCol] =
        previousPlacement[atom];
    const auto& [nearestSlm, nearestRow, nearestCol] =
        architecture_.get().nearestStorageSite(previousSlm, previousRow,
                                               previousCol);
    const auto discreteColumnOfNearestSite =
        discreteTargetColumns.at(std::tie(nearestSlm, nearestCol));
    minDiscreteColumnOfNearestSite =
        std::min(minDiscreteColumnOfNearestSite, discreteColumnOfNearestSite);
    maxDiscreteColumnOfNearestSite =
        std::max(maxDiscreteColumnOfNearestSite, discreteColumnOfNearestSite);
    auto& job = atomJobs.emplace_back();
    job.qubit = atom;
    job.currentSite =
        std::array{discreteRows.at(std::tie(previousSlm, previousRow)),
                   discreteColumns.at(std::tie(previousSlm, previousCol))};
    size_t rLow = 0;
    size_t rHigh = nearestSlm.get().nRows;
    size_t cLow = 0;
    size_t cHigh = nearestSlm.get().nCols;
    if (useWindow_) {
      rLow = nearestRow > windowMinHeight_ / 2
                 ? nearestRow - (windowMinHeight_ / 2)
                 : 0;
      rHigh = std::min(nearestRow + (windowMinHeight_ / 2) + 1,
                       nearestSlm.get().nRows);
      cLow = nearestCol > windowMinWidth_ / 2
                 ? nearestCol - (windowMinWidth_ / 2)
                 : 0;
      cHigh = std::min(nearestCol + (windowMinWidth_ / 2) + 1,
                       nearestSlm.get().nCols);
    }
    for (size_t r = rLow; r < rHigh; ++r) {
      for (size_t c = cLow; c < cHigh; ++c) {
        if (occupiedStorageSites.find(std::tie(nearestSlm, r, c)) ==
            occupiedStorageSites.end()) {
          const auto distance = static_cast<float>(architecture_.get().distance(
              previousSlm, previousRow, previousCol, nearestSlm, r, c));
          job.options.emplace_back(AtomJob::Option{
              {discreteTargetRows.at(std::tie(nearestSlm, r)),
               discreteTargetColumns.at(std::tie(nearestSlm, c))},
              distance});
        }
      }
    }
    size_t expansion = 0;
    while (useWindow_ && static_cast<double>(job.options.size()) <
                             windowShare_ * static_cast<double>(nJobs)) {
      // window does not contain enough options, so expand it
      ++expansion;
      size_t windowWidth = 0;
      size_t windowHeight = 0;
      if (windowRatio_ < 1.0) {
        // landscpe ==> expand width and adjust height
        // the overall width and height is divided by 2 later, hence an
        // expansion of 2 is needed to actually increase the window size
        windowWidth = windowMinWidth_ + 2 * expansion;
        windowHeight = static_cast<size_t>(
            std::round(windowRatio_ * static_cast<double>(windowWidth)));
      } else {
        // portrait ==> expand height and adjust width
        // the overall width and height is divided by 2 later, hence an
        // expansion of 2 is needed to actually increase the window size
        windowHeight = windowMinHeight_ + 2 * expansion;
        windowWidth = static_cast<size_t>(
            std::round(static_cast<double>(windowHeight) / windowRatio_));
      }
      auto rLowNew =
          nearestRow > windowHeight / 2 ? nearestRow - (windowHeight / 2) : 0;
      auto rHighNew =
          std::min(nearestRow + (windowHeight / 2) + 1, nearestSlm.get().nRows);
      auto cLowNew =
          nearestCol > windowWidth / 2 ? nearestCol - (windowWidth / 2) : 0;
      auto cHighNew =
          std::min(nearestCol + (windowWidth / 2) + 1, nearestSlm.get().nCols);
      if (rLowNew < rLow) {
        assert(rLow - rLowNew == 1);
        for (size_t c = cLowNew; c < cHighNew; ++c) {
          if (occupiedStorageSites.find(std::tie(nearestSlm, rLowNew, c)) ==
              occupiedStorageSites.end()) {
            const auto distance =
                static_cast<float>(architecture_.get().distance(
                    previousSlm, previousRow, previousCol, nearestSlm, rLowNew,
                    c));
            job.options.emplace_back(AtomJob::Option{
                {discreteTargetRows.at(std::tie(nearestSlm, rLowNew)),
                 discreteTargetColumns.at(std::tie(nearestSlm, c))},
                distance});
          }
        }
      }
      if (rHighNew > rHigh) {
        assert(rHighNew - rHigh == 1);
        for (size_t c = cLowNew; c < cHighNew; ++c) {
          // NOTE: we have to use rHighNew - 1 here, which is equal to rHigh
          if (occupiedStorageSites.find(std::tie(nearestSlm, rHigh, c)) ==
              occupiedStorageSites.end()) {
            const auto distance =
                static_cast<float>(architecture_.get().distance(
                    previousSlm, previousRow, previousCol, nearestSlm, rHigh,
                    c));
            job.options.emplace_back(AtomJob::Option{
                {discreteTargetRows.at(std::tie(nearestSlm, rHigh)),
                 discreteTargetColumns.at(std::tie(nearestSlm, c))},
                distance});
          }
        }
      }
      if (cLowNew < cLow) {
        assert(cLow - cLowNew == 1);
        for (size_t r = rLow; r < rHigh; ++r) {
          if (occupiedStorageSites.find(std::tie(nearestSlm, r, cLowNew)) ==
              occupiedStorageSites.end()) {
            const auto distance =
                static_cast<float>(architecture_.get().distance(
                    previousSlm, previousRow, previousCol, nearestSlm, r,
                    cLowNew));
            job.options.emplace_back(AtomJob::Option{
                {discreteTargetRows.at(std::tie(nearestSlm, r)),
                 discreteTargetColumns.at(std::tie(nearestSlm, cLowNew))},
                distance});
          }
        }
      }
      if (cHighNew > cHigh) {
        assert(cHighNew - cHigh == 1);
        for (size_t r = rLow; r < rHigh; ++r) {
          // NOTE: we have to use cHighNew - 1 here, which is equal to cHigh
          if (occupiedStorageSites.find(std::tie(nearestSlm, r, cHigh)) ==
              occupiedStorageSites.end()) {
            const auto distance =
                static_cast<float>(architecture_.get().distance(
                    previousSlm, previousRow, previousCol, nearestSlm, r,
                    cHigh));
            job.options.emplace_back(AtomJob::Option{
                {discreteTargetRows.at(std::tie(nearestSlm, r)),
                 discreteTargetColumns.at(std::tie(nearestSlm, cHigh))},
                distance});
          }
        }
      }
      rLow = rLowNew;
      rHigh = rHighNew;
      cLow = cLowNew;
      cHigh = cHighNew;
    }
    std::sort(
        job.options.begin(), job.options.end(),
        [](const AtomJob::Option& lhs, const AtomJob::Option& rhs) -> bool {
          return lhs.distance < rhs.distance;
        });
  }
  //===------------------------------------------------------------------===//
  // Get the extent of discrete source and target
  //===------------------------------------------------------------------===//
  const uint8_t maxDiscreteSourceRow =
      std::max_element(discreteRows.begin(), discreteRows.end(),
                       [](const auto& lhs, const auto& rhs) {
                         return lhs.second < rhs.second;
                       })
          ->second;
  const uint8_t maxDiscreteSourceColumn =
      std::max_element(discreteColumns.begin(), discreteColumns.end(),
                       [](const auto& lhs, const auto& rhs) {
                         return lhs.second < rhs.second;
                       })
          ->second;
  const uint8_t maxDiscreteTargetRow =
      std::max_element(discreteTargetRows.begin(), discreteTargetRows.end(),
                       [](const auto& lhs, const auto& rhs) {
                         return lhs.second < rhs.second;
                       })
          ->second;
  const uint8_t maxDiscreteTargetColumn =
      std::max_element(discreteTargetColumns.begin(),
                       discreteTargetColumns.end(),
                       [](const auto& lhs, const auto& rhs) {
                         return lhs.second < rhs.second;
                       })
          ->second;
  const std::array<float, 2> scaleFactors{
      std::min(1.F, static_cast<float>(1 + maxDiscreteTargetRow) /
                        static_cast<float>(1 + maxDiscreteSourceRow)),
      std::min(
          std::max(1.F, static_cast<float>(maxDiscreteColumnOfNearestSite -
                                           minDiscreteColumnOfNearestSite) /
                            static_cast<float>(maxDiscreteSourceColumn)),
          static_cast<float>(1 + maxDiscreteTargetColumn) /
              static_cast<float>(1 + maxDiscreteSourceColumn))};
  //===------------------------------------------------------------------===//
  // Run the A* algorithm
  //===------------------------------------------------------------------===//
  /// A list of all nodes that have been created so far.
  /// This list is dynamically extended when new nodes are created.
  /// This happens when a node is expanded by calling getNeighbors.
  std::deque<std::unique_ptr<AtomNode>> nodes;
  nodes.emplace_back(std::make_unique<AtomNode>());
  const auto deepeningFactor = deepeningFactor_;
  const auto& path = aStarTreeSearch<AtomNode>(
      *nodes.front(),
      [&nodes, &atomJobs](const auto& node) {
        return getAtomPlacementNeighbors(nodes, atomJobs, std::move(node));
      },
      [nJobs](const auto& node) { return isGoal(nJobs, std::move(node)); },
      [](const auto& node) { return getCost(std::move(node)); },
      [&atomJobs, deepeningFactor, &scaleFactors](const auto& node) {
        return getAtomPlacementHeuristic(atomJobs, deepeningFactor,
                                         scaleFactors, std::move(node));
      });
  //===------------------------------------------------------------------===//
  // Extract the final mapping
  //===------------------------------------------------------------------===//
  std::unordered_map<
      uint8_t,
      std::unordered_map<uint8_t, std::tuple<std::reference_wrapper<const SLM>,
                                             size_t, size_t>>>
      targetSites;
  for (const auto& [row, r] : discreteTargetRows) {
    const SLM& slm = row.first.get();
    auto& targetSitesForThisRow = targetSites.try_emplace(r).first->second;
    for (const auto& [column, c] : discreteTargetColumns) {
      if (&slm == &column.first.get()) {
        targetSitesForThisRow.emplace(c,
                                      std::tie(slm, row.second, column.second));
      }
    }
  }
  assert(!targetSites.empty());
  assert(path.size() == nJobs + 1);
  for (size_t i = 0; i < nJobs; ++i) {
    const auto& job = atomJobs[i];
    const auto& option = *path[i + 1].get().option;
    const auto atom = job.qubit;
    const auto& [row, col] = option.site;
    currentPlacement[atom] = targetSites.at(row).at(col);
  }
  return currentPlacement;
}

template <class Node> auto AStarPlacer::getCost(const Node& node) -> float {
  float cost = 0.0;
  for (const auto d : node.maxDistancesOfPlacedAtomsPerGroup) {
    cost += std::sqrt(d);
  }
  return cost;
}

auto AStarPlacer::sumStdDeviationForGroups(
    const std::array<float, 2>& scaleFactors,
    const std::vector<std::array<std::map<uint8_t, uint8_t>, 2>>& groups)
    -> float {
  float sumStdDev = 0.F;
  for (const auto& groupPair : groups) {
    for (size_t i = 0; i < 2; ++i) {
      const auto& group = groupPair[i];
      std::vector<float> diffs;
      diffs.reserve(groups.size());
      for (const auto& [key, value] : group) {
        diffs.emplace_back(static_cast<float>(value) -
                           (scaleFactors[i] * static_cast<float>(key)));
      }
      const auto n = static_cast<float>(group.size());
      const auto mean = std::accumulate(diffs.cbegin(), diffs.cend(), 0.F) / n;
      const auto variance =
          std::accumulate(diffs.cbegin(), diffs.cend(), 0.F,
                          [mean](const float acc, const float diff) -> float {
                            return acc + ((diff - mean) * (diff - mean));
                          }) /
          n;
      sumStdDev += std::sqrt(variance);
    }
  }
  return sumStdDev;
}

auto AStarPlacer::getAtomPlacementHeuristic(
    const std::vector<AtomJob>& atomJobs, const float deepeningFactor,
    const std::array<float, 2>& scaleFactors, const AtomNode& node) -> float {
  const auto nAtomJobs = atomJobs.size();
  const auto nUnplacedAtoms =
      static_cast<float>(nAtomJobs - node.consumedFreeSites.size());
  float maxDistanceOfUnplacedAtom = 0.0;
  for (size_t i = node.consumedFreeSites.size(); i < nAtomJobs; ++i) {
    for (const auto& option : atomJobs[i].options) {
      if (node.consumedFreeSites.find(option.site) ==
          node.consumedFreeSites.end()) {
        // this assumes that the first found free site is the nearest free site
        // for that atom. This requires that the job options are sorted by
        // distance.
        maxDistanceOfUnplacedAtom =
            std::max(maxDistanceOfUnplacedAtom, option.distance);
        break; // exit when the first free site is found
      }
    }
  }
  float heuristic = maxDistanceOfUnplacedAtom <= node.maxDistanceOfPlacedAtom
                        ? 0.F
                        : std::sqrt(maxDistanceOfUnplacedAtom) -
                              std::sqrt(node.maxDistanceOfPlacedAtom);
  heuristic += deepeningFactor *
               (sumStdDeviationForGroups(scaleFactors, node.groups) + 0.2F) *
               nUnplacedAtoms;
  return heuristic;
}

auto AStarPlacer::getGatePlacementHeuristic(
    const std::vector<GateJob>& gateJobs, const float deepeningFactor,
    const std::array<float, 2>& scaleFactors, const GateNode& node) -> float {
  const auto nGateJobs = gateJobs.size();
  assert(node.consumedFreeSites.size() % 2 == 0);
  const auto nUnplacedGates = // NOLINTNEXTLINE(bugprone-integer-division)
      static_cast<float>(nGateJobs - (node.consumedFreeSites.size() / 2));
  float maxDistanceOfUnplacedAtom = 0.0;
  for (size_t i = node.consumedFreeSites.size() / 2; i < nGateJobs; ++i) {
    for (const auto& option : gateJobs[i].options) {
      // this assumes that the first found pair of free sites is the nearest
      // pair of free sites for that gate. This requires that the job options
      // are sorted by distance.
      if (std::all_of(option.sites.cbegin(), option.sites.cend(),
                      [&node](const DiscreteSite& site) -> bool {
                        return node.consumedFreeSites.find(site) ==
                               node.consumedFreeSites.end();
                      })) {
        maxDistanceOfUnplacedAtom =
            std::max(maxDistanceOfUnplacedAtom,
                     *std::max_element(option.distance.cbegin(),
                                       option.distance.cend()));
        break; // exit when the first free site pair is found
      }
    }
  }
  float heuristic = maxDistanceOfUnplacedAtom <= node.maxDistanceOfPlacedAtom
                        ? 0.F
                        : std::sqrt(maxDistanceOfUnplacedAtom) -
                              std::sqrt(node.maxDistanceOfPlacedAtom);
  heuristic += deepeningFactor *
               (sumStdDeviationForGroups(scaleFactors, node.groups) + 0.2F) *
               nUnplacedGates;
  return heuristic;
}

auto AStarPlacer::getAtomPlacementNeighbors(
    std::deque<std::unique_ptr<AtomNode>>& nodes,
    const std::vector<AtomJob>& atomJobs, const AtomNode& node)
    -> std::vector<std::reference_wrapper<const AtomNode>> {
  const size_t atomToBePlacedNext = node.consumedFreeSites.size();
  const auto& atomJob = atomJobs[atomToBePlacedNext];
  std::vector<std::reference_wrapper<const AtomNode>> neighbors;
  for (const auto& option : atomJob.options) {
    const auto& [site, distance] = option;
    // skip the site is already consumed
    if (node.consumedFreeSites.find(site) != node.consumedFreeSites.end()) {
      continue;
    }
    // make a copy of node, the parent of neighbor
    AtomNode& neighbor = *nodes.emplace_back(std::make_unique<AtomNode>(node));
    neighbor.option = &option;
    neighbor.maxDistanceOfPlacedAtom =
        std::max(node.maxDistanceOfPlacedAtom, distance);
    neighbor.consumedFreeSites.emplace(site);
    // check whether the current placement is compatible with any existing group
    checkCompatibilityAndAddPlacement(
        atomJob.currentSite.front(), site.front(), atomJob.currentSite.back(),
        site.back(), distance, neighbor.groups,
        neighbor.maxDistancesOfPlacedAtomsPerGroup);
    // add the neighbor to the list of neighbors to be returned
    neighbors.emplace_back(neighbor);
  }
  return neighbors;
}

auto AStarPlacer::getGatePlacementNeighbors(
    std::deque<std::unique_ptr<GateNode>>& nodes,
    const std::vector<GateJob>& gateJobs, const GateNode& node)
    -> std::vector<std::reference_wrapper<const GateNode>> {
  const size_t gateToBePlacedNext = node.consumedFreeSites.size() / 2;
  const auto& gateJob = gateJobs[gateToBePlacedNext];
  std::vector<std::reference_wrapper<const GateNode>> neighbors;
  // Get the current placement of the atoms that must be placed next
  const auto& [currentSiteOfLeftAtom, currentSiteOfRightAtom] =
      gateJob.currentSites;
  for (const auto& option : gateJob.options) {
    const auto& [sites, distances] = option;
    const auto& [leftSite, rightSite] = sites;
    // skip if one of the sites is already consumed
    if (node.consumedFreeSites.find(leftSite) != node.consumedFreeSites.end() ||
        node.consumedFreeSites.find(rightSite) !=
            node.consumedFreeSites.end()) {
      continue;
    }
    // make a copy of node, the parent of neighbor as use this as a starting
    // point for the new node
    GateNode& neighbor = *nodes.emplace_back(std::make_unique<GateNode>(node));
    neighbor.option = &option;
    neighbor.maxDistanceOfPlacedAtom =
        std::max(node.maxDistanceOfPlacedAtom,
                 *std::max_element(distances.cbegin(), distances.cend()));
    neighbor.consumedFreeSites.emplace(leftSite);
    neighbor.consumedFreeSites.emplace(rightSite);
    // check whether the current placement is compatible with any existing
    // horizontal group
    checkCompatibilityAndAddPlacement(
        currentSiteOfLeftAtom.front(), leftSite.front(),
        currentSiteOfLeftAtom.back(), leftSite.back(), distances.front(),
        neighbor.groups, neighbor.maxDistancesOfPlacedAtomsPerGroup);
    checkCompatibilityAndAddPlacement(
        currentSiteOfRightAtom.front(), rightSite.front(),
        currentSiteOfRightAtom.back(), rightSite.back(), distances.back(),
        neighbor.groups, neighbor.maxDistancesOfPlacedAtomsPerGroup);
    // add the final neighbor to the list of neighbors to be returned
    neighbors.emplace_back(neighbor);
  }
  return neighbors;
}

auto AStarPlacer::checkCompatibilityWithGroup(
    const uint8_t key, const uint8_t value,
    const std::map<uint8_t, uint8_t>& group)
    -> std::optional<
        std::pair<std::map<uint8_t, uint8_t>::const_iterator, bool>> {
  if (auto it = group.lower_bound(key); it != group.end()) {
    // an assignment for this key already exists in this group
    if (const auto& [upperKey, upperValue] = *it; upperKey == key) {
      if (upperValue == value) {
        // new placement is compatible with this group and key already exists
        return std::pair{it, true};
      }
    } else {
      // if (upperKey > key)
      if (it != group.begin()) {
        // it can be safely decremented
        if (const auto& [_, lowerValue] = *std::prev(it);
            lowerValue < value && value < upperValue) {
          // new placement is compatible with this group
          return std::pair{it, false};
        }
      } else {
        // if (it == hGroup.begin())
        if (value < upperValue) {
          // new placement is compatible with this group
          return std::pair{it, false};
        }
      }
    }
  } else {
    // if (it == hGroup.end())
    // it can be safely decremented because the group must contain
    // at least one element
    if (const auto& [_, lowerValue] = *std::prev(it); lowerValue < value) {
      // new placement is compatible with this group
      return std::pair{it, false};
    }
  }
  return std::nullopt;
}

auto AStarPlacer::checkCompatibilityAndAddPlacement(
    const uint8_t hKey, const uint8_t hValue, const uint8_t vKey,
    const uint8_t vValue, const float distance,
    std::vector<std::array<std::map<uint8_t, uint8_t>, 2>>& groups,
    std::vector<float>& maxDistances) -> bool {
  size_t i = 0;
  for (auto& group : groups) {
    auto& [hGroup, vGroup] = group;
    if (const auto& hCompatible =
            checkCompatibilityWithGroup(hKey, hValue, hGroup);
        hCompatible) {
      if (const auto& vCompatible =
              checkCompatibilityWithGroup(vKey, vValue, vGroup)) {
        const auto& [hIt, hExists] = *hCompatible;
        const auto& [vIt, vExists] = *vCompatible;
        // new placement is compatible with this group
        if (!hExists) {
          hGroup.emplace_hint(hIt, hKey, hValue);
        }
        if (!vExists) {
          vGroup.emplace_hint(vIt, vKey, vValue);
        }
        maxDistances[i] = std::max(maxDistances[i], distance);
        return true;
      }
    }
    ++i;
  }
  // no compatible group could be found and a new group is created
  auto& [hGroup, vGroup] = groups.emplace_back();
  hGroup.emplace(hKey, hValue);
  vGroup.emplace(vKey, vValue);
  maxDistances.emplace_back(distance);
  return false;
}

AStarPlacer::AStarPlacer(const Architecture& architecture,
                         const nlohmann::json& config)
    : architecture_(architecture) {
  // get first storage SLM and first entanglement SLM
  const auto& firstStorageSlm = *architecture_.get().storageZones.front();
  const auto& firstEntanglementSlm =
      architecture_.get().entanglementZones.front()->front();
  // check which side of the first storage SLM is closer to the entanglement
  // SLM
  if (firstStorageSlm.location.second < firstEntanglementSlm.location.second) {
    // if the entanglement SLM is closer to the last row of the storage SLM
    // start initial placement of the atoms in the last row instead of the
    // first and hence revert initial placement
    reverseInitialPlacement_ = true;
  }
  // check whether the config contains information on the windowed placement
  if (const auto& configIt = config.find("a_star_placer");
      configIt != config.end() && configIt->is_object()) {
    bool useWindowSet = false;
    bool windowMinWidthSet = false;
    bool windowRatioSet = false;
    bool windowShareSet = false;
    bool deepeningFactorSet = false;
    bool lookaheadFactorSet = false;
    for (const auto& [key, value] : configIt.value().items()) {
      if (key == "use_window") {
        if (value.is_boolean()) {
          useWindow_ = value;
          useWindowSet = true;
        } else {
          std::ostringstream oss;
          oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer "
                 "contains an invalid "
                 "value for use_window. Using default.\n";
          std::cout << oss.str();
        }
      } else if (key == "window_min_width") {
        if (value.is_number_unsigned()) {
          windowMinWidth_ = value;
          windowMinWidthSet = true;
        } else {
          std::ostringstream oss;
          oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer "
                 "contains an invalid value for window_min_width. Using "
                 "default.\n";
          std::cout << oss.str();
        }
      } else if (key == "window_ratio") {
        if (value.is_number()) {
          windowRatio_ = value;
          windowRatioSet = true;
        } else {
          std::ostringstream oss;
          oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer "
                 "contains an invalid value for window_ratio. Using default.\n";
          std::cout << oss.str();
        }
      } else if (key == "window_share") {
        if (value.is_number()) {
          windowShare_ = value;
          windowShareSet = true;
        } else {
          std::ostringstream oss;
          oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer "
                 "contains an invalid value for window_share. Using default.\n";
          std::cout << oss.str();
        }
      } else if (key == "deepening_factor") {
        if (value.is_number()) {
          deepeningFactor_ = value;
          deepeningFactorSet = true;
        } else {
          std::ostringstream oss;
          oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer "
                 "contains an invalid value for deepening_factor. Using "
                 "default.\n";
          std::cout << oss.str();
        }
      } else if (key == "lookahead_factor") {
        if (value.is_number()) {
          lookaheadFactor_ = value;
          lookaheadFactorSet = true;
        } else {
          std::ostringstream oss;
          oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer "
                 "contains an invalid value for lookahead_factor. Using "
                 "default.\n";
          std::cout << oss.str();
        }
      } else {
        std::ostringstream oss;
        oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer contains "
               "an unknown key: "
            << key << ". Ignoring.\n";
        std::cout << oss.str();
      }
    }
    if (!useWindowSet) {
      std::cout << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer does "
                   "not contain a setting for use_window. Using default.\n";
    }
    if (useWindow_) {
      if (!windowMinWidthSet) {
        std::cout << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer "
                     "does not contain a setting for window_min_width. Using "
                     "default.\n";
      }
      if (!windowRatioSet) {
        std::cout
            << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer "
               "does not contain a setting for window_ratio. Using default.\n";
      }
      if (!windowShareSet) {
        std::cout
            << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer "
               "does not contain a setting for window_share. Using default.\n";
      }
    }
    if (windowMinWidthSet) {
      windowMinHeight_ = static_cast<size_t>(
          std::round(windowRatio_ * static_cast<double>(windowMinWidth_)));
    }
    if (!deepeningFactorSet) {
      std::cout
          << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer does "
             "not contain a setting for deepening_factor. Using default.\n";
    }
    if (!lookaheadFactorSet) {
      std::cout
          << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer does "
             "not contain a setting for lookahead_factor. Using default.\n";
    }
  } else {
    std::cout << "\033[1;35m[WARN]\033[0m Configuration does not contain "
                 "settings for "
                 "AStarPlacer or is malformed. Using default settings.\n";
  }
}

auto AStarPlacer::place(
    const size_t nQubits,
    const std::vector<std::vector<std::array<qc::Qubit, 2>>>&
        twoQubitGateLayers,
    const std::vector<std::unordered_set<qc::Qubit>>& reuseQubits)
    -> std::vector<std::vector<
        std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>> {
  std::vector<std::vector<
      std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>>
      placement;
  placement.reserve((2 * twoQubitGateLayers.size()) + 1);
  placement.emplace_back(makeInitialPlacement(nQubits));
  for (size_t layer = 0; layer < twoQubitGateLayers.size(); ++layer) {
    const auto& [gatePlacement, qubitPlacement] = makeIntermediatePlacement(
        placement.back(),
        layer == 0 ? std::unordered_set<qc::Qubit>{} : reuseQubits[layer - 1],
        layer == reuseQubits.size() ? std::unordered_set<qc::Qubit>{}
                                    : reuseQubits[layer],
        twoQubitGateLayers[layer]);
    placement.emplace_back(gatePlacement);
    placement.emplace_back(qubitPlacement);
  }
  return placement;
}
} // namespace na::azac
