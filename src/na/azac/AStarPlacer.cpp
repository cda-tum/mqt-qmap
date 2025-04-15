#include "na/azac/AStarPlacer.hpp"

#include "ir/Definitions.hpp"
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
    double priority_;  //< sum of cost and heuristic
    const Node* node_; //< pointer to the node
    // pointer to the parent item to reconstruct the path in the end
    Item* parent_;

    Item(const double priority, const Node& node, Item* parent)
        : priority_(priority), node_(&node), parent_(parent) {
      assert(!std::isnan(priority));
    }
  };
  // compare function for the open set
  struct ItemCompare {
    auto operator()(const Item* a, const Item* b) const -> bool {
      // this way, the item with the lowest priority is on top of the heap
      return a->priority_ > b->priority_;
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
    if (isGoal(*itm->node_)) {
      // reconstruct the path from the goal to the start and then reverse it
      std::vector<std::reference_wrapper<const Node>> path;
      for (; itm != nullptr; itm = itm->parent_) {
        path.emplace_back(*itm->node_);
      }
      std::reverse(path.begin(), path.end());
      return path;
    }
    // expand the current node by adding all neighbors to the open set
    const auto& neighbors = getNeighbors(*itm->node_);
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
auto AStarPlacer::isGoal(const size_t nGates, const GateNode& node) -> bool {
  return node.level == nGates;
}
auto AStarPlacer::isGoal(const size_t nAtoms, const AtomNode& node) -> bool {
  return node.level == nAtoms;
}
auto AStarPlacer::discretizePlacementOfAtoms(
    const Placement& placement, const std::vector<qc::Qubit>& atoms) const
    -> std::pair<RowColumnMap<uint8_t>, RowColumnMap<uint8_t>> {
  std::map<size_t, RowColumnSet> rows;
  std::map<size_t, RowColumnSet> columns;
  for (const auto atom : atoms) {
    const auto& [slm, r, c] = placement[atom];
    const auto& [x, y] = architecture_.get().exactSlmLocation(slm, r, c);
    rows.try_emplace(y).first->second.emplace(std::pair{std::cref(slm), r});
    columns.try_emplace(x).first->second.emplace(std::pair{std::cref(slm), c});
  }
  RowColumnMap<uint8_t> rowIndices;
  uint8_t rowIndex = 0;
  for (const auto& [_, sites] : rows) {
    for (const auto& site : sites) {
      rowIndices.emplace(site, rowIndex);
    }
    ++rowIndex;
  }
  RowColumnMap<uint8_t> columnIndices;
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
    const SiteSet& occupiedSites) const
    -> std::pair<RowColumnMap<uint8_t>, RowColumnMap<uint8_t>> {
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
                       std::pair{std::cref(*slm), r});
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
                          std::pair{std::cref(*slm), c});
          break;
        }
      }
    }
  }
  RowColumnMap<uint8_t> rowIndices;
  uint8_t rowIndex = 0;
  for (const auto& [_, site] : rows) {
    rowIndices.emplace(site, rowIndex++);
  }
  RowColumnMap<uint8_t> columnIndices;
  uint8_t columnIndex = 0;
  for (const auto& [_, site] : columns) {
    columnIndices.emplace(site, columnIndex++);
  }
  return std::pair{rowIndices, columnIndices};
}

auto AStarPlacer::discretizeNonOccupiedEntanglementSites(
    const SiteSet& occupiedSites) const
    -> std::pair<RowColumnMap<uint8_t>, RowColumnMap<uint8_t>> {
  std::map<size_t, RowColumnSet> rows;
  std::map<size_t, RowColumnSet> columns;
  for (const auto& zone : architecture_.get().entanglementZones) {
    for (const auto& slm : *zone) {
      // find rows with free sites
      for (size_t r = 0; r < slm.nRows; ++r) {
        for (size_t c = 0; c < slm.nCols; ++c) {
          if (occupiedSites.find(std::tie(slm, r, c)) == occupiedSites.end()) {
            // free site in row r found at column c
            rows.try_emplace(slm.location.second +
                             (slm.siteSeparation.second * r))
                .first->second.emplace(std::pair{std::cref(slm), r});
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
                .first->second.emplace(std::pair{std::cref(slm), c});
            break;
          }
        }
      }
    }
  }
  RowColumnMap<uint8_t> rowIndices;
  uint8_t rowIndex = 0;
  for (const auto& [_, sites] : rows) {
    for (const auto& site : sites) {
      rowIndices.emplace(site, rowIndex);
    }
    ++rowIndex;
  }
  RowColumnMap<uint8_t> columnIndices;
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
    -> Placement {
  auto slmIt = architecture_.get().storageZones.cbegin();
  std::size_t c = 0;
  std::int64_t r = reverseInitialPlacement_
                       ? static_cast<std::int64_t>((*slmIt)->nRows) - 1
                       : 0;
  const std::int64_t step = reverseInitialPlacement_ ? -1 : 1;
  Placement initialPlacement;
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
    const Placement& previousPlacement,
    const std::unordered_set<qc::Qubit>& previousReuseQubits,
    const std::unordered_set<qc::Qubit>& reuseQubits,
    const TwoQubitGateLayer& twoQubitGates,
    const TwoQubitGateLayer& nextTwoQubitGates)
    -> std::pair<Placement, Placement> {
  const auto& gatePlacement = placeGatesInEntanglementZone(
      previousPlacement, previousReuseQubits, twoQubitGates, reuseQubits,
      nextTwoQubitGates);
  return {gatePlacement,
          placeAtomsInStorageZone(gatePlacement, reuseQubits, twoQubitGates,
                                  nextTwoQubitGates)};
}

auto AStarPlacer::addGateOption(
    const RowColumnMap<uint8_t>& discreteTargetRows,
    const RowColumnMap<uint8_t>& discreteTargetColumns, const SLM& leftSlm,
    const size_t leftRow, const size_t leftCol, const SLM& rightSlm,
    const size_t rightRow, const size_t rightCol, const SLM& nearestSlm,
    const size_t r, const size_t c, GateJob& job) const -> void {
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
            std::array{
                discreteTargetRows.at(std::pair{std::cref(nearestSlm), r}),
                discreteTargetColumns.at(std::pair{std::cref(nearestSlm), c})},
            std::array{
                discreteTargetRows.at(std::pair{std::cref(otherSlm), otherRow}),
                discreteTargetColumns.at(
                    std::pair{std::cref(otherSlm), otherCol})}},
        std::array{dis1, dis4}});
  } else {
    job.options.emplace_back(GateJob::Option{
        std::array{
            std::array{
                discreteTargetRows.at(std::pair{std::cref(otherSlm), otherRow}),
                discreteTargetColumns.at(
                    std::pair{std::cref(otherSlm), otherCol})},
            std::array{
                discreteTargetRows.at(std::pair{std::cref(nearestSlm), r}),
                discreteTargetColumns.at(std::pair{std::cref(nearestSlm), c})}},
        std::array{dis2, dis3}});
  }
}

auto AStarPlacer::placeGatesInEntanglementZone(
    const Placement& previousPlacement,
    const std::unordered_set<qc::Qubit>& reuseQubits,
    const TwoQubitGateLayer& twoQubitGates,
    const std::unordered_set<qc::Qubit>& nextReuseQubits,
    const TwoQubitGateLayer& nextTwoQubitGates) -> Placement {
  // Duplicate the previous placement as a starting point for the current
  auto currentPlacement = previousPlacement;
  //===------------------------------------------------------------------===//
  // Find gates and atoms that must be placed
  //===------------------------------------------------------------------===//
  std::set<std::pair<double, QubitPair>, std::greater<>> gatesToPlace;
  std::vector<qc::Qubit> atomsToPlace;
  for (const auto& gate : twoQubitGates) {
    const auto& [first, second] = gate;
    if (const auto firstQubitReuse =
            reuseQubits.find(first) != reuseQubits.end() &&
            std::get<0>(previousPlacement[first]).get().isEntanglement();
        !firstQubitReuse &&
        (reuseQubits.find(second) == reuseQubits.end() ||
         std::get<0>(previousPlacement[second]).get().isStorage())) {
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
  SiteSet occupiedEntanglementSites{};
  for (const auto qubit : reuseQubits) {
    const auto& [slm, r, c] = previousPlacement[qubit];
    if (slm.get().isEntanglement()) {
      occupiedEntanglementSites.emplace(slm, r, c);
      // also add the interaction partner's site to the set
      const auto& [otherSlm, otherRow, otherCol] =
          architecture_.get().otherEntanglementSite(slm, r, c);
      occupiedEntanglementSites.emplace(otherSlm, otherRow, otherCol);
    }
  }
  //===------------------------------------------------------------------===//
  // Discretize the free sites for the atoms to be placed
  //===------------------------------------------------------------------===//
  const auto& [discreteTargetRows, discreteTargetColumns] =
      discretizeNonOccupiedEntanglementSites(occupiedEntanglementSites);
  std::unordered_map<uint8_t, std::unordered_map<uint8_t, Site>> targetSites;
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
        std::array{discreteRows.at(std::pair{std::cref(leftSlm), leftRow}),
                   discreteColumns.at(std::pair{std::cref(leftSlm), leftCol})},
        std::array{
            discreteRows.at(std::pair{std::cref(rightSlm), rightRow}),
            discreteColumns.at(std::pair{std::cref(rightSlm), rightCol})}};
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
        // landscape ==> expand width and adjust height
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
    // Determine whether lookahead for the gate should be considered.
    // That is the case it the gate to be placed contains a reuse qubit
    // because then we do not only decide the position of the gate in this layer
    // but also of the gate in the next layer.
    bool leftReuse = nextReuseQubits.find(leftAtom) != nextReuseQubits.end();
    bool rightReuse = nextReuseQubits.find(rightAtom) != nextReuseQubits.end();
    qc::Qubit nextInteractionPartner = 0;
    if (leftReuse || rightReuse) {
      for (const auto& nextGate : nextTwoQubitGates) {
        const auto& [nextLeftAtom, nextRightAtom] = nextGate;
        if (leftReuse) {
          if (nextLeftAtom == leftAtom) {
            nextInteractionPartner = nextRightAtom;
            break;
          }
          if (nextRightAtom == leftAtom) {
            nextInteractionPartner = nextLeftAtom;
            break;
          }
        }
        if (rightReuse) {
          if (nextLeftAtom == rightAtom) {
            nextInteractionPartner = nextRightAtom;
            break;
          }
          if (nextRightAtom == rightAtom) {
            nextInteractionPartner = nextLeftAtom;
            break;
          }
        }
      }
      job.meanLookaheadCost = 0.0F;
      const auto& [nextSlm, nextRow, nextCol] =
          previousPlacement[nextInteractionPartner];
      for (auto& option : job.options) {
        const auto& [row, col] = option.sites[leftReuse ? 0 : 1];
        const auto& [targetSlm, targetRow, targetCol] =
            targetSites.at(row).at(col);
        const auto distance = static_cast<float>(architecture_.get().distance(
            nextSlm, nextRow, nextCol, targetSlm, targetRow, targetCol));
        option.lookaheadCost = lookaheadFactor_ * std::sqrt(distance);
        job.meanLookaheadCost += option.lookaheadCost;
      }
      job.meanLookaheadCost /= static_cast<float>(job.options.size());
    }
  }
  //===------------------------------------------------------------------===//
  // Get the extent of discrete source and target
  //===------------------------------------------------------------------===//
  assert(!discreteRows.empty()); // ==> the following std::max_element does not
                                 // return a nullptr
  const uint8_t maxDiscreteSourceRow =
      std::max_element(discreteRows.begin(), discreteRows.end(),
                       [](const auto& lhs, const auto& rhs) {
                         return lhs.second < rhs.second;
                       })
          ->second;
  assert(
      !discreteColumns.empty()); // ==> the following std::max_element does not
  // return a nullptr
  const uint8_t maxDiscreteSourceColumn =
      std::max_element(discreteColumns.begin(), discreteColumns.end(),
                       [](const auto& lhs, const auto& rhs) {
                         return lhs.second < rhs.second;
                       })
          ->second;
  assert(!discreteTargetRows
              .empty()); // ==> the following std::max_element does not
  // return a nullptr
  const uint8_t maxDiscreteTargetRow =
      std::max_element(discreteTargetRows.begin(), discreteTargetRows.end(),
                       [](const auto& lhs, const auto& rhs) {
                         return lhs.second < rhs.second;
                       })
          ->second;
  assert(!discreteTargetColumns
              .empty()); // ==> the following std::max_element does not
  // return a nullptr
  const uint8_t maxDiscreteTargetColumn =
      std::max_element(discreteTargetColumns.begin(),
                       discreteTargetColumns.end(),
                       [](const auto& lhs, const auto& rhs) {
                         return lhs.second < rhs.second;
                       })
          ->second;
  // return a nullptr
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
  const auto deepeningValue = deepeningValue_;
  const auto& path = aStarTreeSearch<GateNode>(
      *nodes.front(),
      [&nodes, &gateJobs](const auto& node) {
        return getNeighbors(nodes, gateJobs, std::move(node));
      },
      [nJobs](const auto& node) { return isGoal(nJobs, std::move(node)); },
      [](const auto& node) { return getCost(std::move(node)); },
      [&gateJobs, deepeningFactor, deepeningValue,
       &scaleFactors](const auto& node) {
        return getHeuristic(gateJobs, deepeningFactor, deepeningValue,
                            scaleFactors, std::move(node));
      });
  //===------------------------------------------------------------------===//
  // Extract the final mapping
  //===------------------------------------------------------------------===//
  assert(path.size() == nJobs + 1);
  for (size_t i = 0; i < nJobs; ++i) {
    const auto& job = gateJobs[i];
    const auto& option = job.options[path[i + 1].get().option];
    for (size_t j = 0; j < 2; ++j) {
      const auto atom = job.qubits[j];
      const auto& [row, col] = option.sites[j];
      currentPlacement[atom] = targetSites.at(row).at(col);
    }
  }
  return currentPlacement;
}

auto AStarPlacer::placeAtomsInStorageZone(
    const Placement& previousPlacement,
    const std::unordered_set<qc::Qubit>& reuseQubits,
    const TwoQubitGateLayer& twoQubitGates,
    const TwoQubitGateLayer& nextTwoQubitGates) -> Placement {
  // Duplicate the previous placement as a starting point for the current
  auto currentPlacement = previousPlacement;
  if (twoQubitGates.empty()) {
    return currentPlacement;
  }
  //===------------------------------------------------------------------===//
  // Find atoms that must be placed
  //===------------------------------------------------------------------===//
  std::vector<qc::Qubit> atomsToPlace;
  double maxDistance = 0.0;
  size_t atomIndexWithMaxDistance = 0;
  for (const auto& gate : twoQubitGates) {
    for (const auto qubit : gate) {
      const auto& [slm, r, c] = previousPlacement[qubit];
      const auto& [nearestSlm, nearestRow, nearestCol] =
          architecture_.get().nearestStorageSite(slm, r, c);
      const auto distance = architecture_.get().distance(
          slm, r, c, nearestSlm, nearestRow, nearestCol);
      if (distance > maxDistance) {
        maxDistance = distance;
        atomIndexWithMaxDistance = atomsToPlace.size();
      }
      atomsToPlace.emplace_back(qubit);
    }
  }
  //===------------------------------------------------------------------===//
  // Discretize the previous placement of the atoms to be placed
  //===------------------------------------------------------------------===//
  // Place the atoms with the longest distance first, but the place atoms in
  // increasing order of distance to the first atom
  // swap atoms at index 0 and atomIndexWithMaxDistance
  std::swap(atomsToPlace.front(), atomsToPlace.at(atomIndexWithMaxDistance));
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
  SiteSet occupiedStorageSites{};
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
  std::unordered_map<uint8_t, std::unordered_map<uint8_t, Site>> targetSites;
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
        discreteTargetColumns.at(std::pair{std::cref(nearestSlm), nearestCol});
    minDiscreteColumnOfNearestSite =
        std::min(minDiscreteColumnOfNearestSite, discreteColumnOfNearestSite);
    maxDiscreteColumnOfNearestSite =
        std::max(maxDiscreteColumnOfNearestSite, discreteColumnOfNearestSite);
    auto& job = atomJobs.emplace_back();
    job.atom = atom;
    job.currentSite = std::array{
        discreteRows.at(std::pair{std::cref(previousSlm), previousRow}),
        discreteColumns.at(std::pair{std::cref(previousSlm), previousCol})};
    if (reuseQubits.find(atom) != reuseQubits.end()) {
      // atom can be reused, so we add an option for the atom to stay at the
      // current site
      job.options.emplace_back(AtomJob::Option{{0, 0}, 0.0F, true});
    }
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
              {discreteTargetRows.at(std::pair{std::cref(nearestSlm), r}),
               discreteTargetColumns.at(std::pair{std::cref(nearestSlm), c})},
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
                {discreteTargetRows.at(
                     std::pair{std::cref(nearestSlm), rLowNew}),
                 discreteTargetColumns.at(std::pair{std::cref(nearestSlm), c})},
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
                {discreteTargetRows.at(std::pair{std::cref(nearestSlm), rHigh}),
                 discreteTargetColumns.at(std::pair{std::cref(nearestSlm), c})},
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
                {discreteTargetRows.at(std::pair{std::cref(nearestSlm), r}),
                 discreteTargetColumns.at(
                     std::pair{std::cref(nearestSlm), cLowNew})},
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
                {discreteTargetRows.at(std::pair{std::cref(nearestSlm), r}),
                 discreteTargetColumns.at(
                     std::pair{std::cref(nearestSlm), cHigh})},
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
    // Determine whether lookahead for the gate should be considered.
    // That is the case it the gate to be placed contains a reuse qubit
    // because then we do not only decide the position of the gate in this layer
    // but also of the gate in the next layer.
    for (const auto& nextGate : nextTwoQubitGates) {
      const auto& [nextLeftAtom, nextRightAtom] = nextGate;
      if (nextLeftAtom == atom || nextRightAtom == atom) {
        const qc::Qubit nextInteractionPartner =
            nextLeftAtom == atom ? nextRightAtom : nextLeftAtom;
        const auto& [nextSlm, nextRow, nextCol] =
            previousPlacement[nextInteractionPartner];
        // use this attribute to sum up the costs and calculate the mean by
        // dividing by the number of options
        job.meanLookaheadCost = 0.0F;
        for (auto& option : job.options) {
          if (option.reuse) {
            const auto distance =
                static_cast<float>(architecture_.get().distance(
                    nextSlm, nextRow, nextCol, previousSlm, previousRow,
                    previousCol));
            // NOTE: the multiplication with the lookahead factor is missing
            // here on purpose as this is the distance the "reuse" costs by
            // taking the distance of the next interaction partner as cost here
            option.lookaheadCost =
                std::max(0.0F, std::sqrt(distance) - reuseLevel_);
          } else {
            const auto& [row, col] = option.site;
            const auto& [targetSlm, targetRow, targetCol] =
                targetSites.at(row).at(col);
            const auto distance = static_cast<float>(
                architecture_.get().distance(nextSlm, nextRow, nextCol,
                                             targetSlm, targetRow, targetCol));
            option.lookaheadCost = lookaheadFactor_ * std::sqrt(distance);
          }
          job.meanLookaheadCost += option.lookaheadCost;
        }
        job.meanLookaheadCost /= static_cast<float>(job.options.size());
        break;
      }
    }
  }
  //===------------------------------------------------------------------===//
  // Get the extent of discrete source and target
  //===------------------------------------------------------------------===//
  assert(!discreteRows.empty()); // ==> the following std::max_element does not
                                 // return a nullptr
  const uint8_t maxDiscreteSourceRow =
      std::max_element(discreteRows.begin(), discreteRows.end(),
                       [](const auto& lhs, const auto& rhs) {
                         return lhs.second < rhs.second;
                       })
          ->second;
  assert(
      !discreteColumns.empty()); // ==> the following std::max_element does not
  // return a nullptr
  const uint8_t maxDiscreteSourceColumn =
      std::max_element(discreteColumns.begin(), discreteColumns.end(),
                       [](const auto& lhs, const auto& rhs) {
                         return lhs.second < rhs.second;
                       })
          ->second;
  assert(!discreteTargetRows
              .empty()); // ==> the following std::max_element does not
  // return a nullptr
  const uint8_t maxDiscreteTargetRow =
      std::max_element(discreteTargetRows.begin(), discreteTargetRows.end(),
                       [](const auto& lhs, const auto& rhs) {
                         return lhs.second < rhs.second;
                       })
          ->second;
  assert(!discreteTargetColumns
              .empty()); // ==> the following std::max_element does not
  // return a nullptr
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
  const auto deepeningValue = deepeningValue_;
  const auto& path = aStarTreeSearch<AtomNode>(
      *nodes.front(),
      [&nodes, &atomJobs](const auto& node) {
        return getNeighbors(nodes, atomJobs, std::move(node));
      },
      [nJobs](const auto& node) { return isGoal(nJobs, std::move(node)); },
      [](const auto& node) { return getCost(std::move(node)); },
      [&atomJobs, deepeningFactor, deepeningValue,
       &scaleFactors](const auto& node) {
        return getHeuristic(atomJobs, deepeningFactor, deepeningValue,
                            scaleFactors, std::move(node));
      });
  //===------------------------------------------------------------------===//
  // Extract the final mapping
  //===------------------------------------------------------------------===//
  assert(path.size() == nJobs + 1);
  for (size_t i = 0; i < nJobs; ++i) {
    const auto& job = atomJobs[i];
    const auto& option = job.options[path[i + 1].get().option];
    if (!option.reuse) {
      const auto atom = job.atom;
      const auto& [row, col] = option.site;
      currentPlacement[atom] = targetSites.at(row).at(col);
    }
  }
  return currentPlacement;
}

auto AStarPlacer::getCost(const GateNode& node) -> float {
  float cost = node.lookaheadCost;
  for (const auto d : node.maxDistancesOfPlacedAtomsPerGroup) {
    cost += std::sqrt(d);
  }
  return cost;
}

auto AStarPlacer::getCost(const AtomNode& node) -> float {
  float cost = node.lookaheadCost;
  for (const auto d : node.maxDistancesOfPlacedAtomsPerGroup) {
    cost += std::sqrt(d);
  }
  return cost;
}

auto AStarPlacer::sumStdDeviationForGroups(
    const std::array<float, 2>& scaleFactors,
    const std::vector<CompatibilityGroup>& groups) -> float {
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

auto AStarPlacer::getHeuristic(const std::vector<AtomJob>& atomJobs,
                               const float deepeningFactor,
                               const float deepeningValue,
                               const std::array<float, 2>& scaleFactors,
                               const AtomNode& node) -> float {
  const auto nAtomJobs = atomJobs.size();
  const auto nUnplacedAtoms = static_cast<float>(nAtomJobs - node.level);
  float maxDistanceOfUnplacedAtom = 0.0F;
  float accMinLookaheadCost = 0.0F;
  for (size_t i = node.level; i < nAtomJobs; ++i) {
    const auto& job = atomJobs[i];
    accMinLookaheadCost += job.meanLookaheadCost;
    for (const auto& option : job.options) {
      if (option.reuse) {
        // since them jobs are sorted by distance and this has distance 0 as the
        // atom does not need to be moved at all, this job will always be the
        // first one for atoms that may be reused
        break;
      }
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
  float maxDistanceOfPlacedAtom = 0.0F;
  for (const auto d : node.maxDistancesOfPlacedAtomsPerGroup) {
    maxDistanceOfPlacedAtom = std::max(maxDistanceOfPlacedAtom, d);
  }
  float heuristic = maxDistanceOfUnplacedAtom <= maxDistanceOfPlacedAtom
                        ? 0.F
                        : std::sqrt(maxDistanceOfUnplacedAtom) -
                              std::sqrt(maxDistanceOfPlacedAtom);
  heuristic += accMinLookaheadCost;
  heuristic +=
      deepeningFactor *
      (sumStdDeviationForGroups(scaleFactors, node.groups) + deepeningValue) *
      nUnplacedAtoms;
  return heuristic;
}

auto AStarPlacer::getHeuristic(const std::vector<GateJob>& gateJobs,
                               const float deepeningFactor,
                               const float deepeningValue,
                               const std::array<float, 2>& scaleFactors,
                               const GateNode& node) -> float {
  const auto nGateJobs = gateJobs.size();
  const auto nUnplacedGates = static_cast<float>(nGateJobs - node.level);
  float maxDistanceOfUnplacedAtom = 0.0F;
  float accMeanLookaheadCost = 0.0F;
  for (size_t i = node.level; i < nGateJobs; ++i) {
    const auto& job = gateJobs[i];
    accMeanLookaheadCost += job.meanLookaheadCost;
    for (const auto& option : job.options) {
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
  float maxDistanceOfPlacedAtom = 0.0F;
  for (const auto d : node.maxDistancesOfPlacedAtomsPerGroup) {
    maxDistanceOfPlacedAtom = std::max(maxDistanceOfPlacedAtom, d);
  }
  float heuristic = maxDistanceOfUnplacedAtom <= maxDistanceOfPlacedAtom
                        ? 0.F
                        : std::sqrt(maxDistanceOfUnplacedAtom) -
                              std::sqrt(maxDistanceOfPlacedAtom);
  heuristic += accMeanLookaheadCost;
  heuristic +=
      deepeningFactor *
      (sumStdDeviationForGroups(scaleFactors, node.groups) + deepeningValue) *
      nUnplacedGates;
  return heuristic;
}

auto AStarPlacer::getNeighbors(std::deque<std::unique_ptr<AtomNode>>& nodes,
                               const std::vector<AtomJob>& atomJobs,
                               const AtomNode& node)
    -> std::vector<std::reference_wrapper<const AtomNode>> {
  const size_t atomToBePlacedNext = node.level;
  const auto& atomJob = atomJobs[atomToBePlacedNext];
  std::vector<std::reference_wrapper<const AtomNode>> neighbors;
  assert(atomJob.options.size() <= std::numeric_limits<uint16_t>::max());
  for (uint16_t i = 0; i < static_cast<uint16_t>(atomJob.options.size()); ++i) {
    const auto& option = atomJob.options[i];
    const auto& [site, distance, reuse, lookaheadCost] = option;
    // skip the sites that are already consumed
    if (!reuse &&
        node.consumedFreeSites.find(site) != node.consumedFreeSites.end()) {
      continue;
    }
    // make a copy of node, the parent of child
    AtomNode& child = *nodes.emplace_back(std::make_unique<AtomNode>(node));
    if (!reuse) {
      child.consumedFreeSites.emplace(site);
      // check whether the current placement is compatible with any existing
      // group
      checkCompatibilityAndAddPlacement(
          atomJob.currentSite.front(), site.front(), atomJob.currentSite.back(),
          site.back(), distance, child.groups,
          child.maxDistancesOfPlacedAtomsPerGroup);
    }
    child.option = i;
    ++child.level;
    child.lookaheadCost += lookaheadCost;
    // add the child to the list of children to be returned
    neighbors.emplace_back(child);
  }
  return neighbors;
}

auto AStarPlacer::getNeighbors(std::deque<std::unique_ptr<GateNode>>& nodes,
                               const std::vector<GateJob>& gateJobs,
                               const GateNode& node)
    -> std::vector<std::reference_wrapper<const GateNode>> {
  const size_t gateToBePlacedNext = node.level;
  const auto& gateJob = gateJobs[gateToBePlacedNext];
  std::vector<std::reference_wrapper<const GateNode>> neighbors;
  // Get the current placement of the atoms that must be placed next
  const auto& [currentSiteOfLeftAtom, currentSiteOfRightAtom] =
      gateJob.currentSites;
  assert(gateJob.options.size() <= std::numeric_limits<uint16_t>::max());
  for (uint16_t i = 0; i < static_cast<uint16_t>(gateJob.options.size()); ++i) {
    const auto& option = gateJob.options[i];
    const auto& [sites, distances, lookaheadCost] = option;
    const auto& [leftSite, rightSite] = sites;
    // skip if one of the sites is already consumed
    if (node.consumedFreeSites.find(leftSite) != node.consumedFreeSites.end() ||
        node.consumedFreeSites.find(rightSite) !=
            node.consumedFreeSites.end()) {
      continue;
    }
    // make a copy of node, the parent of child as use this as a starting
    // point for the new node
    GateNode& child = *nodes.emplace_back(std::make_unique<GateNode>(node));
    ++child.level;
    child.option = i;
    child.consumedFreeSites.emplace(leftSite);
    child.consumedFreeSites.emplace(rightSite);
    // check whether the current placement is compatible with any existing
    // horizontal group
    checkCompatibilityAndAddPlacement(
        currentSiteOfLeftAtom.front(), leftSite.front(),
        currentSiteOfLeftAtom.back(), leftSite.back(), distances.front(),
        child.groups, child.maxDistancesOfPlacedAtomsPerGroup);
    checkCompatibilityAndAddPlacement(
        currentSiteOfRightAtom.front(), rightSite.front(),
        currentSiteOfRightAtom.back(), rightSite.back(), distances.back(),
        child.groups, child.maxDistancesOfPlacedAtomsPerGroup);
    child.lookaheadCost += lookaheadCost;
    // add the final child to the list of children to be returned
    neighbors.emplace_back(child);
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
        if (const auto lowerValue = std::prev(it)->second;
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
    if (const auto lowerValue = std::prev(it)->second; lowerValue < value) {
      // new placement is compatible with this group
      return std::pair{it, false};
    }
  }
  return std::nullopt;
}

auto AStarPlacer::checkCompatibilityAndAddPlacement(
    const uint8_t hKey, const uint8_t hValue, const uint8_t vKey,
    const uint8_t vValue, const float distance,
    std::vector<CompatibilityGroup>& groups, std::vector<float>& maxDistances)
    -> bool {
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
    bool deepeningValueSet = false;
    bool lookaheadFactorSet = false;
    bool reuseLevelSet = false;
    for (const auto& [key, value] : configIt.value().items()) {
      if (key == "use_window") {
        useWindowSet = true;
        if (value.is_boolean()) {
          useWindow_ = value;
        } else {
          std::ostringstream oss;
          oss << std::boolalpha;
          oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer "
                 "contains an invalid "
                 "value for use_window. Using default ("
              << useWindow_ << ").\n";
          std::cout << oss.str();
        }
      } else if (key == "window_min_width") {
        windowMinWidthSet = true;
        if (value.is_number_unsigned()) {
          windowMinWidth_ = value;
        } else {
          std::ostringstream oss;
          oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer "
                 "contains an invalid value for window_min_width. Using "
                 "default ("
              << windowMinWidth_ << ").\n";
          std::cout << oss.str();
        }
      } else if (key == "window_ratio") {
        windowRatioSet = true;
        if (value.is_number()) {
          windowRatio_ = value;
        } else {
          std::ostringstream oss;
          oss << std::setprecision(4);
          oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer "
                 "contains an invalid value for window_ratio. Using default ("
              << windowRatio_ << ").\n";
          std::cout << oss.str();
        }
      } else if (key == "window_share") {
        windowShareSet = true;
        if (value.is_number()) {
          windowShare_ = value;
        } else {
          std::ostringstream oss;
          oss << std::setprecision(4);
          oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer "
                 "contains an invalid value for window_share. Using default ("
              << windowShare_ << ").\n";
          std::cout << oss.str();
        }
      } else if (key == "deepening_factor") {
        deepeningFactorSet = true;
        if (value.is_number()) {
          deepeningFactor_ = value;
        } else {
          std::ostringstream oss;
          oss << std::setprecision(4);
          oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer "
                 "contains an invalid value for deepening_factor. Using "
                 "default ("
              << deepeningFactor_ << ").\n";
          std::cout << oss.str();
        }
      } else if (key == "deepening_value") {
        deepeningValueSet = true;
        if (value.is_number()) {
          deepeningValue_ = value;
        } else {
          std::ostringstream oss;
          oss << std::setprecision(4);
          oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer "
                 "contains an invalid value for deepening_value. Using "
                 "default ("
              << deepeningValue_ << ").\n";
          std::cout << oss.str();
        }
      } else if (key == "lookahead_factor") {
        lookaheadFactorSet = true;
        if (value.is_number()) {
          lookaheadFactor_ = value;
        } else {
          std::ostringstream oss;
          oss << std::setprecision(4);
          oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer "
                 "contains an invalid value for lookahead_factor. Using "
                 "default ("
              << lookaheadFactor_ << ").\n";
          std::cout << oss.str();
        }
      } else if (key == "reuse_level") {
        reuseLevelSet = true;
        if (value.is_number()) {
          reuseLevel_ = value;
        } else {
          std::ostringstream oss;
          oss << std::setprecision(4);
          oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer "
                 "contains an invalid value for reuse_level. Using "
                 "default ("
              << reuseLevel_ << ").\n";
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
      std::ostringstream oss;
      oss << std::boolalpha;
      oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer does "
             "not contain a setting for use_window. Using default ("
          << useWindow_ << ").\n";
      std::cout << oss.str();
    }
    if (useWindow_) {
      if (!windowMinWidthSet) {
        std::ostringstream oss;
        oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer "
               "does not contain a setting for window_min_width. Using "
               "default ("
            << windowMinWidth_ << ").\n";
        std::cout << oss.str();
      }
      if (!windowRatioSet) {
        std::ostringstream oss;
        oss << std::setprecision(4);
        oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer "
               "does not contain a setting for window_ratio. Using default ("
            << windowRatio_ << ").\n";
        std::cout << oss.str();
      }
      if (!windowShareSet) {
        std::ostringstream oss;
        oss << std::setprecision(4);
        oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer "
               "does not contain a setting for window_share. Using default ("
            << windowShare_ << ").\n";
        std::cout << oss.str();
      }
    }
    if (windowMinWidthSet) {
      windowMinHeight_ = static_cast<size_t>(
          std::round(windowRatio_ * static_cast<double>(windowMinWidth_)));
    }
    if (!deepeningFactorSet) {
      std::ostringstream oss;
      oss << std::setprecision(4);
      oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer does "
             "not contain a setting for deepening_factor. Using default ("
          << deepeningFactor_ << ").\n";
      std::cout << oss.str();
    }
    if (!deepeningValueSet) {
      std::ostringstream oss;
      oss << std::setprecision(4);
      oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer does "
             "not contain a setting for deepening_value. Using default ("
          << deepeningValue_ << ").\n";
      std::cout << oss.str();
    }
    if (!lookaheadFactorSet) {
      std::ostringstream oss;
      oss << std::setprecision(4);
      oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer does "
             "not contain a setting for lookahead_factor. Using default ("
          << lookaheadFactor_ << ").\n";
      std::cout << oss.str();
    }
    if (!reuseLevelSet) {
      std::ostringstream oss;
      oss << std::setprecision(4);
      oss << "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer does "
             "not contain a setting for reuse_level. Using default ("
          << reuseLevel_ << ").\n";
      std::cout << oss.str();
    }
  } else {
    std::ostringstream oss;
    oss << std::boolalpha << std::setprecision(4);
    oss << "\033[1;35m[WARN]\033[0m Configuration does not contain "
           "settings for AStarPlacer or is malformed. Using default settings "
           "("
        << "\"use_window\" :" << useWindow_
        << ", \"window_min_width\" :" << windowMinWidth_
        << ", \"window_ratio\" :" << windowRatio_
        << ", \"window_share\" :" << windowShare_
        << ", \"deepening_factor\" :" << deepeningFactor_
        << ", \"deepening_value\" :" << deepeningValue_
        << ", \"lookahead_factor\" :" << lookaheadFactor_
        << ", \"reuse_level\" :" << reuseLevel_ << ").\n";
    std::cout << oss.str();
  }
}

auto AStarPlacer::place(
    const size_t nQubits,
    const std::vector<TwoQubitGateLayer>& twoQubitGateLayers,
    const std::vector<std::unordered_set<qc::Qubit>>& reuseQubits)
    -> std::vector<Placement> {
  std::vector<Placement> placement;
  placement.reserve((2 * twoQubitGateLayers.size()) + 1);
  placement.emplace_back(makeInitialPlacement(nQubits));
  for (size_t layer = 0; layer < twoQubitGateLayers.size(); ++layer) {
    const auto& [gatePlacement, qubitPlacement] = makeIntermediatePlacement(
        placement.back(),
        layer == 0 ? std::unordered_set<qc::Qubit>{} : reuseQubits[layer - 1],
        layer == reuseQubits.size() ? std::unordered_set<qc::Qubit>{}
                                    : reuseQubits[layer],
        twoQubitGateLayers[layer],
        layer == twoQubitGateLayers.size() - 1 ? TwoQubitGateLayer{}
                                               : twoQubitGateLayers[layer + 1]);
    placement.emplace_back(gatePlacement);
    placement.emplace_back(qubitPlacement);
  }
  return placement;
}
} // namespace na::azac
