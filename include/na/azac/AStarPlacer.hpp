#pragma once

#include "Definitions.hpp"
#include "na/azac/Architecture.hpp"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <nlohmann/json_fwd.hpp>
#include <tuple>
#include <unordered_set>
#include <vector>

namespace na {
class AStarPlacer {
  std::reference_wrapper<const Architecture> architecture_;
  bool reverseInitialPlacement = false;

  auto makeInitialPlacement(const size_t nQubits) const
      -> std::vector<std::tuple<const SLM&, size_t, size_t>> {
    auto slmIt = architecture_.get().storageZones.cbegin();
    std::size_t c = 0;
    std::int64_t r = reverseInitialPlacement
                         ? static_cast<std::int64_t>((*slmIt)->nRows) - 1
                         : 0;
    const std::int64_t step = reverseInitialPlacement ? -1 : 1;
    std::vector<std::tuple<const SLM&, size_t, size_t>> initialPlacement;
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

  /// This function places the qubits corresponding to gates in the entanglement
  /// zone. After this placement has been performed, the activation of the
  /// Rydberg beam will execute the gates in the given layer. Afterward, the
  /// next placement for moving (non-reuse) qubits back to the storage zone is
  /// determined by @ref placeQubitsInStorageZone.
  auto makeIntermediatePlacement(
      const std::vector<std::tuple<const SLM&, size_t, size_t>>&
          previousPlacement,
      const std::unordered_set<qc::Qubit>& reuseQubits,
      const std::vector<std::pair<qc::Qubit, qc::Qubit>>& twoQubitGates) const
      -> std::vector<std::tuple<const SLM&, size_t, size_t>> {
    // Duplicate the previous placement as a starting point for the current
    std::vector<std::tuple<const SLM&, size_t, size_t>> currentPlacement =
        previousPlacement;
    //===------------------------------------------------------------------===//
    // Find gates that must be placed
    //===------------------------------------------------------------------===//
    std::vector<std::pair<std::pair<size_t, size_t>, double>> gatesToPlace;
    std::vector<std::pair<size_t, double>> atomsToPlace;
    for (const auto& gate : twoQubitGates) {
      if (reuseQubits.find(gate.first) != reuseQubits.end()) {
        assert(std::get<0>(previousPlacement[gate.first]).isEntanglement());
        currentPlacement[gate.second] =
            std::apply(architecture_.get().otherEntanglementSite,
                       previousPlacement[gate.first]);
      } else if (reuseQubits.find(gate->second) != reuseQubits.end()) {
        assert(std::get<0>(previousPlacement[gate->second])->isEntanglement());
        currentPlacement[gate->first] =
            architecture.otherEntanglementSite(previousPlacement[gate->second]);
      } else {
        const auto& nearestEntanglementSite =
            architecture.nearestEntanglementSite(
                previousPlacement[gate->first],
                previousPlacement[gate->second]);
        const auto& otherNearestEntanglementSite =
            architecture.otherEntanglementSite(nearestEntanglementSite);
        //===------------------------------------------------------------------===//
        // Calculate the various distances to the entanglement sites
        //===------------------------------------------------------------------===//
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
        const auto dis1 = architecture.distance(previousPlacement[gate->first],
                                                nearestEntanglementSite);
        const auto dis2 = architecture.distance(previousPlacement[gate->second],
                                                nearestEntanglementSite);
        const auto dis3 = architecture.distance(previousPlacement[gate->first],
                                                otherNearestEntanglementSite);
        const auto dis4 = architecture.distance(previousPlacement[gate->second],
                                                otherNearestEntanglementSite);
        // If the situation is as depicted in the example above
        if (dis1 + dis4 <= dis2 + dis3) {
          gatesToPlace.emplace_back(
              std::pair{atomsToPlace.size(), atomsToPlace.size() + 1},
              std::max(dis1, dis4));
          atomsToPlace.emplace_back(gate->first, dis1);
          atomsToPlace.emplace_back(gate->second, dis2);
        } else {
          // otherwise, either the entanglement sites or storage sites are
          // flipped
          gatesToPlace.emplace_back(
              std::pair{atomsToPlace.size(), atomsToPlace.size() + 1},
              std::max(dis2, dis3));
          atomsToPlace.emplace_back(gate->first, dis2);
          atomsToPlace.emplace_back(gate->second, dis3);
        }
      }
    }
    nAtoms = atomsToPlace.size();
    //===------------------------------------------------------------------===//
    // Order the gates to be placed by the distance to the nearest entanglement
    // site
    //===------------------------------------------------------------------===//
    std::sort(gatesToPlace.begin(), gatesToPlace.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });
    //===------------------------------------------------------------------===//
    // Extract occupied entanglement sites from the previous placement
    //===------------------------------------------------------------------===//
    // This set will only contain of the first SLM in a pair of entanglement
    // SLMs and represents the occupied pair of sites
    std::set<std::tuple<const std::vector<std::unique_ptr<SLM>>*, std::size_t,
                        std::size_t>>
        occupiedEntanglementSites{};
    for (const auto qubit : reuseQubits) {
      assert(std::get<0>(previousPlacement[qubit])->isEntanglement());
      occupiedEntanglementSites.emplace(
          std::get<0>(previousPlacement[qubit])->entanglementZone,
          std::get<1>(previousPlacement[qubit]),
          std::get<2>(previousPlacement[qubit]));
    }
    //===------------------------------------------------------------------===//
    // Discretize the previous placement of the qubits to be placed
    //===------------------------------------------------------------------===//
    std::vector<std::pair<size_t, std::pair<const SLM*, size_t>>>
        atomsToPlaceWithRow;
    std::vector<std::pair<size_t, std::pair<const SLM*, size_t>>>
        atomsToPlaceWithColumn;
    atomsToPlaceWithRow.reserve(nAtoms);
    atomsToPlaceWithColumn.reserve(nAtoms);
    for (size_t atom = 0; atom < nAtoms; ++atom) {
      const auto& [slm, r, c] = previousPlacement[atomsToPlace[atom].first];
      atomsToPlaceWithRow.emplace_back(atom, std::pair{slm, r});
      atomsToPlaceWithColumn.emplace_back(atom, std::pair{slm, c});
    }
    std::sort(
        atomsToPlaceWithRow.begin(), atomsToPlaceWithRow.end(),
        [](const std::pair<size_t, std::pair<const SLM*, size_t>>& a,
           const std::pair<size_t, std::pair<const SLM*, size_t>>& b) {
          return a.second.first->location.second +
                     a.second.first->siteSeparation.second * a.second.second <
                 b.second.first->location.second +
                     b.second.first->siteSeparation.second * b.second.second;
        });
    std::sort(
        atomsToPlaceWithColumn.begin(), atomsToPlaceWithColumn.end(),
        [](const std::pair<size_t, std::pair<const SLM*, size_t>>& a,
           const std::pair<size_t, std::pair<const SLM*, size_t>>& b) {
          return a.second.first->location.first +
                     a.second.first->siteSeparation.first * a.second.second <
                 b.second.first->location.first +
                     b.second.first->siteSeparation.first * b.second.second;
        });
    //===------------------------------------------------------------------===//
    // Discretize the free sites for the atoms to be placed
    //===------------------------------------------------------------------===//
    std::vector<std::pair<std::pair<const SLM*, size_t>, size_t>>
        discreteTargetRows;
    std::vector<std::pair<std::pair<const SLM*, size_t>, size_t>>
        discreteTargetColumns;
    for (const auto& entangleSlms : architecture.entanglementZones) {
      for (const auto& storageSlm : entangleSlms) {
        // find rows with free sites
        for (size_t r = 0; r < storageSlm->nRows; ++r) {
          for (size_t c = 0; c < storageSlm->nCols; ++c) {
            if (occupiedEntanglementSites.find(
                    {storageSlm->entanglementZone, r, c}) ==
                occupiedEntanglementSites.end()) {
              // free site in row r found at column c
              discreteTargetRows.emplace_back(
                  std::pair{storageSlm.get(), r},
                  storageSlm->location.second +
                      (storageSlm->siteSeparation.second * r));
              break;
            }
          }
        }
        // find columns with free sites
        for (size_t c = 0; c < storageSlm->nCols; ++c) {
          for (size_t r = 0; r < storageSlm->nRows; ++r) {
            if (occupiedEntanglementSites.find(
                    {storageSlm->entanglementZone, r, c}) ==
                occupiedEntanglementSites.end()) {
              // free site in column c found at row r
              discreteTargetColumns.emplace_back(
                  std::pair{storageSlm.get(), c},
                  storageSlm->location.first +
                      (storageSlm->siteSeparation.first * c));
              break;
            }
          }
        }
      }
    }
    std::sort(discreteTargetRows.begin(), discreteTargetRows.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });
    std::sort(discreteTargetColumns.begin(), discreteTargetColumns.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });
    //===------------------------------------------------------------------===//
    // Provide the above mapping to SLMs and row/columns from their discrete
    // indices in the reverse direction
    //===------------------------------------------------------------------===//
    std::unordered_map<const SLM*, std::unordered_map<size_t, size_t>>
        rowsWithFreeSites;
    std::unordered_map<const SLM*, std::unordered_map<size_t, size_t>>
        columnsWithFreeSites;
    for (size_t r = 0; r < discreteTargetRows.size(); ++r) {
      if (rowsWithFreeSites.find(discreteTargetRows[r].first.first) ==
          rowsWithFreeSites.end()) {
        rowsWithFreeSites.emplace(discreteTargetRows[r].first.first,
                                  std::unordered_map<size_t, size_t>{});
      }
      rowsWithFreeSites[discreteTargetRows[r].first.first].emplace(
          discreteTargetRows[r].first.second, r);
    }
    for (size_t c = 0; c < discreteTargetColumns.size(); ++c) {
      if (columnsWithFreeSites.find(discreteTargetColumns[c].first.first) ==
          columnsWithFreeSites.end()) {
        columnsWithFreeSites.emplace(discreteTargetColumns[c].first.first,
                                     std::unordered_map<size_t, size_t>{});
      }
      columnsWithFreeSites[discreteTargetColumns[c].first.first].emplace(
          discreteTargetColumns[c].first.second, c);
    }
    //===------------------------------------------------------------------===//
    // Initialize the nearest free sites for each atom
    //===------------------------------------------------------------------===//
    const bool useWindowedPlacement =
        static_cast<T*>(this)->useWindowedPlacement;
    const size_t windowHeight = static_cast<T*>(this)->windowHeight;
    const size_t windowWidth = static_cast<T*>(this)->windowWidth;
    nearestFreeSitesForEachAtom.clear();
    nearestFreeSitesForEachAtom.reserve(atomsToPlace.size());
    for (const auto& [gate, _] : gatesToPlace) {
      const auto& [atom1, atom2] = gate;
      const auto& [nearestSLM, nearestRow, nearestCol] =
          architecture.nearestEntanglementSite(
              previousPlacement[atomsToPlace[atom1].first],
              previousPlacement[atomsToPlace[atom2].first]);
      auto& nearestFreeSitesForGate =
          nearestFreeSitesForEachGate.emplace_back();
      size_t rLow = 0;
      size_t rHigh = nearestSLM->nRows;
      size_t cLow = 0;
      size_t cHigh = nearestSLM->nCols;
      if (useWindowedPlacement) {
        rLow =
            nearestRow > windowHeight / 2 ? nearestRow - (windowHeight / 2) : 0;
        rHigh =
            std::max(nearestRow + (windowHeight / 2) + 1, nearestSLM->nRows);
        cLow =
            nearestCol > windowWidth / 2 ? nearestCol - (windowWidth / 2) : 0;
        cHigh = std::max(nearestCol + (windowWidth / 2) + 1, nearestSLM->nCols);
      }
      for (size_t r = rLow; r < rHigh; ++r) {
        for (size_t c = cLow; c < cHigh; ++c) {
          if (occupiedEntanglementSites.find(
                  {nearestSLM->entanglementZone, r, c}) ==
              occupiedEntanglementSites.end()) {
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
                architecture.otherEntanglementSite(*nearestSLM, r, c);
            const auto dis1 = architecture.distance(
                previousPlacement[atomsToPlace[atom1].first],
                {nearestSLM, r, c});
            const auto dis2 = architecture.distance(
                previousPlacement[atomsToPlace[atom2].first],
                {nearestSLM, r, c});
            const auto dis3 = architecture.distance(
                previousPlacement[atomsToPlace[atom1].first],
                {otherSlm, otherRow, otherCol});
            const auto dis4 = architecture.distance(
                previousPlacement[atomsToPlace[atom2].first],
                {otherSlm, otherRow, otherCol});
            if (dis1 + dis4 <= dis2 + dis3) {
              nearestFreeSitesForGate.emplace_back(
                  std::pair{
                      std::pair{rowsWithFreeSites[nearestSLM][r],
                                columnsWithFreeSites[nearestSLM][c]},
                      std::pair{rowsWithFreeSites[otherSlm][otherRow],
                                columnsWithFreeSites[otherSlm][otherCol]}},
                  std::max(dis1, dis4));
            } else {
              nearestFreeSitesForGate.emplace_back(
                  std::pair{std::pair{rowsWithFreeSites[otherSlm][otherRow],
                                      columnsWithFreeSites[otherSlm][otherCol]},
                            std::pair{rowsWithFreeSites[nearestSLM][r],
                                      columnsWithFreeSites[nearestSLM][c]}},
                  std::max(dis2, dis3));
            }
          }
        }
      }
    }
    //===------------------------------------------------------------------===//
    // Run the A* algorithm
    //===------------------------------------------------------------------===//
    nodes.clear();
    nodes.emplace_back(0.0, std::unordered_set<std::pair<size_t, size_t>>{},
                       std::vector<std::map<size_t, size_t>>{},
                       std::vector<double>{},
                       std::vector<std::map<size_t, size_t>>{},
                       std::vector<double>{}); // root node
    const Node& finalNode =
        *aStarTreeSearch(*nodes.front(), &getGatePlacementNeighbors, &isGoal,
                         &getCost, &getHeuristic)
             .back();
    //===------------------------------------------------------------------===//
    // Extract the final mapping
    //===------------------------------------------------------------------===//
    // the following two maps combine all hGroups and vGroups of the final node,
    // respectively, and hence serve as a mapping from rows and columns from the
    // previous placement to rows and columns of the next placement,
    // respectively.
    std::unordered_map<size_t, size_t> rowMapping;
    std::unordered_map<size_t, size_t> columnMapping;
    for (const auto& hGroup : finalNode.hGroups) {
      rowMapping.insert(hGroup.begin(), hGroup.end());
    }
    for (const auto& vGroup : finalNode.vGroups) {
      columnMapping.insert(vGroup.begin(), vGroup.end());
    }
    std::vector<const SLM*> finalSLMs;
    std::vector<size_t> finalRows;
    std::vector<size_t> finalColumns;
    for (size_t i = 0; i < nAtoms; ++i) {
      const auto& [startRow, startColumn] = currentSitesForEachAtom[i];
      const auto targetRow = rowMapping.at(startRow);
      const auto targetColumn = columnMapping.at(startColumn);
      finalSLMs.emplace_back(discreteTargetRows.at(targetRow).first.first);
      assert(discreteTargetRows.at(targetRow).first ==
             discreteTargetColumns.at(targetColumn).first);
      finalRows.emplace_back(discreteTargetRows.at(targetRow).second);
      finalColumns.emplace_back(discreteTargetColumns.at(targetColumn).second);
    }
    //===------------------------------------------------------------------===//
    // Store the final mapping
    //===------------------------------------------------------------------===//
    for (size_t i = 0; i < nAtoms; ++i) {
      const auto atom = atomsToPlace[i].first;
      currentPlacement[atom] =
          std::tuple{finalSLMs[i], finalRows[i], finalColumns[i]};
    }
  }

  /// This function places qubits from the entanglement zone in the storage
  /// zone after a rydberg gate has been performed.
  /// It initializes the graph structure for the A* algorithm.
  /// Afterward, the A* algorithm is called to find the optimal mapping.
  auto placeQubitsInStorageZone(const size_t layer) -> void {
    //===------------------------------------------------------------------===//
    // Retrieve references to required data structures
    //===------------------------------------------------------------------===//
    const Architecture& architecture = static_cast<T*>(this)->getArchitecture();
    // placement of atoms in the previous stage, i.e., when the last gates were
    // executed
    const std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>&
        previousPlacement = static_cast<T*>(this)->getQubitMapping().back();
    // gates that were executed in the previous stage
    const std::vector<const std::pair<qc::Qubit, qc::Qubit>*>& gates =
        static_cast<T*>(this)->getGateScheduling().at(layer);
    // qubits that are reused in the next stage and hence remain at their last
    // position
    const std::unordered_set<std::size_t>& reuseQubits =
        static_cast<T*>(this)->getReuseQubits().at(layer);
    //===------------------------------------------------------------------===//
    // Extract occupied storage sites from the previous placement
    //===------------------------------------------------------------------===//
    std::set<std::tuple<const SLM*, std::size_t, std::size_t>>
        occupiedStorageSites{};
    for (const auto& [slm, r, c] : previousPlacement) {
      if (slm->isStorage()) {
        occupiedStorageSites.emplace(slm, r, c);
      }
    }
    //===------------------------------------------------------------------===//
    // Get the atoms that must be placed with the distance to their nearest free
    // site
    //===------------------------------------------------------------------===//
    std::vector<std::pair<size_t, double>> atomsToPlace;
    atomsToPlace.reserve((gates.size() * 2) - reuseQubits.size());
    for (const auto* atoms : gates) {
      if (const auto atom = atoms->first;
          reuseQubits.find(atom) == reuseQubits.end()) {
        const auto& currentSite = previousPlacement[atom];
        const auto& nearestSite = architecture.nearestStorageSite(currentSite);
        atomsToPlace.emplace_back(
            atom, architecture.distance(currentSite, nearestSite));
      }
      if (const auto atom = atoms->second;
          reuseQubits.find(atom) == reuseQubits.end()) {
        const auto& currentSite = previousPlacement[atom];
        const auto& nearestSite = architecture.nearestStorageSite(currentSite);
        atomsToPlace.emplace_back(
            atom, architecture.distance(currentSite, nearestSite));
      }
    }
    nAtoms = atomsToPlace.size();
    //===------------------------------------------------------------------===//
    // Order the atoms to be placed by the distance to the nearest free site
    //===------------------------------------------------------------------===//
    std::sort(atomsToPlace.begin(), atomsToPlace.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });
    //===------------------------------------------------------------------===//
    // Discretize the previous placement of the atoms to be placed
    //===------------------------------------------------------------------===//
    std::vector<std::pair<size_t, std::pair<const SLM*, size_t>>>
        atomsToPlaceWithRow;
    std::vector<std::pair<size_t, std::pair<const SLM*, size_t>>>
        atomsToPlaceWithColumn;
    atomsToPlaceWithRow.reserve(nAtoms);
    atomsToPlaceWithColumn.reserve(nAtoms);
    for (size_t atom = 0; atom < nAtoms; ++atom) {
      const auto& [slm, r, c] = previousPlacement[atomsToPlace[atom].first];
      atomsToPlaceWithRow.emplace_back(atom, std::pair{slm, r});
      atomsToPlaceWithColumn.emplace_back(atom, std::pair{slm, c});
    }
    std::sort(
        atomsToPlaceWithRow.begin(), atomsToPlaceWithRow.end(),
        [](const std::pair<size_t, std::pair<const SLM*, size_t>>& a,
           const std::pair<size_t, std::pair<const SLM*, size_t>>& b) {
          return a.second.first->location.second +
                     a.second.first->siteSeparation.second * a.second.second <
                 b.second.first->location.second +
                     b.second.first->siteSeparation.second * b.second.second;
        });
    std::sort(
        atomsToPlaceWithColumn.begin(), atomsToPlaceWithColumn.end(),
        [](const std::pair<size_t, std::pair<const SLM*, size_t>>& a,
           const std::pair<size_t, std::pair<const SLM*, size_t>>& b) {
          return a.second.first->location.first +
                     a.second.first->siteSeparation.first * a.second.second <
                 b.second.first->location.first +
                     b.second.first->siteSeparation.first * b.second.second;
        });
    std::vector<std::pair<const SLM*, size_t>> startRows;
    std::vector<std::pair<const SLM*, size_t>> startColumns;
    std::vector<size_t> atomsStartRow(nAtoms);
    std::vector<size_t> atomsStartColumn(nAtoms);
    for (const auto& [atom, row] : atomsToPlaceWithRow) {
      if (startRows.empty() || startRows.back().first != row.first ||
          startRows.back().second != row.second) {
        startRows.emplace_back(row);
      }
      atomsStartRow[atom] = startRows.size() - 1;
    }
    for (const auto& [atom, column] : atomsToPlaceWithColumn) {
      if (startColumns.empty() || startColumns.back().first != column.first ||
          startColumns.back().second != column.second) {
        startColumns.emplace_back(column);
      }
      atomsStartColumn[atom] = startColumns.size() - 1;
    }
    currentSitesForEachAtom.clear();
    currentSitesForEachAtom.reserve(nAtoms);
    for (size_t atom = 0; atom < nAtoms; ++atom) {
      currentSitesForEachAtom.emplace_back(atomsStartRow[atom],
                                           atomsStartColumn[atom]);
    }
    //===------------------------------------------------------------------===//
    // Discretize the free sites for the atoms to be placed
    //===------------------------------------------------------------------===//
    std::vector<std::pair<std::pair<const SLM*, size_t>, size_t>>
        discreteTargetRows;
    std::vector<std::pair<std::pair<const SLM*, size_t>, size_t>>
        discreteTargetColumns;
    for (const auto& storageSlm : architecture.storageZones) {
      // find rows with free sites
      for (size_t r = 0; r < storageSlm->nRows; ++r) {
        for (size_t c = 0; c < storageSlm->nCols; ++c) {
          if (occupiedStorageSites.find({storageSlm.get(), r, c}) ==
              occupiedStorageSites.end()) {
            // free site in row r found at column c
            discreteTargetRows.emplace_back(
                std::pair{storageSlm.get(), r},
                storageSlm->location.second +
                    (storageSlm->siteSeparation.second * r));
            break;
          }
        }
      }
      // find columns with free sites
      for (size_t c = 0; c < storageSlm->nCols; ++c) {
        for (size_t r = 0; r < storageSlm->nRows; ++r) {
          if (occupiedStorageSites.find({storageSlm.get(), r, c}) ==
              occupiedStorageSites.end()) {
            // free site in column c found at row r
            discreteTargetColumns.emplace_back(
                std::pair{storageSlm.get(), c},
                storageSlm->location.first +
                    (storageSlm->siteSeparation.first * c));
            break;
          }
        }
      }
    }
    std::sort(discreteTargetRows.begin(), discreteTargetRows.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });
    std::sort(discreteTargetColumns.begin(), discreteTargetColumns.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });
    //===------------------------------------------------------------------===//
    // Provide the above mapping to SLMs and row/columns from their discrete
    // indices in the reverse direction
    //===------------------------------------------------------------------===//
    std::unordered_map<const SLM*, std::unordered_map<size_t, size_t>>
        rowsWithFreeSites;
    std::unordered_map<const SLM*, std::unordered_map<size_t, size_t>>
        columnsWithFreeSites;
    for (size_t r = 0; r < discreteTargetRows.size(); ++r) {
      if (rowsWithFreeSites.find(discreteTargetRows[r].first.first) ==
          rowsWithFreeSites.end()) {
        rowsWithFreeSites.emplace(discreteTargetRows[r].first.first,
                                  std::unordered_map<size_t, size_t>{});
      }
      rowsWithFreeSites[discreteTargetRows[r].first.first].emplace(
          discreteTargetRows[r].first.second, r);
    }
    for (size_t c = 0; c < discreteTargetColumns.size(); ++c) {
      if (columnsWithFreeSites.find(discreteTargetColumns[c].first.first) ==
          columnsWithFreeSites.end()) {
        columnsWithFreeSites.emplace(discreteTargetColumns[c].first.first,
                                     std::unordered_map<size_t, size_t>{});
      }
      columnsWithFreeSites[discreteTargetColumns[c].first.first].emplace(
          discreteTargetColumns[c].first.second, c);
    }
    //===------------------------------------------------------------------===//
    // Initialize the nearest free sites for each atom
    //===------------------------------------------------------------------===//
    const bool useWindowedPlacement =
        static_cast<T*>(this)->useWindowedPlacement;
    const size_t windowHeight = static_cast<T*>(this)->windowHeight;
    const size_t windowWidth = static_cast<T*>(this)->windowWidth;
    nearestFreeSitesForEachAtom.clear();
    nearestFreeSitesForEachAtom.reserve(atomsToPlace.size());
    for (const auto& [atom, _] : atomsToPlace) {
      const auto& [nearestSLM, nearestRow, nearestCol] =
          architecture.nearestStorageSite(previousPlacement[atom]);
      auto& nearestFreeSitesForAtom =
          nearestFreeSitesForEachAtom.emplace_back();
      size_t rLow = 0;
      size_t rHigh = nearestSLM->nRows;
      size_t cLow = 0;
      size_t cHigh = nearestSLM->nCols;
      if (useWindowedPlacement) {
        rLow =
            nearestRow > windowHeight / 2 ? nearestRow - (windowHeight / 2) : 0;
        rHigh =
            std::max(nearestRow + (windowHeight / 2) + 1, nearestSLM->nRows);
        cLow =
            nearestCol > windowWidth / 2 ? nearestCol - (windowWidth / 2) : 0;
        cHigh = std::max(nearestCol + (windowWidth / 2) + 1, nearestSLM->nCols);
      }
      for (size_t r = rLow; r < rHigh; ++r) {
        for (size_t c = cLow; c < cHigh; ++c) {
          if (occupiedStorageSites.find({nearestSLM, r, c}) ==
              occupiedStorageSites.end()) {
            nearestFreeSitesForAtom.emplace_back(
                std::pair{rowsWithFreeSites[nearestSLM][r],
                          columnsWithFreeSites[nearestSLM][c]},
                architecture.distance(previousPlacement[atom],
                                      {nearestSLM, r, c}));
          }
        }
      }
    }
    //===------------------------------------------------------------------===//
    // Run the A* algorithm
    //===------------------------------------------------------------------===//
    nodes.clear();
    nodes.emplace_back(0.0, std::unordered_set<std::pair<size_t, size_t>>{},
                       std::vector<std::map<size_t, size_t>>{},
                       std::vector<double>{},
                       std::vector<std::map<size_t, size_t>>{},
                       std::vector<double>{}); // root node
    const Node& finalNode =
        *aStarTreeSearch(*nodes.front(), &getQubitPlacementNeighbors, &isGoal,
                         &getCost, &getHeuristic)
             .back();
    //===------------------------------------------------------------------===//
    // Extract the final mapping
    //===------------------------------------------------------------------===//
    // the following two maps combine all hGroups and vGroups of the final node,
    // respectively, and hence serve as a mapping from rows and columns from the
    // previous placement to rows and columns of the next placement,
    // respectively.
    std::unordered_map<size_t, size_t> rowMapping;
    std::unordered_map<size_t, size_t> columnMapping;
    for (const auto& hGroup : finalNode.hGroups) {
      rowMapping.insert(hGroup.begin(), hGroup.end());
    }
    for (const auto& vGroup : finalNode.vGroups) {
      columnMapping.insert(vGroup.begin(), vGroup.end());
    }
    std::vector<const SLM*> finalSLMs;
    std::vector<size_t> finalRows;
    std::vector<size_t> finalColumns;
    for (size_t i = 0; i < nAtoms; ++i) {
      const auto& [startRow, startColumn] = currentSitesForEachAtom[i];
      const auto targetRow = rowMapping.at(startRow);
      const auto targetColumn = columnMapping.at(startColumn);
      finalSLMs.emplace_back(discreteTargetRows.at(targetRow).first.first);
      assert(discreteTargetRows.at(targetRow).first ==
             discreteTargetColumns.at(targetColumn).first);
      finalRows.emplace_back(discreteTargetRows.at(targetRow).second);
      finalColumns.emplace_back(discreteTargetColumns.at(targetColumn).second);
    }
    //===------------------------------------------------------------------===//
    // Store the final mapping
    //===------------------------------------------------------------------===//
    std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>&
        currentPlacement =
            static_cast<T*>(this)->getQubitMapping().emplace_back(
                previousPlacement);
    for (size_t i = 0; i < nAtoms; ++i) {
      const auto atom = atomsToPlace[i].first;
      currentPlacement[atom] =
          std::tuple{finalSLMs[i], finalRows[i], finalColumns[i]};
    }
  }

  /// A node representing one stage in the process of placing all atoms
  /// that must be moved for the next stage starting from the last mapping
  /// until a new mapping is found satisfying all constraints of the next
  /// stage
  struct Node {
    /// the maximum distance an already placed atom must travel to its
    /// target location
    float maxDistanceOfPlacedAtom = 0.0;
    /// a set of all sites that are already occupied by an atom due to the
    /// current placement
    std::unordered_set<std::pair<uint16_t, uint16_t>> consumedFreeSites{};
    /// a binary search tree representing the horizontal groups
    /// @see getNeighbors for more details
    std::vector<std::map<uint16_t, uint16_t>> hGroups;
    /// the maximum distance of placed atoms in every group to their
    /// target location
    std::vector<float> maxDistancesOfPlacedAtomsPerHGroup;
    /// @see hGroups
    std::vector<std::map<uint16_t, uint16_t>> vGroups;
    /// @see maxDistancesOfPlacedAtomsPerHGroup
    std::vector<float> maxDistancesOfPlacedAtomsPerVGroup;
  };

  /// A list of all nodes that have been created so far.
  /// This list is dynamically extended when new nodes are created.
  /// This happens when a node is expanded by calling getNeighbors.
  std::vector<std::unique_ptr<Node>> nodes;

  /// The number of atoms that must be placed in this stage.
  size_t nAtoms = 0;

  /// The discrete location of the atoms to be placed. The current location
  /// is taken from the previous stage.
  std::vector<std::pair<size_t, size_t>> currentSitesForEachAtom;

  /// The atoms that must be placed in this stage are numbered from 0 to
  /// nAtoms - 1.
  /// This vector lists all free sites for each atom ordered ascending by the
  /// distance to the atom.
  /// The distance itself is the second element of the pair.
  /// The set of free sites per atom may be limited by a window size that
  /// restricts the sites to be considered to be within the window around the
  /// very nearest site.
  std::vector<std::vector<std::pair<std::pair<size_t, size_t>, float>>>
      nearestFreeSitesForEachAtom;

  /// todo
  std::vector<std::vector<std::pair<
      std::pair<std::pair<size_t, size_t>, std::pair<size_t, size_t>>, float>>>
      nearestFreeSitesForEachGate;

  /// todo
  std::vector<std::pair<size_t, size_t>> atomsOfGatesToBePlaced;

  auto isGoal(const Node& node) -> bool {
    return node.consumedFreeSites.size() == nAtoms;
  }

  /// @brief Returns the cost of a node, i.e., the total cost to reach that node
  /// from the start node.
  /// @details Different groups cannot be rearranged concurrently in one step.
  /// Hence, we add up the time it takes to perform the rearrangements of one
  /// group in one step and sum it up over all groups.
  /// This will not resemble the exact time to rearrange all costs because at
  /// this point it is not clear how the horizontal and vertical groups can be
  /// combined.
  auto getCost(const Node& node) -> double {
    double cost = 0.0;
    for (const auto d : node.maxDistancesOfPlacedAtomsPerHGroup) {
      cost += std::sqrt(d);
    }
    for (const auto d : node.maxDistancesOfPlacedAtomsPerVGroup) {
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
  double getHeuristic(const Node& node) {
    float maxDistanceOfUnplacedAtom = 0.0;
    for (size_t i = node.consumedFreeSites.size() + 1; i < nAtoms; ++i) {
      for (const auto& site : nearestFreeSitesForEachAtom[i]) {
        if (node.consumedFreeSites.find(site.first) ==
            node.consumedFreeSites.end()) {
          maxDistanceOfUnplacedAtom =
              std::max(maxDistanceOfUnplacedAtom, site.second);
          break;
        }
      }
    }
    // We can multiply the difference by 2 because the cost function considers
    // this difference for the horizontal and vertical groups.
    return 2 *
           std::sqrt(maxDistanceOfUnplacedAtom - node.maxDistanceOfPlacedAtom);
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
  auto getQubitPlacementNeighbors(const Node& node)
      -> std::vector<const Node*> {
    size_t atomToBePlacedNext = node.consumedFreeSites.size();
    std::vector<const Node*> neighbors;
    for (const auto [site, distance] :
         nearestFreeSitesForEachAtom[atomToBePlacedNext]) {
      // assume nodes is of type std::vector<std::unique_ptr<Node>>
      // make a copy of node, the parent of neighbor
      Node& neighbor = *nodes.emplace_back(std::make_unique<Node>(node));
      neighbor.maxDistanceOfPlacedAtom =
          std::max(node.maxDistanceOfPlacedAtom, distance);
      neighbor.consumedFreeSites.emplace(site);
      //===------------------------------------------------------------------===//
      // check whether the current placement is compatible with any existing
      // horizontal group
      //===------------------------------------------------------------------===//
      checkCompatibilityAndAddPlacement(
          currentSitesForEachAtom[atomToBePlacedNext].first, site.first,
          distance, neighbor.hGroups,
          neighbor.maxDistancesOfPlacedAtomsPerHGroup);
      //===------------------------------------------------------------------===//
      // do the same for the vertical group
      //===------------------------------------------------------------------===//
      checkCompatibilityAndAddPlacement(
          currentSitesForEachAtom[atomToBePlacedNext].second, site.second,
          distance, neighbor.vGroups,
          neighbor.maxDistancesOfPlacedAtomsPerVGroup);
      //===------------------------------------------------------------------===//
      neighbors.emplace_back(neighbor);
    }
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
  auto getGatePlacementNeighbors(const Node& node) -> std::vector<const Node*> {
    const auto gateToBePlacedNext = node.consumedFreeSites.size();
    const auto& atomsToBePlacedNext =
        atomsOfGatesToBePlaced[gateToBePlacedNext];
    std::vector<const Node*> neighbors;
    for (const auto& [pair, distances] :
         nearestFreeSitesForEachGate[gateToBePlacedNext]) {
      const auto& [leftSite, rightSite] = pair;
      // assume nodes is of type std::vector<std::unique_ptr<Node>>
      // make a copy of node, the parent of neighbor
      Node& neighbor = *nodes.emplace_back(std::make_unique<Node>(node));
      neighbor.maxDistanceOfPlacedAtom =
          std::max(node.maxDistanceOfPlacedAtom,
                   std::max(distances.first, distances.second));
      neighbor.consumedFreeSites.emplace(leftSite);
      neighbor.consumedFreeSites.emplace(rightSite);
      //===------------------------------------------------------------------===//
      // Get the current placement of the atoms that must be placed next
      //===------------------------------------------------------------------===//
      const auto& currentSiteOfLeftAtom =
          currentSitesForEachAtom[atomsToBePlacedNext.first];
      const auto& currentSiteOfRightAtom =
          currentSitesForEachAtom[atomsToBePlacedNext.second];
      //===------------------------------------------------------------------===//
      // check whether the current placement is compatible with any existing
      // horizontal group
      //===------------------------------------------------------------------===//
      checkCompatibilityAndAddPlacement(
          currentSiteOfLeftAtom.first, leftSite.first, distances.first,
          neighbor.hGroups, neighbor.maxDistancesOfPlacedAtomsPerHGroup);
      checkCompatibilityAndAddPlacement(
          currentSiteOfRightAtom.first, rightSite.first, distances.second,
          neighbor.hGroups, neighbor.maxDistancesOfPlacedAtomsPerHGroup);
      //===------------------------------------------------------------------===//
      // do the same for the vertical group
      //===------------------------------------------------------------------===//
      checkCompatibilityAndAddPlacement(
          currentSiteOfLeftAtom.second, leftSite.second, distances.first,
          neighbor.vGroups, neighbor.maxDistancesOfPlacedAtomsPerVGroup);
      checkCompatibilityAndAddPlacement(
          currentSiteOfRightAtom.second, rightSite.second, distances.second,
          neighbor.vGroups, neighbor.maxDistancesOfPlacedAtomsPerVGroup);
      //===------------------------------------------------------------------===//
      neighbors.emplace_back(neighbor);
    }
  }

  /// Checks for the new placement of the atom whether it is compatible with
  /// one of the existing groups. If yes, the new placement is added to the
  /// respective group. Otherwise, a new group is formed with the new placement.
  /// @param key the start index of the new placement
  /// @param value the target index of the new placement
  /// @param groups the groups to which the new placement can be added
  /// @param maxDistances the maximum distances of placed atoms in each group
  /// @return true if the new placement could be added to an existing group
  auto checkCompatibilityAndAddPlacement(
      const uint16_t key, const uint16_t value, const float distance,
      std::vector<std::map<uint16_t, uint16_t>>& groups,
      std::vector<float>& maxDistances) -> bool {
    size_t i = 0;
    for (auto& hGroup : groups) {
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
    if (i == groups.size()) {
      // no compatible group could be found and a new group is created
      groups.emplace_back().emplace(key, value);
      maxDistances.emplace_back(distance);
      return false;
    }
    groups[i].emplace(key, value);
    maxDistances[i] = std::max(maxDistances[i], distance);
    return true;
  }

protected:
  AStarPlacer(const Architecture& architecture, const nlohmann::json& config)
      : architecture_(architecture) {
    // get first storage SLM and first entanglement SLM
    const auto& firstStorageSlm = *architecture_.get().storageZones.front();
    const auto& firstEntanglementSlm =
        *architecture_.get().entanglementZones.front().front();
    // check which side of the first storage SLM is closer to the entanglement
    // SLM
    if (firstStorageSlm.location.second <
        firstEntanglementSlm.location.second) {
      // if the entanglement SLM is closer to the last row of the storage SLM
      // start initial placement of the atoms in the last row instead of the
      // first and hence revert initial placement
      reverseInitialPlacement = true;
    }
  }
  [[nodiscard]] auto
  place(const size_t nQubits,
        const std::vector<std::vector<std::pair<qc::Qubit, qc::Qubit>>>&
            twoQubitGateLayers,
        const std::vector<std::unordered_set<qc::Qubit>>& reuseQubits) const
      -> std::vector<std::vector<std::tuple<const SLM&, size_t, size_t>>> {
    std::vector<std::vector<std::tuple<const SLM&, size_t, size_t>>> placement;
    placement.reserve((2 * twoQubitGateLayers.size()) + 1);
    placement.emplace_back(makeInitialPlacement(nQubits));
    for (size_t layer = 0; layer < twoQubitGateLayers.size(); ++layer) {
      placement.emplace_back(makeIntermediatePlacement(
          placement.back(), reuseQubits[layer], twoQubitGateLayers[layer]));
    }
  }
};
} // namespace na
