#include "na/nalac/NAMapper.hpp"

#include "Definitions.hpp"
#include "datastructures/Layer.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/CompoundOperation.hpp"
#include "ir/operations/OpType.hpp"
#include "ir/operations/Operation.hpp"
#include "na/Architecture.hpp"
#include "na/Configuration.hpp"
#include "na/NAComputation.hpp"
#include "na/entities/Atom.hpp"
#include "na/entities/Zone.hpp"
#include "na/nalac/NAGraphAlgorithms.hpp"
#include "na/operations/GlobalCZOp.hpp"
#include "na/operations/GlobalOp.hpp"
#include "na/operations/GlobalRYOp.hpp"
#include "na/operations/LocalRZOp.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <ratio>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na {

auto NAMapper::validateCircuit() -> void {
  for (const auto& op : initialQc) {
    if (op->isCompoundOperation() && isGlobal(*op, initialQc.getNqubits())) {
      const auto& co = dynamic_cast<qc::CompoundOperation&>(*op);
      if (!arch.isAllowedGlobally(co.at(0)->getType(), 0)) {
        std::stringstream ss;
        ss << "The chosen architecture does not support the operation "
           << op->getType() << " globally.";
        throw std::invalid_argument(ss.str());
      }
    } else if (op->isStandardOperation() && op->isSingleQubitGate()) {
      assert(op->getNcontrols() == 0);
      if (!arch.isAllowedLocally(op->getType(), 0)) {
        std::stringstream ss;
        ss << "The chosen architecture does not support the operation "
           << op->getType() << " locally.";
        throw std::invalid_argument(ss.str());
      }
    } else if (op->isStandardOperation() &&
               op->getNcontrols() + op->getNtargets() == 2) {
      assert(!op->isSingleQubitGate());
      if (!arch.isAllowedLocally(op->getType(), op->getNcontrols())) {
        if (!arch.isAllowedGlobally(op->getType(), op->getNcontrols())) {
          std::stringstream ss;
          ss << "The chosen architecture does not support the operation of "
                "type "
             << op->getType() << " with " << op->getNcontrols() << " controls "
             << " either locally or globally.";
          throw std::invalid_argument(ss.str());
        }
      }
    } else {
      if (op->isCompoundOperation()) {
        throw std::logic_error(
            "Compound operations are only supported when they are global.");
      }
      if (op->isStandardOperation()) {
        throw std::logic_error("Standard operations are only supported when "
                               "they act on one or two qubits.");
      }
      throw std::logic_error(
          "Operation class is not supported. Supported are StandardOperations "
          "and global CompoundOperations.");
    }
  }
}

auto NAMapper::makeLogicalArrays() -> void {
  const auto rows = static_cast<std::int64_t>(config.getPatchRows());
  const auto cols = static_cast<std::int64_t>(config.getPatchCols());
  const NAComputation logicQc = std::move(mappedQc);
  mappedQc = NAComputation();
  std::unordered_map<const na::Atom*, std::vector<std::vector<const na::Atom*>>>
      arrayAtoms{};
  for (const auto& atom : logicQc.getAtoms()) {
    arrayAtoms.emplace(atom.get(), std::vector<std::vector<const na::Atom*>>{});
    arrayAtoms[atom.get()].reserve(rows);
    for (std::int64_t r = 0; r < rows; ++r) {
      arrayAtoms[atom.get()].emplace_back();
      for (std::int64_t c = 0; c < cols; ++c) {
        mappedQc.emplaceInitialLocation(
            arrayAtoms[atom.get()][r].emplace_back(mappedQc.emplaceBackAtom(
                atom->getName() + "_" + std::to_string(r) + "_" +
                std::to_string(c))),
            initialArch.getLocationOffsetBy(
                logicQc.getInitialLocations().at(atom.get()), r, c));
      }
    }
  }
  std::unordered_map<const Zone*, const Zone*> arrayZones;
  for (const auto& zone : logicQc.getZones()) {
    arrayZones.emplace(zone.get(), mappedQc.emplaceBackZone(zone->getName()));
  }
  for (const auto& op : logicQc) {
    if (op->is<GlobalRYOp>()) {
      const auto& ry = op->as<GlobalRYOp>();
      mappedQc.emplaceBack<GlobalRYOp>(ry.getParams().front(),
                                       arrayZones[ry.getZone()]);
    } else if (op->is<GlobalCZOp>()) {
      const auto& cz = op->as<GlobalCZOp>();
      mappedQc.emplaceBack<GlobalCZOp>(arrayZones[cz.getZone()]);
    } else if (op->is<LocalRZOp>()) {
      const auto& rz = op->as<LocalRZOp>();
      const auto& allAtoms = rz.getAtoms();
      const auto& atomsPerRow = std::accumulate(
          allAtoms.cbegin(), allAtoms.cend(),
          std::vector < std::pair<double, std::vector<const na::Atom*>>(),
          [&](auto& acc, const auto* const p) {
            acc[p->y].insert(p->x);
            return acc;
          });
      for (const auto& [y, xs] : atomsPerRow) {
        for (std::int64_t r = 0; r < rows; ++r) {
          std::vector<std::shared_ptr<Location>> positions;
          positions.reserve(xs.size() * static_cast<std::size_t>(cols));
          for (const auto x : xs) {
            for (std::int64_t c = 0; c < cols; ++c) {
              positions.emplace_back(std::make_shared<Location>(
                  initialArch.getLocationOffsetBy({x, y}, r, c)));
            }
          }
          mappedQc.emplaceBack<NALocalOperation>(lop.getType(), lop.getParams(),
                                                 positions);
        }
      }
    } else if (op->isShuttlingOperation()) {
      const auto& sop = dynamic_cast<NAShuttlingOperation&>(*op);
      const auto shuttlingSize = sop.getStart().size();
      const auto pointCount =
          shuttlingSize * static_cast<std::size_t>(rows * cols);
      std::vector<std::shared_ptr<Location>> start;
      start.reserve(pointCount);
      std::vector<std::shared_ptr<Location>> end;
      end.reserve(pointCount);
      for (std::size_t i = 0; i < shuttlingSize; ++i) {
        const auto s = sop.getStart()[i];
        const auto e = sop.getEnd()[i];
        for (std::int64_t r = 0; r < rows; ++r) {
          for (std::int64_t c = 0; c < cols; ++c) {
            start.emplace_back(std::make_shared<Location>(
                initialArch.getLocationOffsetBy(*s, r, c)));
            end.emplace_back(std::make_shared<Location>(
                initialArch.getLocationOffsetBy(*e, r, c)));
          }
        }
      }
      mappedQc.emplaceBack<NAShuttlingOperation>(sop.getType(), start, end);
    } else {
      throw std::logic_error("Operation type is not supported.");
    }
  }
}

/**
 * @brief This function computes the actual movements of atoms.
 * @details This function is called after the circuit is mapped to the
 * architecture. At this point, only the start and end point of each shuttling
 * operations is specified. During shuttling the collisiion with other atoms
 * must be prevented. This is achieved by moving atoms first in y-direction and
 * then in x-direction with possible offsets to avoid collisions. In the example
 * below you can see the shuttling operation before and after this function is
 * called.
 *            Before                                  After
 * ===========================================================================
 *  Atom --> o     o     o <-- End          o     o     o     o <-- End
 *                     ∕                             ┌────────┘
 *     o     o     o ∕   o     o            o     o  ⏐  o     o     o
 *                 ∕                                 ⏐
 *     o     o   ∕ o     o     o            o     o  ⏐  o     o     o
 *             ∕                                     ⏐
 * Start --> o     o     o     o        Start --> o ─┘  o     o     o
 */
auto NAMapper::calculateMovements() -> void {
  const auto prelQC = mappedQc;
  mappedQc.clear(false);
  const auto d = static_cast<std::int64_t>(arch.getMinAtomDistance());
  for (const auto& op : prelQC) {
    if (op->isShuttlingOperation()) {
      const auto& shuttlingOp = dynamic_cast<NAShuttlingOperation&>(*op);
      if (shuttlingOp.getType() == MOVE) {
        std::vector<std::shared_ptr<Location>> hOffsetStart;
        std::vector<std::shared_ptr<Location>> hOffsetEnd;
        std::vector<std::shared_ptr<Location>> vMoveStart;
        std::vector<std::shared_ptr<Location>> vMoveEnd;
        std::vector<std::shared_ptr<Location>> hMoveStart;
        std::vector<std::shared_ptr<Location>> hMoveEnd;
        std::vector<std::shared_ptr<Location>> vOffsetStart;
        std::vector<std::shared_ptr<Location>> vOffsetEnd;
        bool vOffset = false;
        for (std::size_t i = 0; i < shuttlingOp.getStart().size(); ++i) {
          const auto start = *shuttlingOp.getStart()[i];
          const auto end = *shuttlingOp.getEnd()[i];
          const auto dx = end.x - start.x;
          Location const mid = {start.x, end.y};
          if (dx > 0) {
            if (const auto s = arch.getNearestSiteRight(mid, true)) {
              if (arch.getLocationOfSite(*s).x < end.x) {
                // in this case not the entire y-movement is possible in one go
                // we stop earlier to perform the x-movement without collision
                // and then continue with the y-movement (voffset)
                vOffset = true;
              }
            }
          } else if (dx < 0) {
            if (const auto s = arch.getNearestSiteLeft(mid, true)) {
              if (arch.getLocationOfSite(*s).x > end.x) {
                // see comment in the cae above
                vOffset = true;
              }
            }
          }
        }
        for (std::size_t i = 0; i < shuttlingOp.getStart().size(); ++i) {
          auto start = *shuttlingOp.getStart()[i];
          auto end = *shuttlingOp.getEnd()[i];
          const auto dx = end.x - start.x;
          const auto dy = end.y - start.y;
          if (dy > 0) {
            if (const auto s = arch.getNearestSiteDown(start, true)) {
              if (arch.getLocationOfSite(*s).y < end.y) {
                // in this case an atom is on the way
                hOffsetStart.emplace_back(std::make_shared<Location>(start));
                start.x += (dx >= 0 ? d : -d);
                hOffsetEnd.emplace_back(std::make_shared<Location>(start));
              }
            }
          } else if (dy < 0) {
            if (const auto s = arch.getNearestSiteUp(start, true)) {
              if (arch.getLocationOfSite(*s).y > end.y) {
                // in this case an atom is on the way
                hOffsetStart.emplace_back(std::make_shared<Location>(start));
                start.x += (dx >= 0 ? d : -d);
                hOffsetEnd.emplace_back(std::make_shared<Location>(start));
              }
            }
          }
          // at this stage, if an hoffset is necessary, the start position was
          // modified already to reflect the hoffset
          Location mid = {start.x, end.y};
          if (vOffset) {
            mid.y += (dy >= 0 ? -d : d);
          }
          if (start.y != mid.y) {
            vMoveStart.emplace_back(std::make_shared<Location>(start));
            start = mid;
            vMoveEnd.emplace_back(std::make_shared<Location>(start));
          }
          if (start.x != end.x) {
            hMoveStart.emplace_back(std::make_shared<Location>(start));
            start.x = end.x;
            hMoveEnd.emplace_back(std::make_shared<Location>(start));
          }
          if (start.y != end.y) {
            vOffsetStart.emplace_back(std::make_shared<Location>(start));
            start.y = end.y;
            vOffsetEnd.emplace_back(std::make_shared<Location>(start));
          }
        }
        if (!hOffsetStart.empty()) {
          mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, hOffsetStart,
                                                     hOffsetEnd);
        }
        if (!vMoveStart.empty()) {
          mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, vMoveStart,
                                                     vMoveEnd);
        }
        if (!hMoveStart.empty()) {
          mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, hMoveStart,
                                                     hMoveEnd);
        }
        if (!vOffsetStart.empty()) {
          mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, vOffsetStart,
                                                     vOffsetEnd);
        }
      } else {
        mappedQc.emplaceBack(op->clone());
      }
    } else {
      mappedQc.emplaceBack(op->clone());
    }
  }
}

auto NAMapper::checkApplicability(const qc::Operation* op,
                                  const std::vector<Atom>& placement) const
    -> bool {
  if (op->isCompoundOperation()) {
    // is global gate
    return true;
  }
  assert(op->isStandardOperation()); // ensured by preprocess
  if (op->isSingleQubitGate()) {
    assert(arch.isAllowedLocally(op->getType(), op->getNcontrols()));
    assert(op->getNcontrols() == 0);
    // individual gate that can act on one or more atoms
    return std::all_of(
        op->getTargets().cbegin(), op->getTargets().cend(),
        [&](const auto& qubit) {
          const auto& qubitPlacement = placement.at(qubit);
          switch (qubitPlacement.positionStatus) {
          case Atom::LocationStatus::UNDEFINED:
            // check whether the gate is applicable in one of the currently
            // selected zones
            return std::any_of(qubitPlacement.zones.cbegin(),
                               qubitPlacement.zones.cend(), [&](const auto& z) {
                                 return arch.isAllowedLocally(
                                     {op->getType(), 0}, z);
                               });
          case Atom::LocationStatus::DEFINED:
            // check whether the gate is applicable at the current position
            return arch.isAllowedLocallyAt({op->getType(), 0},
                                           *qubitPlacement.currentLocation);
          default:
            qc::unreachable();
          }
        });
  }
  assert(op->getNcontrols() + op->getNtargets() == 2);
  assert(arch.isAllowedGlobally({op->getType(), op->getNcontrols()}));
  // TODO global gate that acts exactly on two atoms
  return false;
}

auto NAMapper::updatePlacement(const qc::Operation* op,
                               std::vector<Atom>& placement) const -> void {
  if (op->isCompoundOperation()) {
    // global gates are represented as compound operations
    return;
  }
  assert(arch.isAllowedLocally(op->getType(), 0));
  assert(op->getNcontrols() == 0);
  // individual gate that can act on one or more atoms
  std::for_each(op->getTargets().cbegin(), op->getTargets().cend(),
                [&](const auto& qubit) {
                  switch (placement.at(qubit).positionStatus) {
                  case Atom::LocationStatus::UNDEFINED:
                    // remove all zones where the gate is not applicable
                    placement.at(qubit).zones.erase(
                        std::remove_if(placement.at(qubit).zones.begin(),
                                       placement.at(qubit).zones.end(),
                                       [&](auto& z) {
                                         return !arch.isAllowedLocally(
                                             {op->getType(), 0}, z);
                                       }),
                        placement.at(qubit).zones.end());
                    break;
                  case Atom::LocationStatus::DEFINED:
                    break;
                  default:
                    qc::unreachable();
                  }
                });
}

auto NAMapper::getMisplacement(const std::vector<Atom>& initial,
                               const std::vector<qc::Qubit>& target,
                               const qc::Qubit& q) -> std::int64_t {
  if (initial.at(q).positionStatus == Atom::LocationStatus::UNDEFINED) {
    return 0;
  }

  std::int64_t misplacement = 0;
  const auto indexOfQ = static_cast<std::size_t>(std::distance(
      target.cbegin(), std::find(target.cbegin(), target.cend(), q)));

  for (std::size_t i = 0; i < target.size(); ++i) {
    if (initial.at(target[i]).positionStatus ==
        Atom::LocationStatus::UNDEFINED) {
      continue;
    }

    if (i < indexOfQ && initial.at(target[i]).currentLocation.first >
                            initial.at(q).currentLocation.first) {
      misplacement += 1;
    }
    if (i > indexOfQ && initial.at(target[i]).currentLocation.first <
                            initial.at(q).currentLocation.first) {
      misplacement -= 1;
    }
  }

  for (const auto& p : target) {
    if (initial.at(p).currentLocation.first <
        initial.at(q).currentLocation.first) {
      misplacement += 1;
    }
  }

  return misplacement + static_cast<std::int64_t>(indexOfQ);
}

auto NAMapper::store(std::vector<bool>& initialFreeSites,
                     std::vector<bool>& currentFreeSites,
                     std::vector<Atom>& placement,
                     std::unordered_set<qc::Qubit>& currentlyShuttling,
                     const std::vector<qc::Qubit>& qubits,
                     const ZoneId destination) -> void {
  // this distance is used for spacing atoms that should interact or pass
  // another atom
  const auto d = static_cast<std::int64_t>(arch.getMinAtomDistance());
  const auto dx = static_cast<std::int64_t>(config.getPatchCols()) *
                  static_cast<std::int64_t>(arch.getNoInteractionRadius());
  { // load atoms that are not already shuttling
    std::vector<std::pair<std::int64_t, std::int64_t>> start;
    std::vector<std::pair<std::int64_t, std::int64_t>> end;
    for (const auto q : qubits) {
      if (currentlyShuttling.find(q) == currentlyShuttling.cend()) {
        start.emplace_back(placement.at(q).currentLocation);
        currentlyShuttling.insert(q);
        currentFreeSites.at(arch.getSiteAt(placement.at(q).currentLocation)) =
            true;
        placement.at(q).currentLocation = {
            placement.at(q).currentLocation.first + d,
            placement.at(q).currentLocation.second};
        end.emplace_back(placement.at(q).currentLocation);
      }
    }
    if (!start.empty()) {
      mappedQc.emplaceBack<NAShuttlingOperation>(LOAD, start, end);
    }
  }
  std::vector<std::tuple<Index, std::size_t>> freeSitesPerRow;
  const auto rowCount = arch.getNrowsInZone(destination);
  freeSitesPerRow.reserve(rowCount);
  for (std::size_t r = 0; r < rowCount; ++r) {
    const auto& sitesInRow = arch.getSitesInRow(destination, r);
    freeSitesPerRow.emplace_back(
        r, std::accumulate(sitesInRow.cbegin(), sitesInRow.cend(), 0UL,
                           [&](const auto& acc, const auto& s) {
                             return acc + (currentFreeSites.at(s) ? 1 : 0);
                           }));
  }
  std::sort(freeSitesPerRow.begin(), freeSitesPerRow.end(),
            [&](const auto& a, const auto& b) {
              return (std::get<1>(a) == std::get<1>(b) &&
                      std::get<0>(a) < std::get<0>(b)) ||
                     std::get<1>(a) > std::get<1>(b);
            });
  std::size_t moveableSpotsNeeded = qubits.size();
  std::vector<std::tuple<Index, std::size_t>> moveableSelectedRows;
  std::size_t firstIWithSameN = 0;
  for (std::size_t i = 0; i < freeSitesPerRow.size(); ++i) {
    const auto [r, n] = freeSitesPerRow[i];
    if (n >= moveableSpotsNeeded) {
      if (n != std::get<1>(freeSitesPerRow.at(firstIWithSameN))) {
        firstIWithSameN = i;
      }
      if (i + 1 == freeSitesPerRow.size() ||
          std::get<1>(freeSitesPerRow[i + 1]) < moveableSpotsNeeded) {
        const auto [rF, nF] = freeSitesPerRow.at(firstIWithSameN);
        moveableSelectedRows.emplace_back(rF, moveableSpotsNeeded);
        freeSitesPerRow.at(firstIWithSameN) = {rF, nF - moveableSpotsNeeded};
        break;
      }
    } else {
      moveableSelectedRows.emplace_back(r, n);
      freeSitesPerRow[i] = {r, 0};
      moveableSpotsNeeded -= n;
      firstIWithSameN = i + 1;
    }
  }
  std::sort(moveableSelectedRows.begin(), moveableSelectedRows.end(),
            [&](const auto& a, const auto& b) {
              return (std::get<0>(a) < std::get<0>(b));
            });
  for (auto& [r, n] : moveableSelectedRows) {
    std::vector<std::shared_ptr<Location>> start;
    std::vector<std::shared_ptr<Location>> end;
    std::vector<std::shared_ptr<Location>> storeStart;
    std::vector<std::shared_ptr<Location>> storeEnd;
    std::size_t notStoredLeft = 0;
    std::size_t j = 0;
    const auto& sitesInRow = arch.getSitesInRow(destination, r);
    const auto y =
        arch.getLocationOfSite(arch.getSitesInRow(destination, r).at(0)).y;
    for (const auto q : qubits) {
      if (currentlyShuttling.find(q) != currentlyShuttling.cend()) {
        start.emplace_back(placement.at(q).currentLocation);
        if (n == currentlyShuttling.size() - notStoredLeft) {
          const auto& site = *std::find_if(
              sitesInRow.cbegin(), sitesInRow.cend(),
              [&](const auto& s) { return currentFreeSites.at(s); });
          const auto& sPos = arch.getLocationOfSite(site);
          placement.at(q).currentLocation =
              std::make_shared<Location>(sPos.x + d, sPos.y);
          end.emplace_back(placement.at(q).currentLocation);
          storeStart.emplace_back(placement.at(q).currentLocation);
          placement.at(q).currentLocation =
              std::make_shared<Location>(sPos.x, sPos.y);
          storeEnd.emplace_back(placement.at(q).currentLocation);
          currentlyShuttling.erase(q);
          currentFreeSites.at(site) = false;
          initialFreeSites.at(site) = false;
          n -= 1;
        } else if (j < sitesInRow.size() &&
                   currentFreeSites.at(sitesInRow.at(j))) {
          const auto& s = sitesInRow.at(j);
          const auto& sPos = arch.getLocationOfSite(s);
          placement.at(q).currentLocation =
              std::make_shared<Location>(sPos.x + d, sPos.y);
          end.emplace_back(placement.at(q).currentLocation);
          storeStart.emplace_back(placement.at(q).currentLocation);
          placement.at(q).currentLocation =
              std::make_shared<Location>(sPos.x, sPos.y);
          storeEnd.emplace_back(placement.at(q).currentLocation);
          currentlyShuttling.erase(q);
          currentFreeSites.at(s) = false;
          initialFreeSites.at(s) = false;
          n -= 1;
        } else if (j < sitesInRow.size()) {
          const auto& s = sitesInRow.at(j);
          const auto& sPos = arch.getLocationOfSite(s);
          placement.at(q).currentLocation =
              std::make_shared<Location>(sPos.x + d, sPos.y);
          end.emplace_back(placement.at(q).currentLocation);
          notStoredLeft += 1;
        } else {
          placement.at(q).currentLocation = std::make_shared<Location>(
              arch.getLocationOfSite(sitesInRow.back()).x +
                  static_cast<std::int64_t>(j - sitesInRow.size() + 1) * dx + d,
              y);
          end.emplace_back(placement.at(q).currentLocation);
        }
        ++j;
      }
    }
    mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, start, end);
    mappedQc.emplaceBack<NAShuttlingOperation>(STORE, storeStart, storeEnd);
  }
}

auto NAMapper::pickUp(std::vector<bool>& initialFreeSites,
                      std::vector<bool>& currentFreeSites,
                      std::vector<Atom>& placement,
                      std::unordered_set<qc::Qubit>& currentlyShuttling,
                      const std::vector<qc::Qubit>& qubitsOrdered) -> void {
  // this distance is used for spacing atoms that should interact or pass
  // another atom
  const auto d = static_cast<std::int64_t>(arch.getMinAtomDistance());
  const auto dx = static_cast<std::int64_t>(config.getPatchCols()) *
                  static_cast<std::int64_t>(arch.getNoInteractionRadius());
  // get a vector of the atoms in the order to pick them up based on
  // their misplacement value
  std::vector pickUpOrder(qubitsOrdered);
  std::sort(pickUpOrder.begin(), pickUpOrder.end(),
            [&](const auto& a, const auto& b) {
              return std::abs(getMisplacement(placement, qubitsOrdered, a)) >
                     std::abs(getMisplacement(placement, qubitsOrdered, b));
            });
  // repeat until all atoms are picked up
  while (!pickUpOrder.empty()) {
    // extract the first atom from the pick up order
    const auto q = pickUpOrder.front();
    pickUpOrder.erase(std::remove(pickUpOrder.begin(), pickUpOrder.end(), q),
                      pickUpOrder.end());
    // get the index of the picked atom among the atoms ordered by their final
    // position
    const auto& i = static_cast<std::uint64_t>(std::distance(
        qubitsOrdered.cbegin(),
        std::find(qubitsOrdered.cbegin(), qubitsOrdered.cend(), q)));
    // if the placement of the atom is undefined, find a good placement for
    // it
    if (placement.at(q).positionStatus == Atom::LocationStatus::UNDEFINED) {
      // calculate not picked up atoms to the left in the resulting order
      // note: all remaining atoms that are not picked up yet have undefined
      // positions as well
      std::size_t notPickedUpLeft = 0U;
      for (std::size_t j = 0; j < i; ++j) {
        const auto& p = qubitsOrdered.at(j);
        // not picked up yet and left of q in the end
        if (currentlyShuttling.find(p) == currentlyShuttling.cend()) {
          ++notPickedUpLeft;
        }
      }
      const auto spotsNeeded = pickUpOrder.size();
      ZoneId zone = 0;
      Index row = 0;
      std::size_t freeSpotsInRow = 0;
      for (const auto& z : placement.at(q).zones) {
        for (std::size_t r = 0; r < arch.getNrowsInZone(z); ++r) {
          const auto& sitesInRow = arch.getSitesInRow(z, r);
          const auto& n =
              std::accumulate(sitesInRow.cbegin(), sitesInRow.cend(), 0UL,
                              [&](const auto& acc, const auto& s) {
                                return acc + (initialFreeSites.at(s) ? 1 : 0);
                              });
          if (freeSpotsInRow <= spotsNeeded && n > freeSpotsInRow) {
            freeSpotsInRow = n;
            zone = z;
            row = r;
          }
        }
      }
      // place q on the i-th free spot where i is the minimum of the number
      // of not picked up atoms to the left and free spots available
      auto possibleSites = arch.getSitesInRow(zone, row);
      possibleSites.erase(
          std::remove_if(possibleSites.begin(), possibleSites.end(),
                         [&](const auto s) { return !initialFreeSites.at(s); }),
          possibleSites.end());
      const auto s =
          possibleSites[std::min(notPickedUpLeft, freeSpotsInRow - 1)];
      placement.at(q).positionStatus = Atom::LocationStatus::DEFINED;
      *placement.at(q).initialLocation = arch.getLocationOfSite(s);
      initialFreeSites.at(s) = false;
      currentFreeSites.at(s) = false;
    }
    // here the position of q is defined
    std::vector<std::shared_ptr<Location>> start;
    std::vector<std::shared_ptr<Location>> end;
    std::vector<std::shared_ptr<Location>> loadStart;
    std::vector<std::shared_ptr<Location>> loadEnd;
    const auto currentX = placement.at(q).currentLocation->x;
    const auto y = placement.at(q).currentLocation->y;
    // pick up q itself
    loadStart.emplace_back(placement.at(q).currentLocation);
    currentFreeSites.at(*arch.getSiteAt(*placement.at(q).currentLocation)) =
        true;
    placement.at(q).currentLocation =
        std::make_shared<Location>(currentX + d, y);
    loadEnd.emplace_back(placement.at(q).currentLocation);
    currentlyShuttling.insert(q);
    // iterate through all not yet picked up atoms to the left and check
    // whether they can be picked up
    const auto nextX =
        arch.getNearestXLeft(currentX, arch.getZoneAt({currentX, y}), true);
    auto x = nextX == currentX ? currentX - dx : nextX + d;
    if (i > 0) {
      for (std::size_t j = i - 1;; --j) {
        const auto p = qubitsOrdered.at(j);
        // if the atom is picked up,
        // move it to the correct row
        if (currentlyShuttling.find(p) != currentlyShuttling.cend()) {
          start.emplace_back(placement.at(p).currentLocation);
          placement.at(p).currentLocation = std::make_shared<Location>(x, y);
          end.emplace_back(placement.at(p).currentLocation);
          const auto xl = arch.getNearestXLeft(x, arch.getZoneAt({x, y}), true);
          const auto nx =
              arch.getNearestXLeft(xl, arch.getZoneAt({xl, y}), true);
          x = nx == xl ? xl - dx : nx + d;
        } else {
          // check whether j can be
          // picked up
          if (placement.at(p).positionStatus == Atom::LocationStatus::DEFINED) {
            if (placement.at(p).currentLocation->y == y &&
                placement.at(p).currentLocation->x <= x - d) {
              // pick up p
              pickUpOrder.erase(
                  std::remove(pickUpOrder.begin(), pickUpOrder.end(), p),
                  pickUpOrder.end());
              x = placement.at(p).currentLocation->x;
              loadStart.emplace_back(placement.at(p).currentLocation);
              currentFreeSites.at(
                  *arch.getSiteAt(*placement.at(p).currentLocation)) = true;
              placement.at(p).currentLocation =
                  std::make_shared<Location>(x + d, y);
              loadEnd.emplace_back(placement.at(p).currentLocation);
              currentlyShuttling.insert(p);
              const auto nx =
                  arch.getNearestXLeft(x, arch.getZoneAt({x, y}), true);
              x = nx == x ? x - dx : nx + d;
            }
          } else {
            // if atom is not placed yet find a site to the left that it can
            // be picked up together with q, i.e. find next free site to the
            // left
            std::int64_t freeX = x;
            bool free = false;
            while (!free) {
              const auto siteOpt = arch.getNearestSiteLeft({freeX, y});
              if (!siteOpt) {
                break;
              }
              const auto site = *siteOpt;
              freeX = arch.getLocationOfSite(site).x;
              if (initialFreeSites.at(site) &&
                  std::find(placement.at(p).zones.cbegin(),
                            placement.at(p).zones.cend(),
                            arch.getZoneOfSite(site)) !=
                      placement.at(p).zones.cend()) {
                // the site is free and satisfies the zone restrictions of
                // the atom
                free = true;
              } else {
                const auto nx = arch.getNearestXLeft(
                    freeX, arch.getZoneAt({freeX - d, y}), true);
                freeX = nx == freeX ? freeX - dx : nx + d;
              }
            }
            if (free) {
              // place p on the free site
              placement.at(p).positionStatus = Atom::LocationStatus::DEFINED;
              *placement.at(p).initialLocation = {freeX, y};
              initialFreeSites.at(
                  *arch.getSiteAt(*placement.at(p).initialLocation)) = false;
              // pick up p
              pickUpOrder.erase(
                  std::remove(pickUpOrder.begin(), pickUpOrder.end(), p),
                  pickUpOrder.end());
              loadStart.emplace_back(placement.at(p).currentLocation);
              placement.at(p).currentLocation =
                  std::make_shared<Location>(freeX + d, y);
              loadEnd.emplace_back(placement.at(p).currentLocation);
              currentlyShuttling.insert(p);
              const auto nx = arch.getNearestXLeft(
                  freeX, arch.getZoneAt({freeX - d, y}), true);
              x = nx == x ? x - dx : nx + d;
            }
          }
        }
        if (j == 0) {
          break;
        }
      }
    }
    // iterate through all not yet picked up atoms to the right and check
    // whether they can be picked up
    x = currentX + d;
    const auto nx1 = arch.getNearestXRight(x, arch.getZoneAt({x, y}));
    x = nx1 == x ? x + dx : nx1 + d;
    for (std::size_t j = i + 1; j < qubitsOrdered.size(); ++j) {
      const auto p = qubitsOrdered[j];
      // if the atom is picked up, move it to the correct row
      if (currentlyShuttling.find(p) != currentlyShuttling.cend()) {
        start.emplace_back(placement.at(p).currentLocation);
        placement.at(p).currentLocation = std::make_shared<Location>(x, y);
        end.emplace_back(placement.at(p).currentLocation);
        const auto nx = arch.getNearestXRight(x, arch.getZoneAt({x, y}));
        x = nx == x ? x + dx : nx + d;
      } else {
        // check whether j can be picked up
        if (placement.at(p).positionStatus == Atom::LocationStatus::DEFINED) {
          if (placement.at(p).currentLocation->y == y &&
              placement.at(p).currentLocation->x >= x - d) {
            // pick up p
            pickUpOrder.erase(
                std::remove(pickUpOrder.begin(), pickUpOrder.end(), p),
                pickUpOrder.end());
            x = placement.at(p).currentLocation->x;
            loadStart.emplace_back(placement.at(p).currentLocation);
            currentFreeSites.at(
                *arch.getSiteAt(*placement.at(p).currentLocation)) = true;
            placement.at(p).currentLocation =
                std::make_shared<Location>(x + d, y);
            loadEnd.emplace_back(placement.at(p).currentLocation);
            currentlyShuttling.insert(p);
            const auto nx =
                arch.getNearestXRight(x + d, arch.getZoneAt({x, y}));
            x = nx == x + d ? nx - d + dx : nx + d;
          }
        } else {
          // if atom is not placed yet find a site to the right that it can
          // be picked up together with q, i.e. find next free site to the
          // right
          std::int64_t freeX = x;
          bool free = false;
          while (!free) {
            const auto siteOpt = arch.getNearestSiteRight({freeX - d, y});
            if (!siteOpt) {
              break;
            }
            const auto site = *siteOpt;
            freeX = arch.getLocationOfSite(site).x;
            if (initialFreeSites.at(site) &&
                std::find(placement.at(p).zones.cbegin(),
                          placement.at(p).zones.cend(),
                          arch.getZoneOfSite(site)) !=
                    placement.at(p).zones.cend()) {
              // the site is free and satisfies the zone restrictions of
              // the atom
              free = true;
            } else {
              freeX =
                  arch.getNearestXRight(freeX + d, arch.getZoneAt({freeX, y})) +
                  d;
            }
          }
          if (free) {
            // place p on the free site
            placement.at(p).positionStatus = Atom::LocationStatus::DEFINED;
            *placement.at(p).initialLocation = {freeX, y};
            initialFreeSites.at(
                *arch.getSiteAt(*placement.at(p).initialLocation)) = false;
            // pick up p
            pickUpOrder.erase(
                std::remove(pickUpOrder.begin(), pickUpOrder.end(), p),
                pickUpOrder.end());
            loadStart.emplace_back(placement.at(p).currentLocation);
            placement.at(p).currentLocation =
                std::make_shared<Location>(freeX + d, y);
            loadEnd.emplace_back(placement.at(p).currentLocation);
            currentlyShuttling.insert(p);
            const auto nx =
                arch.getNearestXRight(freeX + d, arch.getZoneAt({freeX, y}));
            x = nx == freeX + d ? nx - d + dx : nx + d;
          }
        }
      }
    }
    if (!start.empty()) {
      mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, start, end);
    }
    mappedQc.emplaceBack<NAShuttlingOperation>(LOAD, loadStart, loadEnd);
  }
}
/// @returns a pair representing the NA computation in a discrete fashion, which
/// is later transformed into a real NACOmputation using double values for
/// locations. Those cannot be used in maps which we need during the mapping
/// process. The first vector in the pair are the initial positions of atoms,
/// the second vector are the operations that are executed in the NA
/// computation.
auto NAMapper::map(const qc::QuantumComputation& qc) -> std::pair <
    std::vector<std::tuple<Zone, Index, Index>, std::vector<std::tuple<>>> {
  auto startPreprocess = std::chrono::high_resolution_clock::now();
  initialQc = qc;
  const auto nqubits = initialQc.getNqubits();
  std::size_t maxSeqWidth = 0;
  mappedQc = NAComputation();
  preprocess();
  // store the placement of atoms, both the initial one (needed later) and the
  // current one leave atoms unmapped as long as possible. This mighty induce
  // some restrictions on the mapping in which zone the atom may be placed.
  const auto initialZones = arch.getInitialZones();
  std::vector<Atom> placement(nqubits);
  std::for_each(placement.begin(), placement.end(),
                [&](auto& p) { p = Atom(initialZones); });
  std::vector initialFreeSites(arch.getNSites(), true);
  std::vector currentFreeSites(arch.getNSites(), true);
  std::unordered_set<qc::Qubit> currentlyShuttling{};
  // this distance is used for spacing atoms that should interact or pass
  // another atom
  const auto d = static_cast<std::int64_t>(arch.getMinAtomDistance());
  const auto dx = static_cast<std::int64_t>(config.getPatchCols()) *
                  static_cast<std::int64_t>(arch.getNoInteractionRadius());
  // get start time
  auto startMapping = std::chrono::high_resolution_clock::now();
  //============================ START MAPPING ============================
  const qc::Layer layer(initialQc);
  const auto& executableSet = layer.getExecutableSet();
  auto it = executableSet.begin();
  if (config.getMethod() == NAMappingMethod::Naive) {
    for (qc::Qubit q = 0; q < nqubits; ++q) {
      placement.at(q).positionStatus = Atom::LocationStatus::DEFINED;
      const auto s = arch.getSitesInZone(initialZones.front()).at(q);
      *placement.at(q).initialLocation = arch.getLocationOfSite(s);
      initialFreeSites.at(s) = false;
      currentFreeSites.at(s) = false;
    }
  }
  while (it != executableSet.end()) {
    // 1. execute all gates that are directly applicable and do not need
    //    shuttling
    while (it != executableSet.end()) {
      const auto* const op = (*it)->getOperation();
      if (checkApplicability(op, placement)) {
        updatePlacement(op, placement);
        (*it)->execute();
        if (op->isCompoundOperation()) {
          const auto* const co = dynamic_cast<const qc::CompoundOperation*>(op);
          mappedQc.emplaceBack<NAGlobalOperation>(
              FullOpType{co->at(0)->getType(), 0}, co->at(0)->getParameter());
        } else if (isGlobal(*op, nqubits) &&
                   arch.isAllowedGlobally(
                       {op->getType(), op->getNcontrols()})) {
          mappedQc.emplaceBack<NAGlobalOperation>(
              FullOpType{op->getType(), op->getNcontrols()},
              op->getParameter());
        } else {
          // collect executable gates of the same type
          std::vector<std::shared_ptr<Location>> positions = {
              placement.at(op->getTargets().front()).currentLocation};
          for (const auto& v :
               layer.getExecutablesOfType(op->getType(), op->getNcontrols())) {
            const auto* const op2 = v->getOperation();
            if (checkApplicability(op2, placement) &&
                op->getParameter() == op2->getParameter() &&
                std::find(
                    positions.cbegin(), positions.cend(),
                    placement.at(op2->getTargets().front()).currentLocation) ==
                    positions.cend()) {
              updatePlacement(op2, placement);
              v->execute();
              positions.emplace_back(
                  placement.at(op2->getTargets().front()).currentLocation);
            }
          }
          mappedQc.emplaceBack<NALocalOperation>(
              FullOpType{op->getType(), op->getNcontrols()}, op->getParameter(),
              positions);
        }
        it = executableSet.begin();
      } else {
        ++it;
      }
    }
    it = executableSet.begin();
    if (it == executableSet.end()) {
      break;
    }
    // 2. when no such gates are left, extract an interaction graph of gates
    //    of the same type and two targets, i.e. cz gates
    if (config.getMethod() == NAMappingMethod::Naive) {
      const qc::Operation* op = (*it)->getOperation();
      if (op->getType() != qc::OpType::Z || op->getNtargets() != 1 ||
          op->getNcontrols() != 1) {
        throw std::logic_error(
            "Other gates than cz are not supported for mapping yet.");
      }
      (*it)->execute();
      const auto& q1 = op->getTargets().front();
      const auto& q2 = op->getControls().begin()->qubit;
      Location start = *placement.at(q1).currentLocation;
      Location end = start;
      const Location& target = arch.getLocationOfSite(
          arch.getSitesInZone(*arch.getPropertiesOfOperation({op->getType(), 1})
                                   .zones.begin())
              .at(0));
      end.x += d;
      mappedQc.emplaceBack<NAShuttlingOperation>(
          LOAD, std::vector{std::make_shared<Location>(start)},
          std::vector{std::make_shared<Location>(end)});
      start = end;
      end = target;
      end.x += d;
      mappedQc.emplaceBack<NAShuttlingOperation>(
          MOVE, std::vector{std::make_shared<Location>(start)},
          std::vector{std::make_shared<Location>(end)});
      start = end;
      end.x -= d;
      mappedQc.emplaceBack<NAShuttlingOperation>(
          STORE, std::vector{std::make_shared<Location>(start)},
          std::vector{std::make_shared<Location>(end)});
      start = *placement.at(q2).currentLocation;
      end = start;
      end.x += d;
      mappedQc.emplaceBack<NAShuttlingOperation>(
          LOAD, std::vector{std::make_shared<Location>(start)},
          std::vector{std::make_shared<Location>(end)});
      start = end;
      end = target;
      end.y += d;
      mappedQc.emplaceBack<NAShuttlingOperation>(
          MOVE, std::vector{std::make_shared<Location>(start)},
          std::vector{std::make_shared<Location>(end)});
      mappedQc.emplaceBack<NAGlobalOperation>(FullOpType{qc::OpType::Z, 1});
      maxSeqWidth = 1UL;
      mappedQc.emplaceBack<NAShuttlingOperation>(
          MOVE, std::vector{std::make_shared<Location>(end)},
          std::vector{std::make_shared<Location>(start)});
      end = *placement.at(q2).currentLocation;
      mappedQc.emplaceBack<NAShuttlingOperation>(
          STORE, std::vector{std::make_shared<Location>(start)},
          std::vector{std::make_shared<Location>(end)});
      start = target;
      end = start;
      end.x += d;
      mappedQc.emplaceBack<NAShuttlingOperation>(
          LOAD, std::vector{std::make_shared<Location>(start)},
          std::vector{std::make_shared<Location>(end)});
      start = end;
      end = *placement.at(q1).currentLocation;
      end.x += d;
      mappedQc.emplaceBack<NAShuttlingOperation>(
          MOVE, std::vector{std::make_shared<Location>(start)},
          std::vector{std::make_shared<Location>(end)});
      start = end;
      end.x -= d;
      mappedQc.emplaceBack<NAShuttlingOperation>(
          STORE, std::vector{std::make_shared<Location>(start)},
          std::vector{std::make_shared<Location>(end)});
    } else if (config.getMethod() ==
               NAMappingMethod::MaximizeParallelismHeuristic) {
      const auto& graph = layer.constructInteractionGraph(qc::OpType::Z, 1);
      if (graph.getNVertices() == 0) {
        throw std::logic_error(
            "Other gates than cz are not supported for mapping yet.");
        // TODO: support other gates than cz
      }
      const ZoneId interactionZone =
          *arch.getPropertiesOfOperation({qc::OpType::Z, 1}).zones.begin();
      const auto sites = arch.getSitesInRow(interactionZone, 0);
      const auto& sequence =
          NAGraphAlgorithms::computeSequence(graph, sites.size());
      const auto& moveable = sequence.first;
      const auto& fixed = sequence.second;
      // 3. move the atoms accordingly and execute the gates
      // pick up the fixed atoms and move them to the interaction zone
      // get a vector of the fixed atoms ordered by their initial position
      // from left to right
      maxSeqWidth = std::accumulate(
          fixed.cbegin(), fixed.cend(), maxSeqWidth,
          [&](const auto max, const auto& q) {
            assert(q.second >= 0);
            return std::max(max, static_cast<std::size_t>(q.second));
          });
      if (sites.size() < maxSeqWidth) {
        std::stringstream ss;
        ss << "Target site in " << arch.getZoneLabel(interactionZone);
        ss << " zone is out of bounds. Possible reason for this error: The "
              "zone is not wide enough.";
        throw std::logic_error(ss.str());
      }
      std::vector<qc::Qubit> fixedOrdered;
      std::transform(fixed.cbegin(), fixed.cend(), back_inserter(fixedOrdered),
                     [](const auto& v) { return v.first; });
      std::sort(fixedOrdered.begin(), fixedOrdered.end(),
                [&](const auto& a, const auto& b) {
                  return fixed.at(a) < fixed.at(b);
                });
      pickUp(initialFreeSites, currentFreeSites, placement, currentlyShuttling,
             fixedOrdered);
      // all atoms are picked up in order, move them to the destination zone and
      // store them there
      std::vector<std::shared_ptr<Location>> start{};
      start.reserve(fixed.size());
      std::vector<std::shared_ptr<Location>> mid{};
      mid.reserve(fixed.size());
      std::vector<std::shared_ptr<Location>> end{};
      end.reserve(fixed.size());
      for (const auto& [q, x] : fixed) {
        if (currentlyShuttling.find(q) == currentlyShuttling.cend()) {
          std::stringstream ss;
          ss << "Atom " << q << " was unexpectedly not picked up.";
          throw std::logic_error(ss.str());
        }
        assert(x >= 0);
        const auto& p =
            arch.getLocationOfSite(sites.at(static_cast<std::size_t>(x)));
        if (!currentFreeSites.at(*arch.getSiteAt(p))) {
          throw std::logic_error(
              "Target site in interaction zone is unexpectedly occupied.");
        }
        start.emplace_back(placement.at(q).currentLocation);
        placement.at(q).currentLocation =
            std::make_shared<Location>(p.x + d, p.y);
        mid.emplace_back(placement.at(q).currentLocation);
        currentFreeSites.at(*arch.getSiteAt(p)) = false;
        placement.at(q).currentLocation = std::make_shared<Location>(p);
        end.emplace_back(placement.at(q).currentLocation);
      }
      currentlyShuttling.clear();
      mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, start, mid);
      mappedQc.emplaceBack<NAShuttlingOperation>(STORE, mid, end);
      // -----------------------------------------------------------------
      std::vector<qc::Qubit> moveableOrdered;
      std::transform(moveable.at(0).cbegin(), moveable.at(0).cend(),
                     back_inserter(moveableOrdered),
                     [](const auto& v) { return v.first; });
      std::sort(moveableOrdered.begin(), moveableOrdered.end(),
                [&](const auto& a, const auto& b) {
                  return moveable.at(0).at(a) < moveable.at(0).at(b);
                });
      pickUp(initialFreeSites, currentFreeSites, placement, currentlyShuttling,
             moveableOrdered);
      // ------------------------------------------------------------------
      // 4. Apply the cz gates
      std::vector<std::shared_ptr<Location>> startMoveable;
      std::vector<std::shared_ptr<Location>> endMoveable;
      std::transform(moveableOrdered.cbegin(), moveableOrdered.cend(),
                     std::back_inserter(startMoveable), [&](const auto& q) {
                       return placement.at(q).currentLocation;
                     });
      for (const auto& timeframe : moveable) {
        for (const auto q : moveableOrdered) {
          const auto x = timeframe.at(q);
          if (currentlyShuttling.find(q) == currentlyShuttling.cend()) {
            std::stringstream ss;
            ss << "Atom " << q << " was unexpectedly not picked up.";
            throw std::logic_error(ss.str());
          }
          if (x >= 0 && static_cast<std::size_t>(x) < sites.size()) {
            const auto pos =
                arch.getLocationOfSite(sites.at(static_cast<std::size_t>(x)));
            placement.at(q).currentLocation =
                std::make_shared<Location>(pos.x, pos.y + d);
          } else if (x < 0) {
            const auto pos = arch.getLocationOfSite(sites.at(0));
            placement.at(q).currentLocation =
                std::make_shared<Location>(pos.x + (x * dx - d), pos.y + d);
          } else { // x >= sites.size()
            const auto pos = arch.getLocationOfSite(sites.at(sites.size() - 1));
            placement.at(q).currentLocation =
                std::make_shared<Location>(pos.x + (x * dx + d), pos.y + d);
          }
          endMoveable.emplace_back(placement.at(q).currentLocation);
        }
        mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, startMoveable,
                                                   endMoveable);
        mappedQc.emplaceBack<NAGlobalOperation>(FullOpType{qc::OpType::Z, 1});
        for (const auto& q : moveableOrdered) {
          for (const auto& [p, _] : fixed) {
            const auto qPos = *placement.at(q).currentLocation;
            const auto pPos = *placement.at(p).currentLocation;
            if ((qPos - pPos).length() <= arch.getInteractionRadius()) {
              graph.getEdge(p, q)->execute();
            }
          }
        }
        startMoveable.clear();
        startMoveable = endMoveable;
        endMoveable.clear();
      }
      // -----------------------------------------------------------------
      // 5. move the atoms back to the storage zone
      if (initialZones.size() != 1) {
        throw std::logic_error("Currently only one storage zone is supported.");
      }
      const auto storageZones = initialZones.front();
      store(initialFreeSites, currentFreeSites, placement, currentlyShuttling,
            moveableOrdered, storageZones);
      // -------------------------------------------------------------
      std::vector<qc::Qubit> fixedVector;
      fixedVector.reserve(fixed.size());
      std::transform(fixed.cbegin(), fixed.cend(),
                     std::back_inserter(fixedVector),
                     [](const auto& v) { return v.first; });
      std::sort(fixedVector.begin(), fixedVector.end(),
                [&fixed](const auto& a, const auto& b) {
                  return fixed.at(a) < fixed.at(b);
                });
      store(initialFreeSites, currentFreeSites, placement, currentlyShuttling,
            fixedVector, storageZones);
    } else {
      throw std::logic_error("NA mapping method not implemented.");
    }
    // -------------------------------------------------------------
    it = executableSet.begin();
  }
  for (auto& p : placement) {
    if (p.positionStatus == Atom::LocationStatus::UNDEFINED) {
      // find next free site in first zone
      std::vector<Index> possibleSites;
      for (const auto& z : p.zones) {
        for (const auto& s : arch.getSitesInZone(z)) {
          possibleSites.emplace_back(s);
        }
      }
      const auto& freeSite =
          std::find_if(possibleSites.cbegin(), possibleSites.cend(),
                       [&](const auto& s) { return initialFreeSites.at(s); });
      *p.initialLocation = arch.getLocationOfSite(*freeSite);
      p.positionStatus = Atom::LocationStatus::DEFINED;
      initialFreeSites.at(*freeSite) = false;
      currentFreeSites.at(*freeSite) = false;
    }
    mappedQc.emplaceInitialLocation(p.initialLocation);
  }
  //========================= END MAPPING =========================
  // get end time
  auto startPostprocess = std::chrono::high_resolution_clock::now();
  postprocess();
  auto end = std::chrono::high_resolution_clock::now();
  // build remaining statistics
  stats.numInitialGates = qc.getNops();
  stats.numEntanglingGates = static_cast<std::size_t>(
      std::count_if(qc.cbegin(), qc.cend(), [](const auto& op) {
        return (op->getType() == qc::OpType::Z &&
                op->getType() == qc::OpType::X) ||
               op->getNcontrols() > 0;
      }));
  stats.initialDepth = qc.getDepth();
  stats.numMappedGates = mappedQc.size();
  stats.numQubits = nqubits;
  stats.maxSeqWidth = config.getPatchCols() * maxSeqWidth;
  // get the mapping time in milliseconds
  stats.preprocessTime =
      std::chrono::duration<qc::fp, std::milli>(startMapping - startPreprocess)
          .count();
  stats.mappingTime =
      std::chrono::duration<qc::fp, std::milli>(startPostprocess - startMapping)
          .count();
  stats.postprocessTime =
      std::chrono::duration<qc::fp, std::milli>(end - startPostprocess).count();
  done = true;
}
} // namespace na
