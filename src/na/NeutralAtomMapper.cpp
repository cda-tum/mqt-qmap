//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/mqt-qmap for more
// information.
//

#include "NeutralAtomMapper.hpp"

#include "Architecture.hpp"
#include "Layer.hpp"
#include "NAGraphAlgorithms.hpp"
#include "QuantumComputation.hpp"
#include "na/NADefinitions.hpp"
#include "na/operations/NAGlobalOperation.hpp"
#include "na/operations/NALocalOperation.hpp"
#include "na/operations/NAShuttlingOperation.hpp"
#include "operations/Operation.hpp"

#include <__algorithm/remove.h>
#include <__algorithm/remove_if.h>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <vector>

namespace na {

auto NeutralAtomMapper::validateCircuit() -> void {
  for (const auto& op : initialQc) {
    if (op->isCompoundOperation() and isGlobal(*op, initialQc.getNqubits())) {
      const auto& co = dynamic_cast<qc::CompoundOperation&>(*op);
      if (!arch.isAllowedGlobally({co.at(0)->getType(), 0})) {
        std::stringstream ss;
        ss << "The chosen architecture does not support the operation "
           << OpType{op->getType(), 0} << " globally.";
        throw std::invalid_argument(ss.str());
      }
    } else if (op->isStandardOperation()) {
      if (!op->isIndividual() and op->getNcontrols() + op->getNtargets() > 2) {
        std::stringstream ss;
        ss << "The current implementation does not support the operation "
           << OpType{op->getType(), op->getNcontrols()} << " acting on more "
           << "than two qubits.";
        throw std::logic_error(ss.str());
      }
      if (!arch.isAllowedLocally({op->getType(), op->getNcontrols()})) {
        if (!arch.isAllowedGlobally({op->getType(), op->getNcontrols()})) {
          std::stringstream ss;
          ss << "The chosen architecture does not support the operation "
             << OpType{op->getType(), 0} << " either locally or globally.";
          throw std::invalid_argument(ss.str());
        }
        // the gate is global: if it is a 1Q-gate it must be applied globally
        if (op->isIndividual() and !isGlobal(*op, initialQc.getNqubits())) {
          std::stringstream ss;
          ss << "The chosen architecture does not support the operation "
             << OpType{op->getType(), op->getNcontrols()} << " locally.";
          throw std::invalid_argument(ss.str());
        }
      }
    } else {
      throw std::logic_error("Operation type is not supported.");
    }
  }
}

auto NeutralAtomMapper::makeLogicalArrays() -> void {
  const auto logicQC = mappedQc;
  mappedQc.clear();
  mappedQc.clearInitialPositions();
  const auto rows = static_cast<std::int64_t>(config.getPatchRows());
  const auto cols = static_cast<std::int64_t>(config.getPatchCols());
  for (const auto& p : logicQC.getInitialPositions()) {
    for (std::int64_t r = 0; r < rows; ++r) {
      for (std::int64_t c = 0; c < cols; ++c) {
        mappedQc.emplaceInitialPosition(
            std::make_shared<Point>(initialArch.getPositionOffsetBy(*p, r, c)));
      }
    }
  }
  for (const auto& op : logicQC) {
    if (op->isGlobalOperation()) {
      mappedQc.emplaceBack(op->clone());
    } else if (op->isLocalOperation()) {
      const auto& lop          = dynamic_cast<NALocalOperation&>(*op);
      const auto& allPositions = lop.getPositions();
      const auto& posPerRow =
          std::accumulate(allPositions.cbegin(), allPositions.cend(),
                          std::map<std::int64_t, std::set<std::int64_t>>(),
                          [&](auto& acc, const auto& p) {
                            acc[p->y].insert(p->x);
                            return acc;
                          });
      for (const auto& [y, xs] : posPerRow) {
        for (std::int64_t r = 0; r < rows; ++r) {
          std::vector<std::shared_ptr<Point>> positions;
          for (const auto x : xs) {
            for (std::int64_t c = 0; c < cols; ++c) {
              positions.emplace_back(std::make_shared<Point>(
                  initialArch.getPositionOffsetBy({x, y}, r, c)));
            }
          }
          mappedQc.emplaceBack<NALocalOperation>(lop.getType(), lop.getParams(),
                                                 positions);
        }
      }
    } else if (op->isShuttlingOperation()) {
      const auto& sop = dynamic_cast<NAShuttlingOperation&>(*op);
      std::vector<std::shared_ptr<Point>> start;
      std::vector<std::shared_ptr<Point>> end;
      for (std::size_t i = 0; i < sop.getStart().size(); ++i) {
        const auto s = sop.getStart()[i];
        const auto e = sop.getEnd()[i];
        for (std::int64_t r = 0; r < rows; ++r) {
          for (std::int64_t c = 0; c < cols; ++c) {
            start.emplace_back(std::make_shared<Point>(
                initialArch.getPositionOffsetBy(*s, r, c)));
            end.emplace_back(std::make_shared<Point>(
                initialArch.getPositionOffsetBy(*e, r, c)));
          }
        }
      }
      mappedQc.emplaceBack<NAShuttlingOperation>(sop.getType(), start, end);
    } else {
      throw std::logic_error("Operation type is not supported.");
    }
  }
}

auto NeutralAtomMapper::calculateMovements() -> void {
  const auto prelQC = mappedQc;
  mappedQc.clear();
  const auto d = static_cast<std::int64_t>(arch.getMinAtomDistance());
  for (const auto& op : prelQC) {
    if (op->isShuttlingOperation()) {
      const auto& shuttlingOp = dynamic_cast<NAShuttlingOperation&>(*op);
      if (shuttlingOp.getType() == MOVE) {
        std::vector<std::shared_ptr<Point>> hOffsetStart;
        std::vector<std::shared_ptr<Point>> hOffsetEnd;
        std::vector<std::shared_ptr<Point>> vMoveStart;
        std::vector<std::shared_ptr<Point>> vMoveEnd;
        std::vector<std::shared_ptr<Point>> hMoveStart;
        std::vector<std::shared_ptr<Point>> hMoveEnd;
        std::vector<std::shared_ptr<Point>> vOffsetStart;
        std::vector<std::shared_ptr<Point>> vOffsetEnd;
        bool                                vOffset = false;
        for (std::size_t i = 0; i < shuttlingOp.getStart().size(); ++i) {
          const auto  start = *shuttlingOp.getStart()[i];
          const auto  end   = *shuttlingOp.getEnd()[i];
          const auto  dx    = end.x - start.x;
          Point const mid   = {start.x, end.y};
          if (dx > 0) {
            if (arch.hasSiteRight(mid, true)) {
              const auto& s = arch.getNearestSiteRight(mid, true);
              if (arch.getPositionOfSite(s).x < end.x) {
                vOffset = true;
              }
            }
          } else if (dx < 0) {
            if (arch.hasSiteLeft(mid, true)) {
              const auto& s = arch.getNearestSiteLeft(mid, true);
              if (arch.getPositionOfSite(s).x > end.x) {
                vOffset = true;
              }
            }
          }
        }
        for (std::size_t i = 0; i < shuttlingOp.getStart().size(); ++i) {
          auto       start = *shuttlingOp.getStart()[i];
          auto       end   = *shuttlingOp.getEnd()[i];
          const auto dx    = end.x - start.x;
          const auto dy    = end.y - start.y;
          if (dy > 0) {
            if (arch.hasSiteDown(start, true)) {
              const auto& s = arch.getNearestSiteDown(start, true);
              if (arch.getPositionOfSite(s).y < end.y) {
                // in this case an atom is on the way
                hOffsetStart.emplace_back(std::make_shared<Point>(start));
                start.x += (dx >= 0 ? d : -d);
                hOffsetEnd.emplace_back(std::make_shared<Point>(start));
              }
            }
          } else if (dy < 0) {
            if (arch.hasSiteUp(start, true)) {
              const auto& s = arch.getNearestSiteUp(start, true);
              if (arch.getPositionOfSite(s).y > end.y) {
                // in this case an atom is on the way
                hOffsetStart.emplace_back(std::make_shared<Point>(start));
                start.x += (dx >= 0 ? d : -d);
                hOffsetEnd.emplace_back(std::make_shared<Point>(start));
              }
            }
          }
          Point mid = {start.x, end.y};
          if (vOffset) {
            mid.y += (dy >= 0 ? -d : d);
          }
          if (start.y != mid.y) {
            vMoveStart.emplace_back(std::make_shared<Point>(start));
            start = mid;
            vMoveEnd.emplace_back(std::make_shared<Point>(start));
          }
          if (start.x != end.x) {
            hMoveStart.emplace_back(std::make_shared<Point>(start));
            start.x = end.x;
            hMoveEnd.emplace_back(std::make_shared<Point>(start));
          }
          if (start.y != end.y) {
            vOffsetStart.emplace_back(std::make_shared<Point>(start));
            start.y = end.y;
            vOffsetEnd.emplace_back(std::make_shared<Point>(start));
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

auto NeutralAtomMapper::checkApplicability(
    const std::unique_ptr<qc::Operation>& op,
    const std::vector<Atom>&              placement) const -> bool {
  if (op->isCompoundOperation()) {
    return true;
  }
  assert(op->isStandardOperation()); // ensured by preprocess
  if (arch.isAllowedLocally({op->getType(), op->getNcontrols()})) {
    if (op->getNcontrols() == 0) {
      // individual gate that can act on one or more atoms
      return std::all_of(
          op->getTargets().cbegin(), op->getTargets().cend(),
          [&](const auto& qubit) {
            switch (placement[qubit].positionStatus) {
            case Atom::PositionStatus::UNDEFINED:
              // check whether the gate is applicable in one of the currently
              // selected zones
              return std::any_of(
                  placement[qubit].zones.cbegin(),
                  placement[qubit].zones.cend(), [&](const auto& z) {
                    return arch.isAllowedLocally({op->getType(), 0}, z);
                  });
            case Atom::PositionStatus::DEFINED:
              // check whether the gate is applicable at the current position
              return arch.isAllowedLocallyAt({op->getType(), 0},
                                             *placement[qubit].currentPosition);
            }
          });
    }
    // TODO single controlled local gate that acts exactly on two atoms
    return false;
  }
  assert(arch.isAllowedGlobally({op->getType(), op->getNcontrols()}));
  if (op->isIndividual()) {
    assert(isGlobal(*op, initialQc.getNqubits()));
    // purely globally gate can always be applied
    return true;
  }
  assert(op->getNcontrols() + op->getNtargets() == 2);
  // TODO single controlled global gate that acts exactly on two atoms
  return false;
}

auto NeutralAtomMapper::updatePlacement(
    const std::unique_ptr<qc::Operation>& op,
    std::vector<Atom>&                    placement) const -> void {
  if (op->isCompoundOperation()) {
    return;
  }
  if (arch.isAllowedLocally({op->getType(), op->getNcontrols()})) {
    if (op->getNcontrols() == 0) {
      // individual gate that can act on one or more atoms
      std::for_each(op->getTargets().cbegin(), op->getTargets().cend(),
                    [&](const auto& qubit) {
                      switch (placement[qubit].positionStatus) {
                      case Atom::PositionStatus::UNDEFINED:
                        // remove all zones where the gate is not applicable
                        placement[qubit].zones.erase(
                            std::remove_if(placement[qubit].zones.begin(),
                                           placement[qubit].zones.end(),
                                           [&](auto& z) {
                                             return !arch.isAllowedLocally(
                                                 {op->getType(), 0}, z);
                                           }),
                            placement[qubit].zones.end());
                        break;
                      case Atom::PositionStatus::DEFINED:
                        break;
                      }
                    });
    } else {
      // TODO single controlled local gate that acts exactly on two atoms
    }
  } else {
    assert(arch.isAllowedGlobally({op->getType(), op->getNcontrols()}));
    if (op->isIndividual()) {
      assert(isGlobal(*op, initialQc.getNqubits()));
      // purely globally gate can always be applied
      std::for_each(op->getTargets().cbegin(), op->getTargets().cend(),
                    [&](const auto qubit) {
                      switch (placement[qubit].positionStatus) {
                      case Atom::PositionStatus::UNDEFINED:
                        placement[qubit].zones.erase(
                            std::remove_if(placement[qubit].zones.begin(),
                                           placement[qubit].zones.end(),
                                           [&](auto& z) {
                                             return !arch.isAllowedGlobally(
                                                 {op->getType(), 0}, z);
                                           }),
                            placement[qubit].zones.end());
                        break;
                      case Atom::PositionStatus::DEFINED:
                        break;
                      }
                    });
    } else {
      assert(op->getNcontrols() + op->getNtargets() == 2);
      // TODO single controlled global gate that acts exactly on two atoms
    }
  }
}

auto NeutralAtomMapper::getMisplacement(const std::vector<Atom>&      initial,
                                        const std::vector<qc::Qubit>& target,
                                        const qc::Qubit& q) -> std::int64_t {
  if (initial.at(q).positionStatus == Atom::PositionStatus::DEFINED) {
    std::vector<std::size_t> enumerated(target.size());
    std::iota(enumerated.begin(), enumerated.end(), 0);
    const auto indexOfQ = static_cast<std::size_t>(std::distance(
        target.cbegin(), std::find(target.cbegin(), target.cend(), q)));
    return std::accumulate(enumerated.cbegin(), enumerated.cend(), 0,
                           [&](const auto& acc, const auto& i) {
                             if (initial.at(target[i]).positionStatus ==
                                 Atom::PositionStatus::DEFINED) {
                               if (initial.at(target[i]).currentPosition->x >
                                       initial.at(q).currentPosition->x and
                                   i < indexOfQ) {
                                 return acc + 1;
                               }
                               if (initial.at(target[i]).currentPosition->x <
                                       initial.at(q).currentPosition->x and
                                   i > indexOfQ) {
                                 return acc - 1;
                               }
                             }
                             return acc;
                           }) +
           // count the spots the atom must be moved relatively
           static_cast<std::int64_t>(indexOfQ) -
           std::count_if(target.cbegin(), target.cend(), [&](const auto& p) {
             return initial[p].currentPosition->x <
                    initial[q].currentPosition->x;
           });
  }
  return 0;
}

auto NeutralAtomMapper::shuttle(
    std::vector<bool>& initialFreeSites, std::vector<bool>& currentFreeSites,
    std::vector<Atom>&             placement,
    std::unordered_set<qc::Qubit>& currentlyShuttling,
    const std::vector<qc::Qubit>& qubits, const Zone destination) -> void {
  // this distance is used for spacing atoms that should interact or pass
  // another atom
  const auto d  = static_cast<std::int64_t>(arch.getMinAtomDistance());
  const auto dx = static_cast<std::int64_t>(config.getPatchCols()) *
                  static_cast<std::int64_t>(arch.getNoInteractionRadius());
  { // load atoms that are not already shuttling
    std::vector<std::shared_ptr<Point>> start;
    std::vector<std::shared_ptr<Point>> end;
    for (const auto q : qubits) {
      if (currentlyShuttling.find(q) == currentlyShuttling.cend()) {
        start.emplace_back(placement[q].currentPosition);
        const auto pos = *placement[q].currentPosition;
        currentlyShuttling.insert(q);
        currentFreeSites[arch.getSiteAt(*placement[q].currentPosition)] = true;
        placement[q].currentPosition =
            std::make_shared<Point>(pos.x + d, pos.y);
        end.emplace_back(placement[q].currentPosition);
      }
    }
  }
  std::vector<std::tuple<Index, std::size_t>> freeSitesPerRow;
  for (std::size_t r = 0; r < arch.getNrowsInZone(destination); ++r) {
    const auto& sitesInRow = arch.getSitesInRow(destination, r);
    freeSitesPerRow.emplace_back(
        r, std::accumulate(sitesInRow.cbegin(), sitesInRow.cend(), 0UL,
                           [&](const auto& acc, const auto& s) {
                             return acc + (currentFreeSites[s] ? 1 : 0);
                           }));
  }
  std::sort(freeSitesPerRow.begin(), freeSitesPerRow.end(),
            [&](const auto& a, const auto& b) {
              return (std::get<1>(a) == std::get<1>(b) and
                      std::get<0>(a) < std::get<0>(b)) or
                     std::get<1>(a) > std::get<1>(b);
            });
  std::size_t moveableSpotsNeeded = qubits.size();
  std::vector<std::tuple<Index, std::size_t>> moveableSelectedRows;
  std::size_t                                 firstIWithSameN = 0;
  for (std::size_t i = 0; i < freeSitesPerRow.size(); ++i) {
    const auto [r, n] = freeSitesPerRow[i];
    if (n >= moveableSpotsNeeded) {
      if (n != std::get<1>(freeSitesPerRow[firstIWithSameN])) {
        firstIWithSameN = i;
      }
      if (i + 1 == freeSitesPerRow.size() or
          std::get<1>(freeSitesPerRow[i + 1]) < moveableSpotsNeeded) {
        const auto [rF, nF] = freeSitesPerRow[firstIWithSameN];
        moveableSelectedRows.emplace_back(rF, moveableSpotsNeeded);
        freeSitesPerRow[firstIWithSameN] = {rF, nF - moveableSpotsNeeded};
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
    std::vector<std::shared_ptr<Point>> start;
    std::vector<std::shared_ptr<Point>> end;
    std::vector<std::shared_ptr<Point>> storeStart;
    std::vector<std::shared_ptr<Point>> storeEnd;
    std::size_t                         notStoredLeft = 0;
    std::size_t                         j             = 0;
    const auto& sitesInRow = arch.getSitesInRow(destination, r);
    const auto  y =
        arch.getPositionOfSite(arch.getSitesInRow(destination, r)[0]).y;
    for (const auto q : qubits) {
      if (currentlyShuttling.find(q) != currentlyShuttling.cend()) {
        start.emplace_back(placement[q].currentPosition);
        if (n == currentlyShuttling.size() - notStoredLeft) {
          const auto& site =
              *std::find_if(sitesInRow.cbegin(), sitesInRow.cend(),
                            [&](const auto& s) { return currentFreeSites[s]; });
          const auto& sPos = arch.getPositionOfSite(site);
          placement[q].currentPosition =
              std::make_shared<Point>(sPos.x + d, sPos.y);
          end.emplace_back(placement[q].currentPosition);
          storeStart.emplace_back(placement[q].currentPosition);
          placement[q].currentPosition =
              std::make_shared<Point>(sPos.x, sPos.y);
          storeEnd.emplace_back(placement[q].currentPosition);
          currentlyShuttling.erase(q);
          currentFreeSites[site] = false;
          initialFreeSites[site] = false;
          n -= 1;
        } else if (j < sitesInRow.size() and currentFreeSites[sitesInRow[j]]) {
          const auto& s    = sitesInRow[j];
          const auto& sPos = arch.getPositionOfSite(s);
          placement[q].currentPosition =
              std::make_shared<Point>(sPos.x + d, sPos.y);
          end.emplace_back(placement[q].currentPosition);
          storeStart.emplace_back(placement[q].currentPosition);
          placement[q].currentPosition =
              std::make_shared<Point>(sPos.x, sPos.y);
          storeEnd.emplace_back(placement[q].currentPosition);
          currentlyShuttling.erase(q);
          currentFreeSites[s] = false;
          initialFreeSites[s] = false;
          n -= 1;
        } else if (j < sitesInRow.size()) {
          const auto& s    = sitesInRow[j];
          const auto& sPos = arch.getPositionOfSite(s);
          placement[q].currentPosition =
              std::make_shared<Point>(sPos.x + d, sPos.y);
          end.emplace_back(placement[q].currentPosition);
          notStoredLeft += 1;
        } else {
          placement[q].currentPosition = std::make_shared<Point>(
              arch.getPositionOfSite(sitesInRow.back()).x +
                  static_cast<std::int64_t>(j - sitesInRow.size() + 1) * dx + d,
              y);
          end.emplace_back(placement[q].currentPosition);
        }
        ++j;
      }
    }
    mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, start, end);
    mappedQc.emplaceBack<NAShuttlingOperation>(STORE, storeStart, storeEnd);
  }
}

auto NeutralAtomMapper::shuttle(
    std::vector<bool>& initialFreeSites, std::vector<bool>& currentFreeSites,
    std::vector<Atom>&                          placement,
    std::unordered_set<qc::Qubit>&              currentlyShuttling,
    const std::unordered_map<qc::Qubit, Point>& qubits, bool store) -> void {
  // this distance is used for spacing atoms that should interact or pass
  // another atom
  const auto d  = static_cast<std::int64_t>(arch.getMinAtomDistance());
  const auto dx = static_cast<std::int64_t>(config.getPatchCols()) *
                  static_cast<std::int64_t>(arch.getNoInteractionRadius());
  std::vector<qc::Qubit> qubitsOrdered;
  std::transform(qubits.cbegin(), qubits.cend(), back_inserter(qubitsOrdered),
                 [](const auto& v) { return v.first; });
  std::sort(qubitsOrdered.begin(), qubitsOrdered.end(),
            [&](const auto& a, const auto& b) {
              return qubits.at(a).x < qubits.at(b).x;
            });
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
    if (placement[q].positionStatus == Atom::PositionStatus::UNDEFINED) {
      // calculate not picked up atoms to the left in the resulting order
      // note: all remaining atoms that are not picked up yet have undefined
      // positions as well
      auto notPickedUpLeft = 0UL;
      for (std::size_t j = 0; j < i; ++j) {
        const auto& p = qubitsOrdered[j];
        // not picked up yet and left
        // of q in the end
        if (currentlyShuttling.find(p) == currentlyShuttling.cend() and
            qubits.at(p).x < qubits.at(q).x) {
          ++notPickedUpLeft;
        }
      }
      const auto  spotsNeeded    = pickUpOrder.size();
      Zone        zone           = 0;
      Index       row            = 0;
      std::size_t freeSpotsInRow = 0;
      for (const auto& z : placement[q].zones) {
        for (std::size_t r = 0; r < arch.getNrowsInZone(z); ++r) {
          const auto& sitesInRow = arch.getSitesInRow(z, r);
          const auto& n =
              std::accumulate(sitesInRow.cbegin(), sitesInRow.cend(), 0UL,
                              [&](const auto& acc, const auto& s) {
                                return acc + (initialFreeSites[s] ? 1 : 0);
                              });
          if (freeSpotsInRow <= spotsNeeded and n > freeSpotsInRow) {
            freeSpotsInRow = n;
            zone           = z;
            row            = r;
          }
        }
      }
      // place q on the i-th free spot where i is the minimum of the number
      // of not picked up atoms to the left and free spots available
      auto possibleSites = arch.getSitesInRow(zone, row);
      possibleSites.erase(
          std::remove_if(possibleSites.begin(), possibleSites.end(),
                         [&](const auto s) { return !initialFreeSites[s]; }),
          possibleSites.end());
      const auto s =
          possibleSites[std::min(notPickedUpLeft, freeSpotsInRow - 1)];
      placement[q].positionStatus   = Atom::PositionStatus::DEFINED;
      *placement[q].initialPosition = arch.getPositionOfSite(s);
      initialFreeSites[s]           = false;
      currentFreeSites[s]           = false;
    }
    // here the position of q is defined
    std::vector<std::shared_ptr<Point>> start;
    std::vector<std::shared_ptr<Point>> end;
    std::vector<std::shared_ptr<Point>> loadStart;
    std::vector<std::shared_ptr<Point>> loadEnd;
    const auto currentX = placement[q].currentPosition->x;
    const auto y        = placement[q].currentPosition->y;
    // pick up q itself
    loadStart.emplace_back(placement[q].currentPosition);
    currentFreeSites[arch.getSiteAt(*placement[q].currentPosition)] = true;
    placement[q].currentPosition = std::make_shared<Point>(currentX + d, y);
    loadEnd.emplace_back(placement[q].currentPosition);
    currentlyShuttling.insert(q);
    // iterate through all not yet picked up atoms to the left and check
    // whether they can be picked up
    const auto nextX =
        arch.getNearestXLeft(currentX, arch.getZoneAt({currentX, y}), true);
    auto x = nextX == currentX ? currentX - dx : nextX + d;
    if (i > 0) {
      for (std::size_t j = i - 1;; --j) {
        const auto p = qubitsOrdered[j];
        // if the atom is picked up,
        // move it to the correct row
        if (currentlyShuttling.find(p) != currentlyShuttling.cend()) {
          start.emplace_back(placement[p].currentPosition);
          placement[p].currentPosition = std::make_shared<Point>(x, y);
          end.emplace_back(placement[p].currentPosition);
          const auto xl = arch.getNearestXLeft(x, arch.getZoneAt({x, y}), true);
          const auto nx =
              arch.getNearestXLeft(xl, arch.getZoneAt({xl, y}), true);
          x = nx == xl ? xl - dx : nx + d;
        } else {
          // check whether j can be
          // picked up
          if (placement[p].positionStatus == Atom::PositionStatus::DEFINED) {
            if (placement[p].currentPosition->y == y and
                placement[p].currentPosition->x <= x - d) {
              // pick up p
              pickUpOrder.erase(
                  std::remove(pickUpOrder.begin(), pickUpOrder.end(), p),
                  pickUpOrder.end());
              x = placement[p].currentPosition->x;
              loadStart.emplace_back(placement[p].currentPosition);
              currentFreeSites[arch.getSiteAt(*placement[p].currentPosition)] =
                  true;
              placement[p].currentPosition = std::make_shared<Point>(x + d, y);
              loadEnd.emplace_back(placement[p].currentPosition);
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
            bool         free  = false;
            while (!free and arch.hasSiteLeft({freeX, y})) {
              const auto& site = arch.getNearestSiteLeft({freeX, y});
              freeX            = arch.getPositionOfSite(site).x;
              if (initialFreeSites[site] and
                  std::find(
                      placement[p].zones.cbegin(), placement[p].zones.cend(),
                      arch.getZoneOfSite(site)) != placement[p].zones.cend()) {
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
              placement[p].positionStatus   = Atom::PositionStatus::DEFINED;
              *placement[p].initialPosition = {freeX, y};
              initialFreeSites[arch.getSiteAt(*placement[p].initialPosition)] =
                  false;
              // pick up p
              pickUpOrder.erase(
                  std::remove(pickUpOrder.begin(), pickUpOrder.end(), p),
                  pickUpOrder.end());
              loadStart.emplace_back(placement[p].currentPosition);
              placement[p].currentPosition =
                  std::make_shared<Point>(freeX + d, y);
              loadEnd.emplace_back(placement[p].currentPosition);
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
    x              = currentX + d;
    const auto nx1 = arch.getNearestXRight(x, arch.getZoneAt({x, y}));
    x              = nx1 == x ? x + dx : nx1 + d;
    for (std::size_t j = i + 1; j < qubitsOrdered.size(); ++j) {
      const auto p = qubitsOrdered[j];
      // if the atom is picked up, move it to the correct row
      if (currentlyShuttling.find(p) != currentlyShuttling.cend()) {
        start.emplace_back(placement[p].currentPosition);
        placement[p].currentPosition = std::make_shared<Point>(x, y);
        end.emplace_back(placement[p].currentPosition);
        const auto nx = arch.getNearestXRight(x, arch.getZoneAt({x, y}));
        x             = nx == x ? x + dx : nx + d;
      } else {
        // check whether j can be picked up
        if (placement[p].positionStatus == Atom::PositionStatus::DEFINED) {
          if (placement[p].currentPosition->y == y and
              placement[p].currentPosition->x >= x - d) {
            // pick up p
            pickUpOrder.erase(
                std::remove(pickUpOrder.begin(), pickUpOrder.end(), p),
                pickUpOrder.end());
            x = placement[p].currentPosition->x;
            loadStart.emplace_back(placement[p].currentPosition);
            currentFreeSites[arch.getSiteAt(*placement[p].currentPosition)] =
                true;
            placement[p].currentPosition = std::make_shared<Point>(x + d, y);
            loadEnd.emplace_back(placement[p].currentPosition);
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
          bool         free  = false;
          while (!free and arch.hasSiteRight({freeX - d, y})) {
            const auto& site = arch.getNearestSiteRight({freeX - d, y});
            freeX            = arch.getPositionOfSite(site).x;
            if (initialFreeSites[site] and
                std::find(
                    placement[p].zones.cbegin(), placement[p].zones.cend(),
                    arch.getZoneOfSite(site)) != placement[p].zones.cend()) {
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
            placement[p].positionStatus   = Atom::PositionStatus::DEFINED;
            *placement[p].initialPosition = {freeX, y};
            initialFreeSites[arch.getSiteAt(*placement[p].initialPosition)] =
                false;
            // pick up p
            pickUpOrder.erase(
                std::remove(pickUpOrder.begin(), pickUpOrder.end(), p),
                pickUpOrder.end());
            loadStart.emplace_back(placement[p].currentPosition);
            placement[p].currentPosition =
                std::make_shared<Point>(freeX + d, y);
            loadEnd.emplace_back(placement[p].currentPosition);
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
  if (store) {
    // all atoms are picked up in order, move them to the destination zone and
    // store them there
    std::vector<std::shared_ptr<Point>> start;
    std::vector<std::shared_ptr<Point>> mid;
    std::vector<std::shared_ptr<Point>> end;
    for (const auto& [q, p] : qubits) {
      if (currentlyShuttling.find(q) == currentlyShuttling.cend()) {
        std::stringstream ss;
        ss << "Atom " << q << " was unexpectedly not picked up.";
        throw std::logic_error(ss.str());
      }
      if (!currentFreeSites[arch.getSiteAt(p)]) {
        throw std::logic_error(
            "Target site in interaction zone is unexpectedly occupied.");
      }
      start.emplace_back(placement[q].currentPosition);
      placement[q].currentPosition = std::make_shared<Point>(p.x + d, p.y);
      mid.emplace_back(placement[q].currentPosition);
      currentFreeSites[arch.getSiteAt(p)] = false;
      placement[q].currentPosition        = std::make_shared<Point>(p);
      end.emplace_back(placement[q].currentPosition);
    }
    currentlyShuttling.clear();
    mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, start, mid);
    mappedQc.emplaceBack<NAShuttlingOperation>(STORE, mid, end);
  }
}

auto NeutralAtomMapper::map(const qc::QuantumComputation& qc) -> void {
  auto startPreprocess    = std::chrono::high_resolution_clock::now();
  initialQc               = qc;
  const auto  nqubits     = initialQc.getNqubits();
  std::size_t maxSeqWidth = 0;
  mappedQc                = NAComputation();
  preprocess();
  // store the placement of atoms, both the initial one (needed later) and the
  // current one leave atoms unmapped as long as possible. This mighty induce
  // some restrictions on the mapping in which zone the atom may be placed.
  const auto        initialZones = arch.getInitialZones();
  std::vector<Atom> placement(nqubits);
  std::for_each(placement.begin(), placement.end(),
                [&](auto& p) { p = Atom(initialZones); });
  std::vector<bool>             initialFreeSites(arch.getNSites(), true);
  std::vector<bool>             currentFreeSites(arch.getNSites(), true);
  std::unordered_set<qc::Qubit> currentlyShuttling{};
  // this distance is used for spacing atoms that should interact or pass
  // another atom
  const auto d  = static_cast<std::int64_t>(arch.getMinAtomDistance());
  const auto dx = static_cast<std::int64_t>(config.getPatchCols()) *
                  static_cast<std::int64_t>(arch.getNoInteractionRadius());
  // get start time
  auto startMapping = std::chrono::high_resolution_clock::now();
  //============================ START MAPPING ============================
  const qc::Layer layer(initialQc);
  const auto&     executableSet = layer.getExecutableSet();
  auto            it            = (*executableSet)->begin();
  if (config.getMethod() == NaMappingMethod::NAIVE) {
    for (qc::Qubit q = 0; q < nqubits; ++q) {
      placement[q].positionStatus = Atom::PositionStatus::DEFINED;
      const auto s = arch.getSitesInZone(initialZones.front())[q];
      *placement[q].initialPosition = arch.getPositionOfSite(s);
      initialFreeSites[s]           = false;
      currentFreeSites[s]           = false;
    }
  }
  while (it != (*executableSet)->end()) {
    // 1. execute all gates that are directly applicable and do not need
    //    shuttling
    while (it != (*executableSet)->end()) {
      const std::unique_ptr<qc::Operation>& op = *(*it)->getOperation();
      if (checkApplicability(op, placement)) {
        updatePlacement(op, placement);
        (*it)->execute();
        if (op->isCompoundOperation()) {
          const auto& co = dynamic_cast<qc::CompoundOperation&>(*op);
          mappedQc.emplaceBack<NAGlobalOperation>(
              OpType{co.at(0)->getType(), 0}, co.at(0)->getParameter());
        } else if (isGlobal(*op, nqubits) and
                   arch.isAllowedGlobally(
                       {op->getType(), op->getNcontrols()})) {
          mappedQc.emplaceBack<NAGlobalOperation>(
              OpType{op->getType(), op->getNcontrols()}, op->getParameter());
        } else {
          // collect executable gates of the same type
          std::vector<std::shared_ptr<Point>> positions = {
              placement[op->getTargets().front()].currentPosition};
          for (const auto& v :
               layer.getExecutablesOfType(op->getType(), op->getNcontrols())) {
            const auto& op2 = *v->getOperation();
            if (checkApplicability(op2, placement) and
                op->getParameter() == op2->getParameter()) {
              updatePlacement(op2, placement);
              v->execute();
              positions.emplace_back(
                  placement[op2->getTargets().front()].currentPosition);
            }
          }
          mappedQc.emplaceBack<NALocalOperation>(
              OpType{op->getType(), op->getNcontrols()}, op->getParameter(),
              positions);
        }
        it = (*executableSet)->begin();
      } else {
        ++it;
      }
    }
    it = (*executableSet)->begin();
    if (it == (*executableSet)->end()) {
      break;
    }
    if (config.getMethod() == NaMappingMethod::NAIVE) {
      const std::unique_ptr<qc::Operation>& op = *(*it)->getOperation();
      if (op->getType() != qc::OpType::Z or op->getNtargets() != 1 or
          op->getNcontrols() != 1) {
        throw std::logic_error(
            "Other gates than cz are not supported for mapping yet.");
      }
      (*it)->execute();
      const auto&  q1     = op->getTargets().front();
      const auto&  q2     = op->getControls().begin()->qubit;
      Point        start  = *placement[q1].currentPosition;
      Point        end    = start;
      const Point& target = arch.getPositionOfSite(arch.getSitesInZone(
          *arch.getPropertiesOfOperation({op->getType(), 1}).zones.begin())[0]);
      end.x += d;
      mappedQc.emplaceBack<NAShuttlingOperation>(
          LOAD,
          std::vector<std::shared_ptr<Point>>{std::make_shared<Point>(start)},
          std::vector<std::shared_ptr<Point>>{std::make_shared<Point>(end)});
      start = end;
      end   = target;
      end.x += d;
      mappedQc.emplaceBack<NAShuttlingOperation>(
          MOVE,
          std::vector<std::shared_ptr<Point>>{std::make_shared<Point>(start)},
          std::vector<std::shared_ptr<Point>>{std::make_shared<Point>(end)});
      start = end;
      end.x -= d;
      mappedQc.emplaceBack<NAShuttlingOperation>(
          STORE,
          std::vector<std::shared_ptr<Point>>{std::make_shared<Point>(start)},
          std::vector<std::shared_ptr<Point>>{std::make_shared<Point>(end)});
      start = *placement[q2].currentPosition;
      end   = start;
      end.x += d;
      mappedQc.emplaceBack<NAShuttlingOperation>(
          LOAD,
          std::vector<std::shared_ptr<Point>>{std::make_shared<Point>(start)},
          std::vector<std::shared_ptr<Point>>{std::make_shared<Point>(end)});
      start = end;
      end   = target;
      end.y += d;
      mappedQc.emplaceBack<NAShuttlingOperation>(
          MOVE,
          std::vector<std::shared_ptr<Point>>{std::make_shared<Point>(start)},
          std::vector<std::shared_ptr<Point>>{std::make_shared<Point>(end)});
      mappedQc.emplaceBack<NAGlobalOperation>(OpType{qc::OpType::Z, 1});
      maxSeqWidth = 1UL;
      mappedQc.emplaceBack<NAShuttlingOperation>(
          MOVE,
          std::vector<std::shared_ptr<Point>>{std::make_shared<Point>(end)},
          std::vector<std::shared_ptr<Point>>{std::make_shared<Point>(start)});
      end = *placement[q2].currentPosition;
      mappedQc.emplaceBack<NAShuttlingOperation>(
          STORE,
          std::vector<std::shared_ptr<Point>>{std::make_shared<Point>(start)},
          std::vector<std::shared_ptr<Point>>{std::make_shared<Point>(end)});
      start = target;
      end   = start;
      end.x += d;
      mappedQc.emplaceBack<NAShuttlingOperation>(
          LOAD,
          std::vector<std::shared_ptr<Point>>{std::make_shared<Point>(start)},
          std::vector<std::shared_ptr<Point>>{std::make_shared<Point>(end)});
      start = end;
      end   = *placement[q1].currentPosition;
      end.x += d;
      mappedQc.emplaceBack<NAShuttlingOperation>(
          MOVE,
          std::vector<std::shared_ptr<Point>>{std::make_shared<Point>(start)},
          std::vector<std::shared_ptr<Point>>{std::make_shared<Point>(end)});
      start = end;
      end.x -= d;
      mappedQc.emplaceBack<NAShuttlingOperation>(
          STORE,
          std::vector<std::shared_ptr<Point>>{std::make_shared<Point>(start)},
          std::vector<std::shared_ptr<Point>>{std::make_shared<Point>(end)});
    } else if (config.getMethod() == NaMappingMethod::SMART) {
      // 2. when no such gates are left, extract an interaction graph of gates
      //    of the same type and two targets, i.e. cz gates
      const auto& graph = layer.constructInteractionGraph(qc::OpType::Z, 1);
      if (graph.getNVertices() == 0) {
        throw std::logic_error(
            "Other gates than cz are not supported for mapping yet.");
        // TODO: support other gates than cz
      }
      const auto& sequence = NAGraphAlgorithms::computeSequence(graph);
      const auto& moveable = sequence.first;
      const auto& fixed    = sequence.second;
      // 3. move the atoms accordingly and execute the gates
      // pick up the fixed atoms and move them to the interaction zone
      // get a vector of the fixed atoms ordered by their initial position
      // from left to right
      maxSeqWidth = std::accumulate(
          fixed.cbegin(), fixed.cend(), maxSeqWidth,
          [&](const auto max, const auto& q) {
            return std::max(max, static_cast<std::size_t>(q.second));
          });
      const Zone interactionZone =
          *arch.getPropertiesOfOperation({qc::OpType::Z, 1}).zones.begin();
      const auto sites = arch.getSitesInRow(interactionZone, 0);
      std::unordered_map<qc::Qubit, Point> finalFixedPositions;
      for (const auto& [q, x] : fixed) {
        if (static_cast<std::size_t>(x) >= sites.size()) {
          throw std::logic_error(
              "Target site in interaction zone is out of bounds. Possible "
              "reason "
              "for this error: Entangling zone is not wide enough.");
        }
        finalFixedPositions.emplace(
            q, arch.getPositionOfSite(sites[static_cast<std::size_t>(x)]));
      }
      shuttle(initialFreeSites, currentFreeSites, placement, currentlyShuttling,
              finalFixedPositions, true);
      // -----------------------------------------------------------------
      std::unordered_map<qc::Qubit, Point> finalMoveablePositions;
      for (const auto& [q, x] : moveable[0]) {
        finalMoveablePositions.emplace(
            q, arch.getPositionOfSite(sites[static_cast<std::size_t>(x)]));
      }
      shuttle(initialFreeSites, currentFreeSites, placement, currentlyShuttling,
              finalMoveablePositions, false);
      // ------------------------------------------------------------------
      // 4. Apply the cz gates
      std::vector<qc::Qubit> moveableOrdered;
      std::transform(moveable[0].cbegin(), moveable[0].cend(),
                     back_inserter(moveableOrdered),
                     [](const auto& v) { return v.first; });
      std::sort(moveableOrdered.begin(), moveableOrdered.end(),
                [&](const auto& a, const auto& b) {
                  return moveable[0].at(a) < moveable[0].at(b);
                });
      std::vector<std::shared_ptr<Point>> startMoveable;
      std::vector<std::shared_ptr<Point>> endMoveable;
      std::transform(moveableOrdered.cbegin(), moveableOrdered.cend(),
                     std::back_inserter(startMoveable), [&](const auto& q) {
                       return placement[q].currentPosition;
                     });
      for (const auto& timeframe : moveable) {
        for (const auto q : moveableOrdered) {
          const auto x = timeframe.at(q);
          if (currentlyShuttling.find(q) == currentlyShuttling.cend()) {
            std::stringstream ss;
            ss << "Atom " << q << " was unexpectedly not picked up.";
            throw std::logic_error(ss.str());
          }
          if (x >= 0 and static_cast<std::size_t>(x) < sites.size()) {
            const auto pos =
                arch.getPositionOfSite(sites[static_cast<std::size_t>(x)]);
            placement[q].currentPosition =
                std::make_shared<Point>(pos.x, pos.y + d);
          } else if (x < 0) {
            const auto pos = arch.getPositionOfSite(sites[0]);
            placement[q].currentPosition =
                std::make_shared<Point>(pos.x + (x * dx), pos.y + d);
          } else { // x >= sites.size()
            const auto pos = arch.getPositionOfSite(sites[sites.size() - 1]);
            placement[q].currentPosition =
                std::make_shared<Point>(pos.x + (x * dx), pos.y + d);
          }
          endMoveable.emplace_back(placement[q].currentPosition);
        }
        mappedQc.emplaceBack<NAShuttlingOperation>(MOVE, startMoveable,
                                                   endMoveable);
        mappedQc.emplaceBack<NAGlobalOperation>(OpType{qc::OpType::Z, 1});
        for (const auto& q : moveableOrdered) {
          for (const auto& [p, _] : fixed) {
            const auto qPos = *placement[q].currentPosition;
            const auto pPos = *placement[p].currentPosition;
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
      const auto storageZone = initialZones.front();
      shuttle(initialFreeSites, currentFreeSites, placement, currentlyShuttling,
              moveableOrdered, storageZone);
      // -------------------------------------------------------------
      std::vector<qc::Qubit> fixedVector;
      std::transform(fixed.cbegin(), fixed.cend(), back_inserter(fixedVector),
                     [](const auto& v) { return v.first; });
      shuttle(initialFreeSites, currentFreeSites, placement, currentlyShuttling,
              fixedVector, storageZone);
    } else {
      throw std::logic_error("NA mapping method not implemented.");
    }
    // -------------------------------------------------------------
    it = (*executableSet)->begin();
  }
  for (auto& p : placement) {
    if (p.positionStatus == Atom::PositionStatus::UNDEFINED) {
      // find next free site in first zone
      std::vector<Index> possibleSites;
      for (const auto& z : p.zones) {
        for (const auto& s : arch.getSitesInZone(z)) {
          possibleSites.emplace_back(s);
        }
      }
      const auto& freeSite =
          std::find_if(possibleSites.cbegin(), possibleSites.cend(),
                       [&](const auto& s) { return initialFreeSites[s]; });
      *p.initialPosition          = arch.getPositionOfSite(*freeSite);
      p.positionStatus            = Atom::PositionStatus::DEFINED;
      initialFreeSites[*freeSite] = false;
      currentFreeSites[*freeSite] = false;
    }
    mappedQc.emplaceInitialPosition(p.initialPosition);
  }
  //========================= END MAPPING =========================
  // get end time
  auto startPostprocess = std::chrono::high_resolution_clock::now();
  postprocess();
  auto end = std::chrono::high_resolution_clock::now();
  // build remaining statistics
  stats.numInitialGates    = qc.getNops();
  stats.numEntanglingGates = static_cast<std::size_t>(
      std::count_if(qc.cbegin(), qc.cend(), [](const auto& op) {
        return (op->getType() == qc::OpType::Z and
                op->getType() == qc::OpType::X) or
               op->getNcontrols() > 0;
      }));
  stats.initialDepth   = qc.getDepth();
  stats.numMappedGates = mappedQc.size();
  stats.numQubits      = nqubits;
  stats.maxSeqWidth    = config.getPatchCols() * maxSeqWidth;
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
