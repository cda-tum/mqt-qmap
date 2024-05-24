//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "hybridmap/HardwareQubits.hpp"

#include "hybridmap/Mapping.hpp"

namespace qc {
void HardwareQubits::initTrivialSwapDistances() {
  swapDistances = SymmetricMatrix(arch.getNqubits());
  for (uint32_t i = 0; i < arch.getNqubits(); ++i) {
    for (uint32_t j = 0; j < i; ++j) {
      swapDistances(i, j) =
          arch.getSwapDistance(hwToCoordIdx.at(i), hwToCoordIdx.at(j));
    }
  }
}

void HardwareQubits::initNearbyQubits() {
  for (uint32_t i = 0; i < arch.getNqubits(); ++i) {
    computeNearbyQubits(i);
  }
}

void HardwareQubits::computeSwapDistance(HwQubit q1, HwQubit q2) {
  std::queue<HwQubit>  q;
  std::vector<bool>    visited(swapDistances.getSize(), false);
  std::vector<HwQubit> parent(swapDistances.getSize(), q2);

  q.push(q1);
  visited[q1] = true;
  parent[q1]  = q1;
  bool found  = false;
  while (!q.empty() && !found) {
    auto current = q.front();
    q.pop();
    for (const auto& nearbyQubit : nearbyQubits.at(current)) {
      if (!visited[nearbyQubit]) {
        q.push(nearbyQubit);
        visited[nearbyQubit] = true;
        parent[nearbyQubit]  = current;
        if (nearbyQubit == q2) {
          found = true;
          break;
        }
      }
    }
  }
  if (!found) {
    swapDistances(q1, q2) = std::numeric_limits<fp>::infinity();
    return;
  }
  // recreate path
  std::vector<HwQubit> path;
  auto                 current = q2;
  while (current != q1) {
    path.emplace_back(current);
    current = parent[current];
  }
  path.emplace_back(q1);
  // update swap distances along path
  for (uint32_t start = 0; start < path.size() - 1; ++start) {
    for (uint32_t end = start + 1; end < path.size(); ++end) {
      swapDistances(path[start], path[end]) = end - start - 1;
    }
  }
  //
  //  swapDistances(q1, q2) = static_cast<fp>(path.size() - 2);
}

void HardwareQubits::resetSwapDistances() {
  // TODO Improve to only reset the swap distances necessary (use a breadth
  // first search)
  swapDistances = SymmetricMatrix(arch.getNqubits(), -1);
}

void HardwareQubits::move(HwQubit hwQubit, CoordIndex newCoord) {
  if (newCoord >= arch.getNpositions()) {
    throw std::runtime_error("Invalid coordinate");
  }
  // check if new coordinate is already occupied
  for (const auto& [qubit, coord] : hwToCoordIdx) {
    if (coord == newCoord) {
      throw std::runtime_error("Coordinate already occupied");
    }
  }

  // remove qubit from old nearby qubits
  auto prevNearbyQubits = nearbyQubits.at(hwQubit);
  for (const auto& qubit : prevNearbyQubits) {
    nearbyQubits.at(qubit).erase(std::find(
        nearbyQubits.at(qubit).begin(), nearbyQubits.at(qubit).end(), hwQubit));
  }
  // move qubit and compute new nearby qubits
  hwToCoordIdx.at(hwQubit) = newCoord;
  computeNearbyQubits(hwQubit);

  // add qubit to new nearby qubits
  auto newNearbyQubits = nearbyQubits.at(hwQubit);
  for (const auto& qubit : newNearbyQubits) {
    nearbyQubits.at(qubit).insert(hwQubit);
  }

  // update/reset swap distances
  resetSwapDistances();
}

std::vector<Swap> HardwareQubits::getNearbySwaps(qc::HwQubit q) const {
  std::vector<Swap> swaps;
  swaps.reserve(nearbyQubits.size());
  for (const auto& nearbyQubit : nearbyQubits.at(q)) {
    swaps.emplace_back(q, nearbyQubit);
  }
  return swaps;
}

void HardwareQubits::computeNearbyQubits(qc::HwQubit q) {
  std::set<HwQubit> newNearbyQubits;
  auto              coordQ = hwToCoordIdx.at(q);
  for (const auto& coord : hwToCoordIdx) {
    if (coord.first == q) {
      continue;
    }
    if (arch.getEuclidianDistance(coordQ, coord.second) <=
        arch.getInteractionRadius()) {
      newNearbyQubits.insert(coord.first);
    }
  }
  nearbyQubits.insert_or_assign(q, newNearbyQubits);
}

fp HardwareQubits::getAllToAllSwapDistance(std::set<HwQubit>& qubits) {
  // two qubit gates
  if (qubits.size() == 2) {
    auto it = qubits.begin();
    auto q1 = *it;
    auto q2 = *(++it);
    return getSwapDistance(q1, q2);
  }
  // for n > 2 all qubits need to be within the interaction radius of each other
  fp totalDistance = 0;
  for (auto it1 = qubits.begin(); it1 != qubits.end(); ++it1) {
    for (auto it2 = std::next(it1); it2 != qubits.end(); ++it2) {
      totalDistance += getSwapDistance(*it1, *it2);
    }
  }
  return totalDistance;
}

std::set<HwQubit>
HardwareQubits::getBlockedQubits(const std::set<HwQubit>& qubits) {
  std::set<HwQubit> blockedQubits;
  for (const auto& qubit : qubits) {
    for (uint32_t i = 0; i < arch.getNqubits(); ++i) {
      if (i == qubit) {
        continue;
      }
      // TODO improve by using the nearby coords as a preselection
      auto const distance =
          arch.getEuclidianDistance(hwToCoordIdx.at(qubit), hwToCoordIdx.at(i));
      if (distance <= arch.getBlockingFactor() * arch.getInteractionRadius()) {
        blockedQubits.insert(i);
      }
    }
  }
  return blockedQubits;
}

std::set<CoordIndex> HardwareQubits::getNearbyFreeCoordinates(HwQubit q) {
  return getNearbyFreeCoordinatesByCoord(hwToCoordIdx.at(q));
}

std::set<CoordIndex>
HardwareQubits::getNearbyFreeCoordinatesByCoord(qc::CoordIndex idx) {
  std::set<CoordIndex> nearbyFreeCoordinates;
  for (auto const& coordIndex : this->arch.getNearbyCoordinates(idx)) {
    if (!this->isMapped(coordIndex)) {
      nearbyFreeCoordinates.insert(coordIndex);
    }
  }
  return nearbyFreeCoordinates;
}

std::set<CoordIndex>
HardwareQubits::getNearbyOccupiedCoordinates(qc::HwQubit q) {
  auto nearbyHwQubits = this->getNearbyQubits(q);
  return this->getCoordIndices(nearbyHwQubits);
}

std::set<CoordIndex>
HardwareQubits::getNearbyOccupiedCoordinatesByCoord(qc::CoordIndex idx) {
  auto nearbyHwQubits = this->getNearbyQubits(this->getHwQubit(idx));
  return this->getCoordIndices(nearbyHwQubits);
}

std::vector<CoordIndex>
HardwareQubits::findClosestFreeCoord(CoordIndex coord, Direction direction,
                                     const CoordIndices& excludeCoord) {
  // return the closest free coord in general
  // and the closest free coord in the given direction
  std::vector<CoordIndex> closestFreeCoords;
  std::queue<CoordIndex>  queue;
  queue.push(coord);
  std::set<CoordIndex> visited;
  visited.insert(coord);
  bool foundClosest = false;
  while (!queue.empty()) {
    auto currentCoord = queue.front();
    queue.pop();
    auto nearbyCoords = this->arch.getNN(currentCoord);
    for (const auto& nearbyCoord : nearbyCoords) {
      if (std::find(visited.rbegin(), visited.rend(), nearbyCoord) ==
          visited.rend()) {
        visited.insert(nearbyCoord);
        if (!this->isMapped(nearbyCoord) &&
            std::find(excludeCoord.begin(), excludeCoord.end(), nearbyCoord) ==
                excludeCoord.end()) {
          if (!foundClosest) {
            closestFreeCoords.push_back(nearbyCoord);
          }
          foundClosest = true;
          if (direction == arch.getVector(coord, nearbyCoord).direction) {
            closestFreeCoords.push_back(nearbyCoord);
            return closestFreeCoords;
          }
        } else {
          queue.push(nearbyCoord);
        }
      }
    }
  }
  return closestFreeCoords;
}

[[maybe_unused]] fp HardwareQubits::getSwapDistanceMove(CoordIndex idx,
                                                        HwQubit    target) {
  auto nearbyCoords    = this->arch.getNearbyCoordinates(idx);
  fp   minSwapDistance = std::numeric_limits<fp>::infinity();
  if (nearbyCoords.find(this->getCoordIndex(target)) != nearbyCoords.end()) {
    return 0;
  }
  for (const auto& nearbyCoord : nearbyCoords) {
    if (this->isMapped(nearbyCoord)) {
      auto swapDistance =
          this->getSwapDistance(target, this->getHwQubit(nearbyCoord)) + 1;
      if (swapDistance < minSwapDistance) {
        minSwapDistance = swapDistance;
      }
    }
  }
  return minSwapDistance;
}

} // namespace qc
