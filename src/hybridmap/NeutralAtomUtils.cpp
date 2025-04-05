//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "hybridmap/NeutralAtomUtils.hpp"

#include "ir/Definitions.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace na {

bool MoveVector::overlap(const MoveVector& other) const {
  // do not consider direction for overlap
  const auto firstStartX = std::min(xStart, xEnd);
  const auto firstEndX = std::max(xStart, xEnd);
  const auto secondStartX = std::min(other.xStart, other.xEnd);
  const auto secondEndX = std::max(other.xStart, other.xEnd);
  const auto firstStartY = std::min(yStart, yEnd);
  const auto firstEndY = std::max(yStart, yEnd);
  const auto secondStartY = std::min(other.yStart, other.yEnd);
  const auto secondEndY = std::max(other.yStart, other.yEnd);

  // need to compute all combinations, as sometimes the start and end x/y points
  // are the same
  auto overlapXFirstStart =
      firstStartX >= secondStartX && firstStartX <= secondEndX;
  auto overlapXFirstEnd = firstEndX >= secondStartX && firstEndX <= secondEndX;
  auto overlapXSecondStart =
      secondStartX >= firstStartX && secondStartX <= firstEndX;
  auto overlapXSecondEnd = secondEndX >= firstStartX && secondEndX <= firstEndX;
  auto overlapYFirstStart =
      firstStartY >= secondStartY && firstStartY <= secondEndY;
  auto overlapYFirstEnd = firstEndY >= secondStartY && firstEndY <= secondEndY;
  auto overlapYSecondStart =
      secondStartY >= firstStartY && secondStartY <= firstEndY;
  auto overlapYSecondEnd = secondEndY >= firstStartY && secondEndY <= firstEndY;

  return (overlapXFirstStart || overlapXFirstEnd || overlapXSecondStart ||
          overlapXSecondEnd || overlapYFirstStart || overlapYFirstEnd ||
          overlapYSecondStart || overlapYSecondEnd);
}

bool MoveVector::include(const MoveVector& other) const {
  const auto firstStartX = std::min(xStart, xEnd);
  const auto firstEndX = std::max(xStart, xEnd);
  const auto secondStartX = std::min(other.xStart, other.xEnd);
  const auto secondEndX = std::max(other.xStart, other.xEnd);
  const auto firstStartY = std::min(yStart, yEnd);
  const auto firstEndY = std::max(yStart, yEnd);
  const auto secondStartY = std::min(other.yStart, other.yEnd);
  const auto secondEndY = std::max(other.yStart, other.yEnd);

  const auto includeX =
      (secondStartX < firstStartX) && (firstEndX < secondEndX);
  const auto includeY =
      (secondStartY < firstStartY) && (firstEndY < secondEndY);

  return includeX || includeY;
}

void MoveCombs::addMoveComb(const MoveComb& otherMove) {
  for (auto& comb : moveCombs) {
    if (comb == otherMove) {
      comb.cost = std::numeric_limits<qc::fp>::max();
      return;
    }
  }
  moveCombs.emplace_back(otherMove);
}

void MoveCombs::addMoveCombs(const MoveCombs& otherMoveCombs) {
  for (const auto& otherMove : otherMoveCombs.moveCombs) {
    addMoveComb(otherMove);
  }
}

void MoveCombs::removeLongerMoveCombs() {
  size_t minSize = std::numeric_limits<uint32_t>::max();
  for (const auto& comb : moveCombs) {
    minSize = std::min(minSize, comb.size());
  }
  for (auto it = moveCombs.begin(); it != moveCombs.end();) {
    if (it->size() > minSize) {
      it = moveCombs.erase(it);
    } else {
      ++it;
    }
  }
}

} // namespace na
