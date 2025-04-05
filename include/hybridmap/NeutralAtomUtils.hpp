//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "hybridmap/NeutralAtomDefinitions.hpp"
#include "ir/Definitions.hpp"
#include "ir/operations/AodOperation.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace na {

// Enums for the different initial mappings strategies
enum InitialCoordinateMapping : uint8_t { Trivial, Random };
enum InitialMapping : uint8_t { Identity };

[[maybe_unused]] static InitialCoordinateMapping
initialCoordinateMappingFromString(
    const std::string& initialCoordinateMapping) {
  if (initialCoordinateMapping == "trivial" ||
      initialCoordinateMapping == "0") {
    return InitialCoordinateMapping::Trivial;
  }
  if (initialCoordinateMapping == "random" || initialCoordinateMapping == "1") {
    return InitialCoordinateMapping::Random;
  }
  throw std::invalid_argument("Invalid initial coordinate mapping value: " +
                              initialCoordinateMapping);
}

[[maybe_unused]] static InitialMapping
initialMappingFromString(const std::string& initialMapping) {
  if (initialMapping == "identity" || initialMapping == "0") {
    return InitialMapping::Identity;
  }
  throw std::invalid_argument("Invalid initial mapping value: " +
                              initialMapping);
}

/**
 * @brief Helper class to represent a direction in x and y coordinates.
 * @details The boolean value corresponds to right/left and down/up.
 */
struct Direction {
  bool x;
  bool y;

  [[maybe_unused]] Direction(bool xDir, bool yDir) : x(xDir), y(yDir) {}
  Direction(qc::fp deltaX, qc::fp deltaY) : x(deltaX >= 0), y(deltaY >= 0) {}

  [[nodiscard]] bool operator==(const Direction& other) const {
    return x == other.x && y == other.y;
  }
  [[nodiscard]] bool operator!=(const Direction& other) const {
    return !(*this == other);
  }
  [[nodiscard]] int32_t getSignX() const { return x ? 1 : -1; }
  [[nodiscard]] int32_t getSignY() const { return y ? 1 : -1; }
  [[nodiscard]] int32_t getSign(Dimension dim) const {
    return dim == Dimension::X ? getSignX() : getSignY();
  }
};

/**
 * @brief Helper class to represent a move of an atom from one position to
 * another.
 * @details Each move consists in a start and end coordinate and the direction.
 */
struct MoveVector {
  qc::fp xStart;
  qc::fp yStart;
  qc::fp xEnd;
  qc::fp yEnd;
  Direction direction;

  MoveVector(qc::fp xstart, qc::fp ystart, qc::fp xend, qc::fp yend)
      : xStart(xstart), yStart(ystart), xEnd(xend), yEnd(yend),
        direction(xend - xstart, yend - ystart) {}
  MoveVector(std::int64_t xstart, std::int64_t ystart, std::int64_t xend,
             std::int64_t yend)
      : xStart(static_cast<qc::fp>(xstart)),
        yStart(static_cast<qc::fp>(ystart)), xEnd(static_cast<qc::fp>(xend)),
        yEnd(static_cast<qc::fp>(yend)),
        direction(static_cast<qc::fp>(xend - xstart),
                  static_cast<qc::fp>(yend - ystart)) {}

  [[nodiscard]] [[maybe_unused]] bool
  sameDirection(const MoveVector& other) const {
    return direction == other.direction;
  }
  [[nodiscard]] qc::fp getLength() const {
    return std::sqrt(std::pow(xEnd - xStart, 2) + std::pow(yEnd - yStart, 2));
  }
  [[nodiscard]] bool overlap(const MoveVector& other) const;
  [[nodiscard]] bool include(const MoveVector& other) const;
};

/**
 * @brief Helper class to manage multiple atom moves which belong together.
 * @details E.g. a move-away combined with the actual move. These are combined
 * in a MoveComb to facilitate the cost calculation.
 */
struct MoveComb {
  std::vector<AtomMove> moves;
  qc::fp cost = std::numeric_limits<qc::fp>::max();

  MoveComb(std::vector<AtomMove> mov, const qc::fp c)
      : moves(std::move(mov)), cost(c) {}
  MoveComb(AtomMove mov, const qc::fp c) : moves({std::move(mov)}), cost(c) {}

  MoveComb() = default;
  explicit MoveComb(std::vector<AtomMove> mov) : moves(std::move(mov)) {}
  explicit MoveComb(AtomMove mov) : moves({std::move(mov)}) {}

  /**
   * @brief Get the first move of the combination
   * @return The first move of the combination
   */
  [[nodiscard]] AtomMove getFirstMove() const { return moves.front(); }

  /**
   * @brief Get the last move of the combination
   * @return The last move of the combination
   */
  [[nodiscard]] [[maybe_unused]] AtomMove getLastMove() const {
    return moves.back();
  }

  // implement == operator for AtomMove
  [[nodiscard]] bool operator==(const MoveComb& other) const {
    return moves == other.moves;
  }
  [[nodiscard]] bool operator!=(const MoveComb& other) const {
    return !(*this == other);
  }

  /**
   * @brief Append a single move to the end of the combination.
   * @param addMove The move to append
   */
  void append(AtomMove addMove) {
    moves.emplace_back(addMove);
    cost = std::numeric_limits<qc::fp>::max();
  }
  /**
   * @brief Append all moves of another combination to the end of this one.
   * @param addMoveComb The other combination to append
   */
  void append(const MoveComb& addMoveComb) {
    moves.insert(moves.end(), addMoveComb.moves.begin(),
                 addMoveComb.moves.end());
    cost = std::numeric_limits<qc::fp>::max();
  }
  [[nodiscard]] size_t size() const { return moves.size(); }
  [[nodiscard]] bool empty() const { return moves.empty(); }
};

/**
 * @brief Helper class to manage multiple move combinations.
 */
struct MoveCombs {
  std::vector<MoveComb> moveCombs;

  MoveCombs() = default;
  explicit MoveCombs(std::vector<MoveComb> combs)
      : moveCombs(std::move(combs)) {}

  [[nodiscard]] bool empty() const { return moveCombs.empty(); }
  [[nodiscard]] size_t size() const { return moveCombs.size(); }

  // define iterators that iterate over the moveCombs vector
  using iterator = std::vector<MoveComb>::iterator;
  using const_iterator = std::vector<MoveComb>::const_iterator;
  iterator begin() { return moveCombs.begin(); }
  iterator end() { return moveCombs.end(); }
  [[nodiscard]] const_iterator begin() const { return moveCombs.cbegin(); }
  [[nodiscard]] const_iterator end() const { return moveCombs.cend(); }

  /**
   * @brief Add a move combination to the list of move combinations.
   * @param moveComb The move combination to add.
   */
  void addMoveComb(const MoveComb& moveComb);
  /**
   * @brief Add all move combinations of another MoveCombs object to the list of
   * move combinations.
   * @param otherMoveCombs The other MoveCombs object to add.
   */
  void addMoveCombs(const MoveCombs& otherMoveCombs);
  /**
   * @brief Remove all move combinations that are longer than the shortest move
   * combination.
   */
  void removeLongerMoveCombs();
};

/**
 * @brief Helper struct to store the position of a multi qubit gate and the
 * number of moves needed to execute it.
 */
struct MultiQubitMovePos {
  CoordIndices coords;
  size_t nMoves{0};
};

} // namespace na
