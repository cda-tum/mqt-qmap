//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "Definitions.hpp"
#include "hybridmap/NeutralAtomDefinitions.hpp"
#include "utils.hpp"

namespace qc {

/**
 * @brief Symmetric matrix class with same number of rows and columns that
 * allows access by row and column but uses less memory than a full matrix
 */
class SymmetricMatrix {
private:
  std::vector<std::vector<fp>> data;
  uint32_t                     size = 0;

public:
  // Constructors
  SymmetricMatrix() = default;
  explicit SymmetricMatrix(uint32_t size) : size(size) {
    data.resize(size);
    for (uint32_t i = 0; i < size; ++i) {
      data[i].resize(i + 1);
    }
  }

  SymmetricMatrix(uint32_t size, fp value) : size(size) {
    data.resize(size);
    for (uint32_t i = 0; i < size; ++i) {
      data[i].resize(i + 1, value);
    }
  }

  inline fp& operator()(uint32_t row, uint32_t col) {
    if (row < col) {
      return data[col][row];
    }
    return data[row][col];
  }

  [[nodiscard]] inline fp operator()(uint32_t row, uint32_t col) const {
    if (row < col) {
      return data[col][row];
    }
    return data[row][col];
  }

  [[nodiscard]] inline uint32_t getSize() const { return size; }
};

// Enums for the different initial mappings strategies
enum InitialCoordinateMapping { Trivial, Random };
enum InitialMapping { Identity };

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
 * @brief Helperclass to represent a direction in x and y coordinates.
 * @details The boolean value corresponds to right/left and down/up.
 */
struct Direction {
  bool x;
  bool y;

  Direction(bool x, bool y) : x(x), y(y) {}
  Direction(fp deltaX, fp deltaY) : x(deltaX >= 0), y(deltaY >= 0) {}

  [[nodiscard]] inline bool operator==(const Direction& other) const {
    return x == other.x && y == other.y;
  }
  [[nodiscard]] inline bool operator!=(const Direction& other) const {
    return !(*this == other);
  }
  [[nodiscard]] inline int32_t getSignX() const { return x ? 1 : -1; }
  [[nodiscard]] inline int32_t getSignY() const { return y ? 1 : -1; }
};

/**
 * @brief Helperclass to represent a move of an atom from one position to
 * another.
 * @details Each move consists in a start and end coordinate and the direction.
 */
struct MoveVector {
  fp        xStart;
  fp        yStart;
  fp        xEnd;
  fp        yEnd;
  Direction direction;

  MoveVector(fp xStart, fp yStart, fp xEnd, fp yEnd)
      : xStart(xStart), yStart(yStart), xEnd(xEnd), yEnd(yEnd),
        direction(xEnd - xStart, yEnd - yStart) {}

  [[nodiscard]] inline bool sameDirection(const MoveVector& other) const {
    return direction == other.direction;
  }
  [[nodiscard]] inline fp getLength() const {
    return std::sqrt(std::pow(xEnd - xStart, 2) + std::pow(yEnd - yStart, 2));
  }
  [[nodiscard]] bool overlap(const MoveVector& other) const;
  [[nodiscard]] bool include(const MoveVector& other) const;
};

/**
 * @brief Helperclass to manage multiple atom moves which belong together.
 * @details E.g. a move-away combined with the actual move. These are combined
 * in a MoveComb to facilitate the cost calculation.
 */
struct MoveComb {
  std::vector<AtomMove> moves;
  fp                    cost = std::numeric_limits<fp>::quiet_NaN();

  MoveComb(std::vector<AtomMove> moves, fp cost)
      : moves(std::move(moves)), cost(cost) {}
  MoveComb(AtomMove move, fp cost)
      : moves(std::vector<AtomMove>{std::move(move)}), cost(cost) {}

  MoveComb() = default;
  MoveComb(std::vector<AtomMove> moves) : moves(std::move(moves)) {}
  MoveComb(AtomMove move) : moves(std::vector<AtomMove>{std::move(move)}) {}

  /**
   * @brief Get the first move of the combination
   * @return The first move of the combination
   */
  [[nodiscard]] AtomMove inline getFirstMove() const { return *moves.begin(); }

  /**
   * @brief Get the last move of the combination
   * @return The last move of the combination
   */
  [[nodiscard]] AtomMove inline getLastMove() const { return *moves.rbegin(); }

  // implement == operator for AtomMove
  [[nodiscard]] inline bool operator==(const MoveComb& other) const {
    return moves == other.moves;
  }
  [[nodiscard]] inline bool operator!=(const MoveComb& other) const {
    return !(*this == other);
  }

  /**
   * @brief Append a single move to the end of the combination.
   * @param addMove The move to append
   */
  inline void append(AtomMove& addMove) {
    moves.push_back(addMove);
    cost = std::numeric_limits<fp>::quiet_NaN();
  }
  /**
   * @brief Append all moves of another combination to the end of this one.
   * @param addMoveComb The other combination to append
   */
  inline void append(const MoveComb& addMoveComb) {
    moves.insert(moves.end(), addMoveComb.moves.begin(),
                 addMoveComb.moves.end());
    cost = std::numeric_limits<fp>::quiet_NaN();
  }
  [[nodiscard]] inline size_t size() const { return moves.size(); }
  [[nodiscard]] inline bool   empty() const { return moves.empty(); }
};

/**
 * @brief Helperclass to manage multiple move combinations.
 */
struct MoveCombs {
  std::vector<MoveComb> moveCombs;

  MoveCombs() = default;
  MoveCombs(std::vector<MoveComb> moveCombs) : moveCombs(moveCombs) {}

  [[nodiscard]] bool inline empty() const { return moveCombs.empty(); }
  [[nodiscard]] size_t inline size() const { return moveCombs.size(); }

  // define iterators that iterate over the moveCombs vector
  typedef std::vector<MoveComb>::iterator       iterator;
  typedef std::vector<MoveComb>::const_iterator const_iterator;
  iterator       begin() { return moveCombs.begin(); }
  iterator       end() { return moveCombs.end(); }
  const_iterator begin() const { return moveCombs.begin(); }
  const_iterator end() const { return moveCombs.end(); }

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
 * @brief Helperclass to facilitate the handling of x/y coordinates.
 */
class Coordinate {
protected:
  std::uint32_t x;
  std::uint32_t y;

public:
  Coordinate() = default;
  Coordinate(CoordIndex x, CoordIndex y) : x(x), y(y) {}

  [[nodiscard]] inline CoordIndex getX() const { return x; }
  [[nodiscard]] inline CoordIndex getY() const { return y; }
  [[nodiscard]] inline fp         getXfp() const { return static_cast<fp>(x); }
  [[nodiscard]] inline fp         getYfp() const { return static_cast<fp>(y); }
  [[nodiscard]] inline std::pair<CoordIndex, CoordIndex> getXY() const {
    return std::make_pair(x, y);
  }
  [[nodiscard]] inline fp getEuclidianDistance(const Coordinate& c) const {
    return std::sqrt(std::pow(static_cast<fp>(x) - static_cast<fp>(c.x), 2) +
                     std::pow(static_cast<fp>(y) - static_cast<fp>(c.y), 2));
  }
  [[nodiscard]] static inline bool sameX(const Coordinate& c1,
                                         const Coordinate& c2) {
    return c1.x == c2.x;
  }
  [[nodiscard]] static inline bool sameY(const Coordinate& c1,
                                         const Coordinate& c2) {
    return c1.y == c2.y;
  }
  [[nodiscard]] static inline bool sameXorY(const Coordinate& c1,
                                            const Coordinate& c2) {
    return sameX(c1, c2) || sameY(c1, c2);
  }

  [[nodiscard]] inline CoordIndex
  getManhattanDistanceX(const Coordinate& c) const {
    if (x > c.x) {
      return x - c.x;
    } else {
      return c.x - x;
    }
  }
  [[nodiscard]] inline CoordIndex
  getManhattanDistanceY(const Coordinate& c) const {
    if (y > c.y) {
      return y - c.y;
    } else {
      return c.y - y;
    }
  }
};

} // namespace qc
