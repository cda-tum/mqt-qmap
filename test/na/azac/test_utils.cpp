#include "na/azac/Utils.hpp"

#include <cmath>
#include <cstddef>
#include <gtest/gtest.h>
#include <optional>
#include <utility>
#include <vector>

namespace na {

TEST(TestUtils, Distance) {
  const std::pair<std::size_t, std::size_t> a{0, 0};
  const std::pair<std::size_t, std::size_t> b{0, 1};
  const std::pair<std::size_t, std::size_t> c{1, 1};
  const std::pair<std::size_t, std::size_t> d{1, 0};
  EXPECT_DOUBLE_EQ(distance(a, a), 0);
  EXPECT_DOUBLE_EQ(distance(a, b), 1);
  EXPECT_DOUBLE_EQ(distance(a, c), std::sqrt(2));
  EXPECT_DOUBLE_EQ(distance(a, d), 1);
}

TEST(TestUtils, MaximumBipartiteMatching) {
  // We consider the following bipartite graph, where the nodes in the upper row
  // are the sources, and the nodes in the lower row are the sinks.
  //   ┌───┐ ┌───┐ ┌───┐ ┌───┐
  //   │ 0 │ │ 1 │ │ 2 │ │ 3 │ <-- SOURCES
  //   └─┬─┘ └─┬─┘ └─┬─┘ └─┬─┘
  //     │╲     ╲   ╱│╲   ╱│
  //     │  ╲     ╳  │  ╳  │
  //     │    ╲ ╱   ╲│╱   ╲│
  //   ┌─┴─┐ ┌─┴─┐ ┌─┴─┐ ┌─┴─┐
  //   │ 0 │ │ 1 │ │ 2 │ │ 3 │ <-- SINKS
  //   └───┘ └───┘ └───┘ └───┘
  const std::vector<std::vector<std::size_t>> sparseMatrix{/* 0 -> */ {0, 1},
                                                           /* 1 -> */ {2},
                                                           /* 2 -> */ {1, 2, 3},
                                                           /* 3 -> */ {2, 3}};
  const auto matching = maximumBipartiteMatching(sparseMatrix);
  // The result should be the following (unique) maximum matching:
  //   ┌───┐ ┌───┐ ┌───┐ ┌───┐
  //   │ 0 │ │ 1 │ │ 2 │ │ 3 │ <-- SOURCES
  //   └─┬─┘ └─┬─┘ └─┬─┘ └─┬─┘
  //     │      ╲   ╱      │
  //     │        ╳        │
  //     │      ╱   ╲      │
  //   ┌─┴─┐ ┌─┴─┐ ┌─┴─┐ ┌─┴─┐
  //   │ 0 │ │ 1 │ │ 2 │ │ 3 │ <-- SINKS
  //   └───┘ └───┘ └───┘ └───┘
  ASSERT_EQ(matching.size(), 4);
  EXPECT_EQ(matching[0], 0);
  EXPECT_EQ(matching[1], 2);
  EXPECT_EQ(matching[2], 1);
  EXPECT_EQ(matching[3], 3);
  const auto invMatching = maximumBipartiteMatching(sparseMatrix, true);
  ASSERT_EQ(invMatching.size(), 4);
  EXPECT_EQ(invMatching[0], 0);
  EXPECT_EQ(invMatching[1], 2);
  EXPECT_EQ(invMatching[2], 1);
  EXPECT_EQ(invMatching[3], 3);
  // We also test with the inverted graph, i.e., the sources and sinks are
  // labeled in reverse order, but sources stay sources and sinks stay sinks.
  const std::vector<std::vector<std::size_t>> inverseSparseMatrix{
      /* 0 -> */ {0, 1},
      /* 1 -> */ {0, 1, 2},
      /* 2 -> */ {1},
      /* 3 -> */ {2, 3}};
  const auto matchingOfInverse = maximumBipartiteMatching(inverseSparseMatrix);
  ASSERT_EQ(matchingOfInverse.size(), 4);
  EXPECT_EQ(matchingOfInverse[0], 0);
  EXPECT_EQ(matchingOfInverse[1], 2);
  EXPECT_EQ(matchingOfInverse[2], 1);
  EXPECT_EQ(matchingOfInverse[3], 3);
}

TEST(TestUtils, MinimumWeightFullBipartiteMatching) {
  {
    // We consider the following bipartite graph, where the nodes in the upper row
    // are the sources, and the nodes in the lower row are the sinks.
    //         ┌───┐ ┌───┐ ┌───┐
    //         │ 0 │ │ 1 │ │ 2 │ <-- SOURCES
    //         └─┬─┘ └─┬─┘ └─┬─┘
    //          ╱│╲3  ╱│╲4   │╲
    //       2╱  │  ╳  │4 ╲  │2 ╲3
    //      ╱   1│╱2  ╲│    ╲│    ╲
    //   ┌─┴─┐ ┌─┴─┐ ┌─┴─┐ ┌─┴─┐ ┌─┴─┐
    //   │ 0 │ │ 1 │ │ 2 │ │ 3 │ │ 4 │ <-- SINKS
    //   └───┘ └───┘ └───┘ └───┘ └───┘
    const std::vector<std::vector<std::optional<double>>> costMatrix{
      /* 0 -> */ {2, 1, 3, std::nullopt, std::nullopt},
      /* 1 -> */ {std::nullopt, 2, 4, 4, std::nullopt},
      /* 2 -> */ {std::nullopt, std::nullopt, std::nullopt, 2, 3}};
    // The result should be the following (unique) minimum weight full matching
    // and has weight 2 + 2 + 2 = 6:
    //         ┌───┐ ┌───┐ ┌───┐
    //         │ 0 │ │ 1 │ │ 2 │ <-- SOURCES
    //         └─┬─┘ └─┬─┘ └─┬─┘
    //          ╱     ╱      │
    //       2╱     ╱        │2
    //      ╱     ╱2         │
    //   ┌─┴─┐ ┌─┴─┐ ┌───┐ ┌─┴─┐ ┌───┐
    //   │ 0 │ │ 1 │ │ 2 │ │ 3 │ │ 4 │ <-- SINKS
    //   └───┘ └───┘ └───┘ └───┘ └───┘
    const auto matching = minimumWeightFullBipartiteMatching(costMatrix);
    ASSERT_EQ(matching.size(), 3);
    EXPECT_EQ(matching[0], 0);
    EXPECT_EQ(matching[1], 1);
    EXPECT_EQ(matching[2], 3);
  }
  {
    // We also consider the following bipartite graph that is the same graph as
    // the previous one, but with different weights:
    //         ┌───┐ ┌───┐ ┌───┐
    //         │ 0 │ │ 1 │ │ 2 │ <-- SOURCES
    //         └─┬─┘ └─┬─┘ └─┬─┘
    //          ╱│╲1  ╱│╲1   │╲
    //       3╱  │  ╳  │1 ╲  │1 ╲3
    //      ╱   3│╱2  ╲│    ╲│    ╲
    //   ┌─┴─┐ ┌─┴─┐ ┌─┴─┐ ┌─┴─┐ ┌─┴─┐
    //   │ 0 │ │ 1 │ │ 2 │ │ 3 │ │ 4 │ <-- SINKS
    //   └───┘ └───┘ └───┘ └───┘ └───┘
    const std::vector<std::vector<std::optional<double>>> costMatrix{
      /* 0 -> */ {3, 3, 1, std::nullopt, std::nullopt},
      /* 1 -> */ {std::nullopt, 2, 1, 1, std::nullopt},
      /* 2 -> */ {std::nullopt, std::nullopt, std::nullopt, 1, 3}};
    // The result should be the following (unique) minimum weight full matching
    // and has weight 1 + 2 + 1 = 4:
    //         ┌───┐ ┌───┐ ┌───┐
    //         │ 0 │ │ 1 │ │ 2 │ <-- SOURCES
    //         └─┬─┘ └─┬─┘ └─┬─┘
    //            ╲1  ╱      │
    //              ╳        │1
    //            ╱2  ╲      │
    //   ┌───┐ ┌─┴─┐ ┌─┴─┐ ┌─┴─┐ ┌───┐
    //   │ 0 │ │ 1 │ │ 2 │ │ 3 │ │ 4 │ <-- SINKS
    //   └───┘ └───┘ └───┘ └───┘ └───┘
    const auto matching = minimumWeightFullBipartiteMatching(costMatrix);
    ASSERT_EQ(matching.size(), 3);
    EXPECT_EQ(matching[0], 2);
    EXPECT_EQ(matching[1], 1);
    EXPECT_EQ(matching[2], 3);
  }
}

} // namespace na
