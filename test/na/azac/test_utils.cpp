#include "na/azac/Utils.hpp"

#include <cmath>
#include <gtest/gtest.h>
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
  // we consider the following bipartite graph, where the nodes in the upper row
  // are the sources, and the nodes in the lower row are the sinks.
  //
  //   ┌───┐ ┌───┐ ┌───┐ ┌───┐
  //   │ 0 │ │ 1 │ │ 2 │ │ 3 │ <-- SOURCES
  //   └─┬─┘ └─┬─┘ └─┬─┘ └─┬─┘
  //     │╲     ╲   ╱│╲   ╱│
  //     │  ╲     ╳  │  ╳  │
  //     │    ╲ ╱   ╲│╱   ╲│
  //   ┌─┴─┐ ┌─┴─┐ ┌─┴─┐ ┌─┴─┐
  //   │ 0 │ │ 1 │ │ 2 │ │ 3 │ <-- SINKS
  //   └───┘ └───┘ └───┘ └───┘
  const std::vector<std::vector<std::size_t>> sparseMatrix{
      /* 0 -> */ {0, 1},
      /* 1 -> */ {2},
      /* 2 -> */ {1, 2, 3},
      /* 3 -> */ {2, 3}
  };
  const auto matching = maximumBipartiteMatching(sparseMatrix);
  // The result should be the following (unique) maximum matching:
  //
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
  // we also test with the inverted graph, i.e., the sources and sinks are
  // labelled in reverse order.
  const std::vector<std::vector<std::size_t>> inverseSparseMatrix{
    /* 0 -> */ {0, 1},
    /* 1 -> */ {0, 1, 2},
    /* 2 -> */ {1},
    /* 3 -> */ {2, 3}
  };
  const auto matchingOfInverse = maximumBipartiteMatching(inverseSparseMatrix);
  ASSERT_EQ(matchingOfInverse.size(), 4);
  EXPECT_EQ(matchingOfInverse[0], 0);
  EXPECT_EQ(matchingOfInverse[1], 2);
  EXPECT_EQ(matchingOfInverse[2], 1);
  EXPECT_EQ(matchingOfInverse[3], 3);
}

TEST(TestUtils, MinimumWeightFullBipartiteMatching) {}

} // namespace na
