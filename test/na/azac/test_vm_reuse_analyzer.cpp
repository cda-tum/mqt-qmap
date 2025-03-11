#include "na/azac/VMReuseAnalyzer.hpp"

#include <cmath>
#include <cstddef>
#include <gtest/gtest.h>
#include <optional>
#include <utility>
#include <vector>

namespace na {
TEST(VMReuseAnalyzerTest, MaximumBipartiteMatching) {
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
  const auto matching = VMReuseAnalyzer::maximumBipartiteMatching(sparseMatrix);
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
  const auto invMatching =
      VMReuseAnalyzer::maximumBipartiteMatching(sparseMatrix, true);
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
  const auto matchingOfInverse =
      VMReuseAnalyzer::maximumBipartiteMatching(inverseSparseMatrix);
  ASSERT_EQ(matchingOfInverse.size(), 4);
  EXPECT_EQ(matchingOfInverse[0], 0);
  EXPECT_EQ(matchingOfInverse[1], 2);
  EXPECT_EQ(matchingOfInverse[2], 1);
  EXPECT_EQ(matchingOfInverse[3], 3);
}
} // namespace na
