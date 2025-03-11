#include "na/azac/VMPlacer.hpp"

#include <cstddef>
#include <gtest/gtest.h>
#include <optional>
#include <utility>
#include <vector>

namespace na {
TEST(VMPlacerTest, MinimumWeightFullBipartiteMatching1) {
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
  const auto matching =
      VMPlacer::minimumWeightFullBipartiteMatching(costMatrix);
  ASSERT_EQ(matching.size(), 3);
  EXPECT_EQ(matching[0], 0);
  EXPECT_EQ(matching[1], 1);
  EXPECT_EQ(matching[2], 3);
}
TEST(VMPlacerTest, MinimumWeightFullBipartiteMatching2) {
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
  const auto matching =
      VMPlacer::minimumWeightFullBipartiteMatching(costMatrix);
  ASSERT_EQ(matching.size(), 3);
  EXPECT_EQ(matching[0], 2);
  EXPECT_EQ(matching[1], 1);
  EXPECT_EQ(matching[2], 3);
}
TEST(VMPlacerTest, MinimumWeightFullBipartiteMatchingExceptions) {
  EXPECT_THROW(std::ignore =
                   VMPlacer::minimumWeightFullBipartiteMatching({{0}, {0}}),
               std::invalid_argument);
  EXPECT_THROW(std::ignore = VMPlacer::minimumWeightFullBipartiteMatching(
                   {{std::nullopt}}),
               std::invalid_argument);
  EXPECT_THROW(std::ignore = VMPlacer::minimumWeightFullBipartiteMatching(
                   {{0, 0}, {std::nullopt, std::nullopt}}),
               std::invalid_argument);
}
} // namespace na
