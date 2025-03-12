#include "na/azac/VMPlacer.hpp"

#include <cstddef>
#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>
#include <optional>
#include <utility>
#include <vector>

namespace na {
constexpr std::string_view architectureJson = R"({
  "name": "vm_placer_architecture",
  "storage_zones": [{
    "zone_id": 0,
    "slms": [{"id": 0, "site_separation": [3, 3], "r": 20, "c": 20, "location": [0, 0]}],
    "offset": [0, 0],
    "dimension": [60, 60]
  }],
  "entanglement_zones": [{
    "zone_id": 0,
    "slms": [
      {"id": 1, "site_separation": [12, 10], "r": 4, "c": 4, "location": [5, 70]},
      {"id": 2, "site_separation": [12, 10], "r": 4, "c": 4, "location": [7, 70]}
    ],
    "offset": [5, 70],
    "dimension": [50, 40]
  }],
  "aods":[{"id": 0, "site_separation": 2, "r": 20, "c": 20}],
  "arch_range": [[0, 0], [60, 110]],
  "rydberg_range": [[[5, 70], [55, 110]]]
})";
constexpr std::string_view configJson = R"({
})";
class VMPlacerTest : public ::testing::Test {
protected:
  Architecture architecture;
  nlohmann::json config;
  VMPlacer placer;
  VMPlacerTest()
      : architecture(nlohmann::json::parse(architectureJson)),
        config(configJson), placer(architecture, config) {}
};
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
  EXPECT_THAT(matching, ::testing::ElementsAre(0, 1, 3));
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
  EXPECT_THAT(matching, ::testing::ElementsAre(2, 1, 3));
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
