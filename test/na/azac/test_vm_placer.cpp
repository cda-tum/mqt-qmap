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
  "vm_placer" : {
    "use_window" : true,
    "window_size" : 10,
    "dynamic_placement" : true
  }
})";
class VMPlacerPlaceTest : public ::testing::Test {
protected:
  Architecture architecture;
  nlohmann::json config;
  VMPlacer placer;
  VMPlacerPlaceTest()
      : architecture(nlohmann::json::parse(architectureJson)),
        config(nlohmann::json::parse(configJson)),
        placer(architecture, config) {}
};
TEST_F(VMPlacerPlaceTest, Empty) {
  const size_t nQubits = 1;
  EXPECT_THAT(
      placer.place(nQubits,
                   std::vector<std::vector<std::pair<qc::Qubit, qc::Qubit>>>{},
                   std::vector<std::unordered_set<qc::Qubit>>{}),
      ::testing::ElementsAre(::testing::SizeIs(nQubits)));
}
TEST_F(VMPlacerPlaceTest, OneGate) {
  const size_t nQubits = 2;
  EXPECT_THAT(
      placer.place(
          nQubits,
          std::vector<std::vector<std::pair<qc::Qubit, qc::Qubit>>>{{{0U, 1U}}},
          std::vector<std::unordered_set<qc::Qubit>>{}),
      ::testing::ElementsAre(::testing::SizeIs(nQubits),
                             ::testing::SizeIs(nQubits),
                             ::testing::SizeIs(nQubits)));
}
TEST(VMPlacerTest, NoConfig) {
  Architecture architecture(nlohmann::json::parse(architectureJson));
  nlohmann::json config = R"({})"_json;
  std::stringstream buffer;
  std::streambuf* oldCout = std::cout.rdbuf(buffer.rdbuf());
  std::ignore = VMPlacer(architecture, config);
  EXPECT_EQ(buffer.str(),
            "[WARN] Configuration does not contain settings for VMPlacer or is "
            "malformed. Using default settings.\n");
  std::cout.rdbuf(oldCout);
}
TEST(VMPlacerTest, InvalidConfig) {
  Architecture architecture(nlohmann::json::parse(architectureJson));
  nlohmann::json config = R"({
  "vm_placer": {
    "use_window": "invalid",
    "window_size": 10,
    "unknown_key": 42
  }
})"_json;
  std::stringstream buffer;
  std::streambuf* oldCout = std::cout.rdbuf(buffer.rdbuf());
  std::ignore = VMPlacer(architecture, config);
  EXPECT_THAT(
      buffer.str(),
      ::testing::AllOf(
          ::testing::MatchesRegex(
              "\\[WARN\\].*\n\\[WARN\\].*\n\\[WARN\\].*\n\\[WARN\\].*\n"),
          ::testing::HasSubstr("[WARN] Configuration for VMPlacer "
                               "contains an invalid value for "
                               "use_window. Using default."),
          ::testing::HasSubstr("[WARN] Configuration for VMPlacer does "
                               "not contain a setting for "
                               "use_window. Using default."),
          ::testing::HasSubstr("[WARN] Configuration for VMPlacer does "
                               "not contain a setting for "
                               "dynamic_placement. Using default."),
          ::testing::HasSubstr("[WARN] Configuration for VMPlacer contains an "
                               "unknown key: unknown_key. Ignoring.")));
  std::cout.rdbuf(oldCout);
}
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
