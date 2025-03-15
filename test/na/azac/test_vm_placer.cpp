#include "na/azac/VMPlacer.hpp"

#include <cstddef>
#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>
#include <map>
#include <optional>
#include <unordered_set>
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
  constexpr size_t nQubits = 1;
  EXPECT_THAT(placer.place(nQubits,
                           std::vector<std::vector<std::array<qc::Qubit, 2>>>{},
                           std::vector<std::unordered_set<qc::Qubit>>{}),
              ::testing::ElementsAre(::testing::SizeIs(nQubits)));
}
TEST_F(VMPlacerPlaceTest, OneGate) {
  constexpr size_t nQubits = 2;
  EXPECT_THAT(placer.place(nQubits,
                           std::vector<std::vector<std::array<qc::Qubit, 2>>>{
                               {{0U, 1U}}},
                           std::vector<std::unordered_set<qc::Qubit>>{}),
              ::testing::ElementsAre(::testing::SizeIs(nQubits),
                                     ::testing::SizeIs(nQubits),
                                     ::testing::SizeIs(nQubits)));
}
TEST_F(VMPlacerPlaceTest, TwoGatesCons) {
  constexpr size_t nQubits = 4;
  const auto& placement = placer.place(
      nQubits,
      std::vector<std::vector<std::array<qc::Qubit, 2>>>{{{0U, 1U}, {2U, 3U}}},
      std::vector<std::unordered_set<qc::Qubit>>{});
  EXPECT_THAT(placement, ::testing::SizeIs(3));
  EXPECT_THAT(placement, ::testing::Each(::testing::SizeIs(nQubits)));
  std::map<size_t, qc::Qubit> qubitsInStorageByX;
  std::unordered_set<size_t> qubitsInStorageYs;
  for (qc::Qubit q = 0; q < placement.front().size(); ++q) {
    const auto& [slm, r, c] = placement.front()[q];
    EXPECT_TRUE(slm.get().isStorage());
    const auto& [x, y] = architecture.exactSlmLocation(slm, r, c);
    qubitsInStorageByX.emplace(x, q);
    qubitsInStorageYs.emplace(y);
  }
  std::vector<qc::Qubit> qubitsInStorageAsc;
  for (const auto& [_, q] : qubitsInStorageByX) {
    qubitsInStorageAsc.push_back(q);
  }
  EXPECT_THAT(qubitsInStorageAsc, ::testing::ElementsAre(0U, 1U, 2U, 3U));
  EXPECT_THAT(qubitsInStorageYs, ::testing::UnorderedElementsAre(19UL * 3));
  std::map<size_t, qc::Qubit> qubitsInEntanglementByX;
  std::unordered_set<size_t> qubitsInEntanglementYs;
  for (qc::Qubit q = 0; q < placement[1].size(); ++q) {
    const auto& [slm, r, c] = placement[1][q];
    EXPECT_TRUE(slm.get().isEntanglement());
    const auto& [x, y] = architecture.exactSlmLocation(slm, r, c);
    qubitsInEntanglementByX.emplace(x, q);
    qubitsInEntanglementYs.emplace(y);
  }
  std::vector<qc::Qubit> qubitsInEntanglementAsc;
  for (const auto& [_, q] : qubitsInEntanglementByX) {
    qubitsInEntanglementAsc.push_back(q);
  }
  EXPECT_THAT(qubitsInEntanglementAsc, ::testing::ElementsAre(0U, 1U, 2U, 3U));
  EXPECT_THAT(qubitsInEntanglementYs, ::testing::UnorderedElementsAre(70UL));
}
TEST_F(VMPlacerPlaceTest, OneGateCross) {
  constexpr size_t nQubits = 2;
  const auto& placement = placer.place(
      nQubits, std::vector<std::vector<std::array<qc::Qubit, 2>>>{{{1U, 0U}}},
      std::vector<std::unordered_set<qc::Qubit>>{});
  EXPECT_THAT(placement, ::testing::SizeIs(3));
  EXPECT_THAT(placement, ::testing::Each(::testing::SizeIs(nQubits)));
  std::map<size_t, qc::Qubit> qubitsInEntanglementByX;
  for (qc::Qubit q = 0; q < placement[1].size(); ++q) {
    const auto& [slm, r, c] = placement[1][q];
    EXPECT_TRUE(slm.get().isEntanglement());
    const auto& [x, y] = architecture.exactSlmLocation(slm, r, c);
    qubitsInEntanglementByX.emplace(x, q);
  }
  std::vector<qc::Qubit> qubitsInEntanglementAsc;
  for (const auto& [_, q] : qubitsInEntanglementByX) {
    qubitsInEntanglementAsc.push_back(q);
  }
  EXPECT_THAT(qubitsInEntanglementAsc, ::testing::ElementsAre(0U, 1U));
}
TEST_F(VMPlacerPlaceTest, TwoGatesZip) {
  constexpr size_t nQubits = 4;
  const auto& placement = placer.place(
      nQubits,
      std::vector<std::vector<std::array<qc::Qubit, 2>>>{{{0U, 2U}, {1U, 3U}}},
      std::vector<std::unordered_set<qc::Qubit>>{});
  EXPECT_THAT(placement, ::testing::SizeIs(3));
  EXPECT_THAT(placement, ::testing::Each(::testing::SizeIs(nQubits)));
  std::map<size_t, qc::Qubit> qubitsInEntanglementByX;
  std::unordered_set<size_t> qubitsInEntanglementYs;
  for (qc::Qubit q = 0; q < placement[1].size(); ++q) {
    const auto& [slm, r, c] = placement[1][q];
    EXPECT_TRUE(slm.get().isEntanglement());
    const auto& [x, y] = architecture.exactSlmLocation(slm, r, c);
    qubitsInEntanglementByX.emplace(x, q);
    qubitsInEntanglementYs.emplace(y);
  }
  std::vector<qc::Qubit> qubitsInEntanglementAsc;
  for (const auto& [_, q] : qubitsInEntanglementByX) {
    qubitsInEntanglementAsc.push_back(q);
  }
  EXPECT_THAT(qubitsInEntanglementAsc, ::testing::ElementsAre(0U, 2U, 1U, 3U));
  EXPECT_THAT(qubitsInEntanglementYs, ::testing::UnorderedElementsAre(70UL));
}
TEST_F(VMPlacerPlaceTest, FullEntanglementZone) {
  constexpr size_t nQubits = 32;
  const auto& placement = placer.place(
      nQubits,
      std::vector<std::vector<std::array<qc::Qubit, 2>>>{{{0U, 1U},
                                                          {2U, 3U},
                                                          {4U, 5U},
                                                          {6U, 7U},
                                                          {8U, 9U},
                                                          {10U, 11U},
                                                          {12U, 13U},
                                                          {14U, 15U},
                                                          {16U, 17U},
                                                          {18U, 19U},
                                                          {20U, 21U},
                                                          {22U, 23U},
                                                          {24U, 25U},
                                                          {26U, 27U},
                                                          {28U, 29U},
                                                          {30U, 31U}}},
      std::vector<std::unordered_set<qc::Qubit>>{});
  EXPECT_THAT(placement, ::testing::SizeIs(3));
  EXPECT_THAT(placement, ::testing::Each(::testing::SizeIs(nQubits)));
  std::unordered_set<std::pair<size_t, size_t>> qubitsLocationsInEntanglement;
  for (qc::Qubit q = 0; q < placement[1].size(); ++q) {
    const auto& [slm, r, c] = placement[1][q];
    EXPECT_TRUE(slm.get().isEntanglement());
    const auto& [x, y] = architecture.exactSlmLocation(slm, r, c);
    qubitsLocationsInEntanglement.emplace(x, y);
  }
  EXPECT_THAT(qubitsLocationsInEntanglement, ::testing::SizeIs(nQubits));
}
TEST_F(VMPlacerPlaceTest, TwoTwoQubitLayerReuse) {
  constexpr size_t nQubits = 3;
  const auto& placement =
      placer.place(nQubits,
                   std::vector<std::vector<std::array<qc::Qubit, 2>>>{
                       {{0U, 1U}}, {{1U, 2U}}},
                   std::vector<std::unordered_set<qc::Qubit>>{{1U}});
  EXPECT_THAT(placement, ::testing::SizeIs(5));
  EXPECT_THAT(placement, ::testing::Each(::testing::SizeIs(nQubits)));
  // Check that qubit 1 remains in the entanglement zone while qubits 0 and 2
  // are placed in the storage zone in the intermediate layer
  EXPECT_TRUE(std::get<0>(placement[2][0]).get().isStorage());
  EXPECT_TRUE(std::get<0>(placement[2][1]).get().isEntanglement());
  EXPECT_TRUE(std::get<0>(placement[2][2]).get().isStorage());
  // Check that qubit 1 remains at the same position from layer 1 through 3
  EXPECT_EQ(std::get<0>(placement[1][1]).get(),
            std::get<0>(placement[2][1]).get());
  EXPECT_EQ(std::get<1>(placement[1][1]), std::get<1>(placement[2][1]));
  EXPECT_EQ(std::get<2>(placement[1][1]), std::get<2>(placement[2][1]));
  EXPECT_EQ(std::get<0>(placement[2][1]).get(),
            std::get<0>(placement[3][1]).get());
  EXPECT_EQ(std::get<1>(placement[2][1]), std::get<1>(placement[3][1]));
  EXPECT_EQ(std::get<2>(placement[2][1]), std::get<2>(placement[3][1]));
}
TEST(VMPlacerTest, NoConfig) {
  Architecture architecture(nlohmann::json::parse(architectureJson));
  nlohmann::json config = R"({})"_json;
  std::stringstream buffer;
  std::streambuf* oldCout = std::cout.rdbuf(buffer.rdbuf());
  std::ignore = VMPlacer(architecture, config);
  EXPECT_EQ(buffer.str(), "\033[1;35m[WARN]\033[0m Configuration does not "
                          "contain settings for VMPlacer or is "
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
  std::cout.rdbuf(oldCout);
  EXPECT_THAT(
      buffer.str(),
      ::testing::AllOf(
          ::testing::MatchesRegex(".*\\[WARN\\].*\n.*\\[WARN\\].*\n.*\\[WARN\\]"
                                  ".*\n.*\\[WARN\\].*\n"),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for VMPlacer "
              "contains an invalid value for "
              "use_window. Using default."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for VMPlacer does "
              "not contain a setting for "
              "use_window. Using default."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for VMPlacer does "
              "not contain a setting for "
              "dynamic_placement. Using default."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for VMPlacer contains an "
              "unknown key: unknown_key. Ignoring.")));
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
