/*
 * Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
 * Copyright (c) 2025 Munich Quantum Software Company GmbH
 * All rights reserved.
 *
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License
 */

#include "na/zoned/AStarPlacer.hpp"

#include <cstddef>
#include <gmock/gmock-function-mocker.h>
#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>
#include <src/gtest-internal-inl.h>
#include <string>
#include <utility>
#include <vector>

namespace na::zoned {
constexpr std::string_view architectureJson = R"({
  "name": "a_star_placer_architecture",
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
  "rydberg_range": [[[5, 70], [55, 110]]]
})";
constexpr std::string_view configJson = R"({
  "a_star_placer": {
    "use_window": true,
    "window_min_width": 4,
    "window_ratio": 1.5,
    "window_share": 0.6,
    "deepening_factor": 0.6,
    "deepening_value": 0.2,
    "lookahead_factor": 0.2,
    "reuse_level": 5.0
  }
})";
class AStarPlacerPlaceTest : public ::testing::Test {
protected:
  Architecture architecture;
  nlohmann::json config;
  AStarPlacer placer;
  AStarPlacerPlaceTest()
      : architecture(nlohmann::json::parse(architectureJson)),
        config(nlohmann::json::parse(configJson)),
        placer(architecture, config) {}
};
TEST_F(AStarPlacerPlaceTest, Empty) {
  constexpr size_t nQubits = 1;
  EXPECT_THAT(placer.place(nQubits,
                           std::vector<std::vector<std::array<qc::Qubit, 2>>>{},
                           std::vector<std::unordered_set<qc::Qubit>>{}),
              ::testing::ElementsAre(::testing::SizeIs(nQubits)));
}
TEST_F(AStarPlacerPlaceTest, OneGate) {
  constexpr size_t nQubits = 2;
  EXPECT_THAT(placer.place(nQubits,
                           std::vector<std::vector<std::array<qc::Qubit, 2>>>{
                               {{0U, 1U}}},
                           std::vector<std::unordered_set<qc::Qubit>>{}),
              ::testing::ElementsAre(::testing::SizeIs(nQubits),
                                     ::testing::SizeIs(nQubits),
                                     ::testing::SizeIs(nQubits)));
}
TEST_F(AStarPlacerPlaceTest, TwoGatesCons) {
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
TEST_F(AStarPlacerPlaceTest, OneGateCross) {
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
    const auto x = architecture.exactSlmLocation(slm, r, c).first;
    qubitsInEntanglementByX.emplace(x, q);
  }
  std::vector<qc::Qubit> qubitsInEntanglementAsc;
  for (const auto& [_, q] : qubitsInEntanglementByX) {
    qubitsInEntanglementAsc.push_back(q);
  }
  EXPECT_THAT(qubitsInEntanglementAsc, ::testing::ElementsAre(0U, 1U));
}
TEST_F(AStarPlacerPlaceTest, TwoGatesZip) {
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
  EXPECT_THAT(qubitsInEntanglementAsc,
              ::testing::AnyOf(::testing::ElementsAre(0U, 2U, 1U, 3U),
                               ::testing::ElementsAre(1U, 3U, 0U, 2U)));
  EXPECT_THAT(qubitsInEntanglementYs, ::testing::UnorderedElementsAre(70UL));
}
TEST_F(AStarPlacerPlaceTest, FullEntanglementZone) {
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
TEST_F(AStarPlacerPlaceTest, TwoTwoQubitLayerReuse) {
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
TEST(AStarPlacerTest, NoSolution) {
  Architecture architecture(nlohmann::json::parse(architectureJson));
  AStarPlacer placer(architecture, R"({
  "a_star_placer": {
    "use_window": true,
    "window_min_width": 0,
    "window_ratio": 1.5,
    "window_share": 0.0,
    "deepening_factor": 0.6,
    "deepening_value": 0.2,
    "lookahead_factor": 0.2,
    "reuse_level": 5.0
  }
})"_json);
  constexpr size_t nQubits = 2;
  EXPECT_THROW(
      std::ignore = placer.place(
          nQubits,
          std::vector<std::vector<std::array<qc::Qubit, 2>>>{{{0U, 1U}}},
          std::vector<std::unordered_set<qc::Qubit>>{}),
      std::runtime_error);
}
TEST(AStarPlacerTest, LimitSpace) {
  Architecture architecture(nlohmann::json::parse(architectureJson));
  AStarPlacer placer(architecture, R"({
  "a_star_placer": {
    "use_window": true,
    "window_min_width": 4,
    "window_ratio": 1.5,
    "window_share": 0.6,
    "deepening_factor": 0.6,
    "deepening_value": 0.2,
    "lookahead_factor": 0.2,
    "reuse_level": 5.0,
    "max_nodes": 2
  }
})"_json);
  constexpr size_t nQubits = 4;
  EXPECT_THROW(std::ignore = placer.place(
                   nQubits,
                   std::vector<std::vector<std::array<qc::Qubit, 2>>>{
                       {{0U, 1U}, {2U, 3U}}},
                   std::vector<std::unordered_set<qc::Qubit>>{}),
               std::runtime_error);
}
TEST(AStarPlacerTest, WindowExpansion) {
  Architecture architecture(nlohmann::json::parse(architectureJson));
  AStarPlacer placer(architecture, R"({
  "a_star_placer": {
    "use_window": true,
    "window_min_width": 1,
    "window_ratio": 1.0,
    "window_share": 1.0,
    "deepening_factor": 0.6,
    "deepening_value": 0.2,
    "lookahead_factor": 0.2,
    "reuse_level": 5.0
  }
})"_json);
  constexpr size_t nQubits = 4;
  EXPECT_NO_THROW(std::ignore = placer.place(
                      nQubits,
                      std::vector<std::vector<std::array<qc::Qubit, 2>>>{
                          {{0U, 3U}, {1U, 2U}}},
                      std::vector<std::unordered_set<qc::Qubit>>{}));
}
TEST(AStarPlacerTest, InitialPlacementForTwoSlms) {
  Architecture architecture(R"({
  "name": "a_star_placer_architecture",
  "storage_zones": [{
    "zone_id": 0,
    "slms": [
      {"id": 0, "site_separation": [3, 3], "r": 2, "c": 20, "location": [0, 0]},
      {"id": 1, "site_separation": [3, 3], "r": 18, "c": 20, "location": [0, 6]}],
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
  "rydberg_range": [[[5, 70], [55, 110]]]
})"_json);
  AStarPlacer placer(architecture, nlohmann::json::parse(configJson));
  constexpr size_t nQubits = 50;
  const auto& placement = placer.place(
      nQubits, std::vector<std::vector<std::array<qc::Qubit, 2>>>{},
      std::vector<std::unordered_set<qc::Qubit>>{});
  EXPECT_THAT(placement, ::testing::ElementsAre(::testing::SizeIs(nQubits)));
  // check that there exists a qubit that is placed in SLM with ID 1
  EXPECT_THAT(placement,
              ::testing::ElementsAre(::testing::Contains(::testing::FieldsAre(
                  ::testing::Field(&SLM::id, ::testing::Eq(1)),
                  ::testing::Lt(18), ::testing::Lt(20)))));
}
TEST(AStarPlacerTest, NoConfig) {
  Architecture architecture(nlohmann::json::parse(architectureJson));
  nlohmann::json config = R"({})"_json;
  std::stringstream buffer;
  std::streambuf* oldCout = std::cout.rdbuf(buffer.rdbuf());
  std::ignore = AStarPlacer(architecture, config);
  EXPECT_EQ(
      buffer.str(),
      "\033[1;35m[WARN]\033[0m Configuration does not contain settings for "
      "AStarPlacer or is malformed. Using default settings ("
      "\"use_window\": true, "
      "\"window_min_width\": 8, "
      "\"window_ratio\": 1, "
      "\"window_share\": 0.6, "
      "\"deepening_factor\": 0.8, "
      "\"deepening_value\": 0.2, "
      "\"lookahead_factor\": 0.2, "
      "\"reuse_level\": 5, "
      "\"max_nodes\": 50000000).\n");
  std::cout.rdbuf(oldCout);
}
TEST(AStarPlacerTest, InvalidConfig) {
  Architecture architecture(nlohmann::json::parse(architectureJson));
  nlohmann::json config = R"({
  "a_star_placer": {
    "use_window": "invalid",
    "window_min_width": "invalid",
    "window_ratio": "invalid",
    "window_share": "invalid",
    "deepening_factor": "invalid",
    "deepening_value": "invalid",
    "lookahead_factor": "invalid",
    "reuse_level": "invalid",
    "max_nodes": "invalid",
    "unknown_key": 42
  }
})"_json;
  std::stringstream buffer;
  std::streambuf* oldCout = std::cout.rdbuf(buffer.rdbuf());
  std::ignore = AStarPlacer(architecture, config);
  std::cout.rdbuf(oldCout);
  EXPECT_THAT(
      buffer.str(),
      ::testing::AllOf(
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer contains "
              "an invalid value for use_window. Using default (true)."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer contains "
              "an invalid value for window_min_width. Using default (8)."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer contains "
              "an invalid value for window_ratio. Using default (1)."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer contains "
              "an invalid value for window_share. Using default (0.6)."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer contains "
              "an invalid value for deepening_factor. Using default (0.8)."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer contains "
              "an invalid value for deepening_value. Using default (0.2)."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer contains "
              "an invalid value for lookahead_factor. Using default (0.2)."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer contains "
              "an invalid value for reuse_level. Using default (5)."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer contains "
              "an invalid value for max_nodes. Using default (50000000)."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer contains "
              "an unknown key: unknown_key. Ignoring.")));
  size_t warnings = 0;
  size_t pos = 0;
  std::string target = "\033[1;35m[WARN]\033[0m";
  while ((pos = buffer.str().find(target, pos)) != std::string::npos) {
    ++warnings;
    pos += target.length();
  }
  EXPECT_EQ(warnings, 10);
}
TEST(AStarPlacerTest, EmptyConfig) {
  Architecture architecture(nlohmann::json::parse(architectureJson));
  nlohmann::json config = R"({
  "a_star_placer": {}
})"_json;
  std::stringstream buffer;
  std::streambuf* oldCout = std::cout.rdbuf(buffer.rdbuf());
  std::ignore = AStarPlacer(architecture, config);
  std::cout.rdbuf(oldCout);
  EXPECT_THAT(
      buffer.str(),
      ::testing::AllOf(
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer does not "
              "contain a value for use_window. Using default (true)."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer does not "
              "contain a value for window_min_width. Using default (8)."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer does not "
              "contain a value for window_ratio. Using default (1)."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer does not "
              "contain a value for window_share. Using default (0.6)."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer does not "
              "contain a value for deepening_factor. Using default (0.8)."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer does not "
              "contain a value for deepening_value. Using default (0.2)."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer does not "
              "contain a value for lookahead_factor. Using default (0.2)."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer does not "
              "contain a value for reuse_level. Using default (5)."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for AStarPlacer does not "
              "contain a value for max_nodes. Using default (50000000).")));
  size_t warnings = 0;
  size_t pos = 0;
  std::string target = "\033[1;35m[WARN]\033[0m";
  while ((pos = buffer.str().find(target, pos)) != std::string::npos) {
    ++warnings;
    pos += target.length();
  }
  EXPECT_EQ(warnings, 9);
}
TEST(AStarPlacerTest, AStarSearch) {
  // for testing purposes, we do not use the structure of nodes and just use
  // their respective address to identify a location in a 4x4 grid that looks
  // like the following, where the cost of each edge is 1:
  // ┌Start┐        ┌─────┐        ┌─────┐        ┌─────┐
  // │  0  ├─────→  │  1  ├─────→  │  2  ├─────→  │  3  │
  // └─────┘        └──┬──┘        └──┬──┘        └──┬──┘
  //    │              │              │              │
  //    ↓              ↓              ↓              ↓
  // ┌─────┐        ┌─────┐        ┌─────┐        ┌─────┐
  // │  4  ├─────→  │  5  ├─────→  │  6  ├─────→  │  7  │
  // └──┬──┘        └──┬──┘        └──┬──┘        └──┬──┘
  //    │              │              │              │
  //    ↓              ↓              ↓              ↓
  // ┌─────┐        ┌─────┐        ┌─────┐        ┌─────┐
  // │  8  ├─────→  │  9  ├─────→  │  10 ├─────→  │  11 │
  // └──┬──┘        └──┬──┘        └──┬──┘        └──┬──┘
  //    │              │              │              │
  //    ↓              ↓              ↓              ↓
  // ┌─────┐        ┌─────┐        ┌Goal=┐        ┌─────┐
  // │  12 ├─────→  │  13 ├─────→  │  14 ├─────→  │  15 │
  // └─────┘        └─────┘        └=====┘        └─────┘
  const std::vector<AStarPlacer::AtomNode> nodes(16);
  std::unordered_map<
      const AStarPlacer::AtomNode*,
      std::vector<std::reference_wrapper<const AStarPlacer::AtomNode>>>
      neighbors{{nodes.data(), {std::cref(nodes[1]), std::cref(nodes[4])}},
                {&nodes[1], {std::cref(nodes[2]), std::cref(nodes[5])}},
                {&nodes[2], {std::cref(nodes[3]), std::cref(nodes[6])}},
                {&nodes[3], {std::cref(nodes[7])}},
                {&nodes[4], {std::cref(nodes[5]), std::cref(nodes[8])}},
                {&nodes[5], {std::cref(nodes[6]), std::cref(nodes[9])}},
                {&nodes[6], {std::cref(nodes[7]), std::cref(nodes[10])}},
                {&nodes[7], {std::cref(nodes[11])}},
                {&nodes[8], {std::cref(nodes[9]), std::cref(nodes[12])}},
                {&nodes[9], {std::cref(nodes[10]), std::cref(nodes[13])}},
                {&nodes[10], {std::cref(nodes[11]), std::cref(nodes[14])}},
                {&nodes[11], {std::cref(nodes[15])}},
                {&nodes[12], {std::cref(nodes[13])}},
                {&nodes[13], {std::cref(nodes[14])}},
                {&nodes[14], {std::cref(nodes[15])}},
                {&nodes[15], {}}};
  const auto path = (AStarPlacer::aStarTreeSearch<AStarPlacer::AtomNode>(
      /* start: */
      nodes[0],
      /* getNeighbors: */
      [&neighbors](const AStarPlacer::AtomNode& node)
          -> std::vector<std::reference_wrapper<const AStarPlacer::AtomNode>> {
        return neighbors.at(&node);
      },
      /* isGoal: */
      [&nodes](const AStarPlacer::AtomNode& node) -> bool {
        return &node == &nodes[14];
      },
      /* getCost: */
      [](const AStarPlacer::AtomNode& /* unused */) -> double { return 1.0; },
      /* getHeuristic: */
      [&nodes](const AStarPlacer::AtomNode& node) -> double {
        const auto* head = nodes.data();
        const auto i = std::distance(head, &node);
        const long x = i % 4;
        const long y = i / 4;
        return std::hypot(x, y);
      },
      1'000'000));
  // convert to const Node* for easier comparison
  std::vector<const AStarPlacer::AtomNode*> pathNodes;
  for (const auto& node : path) {
    pathNodes.emplace_back(&node.get());
  }
  EXPECT_THAT(pathNodes,
              ::testing::ElementsAre(&nodes[0], ::testing::_, ::testing::_,
                                     ::testing::_, ::testing::_, &nodes[14]));
}
} // namespace na::zoned
