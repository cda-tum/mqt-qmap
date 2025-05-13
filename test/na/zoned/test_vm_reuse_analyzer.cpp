#include "na/zoned/VMReuseAnalyzer.hpp"

#include <cmath>
#include <cstddef>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>
#include <gtest/gtest.h>
#include <optional>
#include <utility>
#include <vector>

namespace na::zoned {
constexpr std::string_view architectureJson = R"({
  "name": "asap_scheduler_architecture",
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
class VMReuseAnalyzerAnalyzeTest : public ::testing::Test {
protected:
  Architecture architecture;
  nlohmann::json config;
  VMReuseAnalyzer analyzer;
  VMReuseAnalyzerAnalyzeTest()
      : architecture(nlohmann::json::parse(architectureJson)),
        analyzer{architecture, config} {}
};
TEST_F(VMReuseAnalyzerAnalyzeTest, NoGates) {
  std::vector<std::vector<std::array<qc::Qubit, 2>>> twoQubitGateLayers;
  EXPECT_THAT(analyzer.analyzeReuse(twoQubitGateLayers), ::testing::IsEmpty());
}
TEST_F(VMReuseAnalyzerAnalyzeTest, OneLayer) {
  std::vector<std::vector<std::array<qc::Qubit, 2>>> twoQubitGateLayers{
      {{0, 1}}};
  EXPECT_THAT(analyzer.analyzeReuse(twoQubitGateLayers), ::testing::IsEmpty());
}
TEST_F(VMReuseAnalyzerAnalyzeTest, NoChoice) {
  std::vector<std::vector<std::array<qc::Qubit, 2>>> twoQubitGateLayers{
      {{0, 1}}, {{1, 2}}};
  EXPECT_THAT(analyzer.analyzeReuse(twoQubitGateLayers),
              ::testing::ElementsAre(::testing::UnorderedElementsAre(1U)));
}
TEST_F(VMReuseAnalyzerAnalyzeTest, Unique) {
  std::vector<std::vector<std::array<qc::Qubit, 2>>> twoQubitGateLayers{
      {{0, 1}, {2, 3}, {4, 5}}, {{1, 2}, {3, 4}, {5, 7}}};
  EXPECT_THAT(
      analyzer.analyzeReuse(twoQubitGateLayers),
      ::testing::ElementsAre(::testing::UnorderedElementsAre(1U, 3U, 5U)));
}
TEST_F(VMReuseAnalyzerAnalyzeTest, UniqueUnbalanced) {
  std::vector<std::vector<std::array<qc::Qubit, 2>>> twoQubitGateLayers{
      {{0, 1}, {2, 3}, {4, 5}, {6, 7}}, {{1, 6}, {7, 8}}};
  EXPECT_THAT(analyzer.analyzeReuse(twoQubitGateLayers),
              ::testing::ElementsAre(::testing::UnorderedElementsAre(1U, 7U)));
}
TEST(VMReuseAnalyzerTest, Config) {
  Architecture architecture(nlohmann::json::parse(architectureJson));
  nlohmann::json config = R"({
  "vm_reuse_analyzer": {
    "unknown_key": 42
  }
})"_json;
  std::stringstream buffer;
  std::streambuf* oldCout = std::cout.rdbuf(buffer.rdbuf());
  std::ignore = VMReuseAnalyzer(architecture, config);
  std::cout.rdbuf(oldCout);
  EXPECT_EQ(
      buffer.str(),
      "\033[1;35m[WARN]\033[0m Configuration for VMReuseAnalyzer contains an "
      "unknown key: unknown_key. Ignoring.\n");
}
class VMReuseAnalyzerMaximumBipartiteMatchingTest : public ::testing::Test {
protected:
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
};
TEST_F(VMReuseAnalyzerMaximumBipartiteMatchingTest, Direct) {
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
  EXPECT_THAT(VMReuseAnalyzer::maximumBipartiteMatching(sparseMatrix),
              ::testing::ElementsAre(0, 2, 1, 3));
}
TEST_F(VMReuseAnalyzerMaximumBipartiteMatchingTest, Inverse) {
  EXPECT_THAT(VMReuseAnalyzer::maximumBipartiteMatching(sparseMatrix, true),
              ::testing::ElementsAre(0, 2, 1, 3));
}
TEST(VMReuseAnalyzerMaximumBipartiteMatchingInvertedTest, Direct) {
  // We also test with the inverted graph, i.e., the sources and sinks are
  // labeled in reverse order, but sources stay sources and sinks stay sinks.
  const std::vector<std::vector<std::size_t>> inverseSparseMatrix{
      /* 0 -> */ {0, 1},
      /* 1 -> */ {0, 1, 2},
      /* 2 -> */ {1},
      /* 3 -> */ {2, 3}};
  EXPECT_THAT(VMReuseAnalyzer::maximumBipartiteMatching(inverseSparseMatrix),
              ::testing::ElementsAre(0, 2, 1, 3));
}
} // namespace na::zoned
