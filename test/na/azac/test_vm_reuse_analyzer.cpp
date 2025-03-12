#include "na/azac/VMReuseAnalyzer.hpp"

#include <cmath>
#include <cstddef>
#include <gtest/gtest.h>
#include <optional>
#include <utility>
#include <vector>

namespace na {
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
}
TEST_F(VMReuseAnalyzerMaximumBipartiteMatchingTest, Inverse) {
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
class VMReuseAnalyzerAnalyzeTest : public ::testing::Test {
protected:
  Architecture architecture;
  nlohmann::json config;
  VMReuseAnalyzer analyzer;

public:
  VMReuseAnalyzerAnalyzeTest() : analyzer{architecture, config} {}
};
TEST_F(VMReuseAnalyzerAnalyzeTest, NoGates) {
  std::vector<std::vector<std::pair<qc::Qubit, qc::Qubit>>> twoQubitGateLayers;
  EXPECT_NO_THROW(
      { EXPECT_TRUE(analyzer.analyzeReuse(twoQubitGateLayers).empty()); });
}
TEST_F(VMReuseAnalyzerAnalyzeTest, OneLayer) {
  std::vector<std::vector<std::pair<qc::Qubit, qc::Qubit>>> twoQubitGateLayers{
      {{0, 1}}};
  EXPECT_NO_THROW(
      { EXPECT_TRUE(analyzer.analyzeReuse(twoQubitGateLayers).empty()); });
}
TEST_F(VMReuseAnalyzerAnalyzeTest, NoChoice) {
  std::vector<std::vector<std::pair<qc::Qubit, qc::Qubit>>> twoQubitGateLayers{
      {{0, 1}}, {{1, 2}}};
  EXPECT_NO_THROW({
    const auto& reuseQubits = analyzer.analyzeReuse(twoQubitGateLayers);
    ASSERT_EQ(reuseQubits.size(), 1);
    const auto& reuseQubitsInLayer = reuseQubits.front();
    EXPECT_EQ(reuseQubitsInLayer, (std::unordered_set<qc::Qubit>{1}));
  });
}
TEST_F(VMReuseAnalyzerAnalyzeTest, Unique) {
  std::vector<std::vector<std::pair<qc::Qubit, qc::Qubit>>> twoQubitGateLayers{
      {{0, 1}, {2, 3}, {4, 5}}, {{1, 2}, {3, 4}, {5, 7}}};
  EXPECT_NO_THROW({
    const auto& reuseQubits = analyzer.analyzeReuse(twoQubitGateLayers);
    ASSERT_EQ(reuseQubits.size(), 1);
    const auto& reuseQubitsInLayer = reuseQubits.front();
    EXPECT_EQ(reuseQubitsInLayer, (std::unordered_set<qc::Qubit>{1, 3, 5}));
  });
}
TEST(VMReuseAnalyzerTest, Config) {
  std::stringstream buffer;
  std::streambuf* oldCout = std::cout.rdbuf(buffer.rdbuf());
  Architecture architecture;
  nlohmann::json config;
  std::istringstream iss(R"({
  "vm_reuse_analyzer": {
    "unknown_key": 42
  }
})");
  iss >> config;
  VMReuseAnalyzer analyzer{architecture, config};
  std::cout.rdbuf(oldCout);
  EXPECT_EQ(buffer.str(), "[WARN] Configuration for Placer contains an unknown "
                          "key: unknown_key. Ignoring.\n");
  // silence unused variable warning
  (void)analyzer;
}
} // namespace na
