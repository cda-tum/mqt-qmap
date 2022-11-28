//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#include "gtest/gtest.h"

namespace cs {

struct SynthesisTest {
  std::string description;
  std::string initialTableau;
  std::string targetTableau;
  std::size_t expectedMinimalGates{};
  std::size_t expectedMinimalDepth{};
  std::size_t expectedMinimalTwoQubitGates{};
};

inline void from_json(const nlohmann::json& j, SynthesisTest& test) {
  test.description          = j.at("description").get<std::string>();
  test.initialTableau       = j.at("initial_tableau").get<std::string>();
  test.targetTableau        = j.at("target_tableau").get<std::string>();
  test.expectedMinimalGates = j.at("expected_minimal_gates").get<std::size_t>();
  test.expectedMinimalDepth = j.at("expected_minimal_depth").get<std::size_t>();
  test.expectedMinimalTwoQubitGates =
      j.at("expected_minimal_two_qubit_gates").get<std::size_t>();
}

static std::vector<SynthesisTest> getTests(const std::string& path) {
  std::ifstream  input(path);
  nlohmann::json j;
  input >> j;
  return j;
}

class TableauSynthesisTest : public ::testing::TestWithParam<SynthesisTest> {
protected:
  void SetUp() override {
    if (plog::get() == nullptr) {
      util::init();
    }

    const auto& test = GetParam();

    initialTableau = Tableau(test.initialTableau);
    std::cout << "Initial tableau:\n" << initialTableau;

    targetTableau = Tableau(test.targetTableau);
    std::cout << "Target tableau:\n" << targetTableau;

    synthesizer = CliffordSynthesizer(initialTableau, targetTableau);

    config = Configuration();

    expectedMinimalGates         = test.expectedMinimalGates;
    expectedMinimalDepth         = test.expectedMinimalDepth;
    expectedMinimalTwoQubitGates = test.expectedMinimalTwoQubitGates;
  }

  void TearDown() override {
    std::cout << "Results:\n" << results << "\n";

    resultTableau = synthesizer.getResultTableau();
    std::cout << "Resulting tableau:\n" << resultTableau;
    EXPECT_EQ(resultTableau, targetTableau);

    const auto& resultCircuit = synthesizer.getResultCircuit();
    std::cout << "Resulting Circuit:\n" << resultCircuit;
    consistencyCheck(resultCircuit);
  }

  void consistencyCheck(const qc::QuantumComputation& qc) const {
    auto circuitTableau = initialTableau;
    for (const auto& gate : qc) {
      circuitTableau.applyGate(gate.get());
    }
    EXPECT_EQ(resultTableau, circuitTableau);
  }

  Tableau             initialTableau;
  Tableau             targetTableau;
  Configuration       config;
  CliffordSynthesizer synthesizer;
  Results             results;
  Tableau             resultTableau;
  std::size_t         expectedMinimalGates{};
  std::size_t         expectedMinimalDepth{};
  std::size_t         expectedMinimalTwoQubitGates{};
};

INSTANTIATE_TEST_SUITE_P(
    Json, TableauSynthesisTest,
    testing::ValuesIn(getTests("cliffordsynthesis/tableaus.json")),
    [](const testing::TestParamInfo<TableauSynthesisTest::ParamType>& info) {
      return info.param.description;
    });

TEST_P(TableauSynthesisTest, Gates) {
  config.target = TargetMetric::GATES;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getGates(), expectedMinimalGates);
}

TEST_P(TableauSynthesisTest, Depth) {
  config.target = TargetMetric::DEPTH;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getDepth(), expectedMinimalDepth);
}

TEST_P(TableauSynthesisTest, TwoQubitGates) {
  config.target = TargetMetric::TWO_QUBIT_GATES;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getTwoQubitGates(), expectedMinimalTwoQubitGates);
}

TEST_P(TableauSynthesisTest, GatesMaxSAT) {
  config.target    = TargetMetric::GATES;
  config.useMaxSAT = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getGates(), expectedMinimalGates);
}

TEST_P(TableauSynthesisTest, DepthMaxSAT) {
  config.target    = TargetMetric::DEPTH;
  config.useMaxSAT = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getDepth(), expectedMinimalDepth);
}

TEST_P(TableauSynthesisTest, TwoQubitGatesMaxSAT) {
  config.target    = TargetMetric::TWO_QUBIT_GATES;
  config.useMaxSAT = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getTwoQubitGates(), expectedMinimalTwoQubitGates);
}

} // namespace cs
