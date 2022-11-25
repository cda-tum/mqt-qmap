//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#include "gtest/gtest.h"

namespace cs {

struct OptimizationTest {
  std::string description;
  std::string initialCircuit;
  std::size_t expectedMinimalGates{};
  std::size_t expectedMinimalDepth{};
};

inline void from_json(const nlohmann::json& j, OptimizationTest& test) {
  test.description          = j.at("description").get<std::string>();
  test.initialCircuit       = j.at("initial_circuit").get<std::string>();
  test.expectedMinimalGates = j.at("expected_minimal_gates").get<std::size_t>();
  test.expectedMinimalDepth = j.at("expected_minimal_depth").get<std::size_t>();
}

static std::vector<OptimizationTest> getTests(const std::string& path) {
  std::ifstream  input(path);
  nlohmann::json j;
  input >> j;
  return j;
}

class CircuitOptimizationTest
    : public ::testing::TestWithParam<OptimizationTest> {
protected:
  void SetUp() override {
    if (plog::get() == nullptr) {
      util::init();
    }

    const auto&       test = GetParam();
    std::stringstream ss(test.initialCircuit);
    qc.import(ss, qc::OpenQASM);
    std::cout << "Initial circuit:\n" << qc;

    targetTableau = Tableau(qc);
    std::cout << "Target tableau:\n" << targetTableau;

    synthesizer = CliffordSynthesizer(qc);

    config = Configuration();

    expectedMinimalGates = test.expectedMinimalGates;
    expectedMinimalDepth = test.expectedMinimalDepth;
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

  void consistencyCheck(const qc::QuantumComputation& circ) const {
    const auto circuitTableau = Tableau(circ);
    EXPECT_EQ(resultTableau, circuitTableau);
  }

  qc::QuantumComputation qc{};
  Tableau                targetTableau{};
  Configuration          config;
  CliffordSynthesizer    synthesizer;
  Results                results;
  Tableau                resultTableau;
  std::size_t            expectedMinimalGates{};
  std::size_t            expectedMinimalDepth{};
};

INSTANTIATE_TEST_SUITE_P(
    Json, CircuitOptimizationTest,
    testing::ValuesIn(getTests("cliffordsynthesis/circuits.json")),
    [](const testing::TestParamInfo<CircuitOptimizationTest::ParamType>& info) {
      return info.param.description;
    });

TEST_P(CircuitOptimizationTest, Gates) {
  config.target = TargetMetric::GATES;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getGates(), expectedMinimalGates);
}

TEST_P(CircuitOptimizationTest, Depth) {
  //  util::init("", plog::info);
  config.target = TargetMetric::DEPTH;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getDepth(), expectedMinimalDepth);
}

TEST_P(CircuitOptimizationTest, GatesMaxSAT) {
  config.target    = TargetMetric::GATES;
  config.useMaxSAT = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getGates(), expectedMinimalGates);
}

TEST_P(CircuitOptimizationTest, DepthMaxSAT) {
  config.target    = TargetMetric::DEPTH;
  config.useMaxSAT = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getDepth(), expectedMinimalDepth);
}

} // namespace cs
