//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#include "gtest/gtest.h"

namespace cs {

struct TestConfiguration {
  // given input (either as tableau or as circuit)
  std::string description;
  std::string initialTableau;
  std::string targetTableau;
  std::string initialCircuit;

  // expected output
  std::size_t expectedMinimalGates{};
  std::size_t expectedMinimalDepth{};
  std::size_t expectedMinimalGatesAtMinimalDepth{};
  std::size_t expectedMinimalTwoQubitGates{};
  std::size_t expectedMinimalGatesAtMinimalTwoQubitGates{};
};

// NOLINTNEXTLINE (readability-identifier-naming)
inline void from_json(const nlohmann::json& j, TestConfiguration& test) {
  test.description = j.at("description").get<std::string>();
  if (j.contains("initial_tableau")) {
    test.initialTableau = j.at("initial_tableau").get<std::string>();
  }
  if (j.contains("target_tableau")) {
    test.targetTableau = j.at("target_tableau").get<std::string>();
  }
  if (j.contains("initial_circuit")) {
    test.initialCircuit = j.at("initial_circuit").get<std::string>();
  }

  test.expectedMinimalGates = j.at("expected_minimal_gates").get<std::size_t>();
  test.expectedMinimalDepth = j.at("expected_minimal_depth").get<std::size_t>();
  test.expectedMinimalGatesAtMinimalDepth =
      j.at("expected_minimal_gates_at_minimal_depth").get<std::size_t>();
  test.expectedMinimalTwoQubitGates =
      j.at("expected_minimal_two_qubit_gates").get<std::size_t>();
  test.expectedMinimalGatesAtMinimalTwoQubitGates =
      j.at("expected_minimal_gates_at_minimal_two_qubit_gates")
          .get<std::size_t>();
}

static std::vector<TestConfiguration> getTests(const std::string& path) {
  std::ifstream  input(path);
  nlohmann::json j;
  input >> j;
  return j;
}

class SynthesisTest : public ::testing::TestWithParam<TestConfiguration> {
protected:
  void SetUp() override {
    if (plog::get() == nullptr) {
      util::init();
    }

    test = GetParam();

    if (!test.initialCircuit.empty()) {
      std::stringstream      ss(test.initialCircuit);
      qc::QuantumComputation qc{};
      qc.import(ss, qc::Format::OpenQASM);
      std::cout << "Initial circuit:\n" << qc;
      targetTableau = Tableau(qc);
      if (test.initialTableau.empty()) {
        initialTableau = Tableau(qc.getNqubits());
        synthesizer    = CliffordSynthesizer(qc);
      } else {
        initialTableau = Tableau(test.initialTableau);
        std::cout << "Initial tableau:\n" << initialTableau;
        synthesizer = CliffordSynthesizer(initialTableau, qc);
      }
    } else {
      targetTableau = Tableau(test.targetTableau);
      if (test.initialTableau.empty()) {
        initialTableau = Tableau(targetTableau.getQubitCount());
        synthesizer    = CliffordSynthesizer(targetTableau);
      } else {
        initialTableau = Tableau(test.initialTableau);
        std::cout << "Initial tableau:\n" << initialTableau;
        synthesizer = CliffordSynthesizer(initialTableau, targetTableau);
      }
    }
    std::cout << "Target tableau:\n" << targetTableau;

    config = Configuration();
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
  TestConfiguration   test;
};

INSTANTIATE_TEST_SUITE_P(
    Tableaus, SynthesisTest,
    testing::ValuesIn(getTests("cliffordsynthesis/tableaus.json")),
    [](const testing::TestParamInfo<SynthesisTest::ParamType>& inf) {
      return inf.param.description;
    });

INSTANTIATE_TEST_SUITE_P(
    Circuits, SynthesisTest,
    testing::ValuesIn(getTests("cliffordsynthesis/circuits.json")),
    [](const testing::TestParamInfo<SynthesisTest::ParamType>& inf) {
      return inf.param.description;
    });

TEST_P(SynthesisTest, Gates) {
  config.target = TargetMetric::Gates;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getGates(), test.expectedMinimalGates);
}

TEST_P(SynthesisTest, GatesMaxSAT) {
  config.target    = TargetMetric::Gates;
  config.useMaxSAT = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getGates(), test.expectedMinimalGates);
}

TEST_P(SynthesisTest, Depth) {
  config.target = TargetMetric::Depth;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getDepth(), test.expectedMinimalDepth);
}

TEST_P(SynthesisTest, DepthMaxSAT) {
  config.target    = TargetMetric::Depth;
  config.useMaxSAT = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getDepth(), test.expectedMinimalDepth);
}

TEST_P(SynthesisTest, DepthMinimalGates) {
  config.target                              = TargetMetric::Depth;
  config.minimizeGatesAfterDepthOptimization = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getDepth(), test.expectedMinimalDepth);
  EXPECT_EQ(results.getGates(), test.expectedMinimalGatesAtMinimalDepth);
}

TEST_P(SynthesisTest, DepthMinimalGatesMaxSAT) {
  config.target                              = TargetMetric::Depth;
  config.useMaxSAT                           = true;
  config.minimizeGatesAfterDepthOptimization = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getDepth(), test.expectedMinimalDepth);
  EXPECT_EQ(results.getGates(), test.expectedMinimalGatesAtMinimalDepth);
}

TEST_P(SynthesisTest, TwoQubitGates) {
  config.target = TargetMetric::TwoQubitGates;
  config.tryHigherGateLimitForTwoQubitGateOptimization = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getTwoQubitGates(), test.expectedMinimalTwoQubitGates);
}

TEST_P(SynthesisTest, TwoQubitGatesMaxSAT) {
  config.target = TargetMetric::TwoQubitGates;
  config.tryHigherGateLimitForTwoQubitGateOptimization = true;
  config.useMaxSAT                                     = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getTwoQubitGates(), test.expectedMinimalTwoQubitGates);
}

TEST_P(SynthesisTest, TwoQubitGatesMinimalGates) {
  config.target = TargetMetric::TwoQubitGates;
  config.tryHigherGateLimitForTwoQubitGateOptimization = true;
  config.minimizeGatesAfterTwoQubitGateOptimization    = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getTwoQubitGates(), test.expectedMinimalTwoQubitGates);
  EXPECT_EQ(results.getGates(),
            test.expectedMinimalGatesAtMinimalTwoQubitGates);
}

TEST_P(SynthesisTest, TwoQubitGatesMinimalGatesMaxSAT) {
  config.target = TargetMetric::TwoQubitGates;
  config.tryHigherGateLimitForTwoQubitGateOptimization = true;
  config.minimizeGatesAfterTwoQubitGateOptimization    = true;
  config.useMaxSAT                                     = true;
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getTwoQubitGates(), test.expectedMinimalTwoQubitGates);
  EXPECT_EQ(results.getGates(),
            test.expectedMinimalGatesAtMinimalTwoQubitGates);
}

TEST(SynthesisTest, CompareFullToHalf) {
  Tableau             initialTableau;
  Tableau             targetTableau;
  Configuration       config;
  CliffordSynthesizer synthesizer;
  Results             results;
  Tableau             resultTableau;

  initialTableau = Tableau(1);
  targetTableau  = Tableau("[+X, -Y]");

  config.target    = TargetMetric::Gates;
  config.useMaxSAT = true;

  synthesizer = CliffordSynthesizer(initialTableau, targetTableau);
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getGates(), 1);

  initialTableau = Tableau(1, true);
  targetTableau  = Tableau("[+X, -Y]", "[+X, -Z]");

  config.target    = TargetMetric::Gates;
  config.useMaxSAT = true;

  synthesizer = CliffordSynthesizer(initialTableau, targetTableau);
  synthesizer.synthesize(config);
  results = synthesizer.getResults();

  EXPECT_EQ(results.getGates(), 0);
}

} // namespace cs
