/*
 * This file is part of the MQT QMAP library which is released under the MIT
 * license. See file README.md or go to
 * https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
 */

#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#include "gtest/gtest.h"

using namespace cs;

class TestCliffordSynthesis : public testing::TestWithParam<std::string> {
protected:
  std::string testArchitectureDir = "./architectures/";
  std::string testCalibrationDir  = "./calibration/";
  std::string testExampleDir      = "./examples/cliffordexamples/";

  Architecture ibmqLondon{};

  std::unique_ptr<CliffordSynthesizer> londonOptimizer{};

  void SetUp() override {
    using namespace dd::literals;
    util::init();
    ibmqLondon.loadCouplingMap(testArchitectureDir + "ibmq_london.arch");
    ibmqLondon.loadProperties(testCalibrationDir + "ibmq_london.csv");

    londonOptimizer = std::make_unique<CliffordSynthesizer>();
  }
};

INSTANTIATE_TEST_SUITE_P(CliffordSynthesizer, TestCliffordSynthesis,
                         testing::Values("Stabilizer = ['+IZ', '+XI']",
                                         "Stabilizer = ['-XX', '+XI']",
                                         "Stabilizer = ['+XX', '+XI']",
                                         "Stabilizer = ['+XY', '+XI']",
                                         "Stabilizer = ['+XY', '+ZZ']",
                                         "Stabilizer = ['+XX', '+ZZ']",
                                         "Stabilizer = ['-XX', '+ZZ']",
                                         "Stabilizer = ['+XX', '-ZZ']",
                                         "Stabilizer = ['+ZI', '+IZ']"));

TEST_P(TestCliffordSynthesis, SimpleSynthesis) {
  const auto& line = GetParam();
  Tableau     tableau{};
  tableau.fromString(line);
  CliffordSynthesizer cs{};
  Configuration       configuration{};

  configuration.nqubits         = 2;
  configuration.initialTimestep = 10;
  configuration.initialTableau  = std::make_shared<Tableau>(2);
  configuration.targetTableau   = std::make_shared<Tableau>(tableau);
  cs.synthesize(configuration);

  cs.optimalResults.dump(std::cout);

  EXPECT_EQ(cs.optimalResults.result, logicbase::Result::SAT);
  EXPECT_LE(
      cs.optimalResults.singleQubitGates + cs.optimalResults.twoQubitGates, 5);
  EXPECT_GT(
      cs.optimalResults.singleQubitGates + cs.optimalResults.twoQubitGates, 1);
}

TEST(TestCliffordSynthesis, SanityCheckDepth1) {
  using namespace dd::literals;
  util::init();
  qc::QuantumComputation qc{};
  CliffordSynthesizer    cs{};
  Configuration configuration{false, 2, 10, OptimizationStrategy::UseMinimizer,
                              TargetMetric::DEPTH};
  configuration.verbosity = 5;
  qc.addQubitRegister(2U);
  qc.h(0);
  qc.h(0);
  qc.h(0);
  qc.h(0);
  qc.h(0);

  configuration.targetCircuit =
      std::make_shared<qc::QuantumComputation>(qc.clone());

  cs.synthesize(configuration);

  EXPECT_EQ(cs.optimalResults.depth, 1);
}

TEST(TestCliffordSynthesis, SanityCheckDepth2) {
  using namespace dd::literals;
  util::init();
  qc::QuantumComputation qc{};
  CliffordSynthesizer    cs{};
  Configuration configuration{false, 2, 10, OptimizationStrategy::UseMinimizer,
                              TargetMetric::DEPTH};
  configuration.verbosity = 5;

  qc.addQubitRegister(2U);
  qc.h(0);
  qc.h(1);
  qc.x(0, 1_pc);
  qc.h(1);
  qc.h(0);

  configuration.targetCircuit =
      std::make_shared<qc::QuantumComputation>(qc.clone());

  cs.synthesize(configuration);

  EXPECT_EQ(cs.optimalResults.depth, 1);
}

TEST(TestCliffordSynthesis, SanityCheckGates) {
  util::init();
  qc::QuantumComputation qc{};
  CliffordSynthesizer    cs{};
  Configuration configuration{false, 2, 10, OptimizationStrategy::UseMinimizer,
                              TargetMetric::GATES};
  configuration.verbosity = 5;
  qc.addQubitRegister(2U);
  qc.h(0);
  qc.h(0);
  qc.h(0);
  qc.h(0);
  qc.h(0);

  configuration.targetCircuit =
      std::make_shared<qc::QuantumComputation>(qc.clone());

  cs.synthesize(configuration);

  EXPECT_EQ(cs.optimalResults.singleQubitGates, 1);
}
TEST(TestCliffordSynthesis, SanityCheck2QubitGates) {
  using namespace dd::literals;
  util::init();
  qc::QuantumComputation qc{};
  CliffordSynthesizer    cs{};
  Configuration configuration{false, 2, 10, OptimizationStrategy::UseMinimizer,
                              TargetMetric::TWO_QUBIT_GATES};
  configuration.verbosity = 5;
  qc.addQubitRegister(2U);
  qc.x(0, 1_pc);
  qc.x(0, 1_pc);
  qc.x(0, 1_pc);

  configuration.targetCircuit =
      std::make_shared<qc::QuantumComputation>(qc.clone());

  cs.synthesize(configuration);

  EXPECT_EQ(cs.optimalResults.twoQubitGates, 1);
}

TEST_P(TestCliffordSynthesis, TestDepthOpt) {
  const auto& line = GetParam();
  Tableau     tableau{};
  tableau.fromString(line);
  CliffordSynthesizer cs{};
  Configuration       configuration{};

  configuration.nqubits         = 2;
  configuration.initialTimestep = 4;
  configuration.target          = TargetMetric::DEPTH;
  configuration.initialTableau  = std::make_shared<Tableau>(2);
  configuration.targetTableau   = std::make_shared<Tableau>(tableau);
  cs.synthesize(configuration);

  cs.optimalResults.dump(std::cout);

  EXPECT_EQ(cs.optimalResults.result, logicbase::Result::SAT);
  EXPECT_LE(cs.optimalResults.depth, 4);
  EXPECT_GE(cs.optimalResults.depth, 1);
}

TEST_P(TestCliffordSynthesis, TestFidelityOpt) {
  const auto& line = GetParam();
  Tableau     tableau{};
  tableau.fromString(line);
  Configuration configuration{};

  configuration.nqubits         = 2;
  configuration.initialTimestep = 4;
  configuration.target          = TargetMetric::FIDELITY;
  configuration.strategy        = OptimizationStrategy::UseMinimizer;
  configuration.architecture    = ibmqLondon;
  configuration.initialTableau  = std::make_shared<Tableau>(2);
  configuration.targetTableau   = std::make_shared<Tableau>(tableau);
  londonOptimizer->synthesize(configuration);

  londonOptimizer->optimalResults.dump(std::cout);

  EXPECT_EQ(londonOptimizer->optimalResults.result, logicbase::Result::SAT);
  EXPECT_LE(londonOptimizer->optimalResults.singleQubitGates +
                londonOptimizer->optimalResults.twoQubitGates,
            6);
  EXPECT_GT(londonOptimizer->optimalResults.singleQubitGates +
                londonOptimizer->optimalResults.twoQubitGates,
            3);
}

TEST_P(TestCliffordSynthesis, TestTwoQubitGatesOpt) {
  const auto& line = GetParam();
  Tableau     tableau{};
  tableau.fromString(line);
  CliffordSynthesizer cs{};
  Configuration       configuration{};

  configuration.nqubits         = 2;
  configuration.initialTimestep = 6;
  configuration.target          = TargetMetric::TWO_QUBIT_GATES;
  configuration.initialTableau  = std::make_shared<Tableau>(2);
  configuration.targetTableau   = std::make_shared<Tableau>(tableau);
  cs.synthesize(configuration);

  cs.optimalResults.dump(std::cout);

  EXPECT_EQ(cs.optimalResults.result, logicbase::Result::SAT);
  EXPECT_LE(
      cs.optimalResults.singleQubitGates + cs.optimalResults.twoQubitGates, 6);
  EXPECT_GE(
      cs.optimalResults.singleQubitGates + cs.optimalResults.twoQubitGates, 2);
}

TEST_P(TestCliffordSynthesis, TestStartLow) {
  const auto& line = GetParam();
  Tableau     tableau{};
  tableau.fromString(line);
  CliffordSynthesizer cs{};
  Configuration       configuration{};

  configuration.nqubits         = 2;
  configuration.initialTimestep = 2;
  configuration.target          = TargetMetric::GATES;
  configuration.strategy        = OptimizationStrategy::StartLow;
  configuration.initialTableau  = std::make_shared<Tableau>(2);
  configuration.targetTableau   = std::make_shared<Tableau>(tableau);
  cs.synthesize(configuration);

  cs.optimalResults.dump(std::cout);

  EXPECT_EQ(cs.optimalResults.result, logicbase::Result::SAT);
  EXPECT_GE(cs.optimalResults.initialTimesteps, 2);
}

TEST_P(TestCliffordSynthesis, TestStartHigh) {
  const auto& line = GetParam();
  Tableau     tableau{};
  tableau.fromString(line);
  CliffordSynthesizer cs{};
  Configuration       configuration{};

  configuration.nqubits         = 2;
  configuration.initialTimestep = 100;
  configuration.target          = TargetMetric::GATES;
  configuration.strategy        = OptimizationStrategy::StartHigh;
  configuration.initialTableau  = std::make_shared<Tableau>(2);
  configuration.targetTableau   = std::make_shared<Tableau>(tableau);
  cs.synthesize(configuration);

  cs.optimalResults.dump(std::cout);

  EXPECT_EQ(cs.optimalResults.result, logicbase::Result::SAT);
  EXPECT_LT(cs.optimalResults.initialTimesteps, 100);
}

TEST_P(TestCliffordSynthesis, TestMinMax) {
  const auto& line = GetParam();
  Tableau     tableau{};
  tableau.fromString(line);
  CliffordSynthesizer cs{};
  Configuration       configuration{};

  configuration.nqubits         = 2;
  configuration.initialTimestep = 4;
  configuration.target          = TargetMetric::GATES;
  configuration.strategy        = OptimizationStrategy::MinMax;
  configuration.initialTableau  = std::make_shared<Tableau>(2);
  configuration.targetTableau   = std::make_shared<Tableau>(tableau);
  cs.synthesize(configuration);

  cs.optimalResults.dump(std::cout);

  EXPECT_EQ(cs.optimalResults.result, logicbase::Result::SAT);
  EXPECT_LE(
      cs.optimalResults.singleQubitGates + cs.optimalResults.twoQubitGates, 5);
  EXPECT_GT(
      cs.optimalResults.singleQubitGates + cs.optimalResults.twoQubitGates, 1);
}

TEST(TestCliffordSynthesis, TestSplitIter) {
  util::init();
  using namespace dd::literals;
  qc::QuantumComputation qc{};
  CliffordSynthesizer    optimizer{};
  qc.addQubitRegister(2U);

  qc.x(1);
  qc.z(0);
  qc.y(0);
  qc.x(1, 0_pc);
  qc.h(1);
  qc.s(0);
  qc.x(1);
  qc.sdag(1);
  qc.z(0);
  qc.y(1);
  qc.x(1);
  qc.z(0);
  qc.y(0);
  qc.x(1, 0_pc);
  qc.h(1);
  qc.s(0);
  qc.x(1);
  qc.sdag(1);
  qc.z(0);
  qc.y(1);
  qc.x(1);
  qc.z(0);
  qc.y(0);
  qc.x(1, 0_pc);
  qc.h(1);
  qc.s(0);
  qc.x(1);
  qc.sdag(1);
  qc.z(0);
  qc.y(1);
  qc.x(1);
  qc.z(0);
  qc.y(0);

  Configuration configuration{};
  configuration.nqubits         = 2;
  configuration.initialTimestep = 10;
  configuration.nThreads        = 1;
  configuration.targetCircuit =
      std::make_shared<qc::QuantumComputation>(qc.clone());
  configuration.target   = TargetMetric::GATES;
  configuration.strategy = OptimizationStrategy::SplitIter;

  optimizer.synthesize(configuration);
  EXPECT_EQ(optimizer.optimalResults.result, logicbase::Result::SAT);
}
