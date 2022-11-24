/*
 * This file is part of the MQT QMAP library which is released under the MIT
 * license. See file README.md or go to
 * https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
 */

#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#include "gtest/gtest.h"

using namespace cs;

class TestCliffordSynthesis
    : public testing::TestWithParam<
          std::tuple<OptimizationStrategy, TargetMetric, std::string>> {
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

INSTANTIATE_TEST_SUITE_P(
    CliffordSynthesizer, TestCliffordSynthesis,
    testing::Combine(testing::Values(OptimizationStrategy::UseMinimizer,
                                     OptimizationStrategy::BinarySearch,
                                     OptimizationStrategy::StartLow,
                                     OptimizationStrategy::StartHigh),
                     testing::Values(TargetMetric::TWO_QUBIT_GATES,
                                     TargetMetric::GATES, TargetMetric::DEPTH),
                     testing::Values("['+IZ', '+XI']", "['-XX', '+XI']",
                                     "['+XX', '+XI']", "['+XY', '+XI']",
                                     "['+XY', '+ZZ']", "['+XX', '+ZZ']",
                                     "['-XX', '+ZZ']", "['+XX', '-ZZ']",
                                     "['+ZI', '+IZ']")));

TEST_P(TestCliffordSynthesis, SimpleSynthesis) {
  const auto&         line = GetParam();
  const Tableau       tableau(std::get<2>(line));
  CliffordSynthesizer cs{};
  Configuration       configuration{};

  configuration.nqubits         = 2;
  configuration.initialTimestep = 10;
  configuration.initialTableau  = std::make_shared<Tableau>(2);
  configuration.targetTableau   = std::make_shared<Tableau>(tableau);

  configuration.strategy = std::get<0>(line);
  configuration.target   = std::get<1>(line);

  cs.synthesize(configuration);

  cs.optimalResults.dump(std::cout);

  EXPECT_EQ(cs.optimalResults.result, logicbase::Result::SAT);
}

TEST_P(TestCliffordSynthesis, TestFidelityOpt) {
  const auto&   line = GetParam();
  const Tableau tableau(std::get<2>(line));
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
            5);
  EXPECT_GE(londonOptimizer->optimalResults.singleQubitGates +
                londonOptimizer->optimalResults.twoQubitGates,
            0);
}

TEST(TestCliffordSynthesis, TestSanityCheckFidelity) {
  CouplingMap  cm = getFullyConnectedMap(2);
  Architecture architecture{2, cm, 0.999, 0.99, 0.995};

  qc::QuantumComputation qc{};
  CliffordSynthesizer    optimizer{};
  qc.addQubitRegister(2U);

  qc.h(0);
  qc.h(1);

  Configuration configuration{};
  configuration.nqubits         = 2;
  configuration.initialTimestep = 10;
  configuration.architecture    = architecture;
  configuration.fidelityScaling = 10000;
  configuration.targetCircuit =
      std::make_shared<qc::QuantumComputation>(qc.clone());
  configuration.target   = TargetMetric::FIDELITY;
  configuration.strategy = OptimizationStrategy::UseMinimizer;

  optimizer.synthesize(configuration);

  EXPECT_EQ(optimizer.optimalResults.fidelity, 0.999 * 0.999);
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
