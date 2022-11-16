#include "exact/ExactMapper.hpp"

#include "gtest/gtest.h"

class TestEncodings
    : public testing::TestWithParam<std::pair<Encoding, CommanderGrouping>> {
protected:
  std::string test_example_dir      = "./examples/";
  std::string test_architecture_dir = "./architectures/";
  std::string test_calibration_dir  = "./calibration/";

  qc::QuantumComputation       qc{};
  Configuration                settings{};
  Architecture                 arch{};
  std::unique_ptr<ExactMapper> mapper{};

  void SetUp() override {
    settings.verbose    = true;
    settings.method     = Method::Exact;
    settings.useSubsets = false;
  }
};

INSTANTIATE_TEST_SUITE_P(
    Encodings, TestEncodings,
    testing::Values(std::pair{Encoding::Naive, CommanderGrouping::Halves},
                    std::pair{Encoding::Commander, CommanderGrouping::Halves},
                    std::pair{Encoding::Commander, CommanderGrouping::Fixed2},
                    std::pair{Encoding::Commander, CommanderGrouping::Fixed3},
                    std::pair{Encoding::Bimander, CommanderGrouping::Halves},
                    std::pair{Encoding::Bimander, CommanderGrouping::Fixed2},
                    std::pair{Encoding::Bimander, CommanderGrouping::Fixed3}));

TEST_P(TestEncodings, ThreeToSevenQubits) {
  using namespace dd::literals;

  qc = qc::QuantumComputation(3U);
  qc.x(2, 1_pc);
  qc.x(1, 0_pc);

  arch.loadCouplingMap(AvailableArchitecture::IBMQ_Casablanca);

  mapper = std::make_unique<ExactMapper>(qc, arch);

  const auto& [encoding, grouping] = GetParam();
  settings.encoding                = encoding;
  settings.commanderGrouping       = grouping;

  mapper->map(settings);
  mapper->printResult(std::cout);

  ASSERT_FALSE(mapper->getResults().timeout);
  EXPECT_EQ(mapper->getResults().output.swaps, 0U);
}

TEST_P(TestEncodings, FiveToSevenQubits) {
  using namespace dd::literals;

  qc = qc::QuantumComputation(5U);
  qc.x(1, 0_pc);
  qc.x(2, 0_pc);
  qc.x(3, 0_pc);
  qc.x(4, 0_pc);

  arch.loadCouplingMap(AvailableArchitecture::IBMQ_Casablanca);

  mapper = std::make_unique<ExactMapper>(qc, arch);

  const auto& [encoding, grouping] = GetParam();
  settings.encoding                = encoding;
  settings.commanderGrouping       = grouping;

  mapper->map(settings);
  mapper->printResult(std::cout);

  ASSERT_FALSE(mapper->getResults().timeout);
  EXPECT_EQ(mapper->getResults().output.swaps, 1U);
}
