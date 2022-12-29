//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "exact/ExactMapper.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

class ExactTest : public testing::TestWithParam<std::string> {
protected:
  std::string testExampleDir      = "./examples/";
  std::string testArchitectureDir = "./architectures/";
  std::string testCalibrationDir  = "./calibration/";

  qc::QuantumComputation       qc{};
  Configuration                settings{};
  Architecture                 ibmqYorktown{};
  Architecture                 ibmqLondon{};
  Architecture                 ibmQX4{};
  std::unique_ptr<ExactMapper> ibmqYorktownMapper;
  std::unique_ptr<ExactMapper> ibmqLondonMapper;
  std::unique_ptr<ExactMapper> ibmQX4Mapper;

  void SetUp() override {
    using namespace qc::literals;

    if (::testing::UnitTest::GetInstance()
            ->current_test_info()
            ->value_param() != nullptr) {
      qc.import(testExampleDir + GetParam() + ".qasm");
    } else {
      qc.addQubitRegister(3U);
      qc.x(0, 1_pc);
      qc.x(1, 2_pc);
      qc.x(2, 0_pc);
    }
    ibmqYorktown.loadCouplingMap(AvailableArchitecture::IBMQ_Yorktown);
    ibmqLondon.loadCouplingMap(testArchitectureDir + "ibmq_london.arch");
    ibmqLondon.loadProperties(testCalibrationDir + "ibmq_london.csv");
    ibmQX4.loadCouplingMap(AvailableArchitecture::IBM_QX4);

    ibmqYorktownMapper = std::make_unique<ExactMapper>(qc, ibmqYorktown);
    ibmqLondonMapper   = std::make_unique<ExactMapper>(qc, ibmqLondon);
    ibmQX4Mapper       = std::make_unique<ExactMapper>(qc, ibmQX4);

    settings.verbose = true;
    settings.method  = Method::Exact;
  }
};

INSTANTIATE_TEST_SUITE_P(
    Exact, ExactTest,
    testing::Values("3_17_13", "ex-1_166", "ham3_102", "miller_11", "4gt11_84"),
    [](const testing::TestParamInfo<ExactTest::ParamType>& info) {
      std::string name = info.param;
      std::replace(name.begin(), name.end(), '-', '_');
      std::stringstream ss{};
      ss << name;
      return ss.str();
    });

TEST_P(ExactTest, IndividualGates) {
  settings.layering = Layering::IndividualGates;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() +
                                 "_exact_yorktown_individual.qasm");
  ibmqYorktownMapper->printResult(std::cout);

  ibmqLondonMapper->map(settings);
  ibmqLondonMapper->dumpResult(GetParam() + "_exact_london_individual.qasm");
  ibmqLondonMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, DisjointQubits) {
  settings.layering = Layering::DisjointQubits;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() + "_exact_yorktown_disjoint.qasm");
  ibmqYorktownMapper->printResult(std::cout);

  ibmqLondonMapper->map(settings);
  ibmqLondonMapper->dumpResult(GetParam() + "_exact_london_disjoint.qasm");
  ibmqLondonMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, OddGates) {
  settings.layering = Layering::OddGates;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() + "_exact_yorktown_odd.qasm");
  ibmqYorktownMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, QubitTriangle) {
  settings.layering = Layering::QubitTriangle;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() + "_exact_yorktown_triangle.qasm");
  ibmqYorktownMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, CommanderEncodingfixed3) {
  settings.encoding          = Encoding::Commander;
  settings.commanderGrouping = CommanderGrouping::Fixed3;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() +
                                 "_exact_yorktown_commander_fixed3.qasm");
  ibmqYorktownMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, CommanderEncodingfixed2) {
  settings.encoding          = Encoding::Commander;
  settings.commanderGrouping = CommanderGrouping::Fixed2;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() +
                                 "_exact_yorktown_commander_fixed2.qasm");
  ibmqYorktownMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, CommanderEncodinghalves) {
  settings.encoding          = Encoding::Commander;
  settings.commanderGrouping = CommanderGrouping::Halves;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() +
                                 "_exact_yorktown_commander_halves.qasm");
  ibmqYorktownMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, CommanderEncodinglogarithm) {
  settings.encoding          = Encoding::Commander;
  settings.commanderGrouping = CommanderGrouping::Logarithm;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() +
                                 "_exact_yorktown_commander_log.qasm");
  ibmqYorktownMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, CommanderEncodingUnidirectionalfixed3) {
  settings.encoding          = Encoding::Commander;
  settings.commanderGrouping = CommanderGrouping::Fixed3;
  ibmQX4Mapper->map(settings);
  ibmQX4Mapper->dumpResult(GetParam() + "_exact_QX4_commander_fixed3.qasm");
  ibmQX4Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, CommanderEncodingUnidirectionalfixed2) {
  settings.encoding          = Encoding::Commander;
  settings.commanderGrouping = CommanderGrouping::Fixed2;
  ibmQX4Mapper->map(settings);
  ibmQX4Mapper->dumpResult(GetParam() + "_exact_QX4_commander_fixed2.qasm");
  ibmQX4Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, CommanderEncodingUnidirectionalhalves) {
  settings.encoding          = Encoding::Commander;
  settings.commanderGrouping = CommanderGrouping::Halves;
  ibmQX4Mapper->map(settings);
  ibmQX4Mapper->dumpResult(GetParam() + "_exact_QX4_commander_halves.qasm");
  ibmQX4Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, CommanderEncodingUnidirectionallogarithm) {
  settings.encoding          = Encoding::Commander;
  settings.commanderGrouping = CommanderGrouping::Logarithm;
  ibmQX4Mapper->map(settings);
  ibmQX4Mapper->dumpResult(GetParam() + "_exact_QX4_commander_log.qasm");
  ibmQX4Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, BimanderEncodingfixed3) {
  settings.encoding          = Encoding::Bimander;
  settings.commanderGrouping = CommanderGrouping::Fixed3;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() + "_exact_yorktown_bimander.qasm");
  ibmqYorktownMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, BimanderEncodingfixed2) {
  settings.encoding          = Encoding::Bimander;
  settings.commanderGrouping = CommanderGrouping::Fixed2;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() + "_exact_yorktown_bimander.qasm");
  ibmqYorktownMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, BimanderEncodinghalves) {
  settings.encoding          = Encoding::Bimander;
  settings.commanderGrouping = CommanderGrouping::Halves;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() + "_exact_yorktown_bimander.qasm");
  ibmqYorktownMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, BimanderEncodinglogaritm) {
  settings.encoding          = Encoding::Bimander;
  settings.commanderGrouping = CommanderGrouping::Logarithm;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() + "_exact_yorktown_bimander.qasm");
  ibmqYorktownMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, BimanderEncodingUnidirectionalfixed3) {
  settings.encoding          = Encoding::Bimander;
  settings.commanderGrouping = CommanderGrouping::Fixed3;
  ibmQX4Mapper->map(settings);
  ibmQX4Mapper->dumpResult(GetParam() + "_exact_QX4_bimander.qasm");
  ibmQX4Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, BimanderEncodingUnidirectionalfixed2) {
  settings.encoding          = Encoding::Bimander;
  settings.commanderGrouping = CommanderGrouping::Fixed2;
  ibmQX4Mapper->map(settings);
  ibmQX4Mapper->dumpResult(GetParam() + "_exact_QX4_bimander.qasm");
  ibmQX4Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, BimanderEncodingUnidirectionalhalves) {
  settings.encoding          = Encoding::Bimander;
  settings.commanderGrouping = CommanderGrouping::Halves;
  ibmQX4Mapper->map(settings);
  ibmQX4Mapper->dumpResult(GetParam() + "_exact_QX4_bimander.qasm");
  ibmQX4Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, BimanderEncodingUnidirectionallogarithm) {
  settings.encoding          = Encoding::Bimander;
  settings.commanderGrouping = CommanderGrouping::Logarithm;
  ibmQX4Mapper->map(settings);
  ibmQX4Mapper->dumpResult(GetParam() + "_exact_QX4_bimander.qasm");
  ibmQX4Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, LimitsBidirectional) {
  settings.enableSwapLimits = true;
  settings.useSubsets       = false;
  settings.swapReduction    = SwapReduction::CouplingLimit;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() +
                                 "_exact_yorktown_swapreduct.qasm");
  ibmqYorktownMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, LimitsBidirectionalSubsetSwaps) {
  settings.enableSwapLimits = true;
  settings.useSubsets       = true;
  settings.swapReduction    = SwapReduction::CouplingLimit;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() +
                                 "_exact_yorktown_swapreduct.qasm");
  ibmqYorktownMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, LimitsBidirectionalCustomLimit) {
  settings.enableSwapLimits = true;
  settings.swapReduction    = SwapReduction::Custom;
  settings.swapLimit        = 10;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() +
                                 "_exact_yorktown_swapreduct.qasm");
  ibmqYorktownMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, LimitsUnidirectional) {
  settings.enableSwapLimits = true;
  settings.useSubsets       = false;
  settings.swapReduction    = SwapReduction::CouplingLimit;
  ibmQX4Mapper->map(settings);
  ibmQX4Mapper->dumpResult(GetParam() + "_exact_QX4_swapreduct.qasm");
  ibmQX4Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, LimitsUnidirectionalSubsetSwaps) {
  settings.enableSwapLimits = true;
  settings.useSubsets       = true;
  settings.swapReduction    = SwapReduction::CouplingLimit;
  ibmQX4Mapper->map(settings);
  ibmQX4Mapper->dumpResult(GetParam() + "_exact_QX4_swapreduct.qasm");
  ibmQX4Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, LimitsUnidirectionalCustomLimit) {
  settings.enableSwapLimits = true;
  settings.swapReduction    = SwapReduction::Custom;
  settings.swapLimit        = 10;
  ibmQX4Mapper->map(settings);
  ibmQX4Mapper->dumpResult(GetParam() + "_exact_QX4_swapreduct.qasm");
  ibmQX4Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, IncreasingCustomLimitUnidirectional) {
  settings.enableSwapLimits = true;
  settings.swapReduction    = SwapReduction::Increasing;
  settings.swapLimit        = 3;
  ibmQX4Mapper->map(settings);
  ibmQX4Mapper->dumpResult(GetParam() + "_exact_QX4_swapreduct_inccustom.qasm");
  ibmQX4Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, IncreasingUnidirectional) {
  settings.enableSwapLimits = true;
  settings.swapReduction    = SwapReduction::Increasing;
  settings.swapLimit        = 0;
  ibmQX4Mapper->map(settings);
  ibmQX4Mapper->dumpResult(GetParam() + "_exact_QX4_swapreduct_inc.qasm");
  ibmQX4Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, LimitsBidirectionalBDD) {
  settings.enableSwapLimits = true;
  settings.useSubsets       = false;
  settings.useBDD           = true;
  settings.swapReduction    = SwapReduction::CouplingLimit;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() +
                                 "_exact_yorktown_swapreduct_bdd.qasm");
  ibmqYorktownMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, LimitsBidirectionalSubsetSwapsBDD) {
  settings.enableSwapLimits = true;
  settings.useSubsets       = true;
  settings.useBDD           = true;
  settings.swapReduction    = SwapReduction::CouplingLimit;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() +
                                 "_exact_yorktown_swapreduct_bdd.qasm");
  ibmqYorktownMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, NoSubsets) {
  settings.useSubsets       = false;
  settings.enableSwapLimits = false;
  ibmQX4Mapper->map(settings);
  ibmQX4Mapper->dumpResult(GetParam() + "_exact_QX4_nosubsets.qasm");
  ibmQX4Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, toStringMethods) {
  EXPECT_EQ(toString(InitialLayout::Identity), "identity");
  EXPECT_EQ(toString(InitialLayout::Static), "static");
  EXPECT_EQ(toString(InitialLayout::Dynamic), "dynamic");
  EXPECT_EQ(toString(InitialLayout::None), "none");

  EXPECT_EQ(toString(Layering::IndividualGates), "individual_gates");
  EXPECT_EQ(toString(Layering::DisjointQubits), "disjoint_qubits");
  EXPECT_EQ(toString(Layering::OddGates), "odd_gates");
  EXPECT_EQ(toString(Layering::QubitTriangle), "qubit_triangle");
  EXPECT_EQ(toString(Layering::None), "none");

  EXPECT_EQ(toString(Encoding::Naive), "naive");
  EXPECT_EQ(toString(Encoding::Commander), "commander");
  EXPECT_EQ(toString(Encoding::Bimander), "bimander");

  EXPECT_EQ(toString(CommanderGrouping::Fixed2), "fixed2");
  EXPECT_EQ(toString(CommanderGrouping::Fixed3), "fixed3");
  EXPECT_EQ(toString(CommanderGrouping::Logarithm), "logarithm");
  EXPECT_EQ(toString(CommanderGrouping::Halves), "halves");

  EXPECT_EQ(toString(SwapReduction::CouplingLimit), "coupling_limit");
  EXPECT_EQ(toString(SwapReduction::Custom), "custom");
  EXPECT_EQ(toString(SwapReduction::None), "none");
  EXPECT_EQ(toString(SwapReduction::Increasing), "increasing");

  SUCCEED() << "ToStringMethods working";
}

TEST_F(ExactTest, CircuitWithOnlySingleQubitGates) {
  qc.clear();
  qc.x(0);
  qc.x(1);
  ibmQX4Mapper = std::make_unique<ExactMapper>(qc, ibmQX4);
  ibmQX4Mapper->map(settings);
  ibmQX4Mapper->dumpResult(std::cout, qc::Format::OpenQASM);
  SUCCEED() << "Mapping successful";
}

TEST_F(ExactTest, MapToSubsetNotIncludingQ0) {
  const CouplingMap cm{{0, 1}, {1, 0}, {1, 2}, {2, 1},
                       {2, 3}, {3, 2}, {1, 3}, {3, 1}};
  Architecture      arch(4U, cm);

  auto mapper         = ExactMapper(qc, arch);
  settings.useSubsets = false;
  mapper.map(settings);

  std::ostringstream oss{};
  mapper.dumpResult(oss, qc::Format::OpenQASM);
  auto               qcMapped = qc::QuantumComputation();
  std::istringstream iss{oss.str()};
  qcMapped.import(iss, qc::Format::OpenQASM);
  std::cout << qcMapped << std::endl;
  EXPECT_EQ(qcMapped.initialLayout.size(), 4U);
  EXPECT_EQ(qcMapped.initialLayout[0], 3);
  EXPECT_EQ(qcMapped.outputPermutation.size(), 3U);
  EXPECT_TRUE(qcMapped.garbage.at(3));
}

TEST_F(ExactTest, WCNF) {
  settings.verbose     = false;
  settings.includeWCNF = true;
  ibmqLondonMapper->map(settings);
  ibmqLondonMapper->printResult(std::cout);
  const auto& wcnf = ibmqLondonMapper->getResults().wcnf;
  EXPECT_TRUE(!wcnf.empty());
}

TEST_F(ExactTest, WCNF_not_available) {
  using namespace qc::literals;

  settings.verbose     = false;
  settings.includeWCNF = true;

  auto circ = qc::QuantumComputation(5U);
  circ.h(0);
  circ.x(1, 0_pc);
  circ.x(2, 0_pc);
  circ.x(3, 0_pc);
  circ.x(4, 0_pc);

  auto mapper = ExactMapper(circ, ibmqLondon);

  mapper.map(settings);
  EXPECT_TRUE(mapper.getResults().wcnf.empty());

  auto mapper2      = ExactMapper(circ, ibmqLondon);
  settings.encoding = Encoding::Commander;
  mapper2.map(settings);
  EXPECT_FALSE(mapper2.getResults().wcnf.empty());
}

TEST_F(ExactTest, MapToSubgraph) {
  const auto connectedSubset = std::set<std::uint16_t>{0U, 1U, 2U};

  settings.subgraph = connectedSubset;
  ibmqLondonMapper->map(settings);
  const auto& results = ibmqLondonMapper->getResults();
  EXPECT_FALSE(results.timeout);
}

TEST_F(ExactTest, MapToSubgraphTooSmall) {
  const auto tooSmallSubset = std::set<std::uint16_t>{0U, 1U};

  settings.subgraph = tooSmallSubset;
  ibmqLondonMapper->map(settings);
  const auto& results = ibmqLondonMapper->getResults();
  EXPECT_TRUE(results.timeout);
}

TEST_F(ExactTest, MapToSubgraphNotConnected) {
  const auto nonConnectedSubset = std::set<std::uint16_t>{0U, 2U, 3U};

  settings.subgraph = nonConnectedSubset;
  ibmqLondonMapper->map(settings);
  const auto& results = ibmqLondonMapper->getResults();
  EXPECT_TRUE(results.timeout);
}
TEST_F(ExactTest, CommanderEncodingRigettiArch) {
  Architecture aspen;
  aspen.loadCouplingMap(AvailableArchitecture::Rigetti_Aspen);
  Architecture agave;
  agave.loadCouplingMap(AvailableArchitecture::Rigetti_Agave);

  auto aspenMapper = ExactMapper(qc, aspen);
  auto agaveMapper = ExactMapper(qc, agave);
  aspenMapper.map(settings);
  agaveMapper.map(settings);
  aspenMapper.printResult(std::cout);
  agaveMapper.printResult(std::cout);

  SUCCEED() << "Mapping successful";
}

TEST_F(ExactTest, NoMeasurmentsAdded) {
  // configure to not include measurements after mapping
  settings.addMeasurementsToMappedCircuit = false;

  // perform the mapping
  ibmqLondonMapper->map(settings);

  // get the resulting circuit
  auto              qcMapped = qc::QuantumComputation();
  std::stringstream qasm{};
  ibmqLondonMapper->dumpResult(qasm, qc::Format::OpenQASM);
  qcMapped.import(qasm, qc::Format::OpenQASM);

  // check no measurements were added
  EXPECT_EQ(qcMapped.getNops(), 4U);
  EXPECT_NE(qcMapped.back()->getType(), qc::Measure);
}
