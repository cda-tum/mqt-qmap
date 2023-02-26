//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "heuristic/HeuristicMapper.hpp"

#include "gtest/gtest.h"

class HeuristicTest5Q : public testing::TestWithParam<std::string> {
protected:
  std::string testExampleDir      = "../examples/";
  std::string testArchitectureDir = "../extern/architectures/";
  std::string testCalibrationDir  = "../extern/calibration/";

  qc::QuantumComputation           qc{};
  Architecture                     ibmqYorktown{};
  Architecture                     ibmqLondon{};
  std::unique_ptr<HeuristicMapper> ibmqYorktownMapper;
  std::unique_ptr<HeuristicMapper> ibmqLondonMapper;

  void SetUp() override {
    qc.import(testExampleDir + GetParam() + ".qasm");
    ibmqYorktown.loadCouplingMap(AvailableArchitecture::IbmqYorktown);
    ibmqLondon.loadCouplingMap(testArchitectureDir + "ibmq_london.arch");
    ibmqLondon.loadProperties(testCalibrationDir + "ibmq_london.csv");
    ibmqYorktownMapper = std::make_unique<HeuristicMapper>(qc, ibmqYorktown);
    ibmqLondonMapper   = std::make_unique<HeuristicMapper>(qc, ibmqLondon);
  }
};

TEST(Functionality, EmptyDump) {
  qc::QuantumComputation qc{1};
  qc.x(0);
  Architecture    arch{1, {}};
  HeuristicMapper mapper(qc, arch);
  mapper.dumpResult("test.qasm");
  mapper.map({});
  EXPECT_NO_THROW(mapper.dumpResult("test.qasm"););
  EXPECT_NO_THROW(mapper.dumpResult("test.real"););
  EXPECT_THROW(mapper.dumpResult("test.dummy"), QMAPException);
}

TEST(Functionality, NoMeasurmentsAdded) {
  using namespace qc::literals;
  // construct circuit
  qc::QuantumComputation qc{4U};
  qc.x(1, 0_pc);
  qc.x(1, 2_pc);
  qc.x(1, 3_pc);

  // load architecture
  Architecture arch{};
  arch.loadCouplingMap(AvailableArchitecture::IbmqLondon);

  // create heuristic mapper
  HeuristicMapper mapper(qc, arch);

  // configure to not include measurements after mapping
  auto config                           = Configuration{};
  config.addMeasurementsToMappedCircuit = false;

  // perform the mapping
  mapper.map(config);

  // get the resulting circuit
  auto              qcMapped = qc::QuantumComputation();
  std::stringstream qasm{};
  mapper.dumpResult(qasm, qc::Format::OpenQASM);
  qcMapped.import(qasm, qc::Format::OpenQASM);

  // check no measurements were added
  EXPECT_EQ(qcMapped.getNops(), 3U);
  EXPECT_NE(qcMapped.back()->getType(), qc::Measure);
}

INSTANTIATE_TEST_SUITE_P(
    Heuristic, HeuristicTest5Q,
    testing::Values("3_17_13", "ex-1_166", "ham3_102", "miller_11", "4gt11_84",
                    "4mod5-v0_20", "mod5d1_63"),
    [](const testing::TestParamInfo<HeuristicTest5Q::ParamType>& inf) {
      std::string name = inf.param;
      std::replace(name.begin(), name.end(), '-', '_');
      return name;
    });

TEST_P(HeuristicTest5Q, Identity) {
  Configuration settings{};
  settings.initialLayout = InitialLayout::Identity;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() + "_heuristic_qx4_identity.qasm");
  ibmqYorktownMapper->printResult(std::cout);

  ibmqLondonMapper->map(settings);
  ibmqLondonMapper->dumpResult(GetParam() + "_heuristic_london_identity.qasm");
  ibmqLondonMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTest5Q, Static) {
  Configuration settings{};
  settings.initialLayout = InitialLayout::Static;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() + "_heuristic_qx4_static.qasm");
  ibmqYorktownMapper->printResult(std::cout);
  ibmqLondonMapper->map(settings);
  ibmqLondonMapper->dumpResult(GetParam() + "_heuristic_london_static.qasm");
  ibmqLondonMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTest5Q, Dynamic) {
  Configuration settings{};
  settings.initialLayout = InitialLayout::Dynamic;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() + "_heuristic_qx4_dynamic.qasm");
  ibmqYorktownMapper->printResult(std::cout);
  ibmqLondonMapper->map(settings);
  ibmqLondonMapper->dumpResult(GetParam() + "_heuristic_london_dynamic.qasm");
  ibmqLondonMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

class HeuristicTest16Q : public testing::TestWithParam<std::string> {
protected:
  std::string testExampleDir      = "../examples/";
  std::string testArchitectureDir = "../extern/architectures/";

  qc::QuantumComputation           qc{};
  Architecture                     ibmQX5{};
  std::unique_ptr<HeuristicMapper> ibmQX5Mapper;

  void SetUp() override {
    qc.import(testExampleDir + GetParam() + ".qasm");
    ibmQX5.loadCouplingMap(AvailableArchitecture::IbmQx5);
    ibmQX5Mapper = std::make_unique<HeuristicMapper>(qc, ibmQX5);
  }
};

INSTANTIATE_TEST_SUITE_P(
    Heuristic, HeuristicTest16Q,
    testing::Values("ising_model_10", "rd73_140", "cnt3-5_179", "qft_16"),
    [](const testing::TestParamInfo<HeuristicTest16Q::ParamType>& inf) {
      std::string name = inf.param;
      std::replace(name.begin(), name.end(), '-', '_');
      return name;
    });

TEST_P(HeuristicTest16Q, Dynamic) {
  Configuration settings{};
  settings.initialLayout = InitialLayout::Dynamic;
  ibmQX5Mapper->map(settings);
  ibmQX5Mapper->dumpResult(GetParam() + "_heuristic_qx5_dynamic.qasm");
  ibmQX5Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTest16Q, Disjoint) {
  Configuration settings{};
  settings.layering = Layering::DisjointQubits;
  ibmQX5Mapper->map(settings);
  ibmQX5Mapper->dumpResult(GetParam() + "_heuristic_qx5_disjoint.qasm");
  ibmQX5Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTest16Q, Disjoint2qBlocks) {
  Configuration settings{};
  settings.layering = Layering::Disjoint2qBlocks;
  ibmQX5Mapper->map(settings);
  ibmQX5Mapper->dumpResult(GetParam() + "_heuristic_qx5_disjoint_2q.qasm");
  ibmQX5Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

class HeuristicTest20Q : public testing::TestWithParam<std::string> {
protected:
  std::string testExampleDir      = "../examples/";
  std::string testArchitectureDir = "../extern/architectures/";

  qc::QuantumComputation           qc{};
  Architecture                     arch{};
  std::unique_ptr<HeuristicMapper> tokyoMapper;

  void SetUp() override {
    qc.import(testExampleDir + GetParam() + ".qasm");
    arch.loadCouplingMap(AvailableArchitecture::IbmqTokyo);
    tokyoMapper = std::make_unique<HeuristicMapper>(qc, arch);
  }
};

INSTANTIATE_TEST_SUITE_P(
    Heuristic, HeuristicTest20Q,
    testing::Values("ising_model_10", "rd73_140", "cnt3-5_179", "qft_16",
                    "z4_268"),
    [](const testing::TestParamInfo<HeuristicTest20Q::ParamType>& inf) {
      std::string name = inf.param;
      std::replace(name.begin(), name.end(), '-', '_');
      return name;
    });

TEST_P(HeuristicTest20Q, Dynamic) {
  Configuration settings{};
  settings.initialLayout = InitialLayout::Dynamic;
  tokyoMapper->map(settings);
  tokyoMapper->dumpResult(GetParam() + "_heuristic_tokyo_dynamic.qasm");
  tokyoMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

class HeuristicTest20QTeleport
    : public testing::TestWithParam<std::tuple<std::uint64_t, std::string>> {
protected:
  std::string testExampleDir      = "../examples/";
  std::string testArchitectureDir = "../extern/architectures/";

  qc::QuantumComputation           qc{};
  Architecture                     arch{};
  std::unique_ptr<HeuristicMapper> tokyoMapper;

  void SetUp() override {
    qc.import(testExampleDir + std::get<1>(GetParam()) + ".qasm");
    arch.loadCouplingMap(AvailableArchitecture::IbmqTokyo);
    tokyoMapper = std::make_unique<HeuristicMapper>(qc, arch);
  }
};

INSTANTIATE_TEST_SUITE_P(
    HeuristicTeleport, HeuristicTest20QTeleport,
    testing::Combine(testing::Values(1, 2, 3, 1337, 1338, 3147),
                     testing::Values("ising_model_10", "rd73_140", "cnt3-5_179",
                                     "qft_16", "z4_268")),
    [](const testing::TestParamInfo<HeuristicTest20QTeleport::ParamType>& inf) {
      std::string name = std::get<1>(inf.param);
      std::replace(name.begin(), name.end(), '-', '_');
      std::stringstream ss{};
      ss << name << "_seed" << std::get<0>(inf.param);
      return ss.str();
    });

TEST_P(HeuristicTest20QTeleport, Teleportation) {
  Configuration settings{};
  settings.initialLayout       = InitialLayout::Dynamic;
  settings.teleportationQubits = std::min(
      (arch.getNqubits() - qc.getNqubits()) & ~1U, static_cast<std::size_t>(8));
  settings.teleportationSeed = std::get<0>(GetParam());
  tokyoMapper->map(settings);
  tokyoMapper->dumpResult(std::get<1>(GetParam()) +
                          "_heuristic_tokyo_teleport.qasm");
  tokyoMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

class HeuristicTestFidelity : public testing::TestWithParam<std::string> {
protected:
  std::string testExampleDir      = "../examples/";
  std::string testArchitectureDir = "../extern/architectures/";
  std::string testCalibrationDir  = "../extern/calibration/";

  qc::QuantumComputation           qc{};
  Architecture                     arch{};
  Architecture                     nonFidelityArch{};
  std::unique_ptr<HeuristicMapper> mapper;
  std::unique_ptr<HeuristicMapper> nonFidelityMapper;

  void SetUp() override {
    qc.import(testExampleDir + GetParam() + ".qasm");
    arch.loadCouplingMap(testArchitectureDir + "ibmq_london.arch");
    arch.loadProperties(testCalibrationDir + "ibmq_london.csv");
    mapper = std::make_unique<HeuristicMapper>(qc, arch);
    nonFidelityArch.loadCouplingMap(AvailableArchitecture::IbmqYorktown);
    nonFidelityMapper = std::make_unique<HeuristicMapper>(qc, nonFidelityArch);
  }
};

INSTANTIATE_TEST_SUITE_P(
    Heuristic, HeuristicTestFidelity,
    testing::Values("3_17_13", "ex-1_166", "ham3_102", "miller_11", "4gt11_84",
                    "4mod5-v0_20", "mod5d1_63"),
    [](const testing::TestParamInfo<HeuristicTestFidelity::ParamType>& inf) {
      std::string name = inf.param;
      std::replace(name.begin(), name.end(), '-', '_');
      return name;
    });

TEST_P(HeuristicTestFidelity, Identity) {
  Configuration settings{};
  settings.layering         = Layering::DisjointQubits;
  settings.initialLayout    = InitialLayout::Identity;
  settings.considerFidelity = true;
  mapper->map(settings);
  mapper->dumpResult(GetParam() + "_heuristic_london_fidelity_identity.qasm");
  mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTestFidelity, Static) {
  Configuration settings{};
  settings.layering         = Layering::DisjointQubits;
  settings.initialLayout    = InitialLayout::Static;
  settings.considerFidelity = true;
  mapper->map(settings);
  mapper->dumpResult(GetParam() + "_heuristic_london_fidelity_static.qasm");
  mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTestFidelity, NoFidelity) {
  Configuration settings{};
  settings.layering         = Layering::DisjointQubits;
  settings.initialLayout    = InitialLayout::Static;
  settings.considerFidelity = true;
  nonFidelityMapper->map(settings);
  nonFidelityMapper->dumpResult(GetParam() +
                                "_heuristic_london_nofidelity.qasm");
  nonFidelityMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST(HeuristicTestFidelity, SimpleGrid) {
  Architecture      architecture{};
  const CouplingMap cm = {
      {0, 1}, {1, 0}, {1, 2},  {2, 1},  {2, 3},   {3, 2},

      {0, 4}, {4, 0}, {1, 5},  {5, 1},  {2, 6},   {6, 2},  {3, 7},  {7, 3},

      {4, 5}, {5, 4}, {5, 6},  {6, 5},  {6, 7},   {7, 6},

      {4, 8}, {8, 4}, {5, 9},  {9, 5},  {6, 10},  {10, 6}, {7, 11}, {11, 7},

      {8, 9}, {9, 8}, {9, 10}, {10, 9}, {10, 11}, {11, 10}};
  architecture.loadCouplingMap(12, cm);

  double e5 = 0.99;
  double e4 = 0.9;
  double e3 = 0.5;
  double e2 = 0.4;
  double e1 = 0.1;
  double e0 = 0.01;

  auto props = Architecture::Properties();
  props.setSingleQubitErrorRate(0, "x", e5);
  props.setSingleQubitErrorRate(1, "x", e5);
  props.setSingleQubitErrorRate(2, "x", e5);
  props.setSingleQubitErrorRate(3, "x", e5);
  props.setSingleQubitErrorRate(4, "x", e3);
  props.setSingleQubitErrorRate(5, "x", e3);
  props.setSingleQubitErrorRate(6, "x", e5);
  props.setSingleQubitErrorRate(7, "x", e3);
  props.setSingleQubitErrorRate(8, "x", e2);
  props.setSingleQubitErrorRate(9, "x", e1);
  props.setSingleQubitErrorRate(10, "x", e5);
  props.setSingleQubitErrorRate(11, "x", e2);

  props.setTwoQubitErrorRate(0, 1, e4);
  props.setTwoQubitErrorRate(1, 0, e4);
  props.setTwoQubitErrorRate(1, 2, e4);
  props.setTwoQubitErrorRate(2, 1, e4);
  props.setTwoQubitErrorRate(2, 3, e1);
  props.setTwoQubitErrorRate(3, 2, e1);

  props.setTwoQubitErrorRate(0, 4, e5);
  props.setTwoQubitErrorRate(4, 0, e5);
  props.setTwoQubitErrorRate(1, 5, e5);
  props.setTwoQubitErrorRate(5, 1, e5);
  props.setTwoQubitErrorRate(2, 6, e4);
  props.setTwoQubitErrorRate(6, 2, e4);
  props.setTwoQubitErrorRate(3, 7, e5);
  props.setTwoQubitErrorRate(7, 3, e5);

  props.setTwoQubitErrorRate(4, 5, e3);
  props.setTwoQubitErrorRate(5, 4, e3);
  props.setTwoQubitErrorRate(5, 6, e5);
  props.setTwoQubitErrorRate(6, 5, e5);
  props.setTwoQubitErrorRate(6, 7, e5);
  props.setTwoQubitErrorRate(7, 6, e5);

  props.setTwoQubitErrorRate(4, 8, e0);
  props.setTwoQubitErrorRate(8, 4, e0);
  props.setTwoQubitErrorRate(5, 9, e3);
  props.setTwoQubitErrorRate(9, 5, e3);
  props.setTwoQubitErrorRate(6, 10, e1);
  props.setTwoQubitErrorRate(10, 6, e1);
  props.setTwoQubitErrorRate(7, 11, e3);
  props.setTwoQubitErrorRate(11, 7, e3);

  props.setTwoQubitErrorRate(8, 9, e5);
  props.setTwoQubitErrorRate(9, 8, e5);
  props.setTwoQubitErrorRate(9, 10, e4);
  props.setTwoQubitErrorRate(10, 9, e4);
  props.setTwoQubitErrorRate(10, 11, e5);
  props.setTwoQubitErrorRate(11, 10, e5);

  architecture.loadProperties(props);

  qc::QuantumComputation qc{12};

  for (std::size_t i = 0; i < 50; ++i) {
    qc.x(4);
  }
  qc.x(5);
  qc.x(7);
  for (std::size_t i = 0; i < 5; ++i) {
    qc.x(0, qc::Control{3});
    qc.x(2, qc::Control{9});
  }

  auto mapper = std::make_unique<HeuristicMapper>(qc, architecture);

  Configuration settings{};
  settings.verbose                  = true;
  settings.layering                 = Layering::Disjoint2qBlocks;
  settings.initialLayout            = InitialLayout::Identity;
  settings.considerFidelity         = true;
  settings.preMappingOptimizations  = false;
  settings.postMappingOptimizations = false;
  mapper->map(settings);
  mapper->dumpResult("simple_grid_mapped.qasm");
  mapper->printResult(std::cout);

  auto& result = mapper->getResults();
  EXPECT_EQ(result.input.layers, 1);
  /*
  expected output (order of gates may vary):
  SWAP(2,6)
  SWAP(9,10)
  SWAP(0,1)
  SWAP(1,2)
  SWAP(4,5)
  SWAP(5,9)
  SWAP(4,8)
  X(8)
  X(7)
  X(9) [x50]
  CX(2,3) [x5]
  CX(6,10) [x5]
  */
  EXPECT_EQ(result.output.swaps, 7);

  double c4 = -std::log2(1 - e4);
  double c3 = -std::log2(1 - e3);
  double c2 = -std::log2(1 - e2);
  double c1 = -std::log2(1 - e1);
  double c0 = -std::log2(1 - e0);

  double expectedFidelity = 3 * c4 + 3 * c4 + 3 * c4 + 3 * c4 + 3 * c3 +
                            3 * c3 + 3 * c0 +   // SWAPs
                            c2 + c3 + 50 * c1 + // Xs
                            5 * c1 + 5 * c1;    // CXs
  EXPECT_NEAR(result.output.totalLogFidelity, expectedFidelity, 1e-6);
}

TEST(HeuristicTestFidelity, RemapSingleQubit) {
  Architecture      architecture{};
  const CouplingMap cm = {
      {0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 3},
      {3, 2}, {3, 4}, {4, 3}, {4, 5}, {5, 4},
  };
  architecture.loadCouplingMap(6, cm);

  double e5 = 0.99;
  double e4 = 0.9;
  double e3 = 0.5;
  double e1 = 0.1;
  double e0 = 0.01;

  auto props = Architecture::Properties();
  props.setSingleQubitErrorRate(0, "x", e5);
  props.setSingleQubitErrorRate(1, "x", e5);
  props.setSingleQubitErrorRate(2, "x", e5);
  props.setSingleQubitErrorRate(3, "x", e4);
  props.setSingleQubitErrorRate(4, "x", e4);
  props.setSingleQubitErrorRate(5, "x", e1);

  props.setTwoQubitErrorRate(0, 1, e1);
  props.setTwoQubitErrorRate(1, 0, e1);
  props.setTwoQubitErrorRate(1, 2, e3);
  props.setTwoQubitErrorRate(2, 1, e3);
  props.setTwoQubitErrorRate(2, 3, e5);
  props.setTwoQubitErrorRate(3, 2, e5);
  props.setTwoQubitErrorRate(3, 4, e0);
  props.setTwoQubitErrorRate(4, 3, e0);
  props.setTwoQubitErrorRate(4, 5, e0);
  props.setTwoQubitErrorRate(5, 4, e0);

  architecture.loadProperties(props);

  qc::QuantumComputation qc{12};
  for (std::size_t i = 0; i < 5; ++i) {
    qc.x(0, qc::Control{2});
    qc.x(3);
  }

  auto mapper = std::make_unique<HeuristicMapper>(qc, architecture);

  Configuration settings{};
  settings.layering                 = Layering::Disjoint2qBlocks;
  settings.initialLayout            = InitialLayout::Identity;
  settings.considerFidelity         = true;
  settings.preMappingOptimizations  = false;
  settings.postMappingOptimizations = false;
  mapper->map(settings);
  mapper->dumpResult("remap_single_qubit_mapped.qasm");
  mapper->printResult(std::cout);

  auto& result = mapper->getResults();
  EXPECT_EQ(result.input.layers, 1);
  /*
  expected output (order of gates may vary):
  SWAP(1,2)
  SWAP(3,4)
  SWAP(4,5)
  CX(0,1) [x5]
  X(5) [x5]
  */
  EXPECT_EQ(result.output.swaps, 4);

  double c3 = -std::log2(1 - e3);
  double c1 = -std::log2(1 - e1);
  double c0 = -std::log2(1 - e0);

  double expectedFidelity = 3 * c3 + 3 * c0 + 3 * c0 + // SWAPs
                            5 * c1 +                   // Xs
                            5 * c1;                    // CXs
  EXPECT_NEAR(result.output.totalLogFidelity, expectedFidelity, 1e-6);
}

TEST(HeuristicTestFidelity, QubitRideAlong) {
  Architecture      architecture{};
  const CouplingMap cm = {{0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 3}, {3, 2},
                          {1, 4}, {4, 1}, {2, 5}, {5, 2}, {5, 6}, {6, 5}};
  architecture.loadCouplingMap(7, cm);

  double e5 = 0.99;
  double e4 = 0.9;
  double e3 = 0.5;
  double e1 = 0.1;

  auto props = Architecture::Properties();
  props.setSingleQubitErrorRate(0, "x", e5);
  props.setSingleQubitErrorRate(1, "x", e5);
  props.setSingleQubitErrorRate(2, "x", e5);
  props.setSingleQubitErrorRate(3, "x", e5);
  props.setSingleQubitErrorRate(4, "x", e5);
  props.setSingleQubitErrorRate(5, "x", e5);
  props.setSingleQubitErrorRate(6, "x", e5);

  props.setTwoQubitErrorRate(0, 1, e4);
  props.setTwoQubitErrorRate(1, 0, e4);
  props.setTwoQubitErrorRate(1, 2, e3);
  props.setTwoQubitErrorRate(2, 1, e3);
  props.setTwoQubitErrorRate(2, 3, e1);
  props.setTwoQubitErrorRate(3, 2, e1);
  props.setTwoQubitErrorRate(1, 4, e1);
  props.setTwoQubitErrorRate(4, 1, e1);
  props.setTwoQubitErrorRate(2, 5, e3);
  props.setTwoQubitErrorRate(5, 2, e3);
  props.setTwoQubitErrorRate(5, 6, e4);
  props.setTwoQubitErrorRate(6, 5, e4);

  architecture.loadProperties(props);

  qc::QuantumComputation qc{12};
  for (std::size_t i = 0; i < 5; ++i) {
    qc.x(0, qc::Control{3});
    qc.x(4, qc::Control{6});
  }

  auto mapper = std::make_unique<HeuristicMapper>(qc, architecture);

  Configuration settings{};
  settings.layering                 = Layering::Disjoint2qBlocks;
  settings.initialLayout            = InitialLayout::Identity;
  settings.considerFidelity         = true;
  settings.preMappingOptimizations  = false;
  settings.postMappingOptimizations = false;
  mapper->map(settings);
  mapper->dumpResult("qubit_ride_along_mapped.qasm");
  mapper->printResult(std::cout);

  auto& result = mapper->getResults();
  EXPECT_EQ(result.input.layers, 1);
  /*
  expected output (order of gates may vary):
  SWAP(5,6)
  SWAP(2,5)
  SWAP(0,1)
  SWAP(1,2)
  CX(2,3) [x5]
  CX(4,1) [x5]
  */
  // EXPECT_EQ(result.output.swaps, 4);

  double c4 = -std::log2(1 - e4);
  double c3 = -std::log2(1 - e3);
  double c1 = -std::log2(1 - e1);

  double expectedFidelity = 3 * c4 + 3 * c3 + 3 * c4 + 3 * c3 + // SWAPs
                            5 * c1 + 5 * c1;                    // CXs
  EXPECT_NEAR(result.output.totalLogFidelity, expectedFidelity, 1e-6);
}
