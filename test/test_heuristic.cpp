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

TEST(Functionality, NodeCostCalculation) {
  const double         tolerance = 1e-6;
  const CouplingMap    cm        = {{0, 1}, {1, 2}, {3, 1}, {4, 3}};
  Architecture         arch{5, cm};
  TwoQubitMultiplicity multiplicity = {{{0, 1}, {5, 2}}, {{2, 3}, {0, 1}}};
  std::array<std::int16_t, MAX_DEVICE_QUBITS> qubits    = {4, 3, 1, 2, 0};
  std::array<std::int16_t, MAX_DEVICE_QUBITS> locations = {4, 2, 3, 1, 0};

  std::vector<std::vector<Exchange>> swaps = {
      {Exchange(0, 1, qc::OpType::Teleportation)},
      {Exchange(1, 2, qc::OpType::SWAP)}};

  HeuristicMapper::Node node(qubits, locations, swaps, 5.);
  EXPECT_NEAR(node.costFixed, 5., tolerance);
  node.updateHeuristicCost(arch, multiplicity, true);
  EXPECT_NEAR(node.costHeur,
              COST_UNIDIRECTIONAL_SWAP * 2 + COST_DIRECTION_REVERSE, tolerance);
  node.updateHeuristicCost(arch, multiplicity, false);
  EXPECT_NEAR(node.costHeur,
              COST_UNIDIRECTIONAL_SWAP * 14 + COST_DIRECTION_REVERSE * 3,
              tolerance);
  node.applySWAP({3, 4}, arch);
  node.updateHeuristicCost(arch, multiplicity, true);
  EXPECT_NEAR(node.costFixed, 5. + COST_UNIDIRECTIONAL_SWAP, tolerance);
  EXPECT_NEAR(node.costHeur, COST_UNIDIRECTIONAL_SWAP + COST_DIRECTION_REVERSE,
              tolerance);
  node.lookaheadPenalty = 0.;
  EXPECT_NEAR(node.getTotalCost(),
              5. + COST_UNIDIRECTIONAL_SWAP * 2 + COST_DIRECTION_REVERSE,
              tolerance);
  EXPECT_NEAR(node.getTotalFixedCost(), 5. + COST_UNIDIRECTIONAL_SWAP,
              tolerance);
  node.lookaheadPenalty = 2.;
  EXPECT_NEAR(node.getTotalCost(),
              7. + COST_UNIDIRECTIONAL_SWAP * 2 + COST_DIRECTION_REVERSE,
              tolerance);
  EXPECT_NEAR(node.getTotalFixedCost(), 7. + COST_UNIDIRECTIONAL_SWAP,
              tolerance);
  node.recalculateFixedCost(arch);
  EXPECT_NEAR(node.costFixed, COST_TELEPORTATION + COST_UNIDIRECTIONAL_SWAP * 2,
              tolerance);
  EXPECT_NEAR(node.costHeur, COST_UNIDIRECTIONAL_SWAP + COST_DIRECTION_REVERSE,
              tolerance);
  EXPECT_NEAR(node.getTotalCost(),
              2. + COST_TELEPORTATION + COST_UNIDIRECTIONAL_SWAP * 3 +
                  COST_DIRECTION_REVERSE,
              tolerance);
  EXPECT_NEAR(node.getTotalFixedCost(),
              2. + COST_TELEPORTATION + COST_UNIDIRECTIONAL_SWAP * 2,
              tolerance);
}

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

TEST(Functionality, HeuristicAdmissibility) {
  Architecture      architecture{};
  const CouplingMap cm = {{0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 3},
                          {3, 2}, {3, 4}, {4, 3}, {4, 5}, {5, 4}};
  architecture.loadCouplingMap(6, cm);
  const std::vector<Edge> perms{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}};

  TwoQubitMultiplicity multiplicity = {
      {{0, 4}, {1, 0}}, {{1, 3}, {1, 0}}, {{2, 5}, {1, 0}}};

  // perform depth-limited depth first search
  const std::size_t                  depthLimit = 11;
  std::vector<HeuristicMapper::Node> stack{};
  std::vector<std::size_t>           currentPerm{};

  auto initNode = HeuristicMapper::Node({0, 1, 2, 3, 4, 5}, {0, 1, 2, 3, 4, 5});
  initNode.recalculateFixedCost(architecture);
  initNode.updateHeuristicCost(architecture, multiplicity, true);
  stack.push_back(initNode);
  currentPerm.push_back(perms.size());

  while (!stack.empty()) {
    auto& node = stack.back();
    if (node.done) {
      // check if all nodes in stack have lower or equal cost
      for (auto& prevNode : stack) {
        EXPECT_LE(prevNode.getTotalCost(), node.getTotalCost());
      }
    }
    if (node.done || stack.size() >= depthLimit || currentPerm.back() == 0) {
      stack.pop_back();
      currentPerm.pop_back();
      continue;
    }
    --currentPerm.back();
    auto perm    = perms[currentPerm.back()];
    auto newNode = HeuristicMapper::Node(node.qubits, node.locations,
                                         node.swaps, node.costFixed);
    newNode.applySWAP(perm, architecture);
    newNode.updateHeuristicCost(architecture, multiplicity, true);
    stack.push_back(newNode);
    currentPerm.push_back(perms.size());
  }
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
