/*
 * This file is part of the MQT QMAP library which is released under the MIT
 * license. See file README.md or go to
 * https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
 */

#include "heuristic/HeuristicMapper.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

class HeuristicTest5Q : public testing::TestWithParam<std::string> {
protected:
  std::string test_example_dir      = "./examples/";
  std::string test_architecture_dir = "./architectures/";
  std::string test_calibration_dir  = "./calibration/";

  qc::QuantumComputation           qc{};
  Architecture                     IBMQ_Yorktown{};
  Architecture                     IBMQ_London{};
  std::unique_ptr<HeuristicMapper> IBMQ_Yorktown_mapper;
  std::unique_ptr<HeuristicMapper> IBMQ_London_mapper;

  void SetUp() override {
    qc.import(test_example_dir + GetParam() + ".qasm");
    IBMQ_Yorktown.loadCouplingMap(AvailableArchitecture::IBMQ_Yorktown);
    IBMQ_London.loadCouplingMap(test_architecture_dir + "ibmq_london.arch");
    IBMQ_London.loadProperties(test_calibration_dir + "ibmq_london.csv");
    IBMQ_Yorktown_mapper = std::make_unique<HeuristicMapper>(qc, IBMQ_Yorktown);
    IBMQ_London_mapper   = std::make_unique<HeuristicMapper>(qc, IBMQ_London);
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
  using namespace dd::literals;
  // construct circuit
  qc::QuantumComputation qc{4U};
  qc.x(1, 0_pc);
  qc.x(1, 2_pc);
  qc.x(1, 3_pc);

  // load architecture
  Architecture arch{};
  arch.loadCouplingMap(AvailableArchitecture::IBMQ_London);

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
  mapper.dumpResult(qasm, qc::OpenQASM);
  qcMapped.import(qasm, qc::OpenQASM);

  // check no measurements were added
  EXPECT_EQ(qcMapped.getNops(), 3U);
  EXPECT_NE(qcMapped.back()->getType(), qc::Measure);
}

INSTANTIATE_TEST_SUITE_P(
    Heuristic, HeuristicTest5Q,
    testing::Values("3_17_13", "ex-1_166", "ham3_102", "miller_11", "4gt11_84",
                    "4mod5-v0_20", "mod5d1_63"),
    [](const testing::TestParamInfo<HeuristicTest5Q::ParamType>& info) {
      std::string name = info.param;
      std::replace(name.begin(), name.end(), '-', '_');
      std::stringstream ss{};
      ss << name;
      return ss.str();
    });

TEST_P(HeuristicTest5Q, Identity) {
  Configuration settings{};
  settings.initialLayout = InitialLayout::Identity;
  IBMQ_Yorktown_mapper->map(settings);
  IBMQ_Yorktown_mapper->dumpResult(GetParam() + "_heuristic_qx4_identity.qasm");
  IBMQ_Yorktown_mapper->printResult(std::cout);

  IBMQ_London_mapper->map(settings);
  IBMQ_London_mapper->dumpResult(GetParam() +
                                 "_heuristic_london_identity.qasm");
  IBMQ_London_mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTest5Q, Static) {
  Configuration settings{};
  settings.initialLayout = InitialLayout::Static;
  IBMQ_Yorktown_mapper->map(settings);
  IBMQ_Yorktown_mapper->dumpResult(GetParam() + "_heuristic_qx4_static.qasm");
  IBMQ_Yorktown_mapper->printResult(std::cout);
  IBMQ_London_mapper->map(settings);
  IBMQ_London_mapper->dumpResult(GetParam() + "_heuristic_london_static.qasm");
  IBMQ_London_mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTest5Q, Dynamic) {
  Configuration settings{};
  settings.initialLayout = InitialLayout::Dynamic;
  IBMQ_Yorktown_mapper->map(settings);
  IBMQ_Yorktown_mapper->dumpResult(GetParam() + "_heuristic_qx4_dynamic.qasm");
  IBMQ_Yorktown_mapper->printResult(std::cout);
  IBMQ_London_mapper->map(settings);
  IBMQ_London_mapper->dumpResult(GetParam() + "_heuristic_london_dynamic.qasm");
  IBMQ_London_mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

class HeuristicTest16Q : public testing::TestWithParam<std::string> {
protected:
  std::string test_example_dir      = "../../examples/";
  std::string test_architecture_dir = "../../extern/architectures/";

  qc::QuantumComputation           qc{};
  Architecture                     IBM_QX5{};
  std::unique_ptr<HeuristicMapper> IBM_QX5_mapper;

  void SetUp() override {
    qc.import(test_example_dir + GetParam() + ".qasm");
    IBM_QX5.loadCouplingMap(AvailableArchitecture::IBM_QX5);
    IBM_QX5_mapper = std::make_unique<HeuristicMapper>(qc, IBM_QX5);
  }
};

INSTANTIATE_TEST_SUITE_P(
    Heuristic, HeuristicTest16Q,
    testing::Values("ising_model_10", "rd73_140", "cnt3-5_179", "qft_16"),
    [](const testing::TestParamInfo<HeuristicTest16Q::ParamType>& info) {
      std::string name = info.param;
      std::replace(name.begin(), name.end(), '-', '_');
      std::stringstream ss{};
      ss << name;
      return ss.str();
    });

TEST_P(HeuristicTest16Q, Dynamic) {
  Configuration settings{};
  settings.initialLayout = InitialLayout::Dynamic;
  IBM_QX5_mapper->map(settings);
  IBM_QX5_mapper->dumpResult(GetParam() + "_heuristic_qx5_dynamic.qasm");
  IBM_QX5_mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

class HeuristicTest20Q : public testing::TestWithParam<std::string> {
protected:
  std::string test_example_dir      = "../../examples/";
  std::string test_architecture_dir = "../../extern/architectures/";

  qc::QuantumComputation           qc{};
  Architecture                     arch{};
  std::unique_ptr<HeuristicMapper> tokyo_mapper;

  void SetUp() override {
    qc.import(test_example_dir + GetParam() + ".qasm");
    arch.loadCouplingMap(AvailableArchitecture::IBMQ_Tokyo);
    tokyo_mapper = std::make_unique<HeuristicMapper>(qc, arch);
  }
};

INSTANTIATE_TEST_SUITE_P(
    Heuristic, HeuristicTest20Q,
    testing::Values("ising_model_10", "rd73_140", "cnt3-5_179", "qft_16",
                    "z4_268"),
    [](const testing::TestParamInfo<HeuristicTest20Q::ParamType>& info) {
      std::string name = info.param;
      std::replace(name.begin(), name.end(), '-', '_');
      std::stringstream ss{};
      ss << name;
      return ss.str();
    });

TEST_P(HeuristicTest20Q, Dynamic) {
  Configuration settings{};
  settings.initialLayout = InitialLayout::Dynamic;
  tokyo_mapper->map(settings);
  tokyo_mapper->dumpResult(GetParam() + "_heuristic_tokyo_dynamic.qasm");
  tokyo_mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

class HeuristicTest20QTeleport
    : public testing::TestWithParam<
          std::tuple<unsigned long long, std::string>> {
protected:
  std::string test_example_dir      = "../../examples/";
  std::string test_architecture_dir = "../../extern/architectures/";

  qc::QuantumComputation           qc{};
  Architecture                     arch{};
  std::unique_ptr<HeuristicMapper> tokyo_mapper;

  void SetUp() override {
    qc.import(test_example_dir + std::get<1>(GetParam()) + ".qasm");
    arch.loadCouplingMap(AvailableArchitecture::IBMQ_Tokyo);
    tokyo_mapper = std::make_unique<HeuristicMapper>(qc, arch);
  }
};

INSTANTIATE_TEST_SUITE_P(
    HeuristicTeleport, HeuristicTest20QTeleport,
    testing::Combine(testing::Values(1, 2, 3, 1337, 1338, 3147),
                     testing::Values("ising_model_10", "rd73_140", "cnt3-5_179",
                                     "qft_16", "z4_268")),
    [](const testing::TestParamInfo<HeuristicTest20QTeleport::ParamType>&
           info) {
      std::string name = std::get<1>(info.param);
      std::replace(name.begin(), name.end(), '-', '_');
      std::stringstream ss{};
      ss << name << "_seed" << std::get<0>(info.param);
      return ss.str();
    });

TEST_P(HeuristicTest20QTeleport, Teleportation) {
  Configuration settings{};
  settings.initialLayout = InitialLayout::Dynamic;
  settings.teleportationQubits =
      std::min((arch.getNqubits() - qc.getNqubits()) & ~1u, 8u);
  settings.teleportationSeed = std::get<0>(GetParam());
  tokyo_mapper->map(settings);
  tokyo_mapper->dumpResult(std::get<1>(GetParam()) +
                           "_heuristic_tokyo_teleport.qasm");
  tokyo_mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}
