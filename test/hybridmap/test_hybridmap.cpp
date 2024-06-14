//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "Definitions.hpp"
#include "QuantumComputation.hpp"
#include "hybridmap/HybridNeutralAtomMapper.hpp"
#include "hybridmap/NeutralAtomArchitecture.hpp"
#include "hybridmap/NeutralAtomUtils.hpp"

#include <cstdint>
#include <filesystem>
#include <gtest/gtest.h>
#include <iostream>
#include <string>
#include <tuple>

class NeutralAtomArchitectureTest
    : public ::testing::TestWithParam<std::string> {
protected:
  std::string testArchitecturePath = "architectures/";

  void SetUp() override { testArchitecturePath += GetParam() + ".json"; }
};

TEST_P(NeutralAtomArchitectureTest, LoadArchitecures) {
  std::cout << "wd: " << std::filesystem::current_path() << '\n';
  auto arch = na::NeutralAtomArchitecture(testArchitecturePath);

  // Test get properties
  EXPECT_LE(arch.getNqubits(), arch.getNpositions());
  EXPECT_EQ(arch.getNpositions(), arch.getNrows() * arch.getNcolumns());
  // Test precomputed values
  auto c1 = arch.getCoordinate(0);
  auto c2 = arch.getCoordinate(1);
  EXPECT_GE(arch.getSwapDistance(c1, c2), 0);
  EXPECT_GE(arch.getNAodIntermediateLevels(), 1);
  // Test get parameters
  EXPECT_GE(arch.getGateTime("cz"), 0);
  EXPECT_GE(arch.getGateAverageFidelity("cz"), 0);
  // Test distance functions
  EXPECT_GE(arch.getEuclideanDistance(c1, c2), 0);
  // Test MoveVector functions
  auto mv = arch.getVector(0, 1);
  EXPECT_GE(arch.getVectorShuttlingTime(mv), 0);
}

INSTANTIATE_TEST_SUITE_P(NeutralAtomArchitectureTestSuite,
                         NeutralAtomArchitectureTest,
                         ::testing::Values("rubidium", "rubidium_hybrid",
                                           "rubidium_shuttling"));
class NeutralAtomMapperTestParams
    // parameters are architecture, circuit, gateWeight, shuttlingWeight,
    // lookAheadWeight, initialCoordinateMapping
    : public ::testing::TestWithParam<
          std::tuple<std::string, std::string, qc::fp, qc::fp, qc::fp,
                     na::InitialCoordinateMapping>> {
protected:
  std::string                  testArchitecturePath = "architectures/";
  std::string                  testQcPath           = "circuits/";
  qc::fp                       gateWeight           = 1;
  qc::fp                       shuttlingWeight      = 1;
  qc::fp                       lookAheadWeight      = 1;
  na::InitialCoordinateMapping initialCoordinateMapping =
      na::InitialCoordinateMapping::Trivial;
  // fixed
  qc::fp   decay               = 0.1;
  qc::fp   shuttlingTimeWeight = 0.1;
  uint32_t seed                = 42;

  void SetUp() override {
    auto params = GetParam();
    testArchitecturePath += std::get<0>(params) + ".json";
    testQcPath += std::get<1>(params) + ".qasm";
    gateWeight               = std::get<2>(params);
    shuttlingWeight          = std::get<3>(params);
    lookAheadWeight          = std::get<4>(params);
    initialCoordinateMapping = std::get<5>(params);
  }
};

TEST_P(NeutralAtomMapperTestParams, MapCircuitsIdentity) {
  auto arch = na::NeutralAtomArchitecture(testArchitecturePath);
  na::InitialMapping const initialMapping = na::InitialMapping::Identity;
  na::NeutralAtomMapper    mapper(arch);
  na::MapperParameters     mapperParameters;
  mapperParameters.initialMapping       = initialCoordinateMapping;
  mapperParameters.lookaheadWeightSwaps = lookAheadWeight;
  mapperParameters.lookaheadWeightMoves = lookAheadWeight;
  mapperParameters.decay                = decay;
  mapperParameters.shuttlingTimeWeight  = shuttlingTimeWeight;
  mapperParameters.gateWeight           = gateWeight;
  mapperParameters.shuttlingWeight      = shuttlingWeight;
  mapperParameters.seed                 = seed;
  mapperParameters.verbose              = true;
  mapper.setParameters(mapperParameters);

  qc::QuantumComputation qc(testQcPath);
  auto                   qcMapped = mapper.map(qc, initialMapping);
  mapper.convertToAod();

  auto scheduleResults = mapper.schedule(true, true);

  ASSERT_GT(scheduleResults.totalFidelities, 0);
  ASSERT_GT(scheduleResults.totalIdleTime, 0);
  ASSERT_GT(scheduleResults.totalExecutionTime, 0);
}

INSTANTIATE_TEST_SUITE_P(
    NeutralAtomMapperTestSuite, NeutralAtomMapperTestParams,
    ::testing::Combine(
        ::testing::Values("rubidium", "rubidium_hybrid", "rubidium_shuttling"),
        ::testing::Values("dj_nativegates_rigetti_qiskit_opt3_10", "modulo_2",
                          "multiply_2",
                          "qft_nativegates_rigetti_qiskit_opt3_10",
                          "random_nativegates_rigetti_qiskit_opt3_10"),
        ::testing::Values(1, 0.), ::testing::Values(1, 0.),
        ::testing::Values(0, 0.1),
        ::testing::Values(na::InitialCoordinateMapping::Trivial,
                          na::InitialCoordinateMapping::Random)));

class NeutralAtomMapperTest : public ::testing::Test {
protected:
  std::string testArchitecturePath = "architectures/rubidium_shuttling.json";
  const na::NeutralAtomArchitecture arch =
      na::NeutralAtomArchitecture(testArchitecturePath);
  na::InitialMapping const initialMapping = na::InitialMapping::Identity;
  na::MapperParameters     mapperParameters;
  na::NeutralAtomMapper    mapper;
  qc::QuantumComputation   qc;

  void SetUp() override {
    mapper                          = na::NeutralAtomMapper(arch);
    mapperParameters.initialMapping = na::InitialCoordinateMapping::Trivial;
    mapperParameters.lookaheadWeightSwaps = 0.1;
    mapperParameters.lookaheadWeightMoves = 0.1;
    mapperParameters.decay                = 0;
    mapperParameters.shuttlingTimeWeight  = 0.1;
    mapperParameters.gateWeight           = 1;
    mapperParameters.shuttlingWeight      = 0;
    mapperParameters.seed                 = 43;
    mapperParameters.verbose              = true;
    mapper.setParameters(mapperParameters);
    qc = qc::QuantumComputation(
        "circuits/dj_nativegates_rigetti_qiskit_opt3_10.qasm");
  }
};

TEST_F(NeutralAtomMapperTest, Output) {
  auto qcMapped = mapper.map(qc, initialMapping);

  qcMapped.dumpOpenQASM(std::cout, false);

  auto qcAodMapped = mapper.convertToAod();
  qcAodMapped.dumpOpenQASM(std::cout, false);

  auto scheduleResults = mapper.schedule(true, true);
  std::cout << scheduleResults.toCsv();

  ASSERT_GT(scheduleResults.totalFidelities, 0);
}
