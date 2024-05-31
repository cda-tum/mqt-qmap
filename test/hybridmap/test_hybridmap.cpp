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
#include <fstream>
#include <gtest/gtest.h>
#include <iostream>
#include <string>
#include <tuple>

class NeutralAtomArchitectureTest
    : public ::testing::TestWithParam<std::string> {
protected:
  std::string testArchitecturePath = "hybridmap/architectures/";

  void SetUp() override { testArchitecturePath += GetParam() + ".json"; }
};

TEST_P(NeutralAtomArchitectureTest, LoadArchitecures) {
  std::cout << "wd" << std::filesystem::current_path() << '\n';
  auto arch = qc::NeutralAtomArchitecture(testArchitecturePath);

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
  EXPECT_GE(arch.getEuclidianDistance(c1, c2), 0);
  // Test MoveVector functions
  auto mv = arch.getVector(0, 1);
  EXPECT_GE(arch.getVectorShuttlingTime(mv), 0);
}

INSTANTIATE_TEST_SUITE_P(NeutralAtomArchitectureTestSuite,
                         NeutralAtomArchitectureTest,
                         ::testing::Values("rubidium", "rubidium_hybrid",
                                           "rubidium_shuttling"));
class NeutralAtomMapperTest
    // parameters are architecture, circuit, gateWeight, shuttlingWeight,
    // lookAheadWeight, initialCoordinateMapping
    : public ::testing::TestWithParam<
          std::tuple<std::string, std::string, qc::fp, qc::fp, qc::fp,
                     qc::InitialCoordinateMapping>> {
protected:
  std::string testArchitecturePath = "hybridmap/architectures/";
  std::string testQcPath           = "hybridmap/circuits/";
  qc::fp      gateWeight           = 1;
  qc::fp      shuttlingWeight      = 1;
  qc::fp      lookAheadWeight      = 1;
  qc::InitialCoordinateMapping initialCoordinateMapping =
      qc::InitialCoordinateMapping::Trivial;
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

TEST_P(NeutralAtomMapperTest, MapCircuitsIdentity) {
  auto arch = qc::NeutralAtomArchitecture(testArchitecturePath);
  qc::InitialMapping const initialMapping = qc::InitialMapping::Identity;
  qc::NeutralAtomMapper    mapper(arch, initialCoordinateMapping);
  qc::MapperParameters     mapperParameters;
  mapperParameters.lookaheadWeightSwaps = lookAheadWeight;
  mapperParameters.lookaheadWeightMoves = lookAheadWeight;
  mapperParameters.decay                = decay;
  mapperParameters.shuttlingTimeWeight  = shuttlingTimeWeight;
  mapperParameters.gateWeight           = gateWeight;
  mapperParameters.shuttlingWeight      = shuttlingWeight;
  mapperParameters.seed                 = seed;
  mapper.setParameters(mapperParameters);

  qc::QuantumComputation qc(testQcPath);
  auto                   qcMapped    = mapper.map(qc, initialMapping, true);
  auto                   qcAodMapped = mapper.convertToAod(qcMapped);

  auto scheduleResults = mapper.schedule(true, true);

  ASSERT_GT(scheduleResults.totalFidelities, 0);
  ASSERT_GT(scheduleResults.totalIdleTime, 0);
  ASSERT_GT(scheduleResults.totalExecutionTime, 0);
}

INSTANTIATE_TEST_SUITE_P(
    NeutralAtomMapperTestSuite, NeutralAtomMapperTest,
    ::testing::Combine(
        ::testing::Values("rubidium", "rubidium_hybrid", "rubidium_shuttling"),
        ::testing::Values("dj_nativegates_rigetti_qiskit_opt3_10", "modulo_2",
                          "multiply_2",
                          "qft_nativegates_rigetti_qiskit_opt3_10"),
        ::testing::Values(1, 0.), ::testing::Values(1, 0.),
        ::testing::Values(0, 0.1),
        ::testing::Values(qc::InitialCoordinateMapping::Trivial,
                          qc::InitialCoordinateMapping::Random)));

TEST(NeutralAtomMapperTest, Output) {
  auto arch = qc::NeutralAtomArchitecture(
      "hybridmap/architectures/rubidium_shuttling.json");
  qc::InitialMapping const initialMapping = qc::InitialMapping::Identity;
  qc::InitialCoordinateMapping const initialCoordinateMapping =
      qc::InitialCoordinateMapping::Trivial;
  qc::NeutralAtomMapper mapper(arch, initialCoordinateMapping);
  qc::MapperParameters  mapperParameters;
  mapperParameters.lookaheadWeightSwaps = 0.1;
  mapperParameters.lookaheadWeightMoves = 0.1;
  mapperParameters.decay                = 0;
  mapperParameters.shuttlingTimeWeight  = 0.1;
  mapperParameters.gateWeight           = 1;
  mapperParameters.shuttlingWeight      = 0;
  mapperParameters.seed                 = 42;
  mapper.setParameters(mapperParameters);

  qc::QuantumComputation qc(
      "hybridmap/circuits/dj_nativegates_rigetti_qiskit_opt3_10.qasm");
  auto qcMapped = mapper.map(qc, initialMapping, true);

  std::ofstream dummyFs;
  qcMapped.dumpOpenQASM(dummyFs, false);

  auto qcAodMapped = mapper.convertToAod(qcMapped);
  qcAodMapped.dumpOpenQASM(dummyFs, false);

  auto scheduleResults = mapper.schedule(true, true);
  dummyFs << scheduleResults.toCsv();

  ASSERT_GT(scheduleResults.totalFidelities, 0);
}
