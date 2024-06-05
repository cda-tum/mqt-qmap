//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "QuantumComputation.hpp"
#include "hybridmap/HybridSynthesisMapper.hpp"
#include "hybridmap/NeutralAtomArchitecture.hpp"

#include <cstddef>
#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace na {
class TestParametrizedHybridSynthesisMapper
    : public ::testing::TestWithParam<std::string> {
protected:
  std::string                         testArchitecturePath = "architectures/";
  std::vector<qc::QuantumComputation> circuits;

  void SetUp() override {
    testArchitecturePath += GetParam() + ".json";
    qc::QuantumComputation qc1(3);
    qc1.x(0);
    qc1.cx(0, 1);
    qc1.cx(1, 2);
    circuits.push_back(qc1);

    qc::QuantumComputation qc2(3);
    qc2.move(0, 2);
    qc2.x(0);
    circuits.push_back(qc2);
  }

  // Test the HybridSynthesisMapper class
};

TEST_P(TestParametrizedHybridSynthesisMapper, AdjaencyMatrix) {
  auto arch   = NeutralAtomArchitecture(testArchitecturePath);
  auto mapper = HybridSynthesisMapper(arch);
  mapper.initMapping(3);
  auto adjMatrix = mapper.getCircuitAdjacencyMatrix();
  EXPECT_EQ(adjMatrix.size(), 3);
  EXPECT_TRUE(adjMatrix(0, 2) == 0 || adjMatrix(0, 2) == 1);
}

TEST_P(TestParametrizedHybridSynthesisMapper, EvaluateSynthesisStep) {
  auto arch   = NeutralAtomArchitecture(testArchitecturePath);
  auto mapper = HybridSynthesisMapper(arch);
  mapper.initMapping(3);
  auto best = mapper.evaluateSynthesisSteps(circuits, false);
  EXPECT_GE(best, 0);
}

INSTANTIATE_TEST_SUITE_P(HybridSynthesisMapperTestSuite,
                         TestParametrizedHybridSynthesisMapper,
                         ::testing::Values("rubidium", "rubidium_hybrid",
                                           "rubidium_shuttling"));

class TestHybridSynthesisMapper : public ::testing::Test {
protected:
  NeutralAtomArchitecture arch =
      NeutralAtomArchitecture("architectures/rubidium.json");
  HybridSynthesisMapper  mapper = HybridSynthesisMapper(arch);
  qc::QuantumComputation qc;

  void SetUp() override {
    qc = qc::QuantumComputation(3);
    qc.x(0);
    qc.cx(0, 1);
    qc.cx(1, 2);

    mapper.initMapping(3);
  }
};

TEST_F(TestHybridSynthesisMapper, DirectlyMap) {
  mapper.appendWithoutMapping(qc);
  auto synthesizedQc = mapper.getSynthesizedQc();
  EXPECT_EQ(synthesizedQc.getNqubits(), 3);
  EXPECT_EQ(synthesizedQc.getNops(), 3);
}

TEST_F(TestHybridSynthesisMapper, completelyRemap) {
  mapper.appendWithoutMapping(qc);
  mapper.appendWithoutMapping(qc);
  auto mappedQc = mapper.getMappedQc();
  EXPECT_EQ(mappedQc.getNqubits(), arch.getNpositions());
  EXPECT_GE(mappedQc.getNops(), 3);
  auto mappedQcRemapped = mapper.getMappedQc();
  EXPECT_EQ(mappedQcRemapped.getNqubits(), arch.getNpositions());
  EXPECT_GE(mappedQcRemapped.getNops(), 3);
}

TEST_F(TestHybridSynthesisMapper, MapAppend) {
  mapper.appendWithMapping(qc);
  auto synthesizedQc = mapper.getSynthesizedQc();
  EXPECT_EQ(synthesizedQc.getNqubits(), 3);
  EXPECT_GE(synthesizedQc.getNops(), 3);
}

} // namespace na
