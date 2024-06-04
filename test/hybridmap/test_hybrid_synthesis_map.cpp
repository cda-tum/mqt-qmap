//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "QuantumComputation.hpp"
#include "hybridmap/HybridNeutralAtomMapper.hpp"
#include "hybridmap/HybridSynthesisMapper.hpp"
#include "hybridmap/NeutralAtomArchitecture.hpp"

#include <cstddef>
#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace na {
class TestHybridSynthesisMapper : public ::testing::TestWithParam<std::string> {
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

TEST_P(TestHybridSynthesisMapper, AdjaencyMatrix) {
  auto arch   = NeutralAtomArchitecture(testArchitecturePath);
  auto mapper = HybridSynthesisMapper(arch);
  mapper.initMapping(3);
  auto adjMatrix = mapper.getCircuitAdjacencyMatrix();
  EXPECT_EQ(adjMatrix.size(), 3);
  EXPECT_TRUE(adjMatrix(0, 2) == 0 || adjMatrix(0, 2) == 1);
}

INSTANTIATE_TEST_SUITE_P(HybridSynthesisMapperTestSuite,
                         TestHybridSynthesisMapper,
                         ::testing::Values("rubidium", "rubidium_hybrid",
                                           "rubidium_shuttling"));

} // namespace na
