//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "Architecture.hpp"

#include "gtest/gtest.h"

class TestArchitecture : public testing::TestWithParam<std::string> {
protected:
  std::string test_architecture_dir = "./architectures/";
  std::string test_calibration_dir  = "./calibration/";
};

INSTANTIATE_TEST_SUITE_P(Architecture, TestArchitecture,
                         testing::Values("ibm_qx4.arch", "ibmq_casablanca.arch",
                                         "ibmq_london.arch",
                                         "ibmq_london.csv"));

TEST_P(TestArchitecture, QubitMap) {
  auto&             arch_name = GetParam();
  Architecture      arch{};
  std::stringstream ss{};
  if (arch_name.find(".arch") != std::string::npos) {
    ss << test_architecture_dir << arch_name;
    arch.loadCouplingMap(ss.str());
  } else {
    ss << test_calibration_dir << arch_name;
    arch.loadProperties(ss.str());
  }

  EXPECT_EQ(Architecture::getQubitList(arch.getCouplingMap()).size(),
            arch.getNqubits());
}
TEST_P(TestArchitecture, GetAllConnectedSubsets) {
  auto&             arch_name = GetParam();
  Architecture      arch{};
  std::stringstream ss{};
  if (arch_name.find(".arch") != std::string::npos) {
    ss << test_architecture_dir << arch_name;
    arch.loadCouplingMap(ss.str());
  } else {
    ss << test_calibration_dir << arch_name;
    arch.loadProperties(ss.str());
  }

  EXPECT_EQ(arch.getAllConnectedSubsets(arch.getNqubits()).size(), 1);
  EXPECT_EQ(arch.getAllConnectedSubsets(1).size(), arch.getNqubits());
}
TEST_P(TestArchitecture, GetHighestFidelity) {
  auto&             arch_name = GetParam();
  Architecture      arch{};
  std::stringstream ss{};
  if (arch_name.find(".arch") != std::string::npos) {
    ss << test_architecture_dir << arch_name;
    arch.loadCouplingMap(ss.str());
  } else {
    ss << test_calibration_dir << arch_name;
    arch.loadProperties(ss.str());
  }
  CouplingMap cm{};

  arch.getHighestFidelityCouplingMap(arch.getNqubits(), cm);
  EXPECT_EQ(cm, arch.getCouplingMap());

  arch.getHighestFidelityCouplingMap(1U, cm);
  EXPECT_TRUE(cm.empty());
}
TEST_P(TestArchitecture, ReducedMaps) {
  auto&             arch_name = GetParam();
  Architecture      arch{};
  std::stringstream ss{};
  if (arch_name.find(".arch") != std::string::npos) {
    ss << test_architecture_dir << arch_name;
    arch.loadCouplingMap(ss.str());
  } else {
    ss << test_calibration_dir << arch_name;
    arch.loadProperties(ss.str());
  }

  std::vector<CouplingMap> cms;

  arch.getReducedCouplingMaps(1, cms);

  EXPECT_EQ(cms.size(), arch.getNqubits());
}

TEST(TestArchitecture, ConnectedTest) {
  Architecture architecture{};
  CouplingMap  cm{};

  cm.emplace(std::make_pair(0, 1));
  cm.emplace(std::make_pair(1, 2));
  cm.emplace(std::make_pair(2, 3));
  cm.emplace(std::make_pair(3, 4));
  cm.emplace(std::make_pair(4, 0));

  Architecture::printCouplingMap(cm, std::cout);

  architecture.loadCouplingMap(5, cm);

  std::vector<CouplingMap> cms;

  architecture.getReducedCouplingMaps(2, cms);

  EXPECT_EQ(cms.size(), 5);

  architecture.getReducedCouplingMaps(4, cms);

  EXPECT_EQ(cms.size(), 5);
}

TEST(TestArchitecture, FidelityTest) {
  Architecture architecture{};
  CouplingMap  cm{};

  auto props = Architecture::Properties();
  props.setNqubits(4);
  props.setSingleQubitErrorRate(0, "x", 0.9);
  props.setSingleQubitErrorRate(1, "x", 0.9);
  props.setSingleQubitErrorRate(2, "x", 0.9);
  props.setSingleQubitErrorRate(3, "x", 0.9);

  props.setTwoQubitErrorRate(0, 1, 0.8);
  props.setTwoQubitErrorRate(1, 0, 0.8);
  props.setTwoQubitErrorRate(1, 2, 0.7);
  props.setTwoQubitErrorRate(2, 1, 0.7);
  props.setTwoQubitErrorRate(2, 3, 0.6);
  props.setTwoQubitErrorRate(3, 2, 0.6);

  architecture.loadProperties(props);
  architecture.getHighestFidelityCouplingMap(2, cm);

  const std::vector<std::uint16_t> highestFidelity{2, 3};
  auto                             qubitList = Architecture::getQubitList(cm);

  EXPECT_EQ(qubitList, highestFidelity);
}

TEST(TestArchitecture, FullyConnectedTest) {
  const auto cm = getFullyConnectedMap(3);

  ASSERT_TRUE(cm.size() == 3 * 2);
}

TEST(TestArchitecture, MinimumNumberOfSwapsError) {
  Architecture               architecture{};
  std::vector<std::uint16_t> permutation{1, 1, 2, 3, 4};
  printPi(permutation);
  std::vector<Edge> swaps{};
  EXPECT_THROW(architecture.minimumNumberOfSwaps(permutation, swaps),
               std::runtime_error);
}
