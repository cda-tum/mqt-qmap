//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "Architecture.hpp"

#include "gtest/gtest.h"

class TestArchitecture : public testing::TestWithParam<std::string> {
protected:
  std::string testArchitectureDir = "../extern/architectures/";
  std::string testCalibrationDir  = "../extern/calibration/";
};

INSTANTIATE_TEST_SUITE_P(Architecture, TestArchitecture,
                         testing::Values("ibm_qx4.arch", "ibmq_casablanca.arch",
                                         "ibmq_london.arch",
                                         "ibmq_london.csv"));

TEST_P(TestArchitecture, QubitMap) {
  const auto&       archName = GetParam();
  Architecture      arch{};
  std::stringstream ss{};
  if (archName.find(".arch") != std::string::npos) {
    ss << testArchitectureDir << archName;
    arch.loadCouplingMap(ss.str());
  } else {
    ss << testCalibrationDir << archName;
    arch.loadProperties(ss.str());
  }

  EXPECT_EQ(Architecture::getQubitList(arch.getCouplingMap()).size(),
            arch.getNqubits());
}
TEST_P(TestArchitecture, GetAllConnectedSubsets) {
  const auto&       archName = GetParam();
  Architecture      arch{};
  std::stringstream ss{};
  if (archName.find(".arch") != std::string::npos) {
    ss << testArchitectureDir << archName;
    arch.loadCouplingMap(ss.str());
  } else {
    ss << testCalibrationDir << archName;
    arch.loadProperties(ss.str());
  }

  EXPECT_EQ(arch.getAllConnectedSubsets(arch.getNqubits()).size(), 1);
  EXPECT_EQ(arch.getAllConnectedSubsets(1).size(), arch.getNqubits());
}
TEST_P(TestArchitecture, GetHighestFidelity) {
  const auto&       archName = GetParam();
  Architecture      arch{};
  std::stringstream ss{};
  if (archName.find(".arch") != std::string::npos) {
    ss << testArchitectureDir << archName;
    arch.loadCouplingMap(ss.str());
  } else {
    ss << testCalibrationDir << archName;
    arch.loadProperties(ss.str());
  }
  CouplingMap cm{};

  arch.getHighestFidelityCouplingMap(arch.getNqubits(), cm);
  EXPECT_EQ(cm, arch.getCouplingMap());

  arch.getHighestFidelityCouplingMap(1U, cm);
  EXPECT_TRUE(cm.empty());
}
TEST_P(TestArchitecture, ReducedMaps) {
  const auto&       archName = GetParam();
  Architecture      arch{};
  std::stringstream ss{};
  if (archName.find(".arch") != std::string::npos) {
    ss << testArchitectureDir << archName;
    arch.loadCouplingMap(ss.str());
  } else {
    ss << testCalibrationDir << archName;
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

  ASSERT_TRUE(cm.size() == static_cast<std::size_t>(2 * 3));
}

TEST(TestArchitecture, MinimumNumberOfSwapsError) {
  Architecture               architecture{};
  std::vector<std::uint16_t> permutation{1, 1, 2, 3, 4};
  printPi(permutation);
  std::vector<Edge> swaps{};
  EXPECT_THROW(architecture.minimumNumberOfSwaps(permutation, swaps),
               std::runtime_error);
}

TEST(TestArchitecture, TestCouplingLimitRing) {
  Architecture      architecture{};
  const CouplingMap cm = {{0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 3},
                          {3, 2}, {3, 4}, {4, 3}, {4, 0}, {0, 4}};
  architecture.loadCouplingMap(5, cm);
  EXPECT_EQ(architecture.getCouplingLimit(), 2);
}

TEST(TestArchitecture, FidelityDistanceBidirectionalTest) {
  Architecture      architecture{};
  const CouplingMap cm = {{0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 3}, {3, 2},
                          {1, 4}, {4, 1}, {2, 5}, {5, 2}, {5, 6}, {6, 5}};
  architecture.loadCouplingMap(7, cm);

  auto props = Architecture::Properties();
  props.setSingleQubitErrorRate(0, "x", 0.03);
  props.setSingleQubitErrorRate(1, "x", 0.03);
  props.setSingleQubitErrorRate(2, "x", 0.02);
  props.setSingleQubitErrorRate(3, "x", 0.03);
  props.setSingleQubitErrorRate(4, "x", 0.03);
  props.setSingleQubitErrorRate(5, "x", 0.02);
  props.setSingleQubitErrorRate(6, "x", 0.03);

  props.setTwoQubitErrorRate(0, 1, 0.9);
  props.setTwoQubitErrorRate(1, 0, 0.9);
  props.setTwoQubitErrorRate(1, 2, 0.5);
  props.setTwoQubitErrorRate(2, 1, 0.5);
  props.setTwoQubitErrorRate(2, 3, 0.1);
  props.setTwoQubitErrorRate(3, 2, 0.1);
  props.setTwoQubitErrorRate(1, 4, 0.1);
  props.setTwoQubitErrorRate(4, 1, 0.1);
  props.setTwoQubitErrorRate(2, 5, 0.5);
  props.setTwoQubitErrorRate(5, 2, 0.5);
  props.setTwoQubitErrorRate(5, 6, 0.9);
  props.setTwoQubitErrorRate(6, 5, 0.9);

  architecture.loadProperties(props);

  const Matrix<double> fidDistance = architecture.getFidelityDistanceTable();

  EXPECT_EQ(fidDistance.size(), 7);
  EXPECT_EQ(fidDistance[0].size(), 7);
  EXPECT_NEAR(fidDistance[0][1], -3 * std::log2(1 - 0.9), 1e-6);
  EXPECT_NEAR(fidDistance[0][2], -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.5)),
              1e-6);
  EXPECT_NEAR(
      fidDistance[0][3],
      -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.5) + std::log2(1 - 0.1)),
      1e-6);
  EXPECT_NEAR(fidDistance[0][4], -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.1)),
              1e-6);
  EXPECT_NEAR(
      fidDistance[0][5],
      -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.5) + std::log2(1 - 0.5)),
      1e-6);
  EXPECT_NEAR(fidDistance[0][6],
              -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.5) +
                    std::log2(1 - 0.5) + std::log2(1 - 0.9)),
              1e-6);
  EXPECT_EQ(fidDistance[1].size(), 7);
  EXPECT_NEAR(fidDistance[1][0], -3 * std::log2(1 - 0.9), 1e-6);
  EXPECT_NEAR(fidDistance[1][2], -3 * std::log2(1 - 0.5), 1e-6);
  EXPECT_NEAR(fidDistance[1][3], -3 * (std::log2(1 - 0.5) + std::log2(1 - 0.1)),
              1e-6);
  EXPECT_NEAR(fidDistance[1][4], -3 * std::log2(1 - 0.1), 1e-6);
  EXPECT_NEAR(fidDistance[1][5], -3 * (std::log2(1 - 0.5) + std::log2(1 - 0.5)),
              1e-6);
  EXPECT_NEAR(
      fidDistance[1][6],
      -3 * (std::log2(1 - 0.5) + std::log2(1 - 0.5) + std::log2(1 - 0.9)),
      1e-6);
  EXPECT_EQ(fidDistance[2].size(), 7);
  EXPECT_NEAR(fidDistance[2][0], -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.5)),
              1e-6);
  EXPECT_NEAR(fidDistance[2][1], -3 * std::log2(1 - 0.5), 1e-6);
  EXPECT_NEAR(fidDistance[2][3], -3 * std::log2(1 - 0.1), 1e-6);
  EXPECT_NEAR(fidDistance[2][4], -3 * (std::log2(1 - 0.5) + std::log2(1 - 0.1)),
              1e-6);
  EXPECT_NEAR(fidDistance[2][5], -3 * std::log2(1 - 0.5), 1e-6);
  EXPECT_NEAR(fidDistance[2][6], -3 * (std::log2(1 - 0.5) + std::log2(1 - 0.9)),
              1e-6);
  EXPECT_EQ(fidDistance[3].size(), 7);
  EXPECT_NEAR(
      fidDistance[3][0],
      -3 * (std::log2(1 - 0.1) + std::log2(1 - 0.5) + std::log2(1 - 0.9)),
      1e-6);
  EXPECT_NEAR(fidDistance[3][1], -3 * (std::log2(1 - 0.5) + std::log2(1 - 0.1)),
              1e-6);
  EXPECT_NEAR(fidDistance[3][2], -3 * std::log2(1 - 0.1), 1e-6);
  EXPECT_NEAR(
      fidDistance[3][4],
      -3 * (std::log2(1 - 0.1) + std::log2(1 - 0.5) + std::log2(1 - 0.1)),
      1e-6);
  EXPECT_NEAR(fidDistance[3][5], -3 * (std::log2(1 - 0.1) + std::log2(1 - 0.5)),
              1e-6);
  EXPECT_NEAR(
      fidDistance[3][6],
      -3 * (std::log2(1 - 0.1) + std::log2(1 - 0.5) + std::log2(1 - 0.9)),
      1e-6);
  EXPECT_EQ(fidDistance[4].size(), 7);
  EXPECT_NEAR(fidDistance[4][0], -3 * (std::log2(1 - 0.1) + std::log2(1 - 0.9)),
              1e-6);
  EXPECT_NEAR(fidDistance[4][1], -3 * std::log2(1 - 0.1), 1e-6);
  EXPECT_NEAR(fidDistance[4][2], -3 * (std::log2(1 - 0.1) + std::log2(1 - 0.5)),
              1e-6);
  EXPECT_NEAR(
      fidDistance[4][3],
      -3 * (std::log2(1 - 0.1) + std::log2(1 - 0.5) + std::log2(1 - 0.1)),
      1e-6);
  EXPECT_NEAR(
      fidDistance[4][5],
      -3 * (std::log2(1 - 0.1) + std::log2(1 - 0.5) + std::log2(1 - 0.5)),
      1e-6);
  EXPECT_NEAR(fidDistance[4][6],
              -3 * (std::log2(1 - 0.1) + std::log2(1 - 0.5) +
                    std::log2(1 - 0.5) + std::log2(1 - 0.9)),
              1e-6);
  EXPECT_EQ(fidDistance[5].size(), 7);
  EXPECT_NEAR(
      fidDistance[5][0],
      -3 * (std::log2(1 - 0.5) + std::log2(1 - 0.5) + std::log2(1 - 0.9)),
      1e-6);
  EXPECT_NEAR(fidDistance[5][1], -3 * (std::log2(1 - 0.5) + std::log2(1 - 0.5)),
              1e-6);
  EXPECT_NEAR(fidDistance[5][2], -3 * std::log2(1 - 0.5), 1e-6);
  EXPECT_NEAR(fidDistance[5][3], -3 * (std::log2(1 - 0.5) + std::log2(1 - 0.1)),
              1e-6);
  EXPECT_NEAR(
      fidDistance[5][4],
      -3 * (std::log2(1 - 0.5) + std::log2(1 - 0.5) + std::log2(1 - 0.1)),
      1e-6);
  EXPECT_NEAR(fidDistance[5][6], -3 * std::log2(1 - 0.9), 1e-6);
  EXPECT_EQ(fidDistance[6].size(), 7);
  EXPECT_NEAR(fidDistance[6][0],
              -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.5) +
                    std::log2(1 - 0.5) + std::log2(1 - 0.9)),
              1e-6);
  EXPECT_NEAR(
      fidDistance[6][1],
      -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.5) + std::log2(1 - 0.5)),
      1e-6);
  EXPECT_NEAR(fidDistance[6][2], -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.5)),
              1e-6);
  EXPECT_NEAR(
      fidDistance[6][3],
      -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.5) + std::log2(1 - 0.1)),
      1e-6);
  EXPECT_NEAR(fidDistance[6][4],
              -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.5) +
                    std::log2(1 - 0.5) + std::log2(1 - 0.1)),
              1e-6);
  EXPECT_NEAR(fidDistance[6][5], -3 * std::log2(1 - 0.9), 1e-6);
}

TEST(TestArchitecture, FidelityDistanceSemiBidirectionalTest) {
  Architecture      architecture{};
  const CouplingMap cm = {{0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 3},
                          {3, 2}, {1, 4}, {2, 5}, {5, 2}, {6, 5}};
  architecture.loadCouplingMap(7, cm);

  auto props = Architecture::Properties();
  props.setSingleQubitErrorRate(0, "x", 0.03);
  props.setSingleQubitErrorRate(1, "x", 0.03);
  props.setSingleQubitErrorRate(2, "x", 0.02);
  props.setSingleQubitErrorRate(3, "x", 0.03);
  props.setSingleQubitErrorRate(4, "x", 0.03);
  props.setSingleQubitErrorRate(5, "x", 0.02);
  props.setSingleQubitErrorRate(6, "x", 0.03);

  props.setTwoQubitErrorRate(0, 1, 0.9);
  props.setTwoQubitErrorRate(1, 0, 0.9);
  props.setTwoQubitErrorRate(1, 2, 0.5);
  props.setTwoQubitErrorRate(2, 1, 0.5);
  props.setTwoQubitErrorRate(2, 3, 0.1);
  props.setTwoQubitErrorRate(3, 2, 0.1);
  props.setTwoQubitErrorRate(1, 4, 0.1);
  props.setTwoQubitErrorRate(2, 5, 0.5);
  props.setTwoQubitErrorRate(5, 2, 0.5);
  props.setTwoQubitErrorRate(6, 5, 0.9);

  architecture.loadProperties(props);

  const Matrix<double> fidDistance = architecture.getFidelityDistanceTable();

  EXPECT_EQ(fidDistance.size(), 7);
  EXPECT_EQ(fidDistance[0].size(), 7);
  EXPECT_NEAR(fidDistance[0][1], -3 * std::log2(1 - 0.9), 1e-6);
  EXPECT_NEAR(fidDistance[0][2], -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.5)),
              1e-6);
  EXPECT_NEAR(
      fidDistance[0][3],
      -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.5) + std::log2(1 - 0.1)),
      1e-6);
  EXPECT_NEAR(fidDistance[0][4],
              -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.1)) -
                  2 * (std::log2(1 - 0.03) + std::log2(1 - 0.03)),
              1e-6);
  EXPECT_NEAR(
      fidDistance[0][5],
      -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.5) + std::log2(1 - 0.5)),
      1e-6);
  EXPECT_NEAR(fidDistance[0][6],
              -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.5) +
                    std::log2(1 - 0.5) + std::log2(1 - 0.9)) -
                  2 * (std::log2(1 - 0.02) + std::log2(1 - 0.03)),
              1e-6);
  EXPECT_EQ(fidDistance[1].size(), 7);
  EXPECT_NEAR(fidDistance[1][0], -3 * std::log2(1 - 0.9), 1e-6);
  EXPECT_NEAR(fidDistance[1][2], -3 * std::log2(1 - 0.5), 1e-6);
  EXPECT_NEAR(fidDistance[1][3], -3 * (std::log2(1 - 0.5) + std::log2(1 - 0.1)),
              1e-6);
  EXPECT_NEAR(fidDistance[1][4],
              -3 * std::log2(1 - 0.1) -
                  2 * (std::log2(1 - 0.03) + std::log2(1 - 0.03)),
              1e-6);
  EXPECT_NEAR(fidDistance[1][5], -3 * (std::log2(1 - 0.5) + std::log2(1 - 0.5)),
              1e-6);
  EXPECT_NEAR(
      fidDistance[1][6],
      -3 * (std::log2(1 - 0.5) + std::log2(1 - 0.5) + std::log2(1 - 0.9)) -
          2 * (std::log2(1 - 0.02) + std::log2(1 - 0.03)),
      1e-6);
  EXPECT_EQ(fidDistance[2].size(), 7);
  EXPECT_NEAR(fidDistance[2][0], -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.5)),
              1e-6);
  EXPECT_NEAR(fidDistance[2][1], -3 * std::log2(1 - 0.5), 1e-6);
  EXPECT_NEAR(fidDistance[2][3], -3 * std::log2(1 - 0.1), 1e-6);
  EXPECT_NEAR(fidDistance[2][4],
              -3 * (std::log2(1 - 0.5) + std::log2(1 - 0.1)) -
                  2 * (std::log2(1 - 0.03) + std::log2(1 - 0.03)),
              1e-6);
  EXPECT_NEAR(fidDistance[2][5], -3 * std::log2(1 - 0.5), 1e-6);
  EXPECT_NEAR(fidDistance[2][6],
              -3 * (std::log2(1 - 0.5) + std::log2(1 - 0.9)) -
                  2 * (std::log2(1 - 0.02) + std::log2(1 - 0.03)),
              1e-6);
  EXPECT_EQ(fidDistance[3].size(), 7);
  EXPECT_NEAR(
      fidDistance[3][0],
      -3 * (std::log2(1 - 0.1) + std::log2(1 - 0.5) + std::log2(1 - 0.9)),
      1e-6);
  EXPECT_NEAR(fidDistance[3][1], -3 * (std::log2(1 - 0.5) + std::log2(1 - 0.1)),
              1e-6);
  EXPECT_NEAR(fidDistance[3][2], -3 * std::log2(1 - 0.1), 1e-6);
  EXPECT_NEAR(
      fidDistance[3][4],
      -3 * (std::log2(1 - 0.1) + std::log2(1 - 0.5) + std::log2(1 - 0.1)) -
          2 * (std::log2(1 - 0.03) + std::log2(1 - 0.03)),
      1e-6);
  EXPECT_NEAR(fidDistance[3][5], -3 * (std::log2(1 - 0.1) + std::log2(1 - 0.5)),
              1e-6);
  EXPECT_NEAR(
      fidDistance[3][6],
      -3 * (std::log2(1 - 0.1) + std::log2(1 - 0.5) + std::log2(1 - 0.9)) -
          2 * (std::log2(1 - 0.02) + std::log2(1 - 0.03)),
      1e-6);
  EXPECT_EQ(fidDistance[4].size(), 7);
  EXPECT_NEAR(fidDistance[4][0],
              -3 * (std::log2(1 - 0.1) + std::log2(1 - 0.9)) -
                  2 * (std::log2(1 - 0.03) + std::log2(1 - 0.03)),
              1e-6);
  EXPECT_NEAR(fidDistance[4][1],
              -3 * std::log2(1 - 0.1) -
                  2 * (std::log2(1 - 0.03) + std::log2(1 - 0.03)),
              1e-6);
  EXPECT_NEAR(fidDistance[4][2],
              -3 * (std::log2(1 - 0.1) + std::log2(1 - 0.5)) -
                  2 * (std::log2(1 - 0.03) + std::log2(1 - 0.03)),
              1e-6);
  EXPECT_NEAR(
      fidDistance[4][3],
      -3 * (std::log2(1 - 0.1) + std::log2(1 - 0.5) + std::log2(1 - 0.1)) -
          2 * (std::log2(1 - 0.03) + std::log2(1 - 0.03)),
      1e-6);
  EXPECT_NEAR(
      fidDistance[4][5],
      -3 * (std::log2(1 - 0.1) + std::log2(1 - 0.5) + std::log2(1 - 0.5)) -
          2 * (std::log2(1 - 0.03) + std::log2(1 - 0.03)),
      1e-6);
  EXPECT_NEAR(fidDistance[4][6],
              -3 * (std::log2(1 - 0.1) + std::log2(1 - 0.5) +
                    std::log2(1 - 0.5) + std::log2(1 - 0.9)) -
                  2 * (std::log2(1 - 0.03) + std::log2(1 - 0.03) +
                       std::log2(1 - 0.02) + std::log2(1 - 0.03)),
              1e-6);
  EXPECT_EQ(fidDistance[5].size(), 7);
  EXPECT_NEAR(
      fidDistance[5][0],
      -3 * (std::log2(1 - 0.5) + std::log2(1 - 0.5) + std::log2(1 - 0.9)),
      1e-6);
  EXPECT_NEAR(fidDistance[5][1], -3 * (std::log2(1 - 0.5) + std::log2(1 - 0.5)),
              1e-6);
  EXPECT_NEAR(fidDistance[5][2], -3 * std::log2(1 - 0.5), 1e-6);
  EXPECT_NEAR(fidDistance[5][3], -3 * (std::log2(1 - 0.5) + std::log2(1 - 0.1)),
              1e-6);
  EXPECT_NEAR(
      fidDistance[5][4],
      -3 * (std::log2(1 - 0.5) + std::log2(1 - 0.5) + std::log2(1 - 0.1)) -
          2 * (std::log2(1 - 0.03) + std::log2(1 - 0.03)),
      1e-6);
  EXPECT_NEAR(fidDistance[5][6],
              -3 * std::log2(1 - 0.9) -
                  2 * (std::log2(1 - 0.02) + std::log2(1 - 0.03)),
              1e-6);
  EXPECT_EQ(fidDistance[6].size(), 7);
  EXPECT_NEAR(fidDistance[6][0],
              -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.5) +
                    std::log2(1 - 0.5) + std::log2(1 - 0.9)) -
                  2 * (std::log2(1 - 0.02) + std::log2(1 - 0.03)),
              1e-6);
  EXPECT_NEAR(
      fidDistance[6][1],
      -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.5) + std::log2(1 - 0.5)) -
          2 * (std::log2(1 - 0.02) + std::log2(1 - 0.03)),
      1e-6);
  EXPECT_NEAR(fidDistance[6][2],
              -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.5)) -
                  2 * (std::log2(1 - 0.02) + std::log2(1 - 0.03)),
              1e-6);
  EXPECT_NEAR(
      fidDistance[6][3],
      -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.5) + std::log2(1 - 0.1)) -
          2 * (std::log2(1 - 0.02) + std::log2(1 - 0.03)),
      1e-6);
  EXPECT_NEAR(fidDistance[6][4],
              -3 * (std::log2(1 - 0.9) + std::log2(1 - 0.5) +
                    std::log2(1 - 0.5) + std::log2(1 - 0.1)) -
                  2 * (std::log2(1 - 0.02) + std::log2(1 - 0.03) +
                       std::log2(1 - 0.03) + std::log2(1 - 0.03)),
              1e-6);
  EXPECT_NEAR(fidDistance[6][5],
              -3 * std::log2(1 - 0.9) -
                  2 * (std::log2(1 - 0.02) + std::log2(1 - 0.03)),
              1e-6);
}

TEST(TestArchitecture, FidelitySwapCostTest) {
  Architecture      architecture{};
  const CouplingMap cm = {{0, 1}, {1, 2}, {2, 1}, {2, 3}, {2, 4}, {4, 2}};
  architecture.loadCouplingMap(5, cm);

  auto props = Architecture::Properties();
  props.setSingleQubitErrorRate(0, "x", 0.11);
  props.setSingleQubitErrorRate(1, "x", 0.12);
  props.setSingleQubitErrorRate(2, "x", 0.13);
  props.setSingleQubitErrorRate(3, "x", 0.14);
  props.setSingleQubitErrorRate(4, "x", 0.15);

  props.setTwoQubitErrorRate(0, 1, 0.1);
  props.setTwoQubitErrorRate(1, 2, 0.2);
  props.setTwoQubitErrorRate(2, 1, 0.2);
  props.setTwoQubitErrorRate(2, 3, 0.3);
  props.setTwoQubitErrorRate(2, 4, 0.4);
  props.setTwoQubitErrorRate(4, 2, 0.4);

  architecture.loadProperties(props);

  const Matrix<double> swapFidCost = architecture.getSwapFidelityCost();

  EXPECT_EQ(swapFidCost.size(), 5);
  EXPECT_EQ(swapFidCost[0].size(), 5);
  EXPECT_NEAR(swapFidCost[0][1],
              -3 * std::log2(1 - 0.1) - 2 * std::log2(1 - 0.11) -
                  2 * std::log2(1 - 0.12),
              1e-6);
  EXPECT_GT(swapFidCost[0][2], 1e20);
  EXPECT_GT(swapFidCost[0][3], 1e20);
  EXPECT_GT(swapFidCost[0][4], 1e20);
  EXPECT_EQ(swapFidCost[1].size(), 5);
  EXPECT_NEAR(swapFidCost[1][0],
              -3 * std::log2(1 - 0.1) - 2 * std::log2(1 - 0.11) -
                  2 * std::log2(1 - 0.12),
              1e-6);
  EXPECT_NEAR(swapFidCost[1][2], -3 * std::log2(1 - 0.2), 1e-6);
  EXPECT_GT(swapFidCost[1][3], 1e20);
  EXPECT_GT(swapFidCost[1][4], 1e20);
  EXPECT_EQ(swapFidCost[2].size(), 5);
  EXPECT_GT(swapFidCost[2][0], 1e20);
  EXPECT_NEAR(swapFidCost[2][1], -3 * std::log2(1 - 0.2), 1e-6);
  EXPECT_NEAR(swapFidCost[2][3],
              -3 * std::log2(1 - 0.3) - 2 * std::log2(1 - 0.13) -
                  2 * std::log2(1 - 0.14),
              1e-6);
  EXPECT_NEAR(swapFidCost[2][4], -3 * std::log2(1 - 0.4), 1e-6);
  EXPECT_EQ(swapFidCost[3].size(), 5);
  EXPECT_GT(swapFidCost[3][0], 1e20);
  EXPECT_GT(swapFidCost[3][1], 1e20);
  EXPECT_NEAR(swapFidCost[3][2],
              -3 * std::log2(1 - 0.3) - 2 * std::log2(1 - 0.13) -
                  2 * std::log2(1 - 0.14),
              1e-6);
  EXPECT_GT(swapFidCost[3][4], 1e20);
  EXPECT_EQ(swapFidCost[4].size(), 5);
  EXPECT_GT(swapFidCost[4][0], 1e20);
  EXPECT_GT(swapFidCost[4][1], 1e20);
  EXPECT_NEAR(swapFidCost[4][2], -3 * std::log2(1 - 0.4), 1e-6);
  EXPECT_GT(swapFidCost[4][3], 1e20);
}

TEST(TestArchitecture, FidelityDistanceCheapestPathTest) {
  Architecture      architecture{};
  const CouplingMap cm = {{0, 1}, {1, 0}, {2, 1}, {2, 6}, {6, 2},
                          {0, 5}, {5, 0}, {5, 6}, {6, 5}, {0, 3},
                          {3, 0}, {3, 4}, {4, 3}, {4, 6}, {6, 4}};
  architecture.loadCouplingMap(7, cm);

  auto props = Architecture::Properties();
  props.setSingleQubitErrorRate(0, "x", 0.1);
  props.setSingleQubitErrorRate(1, "x", 0.1);
  props.setSingleQubitErrorRate(2, "x", 0.1);
  props.setSingleQubitErrorRate(3, "x", 0.1);
  props.setSingleQubitErrorRate(4, "x", 0.1);
  props.setSingleQubitErrorRate(5, "x", 0.1);
  props.setSingleQubitErrorRate(6, "x", 0.1);

  props.setTwoQubitErrorRate(0, 1, 0.1);
  props.setTwoQubitErrorRate(1, 0, 0.1);
  props.setTwoQubitErrorRate(2, 1, 0.1);
  props.setTwoQubitErrorRate(2, 6, 0.1);
  props.setTwoQubitErrorRate(6, 2, 0.1);
  props.setTwoQubitErrorRate(0, 5, 0.7);
  props.setTwoQubitErrorRate(5, 0, 0.7);
  props.setTwoQubitErrorRate(5, 6, 0.7);
  props.setTwoQubitErrorRate(6, 5, 0.7);
  props.setTwoQubitErrorRate(0, 3, 0.5);
  props.setTwoQubitErrorRate(3, 0, 0.5);
  props.setTwoQubitErrorRate(3, 4, 0.5);
  props.setTwoQubitErrorRate(4, 3, 0.5);
  props.setTwoQubitErrorRate(4, 6, 0.5);
  props.setTwoQubitErrorRate(6, 4, 0.5);

  architecture.loadProperties(props);

  const Matrix<double> fidDistance = architecture.getFidelityDistanceTable();

  EXPECT_EQ(fidDistance.size(), 7);
  EXPECT_EQ(fidDistance[0].size(), 7);
  EXPECT_NEAR(fidDistance[0][6],
              -3 * 3 * std::log2(1 - 0.1) - 2 * 2 * std::log2(1 - 0.1), 1e-6);
}

TEST(TestArchitecture, DistanceCheapestPathTest) {
  Architecture architecture{};

  // minimum number of backward edges on a path where the same path with forward
  // edges can afford at least 1 more edge and still be cheaper
  std::uint8_t nrEdges =
      1 + static_cast<std::uint8_t>(
              std::ceil(static_cast<double>(COST_BIDIRECTIONAL_SWAP) /
                        (static_cast<double>(COST_UNIDIRECTIONAL_SWAP) -
                         static_cast<double>(COST_BIDIRECTIONAL_SWAP))));

  CouplingMap cm = {};
  for (std::uint8_t i = 0; i < nrEdges; ++i) {
    cm.insert(Edge{i + 1, i});
  }
  for (std::uint8_t i = nrEdges + 1; i < 2 * nrEdges; ++i) {
    cm.insert(Edge{i, i + 1});
    cm.insert(Edge{i + 1, i});
  }
  cm.insert(Edge{0, nrEdges + 1});
  cm.insert(Edge{nrEdges + 1, 0});
  cm.insert(Edge{2 * nrEdges, nrEdges});
  cm.insert(Edge{nrEdges, 2 * nrEdges});
  architecture.loadCouplingMap(static_cast<std::uint16_t>(2 * nrEdges + 1), cm);

  const Matrix<double> distances = architecture.getDistanceTable();

  EXPECT_EQ(distances.size(), 2 * nrEdges + 1);
  EXPECT_EQ(distances[0].size(), 2 * nrEdges + 1);
  EXPECT_NEAR(distances[0][nrEdges], nrEdges * COST_BIDIRECTIONAL_SWAP, 1e-6);
}
