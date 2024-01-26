//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "na/Architecture.hpp"

#include "gtest/gtest.h"
#include <random>
#include <sstream>
#include <string>

class TestArchitecture : public testing::TestWithParam<std::string> {
protected:
  std::string testArchitectureDir = "../examples/na/";
};

INSTANTIATE_TEST_SUITE_P(Architecture, TestArchitecture,
                         testing::Values("nature.csv"));

TEST_P(TestArchitecture, QubitMap) {
  const auto&       archName = GetParam();
  std::stringstream ss;
  ss << testArchitectureDir << archName;
  std::string filename = ss.str();
  na::Architecture arch(filename);

  // EXPECT_EQ();
  // EXPECT_TRUE();
  // EXPECT_THROW();
  // EXPECT_NEAR(,1e-6);
}
