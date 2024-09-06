//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "na/Configuration.hpp"

#include <gtest/gtest.h>
#include <sstream>
#include <stdexcept>
#include <tuple>

TEST(Configuration, MethodOfString) {
  EXPECT_EQ(na::getMethodOfString("naive"), na::NAMappingMethod::Naive);
  EXPECT_EQ(na::getMethodOfString("maximize parallelism"),
            na::NAMappingMethod::MaximizeParallelismHeuristic);
  EXPECT_EQ(na::getMethodOfString("NaIvE"), na::NAMappingMethod::Naive);
  EXPECT_EQ(na::getMethodOfString("mAxImIzE pArAllElIsm"),
            na::NAMappingMethod::MaximizeParallelismHeuristic);
  EXPECT_THROW(std::ignore = na::getMethodOfString("unsupported"),
               std::invalid_argument);
}

TEST(Configuration, Import) {
  EXPECT_THROW(na::Configuration("nonexistent.json"), std::runtime_error);
  std::istringstream configIS(R"(
    {
      "patch": {
        "rows": 2,
        "cols": 3
      },
      "method": "maximize parallelism"
    }
  )");
  const na::Configuration config(configIS);
  EXPECT_EQ(config.getPatchRows(), 2);
  EXPECT_EQ(config.getPatchCols(), 3);
  EXPECT_EQ(config.getMethod(),
            na::NAMappingMethod::MaximizeParallelismHeuristic);
  std::istringstream invalidJson("{name: invalid}");
  EXPECT_THROW(const na::Configuration ignore(invalidJson), std::runtime_error);
}
