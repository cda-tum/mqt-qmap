/*
 * Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
 * Copyright (c) 2025 Munich Quantum Software Company GmbH
 * All rights reserved.
 *
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License
 */

#include "circuit_optimizer/CircuitOptimizer.hpp"
#include "ir/QuantumComputation.hpp"
#include "na/zoned/Compiler.hpp"
#include "qasm3/Importer.hpp"

#include <cstddef>
#include <gtest/gtest.h>
#include <utility>

namespace na::zoned {

constexpr std::string_view architectureSpecification = R"({
  "name": "compiler_architecture",
  "storage_zones": [{
    "zone_id": 0,
    "slms": [{"id": 0, "site_separation": [3, 3], "r": 20, "c": 20, "location": [0, 0]}],
    "offset": [0, 0],
    "dimension": [60, 60]
  }],
  "entanglement_zones": [{
    "zone_id": 0,
    "slms": [
      {"id": 1, "site_separation": [12, 10], "r": 4, "c": 4, "location": [5, 70]},
      {"id": 2, "site_separation": [12, 10], "r": 4, "c": 4, "location": [7, 70]}
    ],
    "offset": [5, 70],
    "dimension": [50, 40]
  }],
  "aods":[{"id": 0, "site_separation": 2, "r": 20, "c": 20}],
  "rydberg_range": [[[5, 70], [55, 110]]]
})";
constexpr std::string_view routingAgnosticConfiguration = R"({
  "placerConfig" : {
    "useWindow" : true,
    "windowSize" : 10,
    "dynamicPlacement" : true
  },
  "codeGeneratorConfig" : {
    "parkingOffset" : 1,
    "warnUnsupportedGates" : false
  }
})";
constexpr std::string_view routingAwareConfiguration = R"({
  "codeGeneratorConfig" : {
    "parkingOffset" : 1,
    "warnUnsupportedGates" : false
  },
  "placerConfig" : {
    "useWindow" : true,
    "windowMinWidth" : 4,
    "windowRatio" : 1.5,
    "windowShare" : 0.6,
    "deepeningFactor" : 0.6,
    "deepeningValue" : 0.2,
    "lookaheadFactor": 0.2,
    "reuseLevel": 5.0
  }
})";
#define COMPILER_TEST(compiler_type, config)                                   \
  TEST(compiler_type##Test, ConstructorWithoutSettings) {                      \
    Architecture architecture(                                                 \
        Architecture::fromJSONString(architectureSpecification));              \
    /* expected not to lead to a segfault */                                   \
    [[maybe_unused]] compiler_type compiler(architecture);                     \
  }                                                                            \
  class compiler_type##Test : public ::testing::TestWithParam<std::string> {   \
    compiler_type::Config settings_;                                           \
    Architecture architecture_;                                                \
                                                                               \
  protected:                                                                   \
    qc::QuantumComputation circ_;                                              \
    compiler_type compiler_;                                                   \
    compiler_type##Test()                                                      \
        : settings_(nlohmann::json(config)),                                   \
          architecture_(                                                       \
              Architecture::fromJSONString(architectureSpecification)),        \
          compiler_(architecture_, settings_) {}                               \
    void SetUp() override { circ_ = qasm3::Importer::importf(GetParam()); }    \
  };                                                                           \
  /*=========================== END TO END TESTS ===========================*/ \
  TEST_P(compiler_type##Test, EndToEnd) {                                      \
    const auto& code = this->compiler_.compile(this->circ_);                   \
    EXPECT_TRUE(code.validate().first);                                        \
    /*===----------------------------------------------------------------===*/ \
    /* write code to a file with extension .naviz in a directory converted */  \
    std::filesystem::path inputFile(GetParam());                               \
    std::filesystem::path outputFile = inputFile.parent_path() / "converted" / \
                                       #compiler_type /                        \
                                       (inputFile.stem().string() + ".naviz"); \
    std::filesystem::create_directories(outputFile.parent_path());             \
    std::ofstream output(outputFile);                                          \
    output << code;                                                            \
    /*===----------------------------------------------------------------===*/ \
    double timeSum = 0;                                                        \
    const nlohmann::json stats = this->compiler_.getStatistics();              \
    for (const auto& [key, value] : stats.items()) {                           \
      if (key != "totalTime") {                                                \
        timeSum += value.get<double>();                                        \
      }                                                                        \
    }                                                                          \
    EXPECT_GE(stats["totalTime"], timeSum);                                    \
  }                                                                            \
  /*========================================================================*/ \
  INSTANTIATE_TEST_SUITE_P(                                                    \
      compiler_type##TestWithCircuits,  /* Custom instantiation name */        \
      compiler_type##Test,              /* Test suite name */                  \
      ::testing::Values(TEST_CIRCUITS), /* Parameters to test with */          \
      [](const ::testing::TestParamInfo<std::string>& pinfo) {                 \
        const auto& path = pinfo.param;                                        \
        const auto& filename = path.substr(path.find_last_of("/") + 1);        \
        return filename.substr(0, filename.find_last_of("."));                 \
      })
/*============================== INSTANTIATIONS ==============================*/
COMPILER_TEST(RoutingAgnosticCompiler, routingAgnosticConfiguration);
COMPILER_TEST(RoutingAwareCompiler, routingAwareConfiguration);
} // namespace na::zoned
