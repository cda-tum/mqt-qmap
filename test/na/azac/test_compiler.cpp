#include "circuit_optimizer/CircuitOptimizer.hpp"
#include "ir/QuantumComputation.hpp"
#include "na/azac/Compiler.hpp"
#include "qasm3/Importer.hpp"

#include <cstddef>
#include <gtest/gtest.h>
#include <utility>

namespace na {

constexpr std::string_view settings = R"({
  "architecture": {
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
    "arch_range": [[0, 0], [60, 110]],
    "rydberg_range": [[[5, 70], [55, 110]]]
  },
  "vm_placer" : {
    "use_window" : true,
    "window_size" : 10,
    "dynamic_placement" : true
  },
  "code_generator" : {
    "parking_offset" : 1
  },
  "a_star_placer" : {
    "use_window" : true,
    "window_height" : 6,
    "window_width" : 4
  }
})";
#define COMPILER_TEST(compiler_type)                                           \
  class compiler_type##Test : public ::testing::TestWithParam<std::string> {   \
    nlohmann::json settings_;                                                  \
    Architecture architecture_;                                                \
                                                                               \
  protected:                                                                   \
    qc::QuantumComputation circ_;                                              \
    compiler_type compiler_;                                                   \
    compiler_type##Test()                                                      \
        : settings_(nlohmann::json::parse(settings)),                          \
          architecture_(settings_["architecture"]),                            \
          compiler_(architecture_, settings_) {}                               \
    void SetUp() override { circ_ = qasm3::Importer::importf(GetParam()); }    \
  };                                                                           \
  /*=========================== END TO END TESTS ===========================*/ \
  TEST_P(compiler_type##Test, EndToEnd) {                                      \
    const auto& code = this->compiler_.compile(this->circ_);                   \
    EXPECT_TRUE(code.validate().first);                                        \
    double timeSum = 0;                                                        \
    const auto& stats = this->compiler_.getStatistics().asJson();              \
    for (const auto& [key, value] : stats.items()) {                           \
      if (key != "total_time") {                                               \
        timeSum += value.get<double>();                                        \
      }                                                                        \
    }                                                                          \
    EXPECT_GE(stats["total_time"], timeSum);                                   \
  }                                                                            \
  /*========================================================================*/ \
  INSTANTIATE_TEST_SUITE_P(                                                    \
      Allcompiler_type##Test,           /* Custom instantiation name */        \
      compiler_type##Test,              /* Test suite name */                  \
      ::testing::Values(TEST_CIRCUITS), /* Parameters to test with */          \
      [](const ::testing::TestParamInfo<std::string>& info) {                  \
        const auto& path = info.param;                                         \
        const auto& filename = path.substr(path.find_last_of("/") + 1);        \
        return filename.substr(0, filename.find_last_of("."));                 \
      })
/*============================== INSTANTIATIONS ==============================*/
COMPILER_TEST(ZACompiler);
COMPILER_TEST(AZACompiler);
} // namespace na
