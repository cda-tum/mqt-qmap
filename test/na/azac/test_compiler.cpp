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
  }
})";
class TestCompiler : public ::testing::TestWithParam<std::string> {
private:
  nlohmann::json settings_;
  Architecture architecture_;

protected:
  qc::QuantumComputation circ_;
  ZACompiler compiler_;
  TestCompiler()
      : settings_(nlohmann::json::parse(settings)),
        architecture_(settings_["architecture"]),
        compiler_(architecture_, settings_) {}
  void SetUp() override {
    const auto& path = GetParam();
    circ_ = qasm3::Importer::importf(path);
    // qc::CircuitOptimizer::flattenOperations(circ_);
  }
};
TEST_P(TestCompiler, EndToEnd) {
  const auto& code = compiler_.compile(circ_);
  EXPECT_TRUE(code.validate().first);
  double timeSum = 0;
  const auto& stats = compiler_.getStatistics().asJson();
  for (const auto& [key, value] : stats.items()) {
    if (key != "total_time") {
      timeSum += value.get<double>();
    }
  }
  EXPECT_GE(stats["total_time"], timeSum);
}
INSTANTIATE_TEST_SUITE_P(TestCompiler, // Custom instantiation name
                         TestCompiler, // Test suite name
                         // Parameters to test with
                         ::testing::Values(TEST_CIRCUITS),
                         [](const ::testing::TestParamInfo<std::string>& info) {
                           const auto& path = info.param;
                           const auto& filename =
                               path.substr(path.find_last_of("/") + 1);
                           return filename.substr(0,
                                                  filename.find_last_of("."));
                         });
} // namespace na
