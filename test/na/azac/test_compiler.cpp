#include "na/azac/Compiler.hpp"

#include <cmath>
#include <cstddef>
#include <gtest/gtest.h>
#include <iostream>
#include <optional>
#include <sstream>
#include <utility>
#include <vector>

namespace na {

class TestAZACompiler : public testing::Test {
protected:
  Compiler compiler;
  void SetUp() override {
    std::istringstream settings(R"({
  "arch_spec": {
    "name": "full_compute_store_architecture",
    "operation_duration": {"rydberg": 0.36, "1qGate": 52, "atom_transfer": 15},
    "operation_fidelity": {
      "two_qubit_gate": 0.995,
      "single_qubit_gate": 0.9997,
      "atom_transfer": 0.999
    },
    "qubit_spec": {"T": 1.5e6},
    "storage_zones": [{
      "zone_id": 0,
      "slms": [{"id": 0, "site_seperation": [3, 3], "r": 100, "c": 100, "location": [0, 0]}],
      "offset": [0, 0],
      "dimenstion": [300, 300]
    }],
    "entanglement_zones": [{
      "zone_id": 0,
      "slms": [
        {"id": 1, "site_seperation": [12, 10], "r": 7, "c": 20, "location": [35, 307]},
        {"id": 2, "site_seperation": [12, 10], "r": 7, "c": 20, "location": [37, 307]}
      ],
      "offset": [35, 307],
      "dimension": [240, 70]
    }],
    "aods":[{"id": 0, "site_seperation": 2, "r": 100, "c": 100}],
    "arch_range": [[0, 0], [297, 402]],
    "rydberg_range": [[[5, 305], [292, 402]]]
  },
  "dependency": true,
  "dir": "result/",
  "routing_strategy": "maximalis_sort",
  "scheduling": "asap",
  "trivial_placement": true,
  "dynamic_placement": true,
  "use_window": true,
  "window_size": 1000,
  "reuse": true,
  "use_verifier": false
})");
    ASSERT_NO_THROW(compiler.loadSettings(settings));
  }
};

TEST_F(TestAZACompiler, LoadSettings) {}

TEST_F(TestAZACompiler, Settings) {
  std::cout << compiler.toString() << '\n';
  EXPECT_FALSE(compiler.toString().empty());
}

} // namespace na
