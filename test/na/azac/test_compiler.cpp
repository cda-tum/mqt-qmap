#include "na/azac/Compiler.hpp"

#include <cmath>
#include <cstddef>
#include <gtest/gtest.h>
#include <optional>
#include <utility>
#include <vector>

namespace na {

class TestCompiler : public testing::Test {
protected:
  Compiler compiler;
  void SetUp() override {
    std::istringstream settings(R"({
  "qasm_list": [
    "benchmark/my/"
  ],
  "zac_setting": [
    {
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
      "trivial_placement": false,
      "dynamic_placement": true,
      "use_window": true,
      "window_size": 1000,
      "reuse": true,
      "use_verifier": true
    }
  ],
  "simulation": false,
  "animation": false
})");
    ASSERT_NO_THROW(compiler.loadSettings(settings));
  }
};

TEST(TestAZAComplier, LoadSettings) {}

} // namespace na
