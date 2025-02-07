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

constexpr std::string_view steaneWithourOneQubitgates = R"(OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
cz q[0],q[3];
cz q[0],q[4];
cz q[1],q[2];
cz q[1],q[5];
cz q[1],q[6];
cz q[2],q[3];
cz q[2],q[4];
cz q[3],q[5];
cz q[4],q[6];
)";

constexpr std::string_view steane = R"(OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
h q;
cz q[0],q[3];
cz q[0],q[4];
cz q[1],q[2];
cz q[1],q[5];
cz q[1],q[6];
cz q[2],q[3];
cz q[2],q[4];
cz q[3],q[5];
cz q[4],q[6];
h q[0];
h q[2];
h q[5];
h q[6];
)";

constexpr std::string_view settings = R"({
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
})";

class TestAZACompiler : public testing::Test {
protected:
  ZACompiler compiler;
  void SetUp() override {
    std::istringstream settingsStream{std::string{settings}};
    ASSERT_NO_THROW(compiler.loadSettings(settingsStream));
  }
};

TEST_F(TestAZACompiler, LoadSettingsNoThrow) {}

TEST_F(TestAZACompiler, PrintSettingsNonEmpty) {
  std::cout << compiler.toString() << '\n';
  EXPECT_FALSE(compiler.toString().empty());
}

TEST_F(TestAZACompiler, SolveNoThrow) {
  EXPECT_NO_THROW(compiler.solve());
}

} // namespace na
