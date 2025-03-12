#include "circuit_optimizer/CircuitOptimizer.hpp"
#include "ir/QuantumComputation.hpp"
#include "na/azac/Compiler.hpp"
#include "qasm3/Importer.hpp"

#include <cstddef>
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>

namespace na {

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
      "slms": [{"id": 0, "site_separation": [3, 3], "r": 100, "c": 100, "location": [0, 0]}],
      "offset": [0, 0],
      "dimension": [300, 300]
    }],
    "entanglement_zones": [{
      "zone_id": 0,
      "slms": [
        {"id": 1, "site_separation": [12, 10], "r": 7, "c": 20, "location": [35, 307]},
        {"id": 2, "site_separation": [12, 10], "r": 7, "c": 20, "location": [37, 307]}
      ],
      "offset": [35, 307],
      "dimension": [240, 70]
    }],
    "aods":[{"id": 0, "site_separation": 2, "r": 100, "c": 100}],
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

constexpr std::string_view steaneWithoutOneQubitGates = R"(OPENQASM 2.0;
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

class TestAZACompiler
    : public ::testing::TestWithParam<std::pair<std::string, std::string>> {
protected:
  qc::QuantumComputation circ;
  void SetUp() override {
    const auto& [name, qasm] = GetParam();
    circ = qasm3::Importer::imports(qasm);
    qc::CircuitOptimizer::flattenOperations(circ);
    std::istringstream settingsStream{std::string{settings}};
  }
};

INSTANTIATE_TEST_SUITE_P(
    SteanStatePreps, // Custom instantiation name
    TestAZACompiler, // Test suite name
    // Parameters to test with
    ::testing::Values(std::pair{"SteaneWithoutOneQubitGates",
                                steaneWithoutOneQubitGates},
                      std::pair{"Steane", steane}),
    [](const ::testing::TestParamInfo<std::pair<std::string, std::string>>&
           info) { return info.param.first; });

} // namespace na
