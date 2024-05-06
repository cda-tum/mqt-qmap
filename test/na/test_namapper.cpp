//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "Architecture.hpp"
#include "Configuration.hpp"
#include "NAMapper.hpp"

#include "gtest/gtest.h"

TEST(NAMapper, Exceptions) {
  std::istringstream archIS(R"({
      "name": "Nature",
      "initialZones": [
          "storage"
      ],
      "zones": [
          {
              "name": "entangling",
              "xmin": -300,
              "xmax": 656,
              "ymin": -10,
              "ymax": 46,
              "fidelity": 0.9959
          },
          {
              "name": "storage",
              "xmin": -300,
              "xmax": 656,
              "ymin": 47,
              "ymax": 121,
              "fidelity": 1
          },
          {
              "name": "readout",
              "xmin": -300,
              "xmax": 656,
              "ymin": 122,
              "ymax": 156,
              "fidelity": 0.99
          }
      ],
      "operations": [
          {
              "name": "rz",
              "type": "local",
              "zones": [
                  "entangling",
                  "storage",
                  "readout"
              ],
              "time": 0.5,
              "fidelity": 0.999
          },
          {
              "name": "ry",
              "type": "global",
              "zones": [
                  "entangling",
                  "storage",
                  "readout"
              ],
              "time": 0.5,
              "fidelity": 0.999
          },
          {
              "name": "cz",
              "type": "global",
              "zones": [
                  "entangling"
              ],
              "time": 0.2,
              "fidelity": 0.9959
          },
          {
              "name": "measure",
              "type": "global",
              "zones": [
                  "readout"
              ],
              "time": 0.2,
              "fidelity": 0.95
          }
      ],
      "decoherence": {
          "t1": 100000000,
          "t2": 1500000
      },
      "interactionRadius": 2,
      "noInteractionRadius": 5,
      "minAtomDistance": 1,
      "shuttling": [
          {
              "rows": 5,
              "columns": 5,
              "xmin": -2.5,
              "xmax": 2.5,
              "ymin": -2.5,
              "ymax": 2.5,
              "move": {
                  "speed": 0.55,
                  "fidelity": 1
              },
              "load": {
                  "time": 20,
                  "fidelity": 1
              },
              "store": {
                  "time": 20,
                  "fidelity": 1
              }
          }
      ]
  })");
  std::stringstream  gridSS;
  gridSS << "x,y\n";
  // entangling zone (4 x 36 = 144 sites)
  for (std::size_t y = 0; y <= 36; y += 12) {
    for (std::size_t x = 3; x <= 353; x += 10) {
      gridSS << x << "," << y << "\n";
    }
  }
  // storage zone (12 x 72 = 864 sites)
  for (std::size_t y = 56; y <= 111; y += 5) {
    for (std::size_t x = 0; x <= 355; x += 5) {
      gridSS << x << "," << y << "\n";
    }
  }
  // readout zone (4 x 72 = 288 sites)
  for (std::size_t y = 131; y <= 146; y += 5) {
    for (std::size_t x = 0; x <= 355; x += 5) {
      gridSS << x << "," << y << "\n";
    }
  }
  // total: 1296 sites
  const auto&       arch = na::Architecture(archIS, gridSS);
  // ---------------------------------------------------------------------
  na::NAMapper mapper(
      arch, na::Configuration(3, 3, na::NAMappingMethod::MaximizeParallelism));
  EXPECT_THROW(std::ignore = mapper.getResult(), std::logic_error);
  EXPECT_THROW(std::ignore = mapper.getStats(), std::logic_error);
  EXPECT_THROW(
      mapper.map(qc::QuantumComputation::fromQASM(
          "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[5];\nx q[0];\n")),
      std::invalid_argument);
  EXPECT_THROW(mapper.map(qc::QuantumComputation::fromQASM(
                   "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg "
                   "q[5];\nry(pi/2) q[0];\n")),
               std::invalid_argument);
  EXPECT_THROW(
      mapper.map(qc::QuantumComputation::fromQASM(
          "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[5];\nrz(pi/2) q;\n")),
      std::invalid_argument);
  EXPECT_THROW(mapper.map(qc::QuantumComputation::fromQASM(
                   "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[5];\nccz "
                   "q[0], q[1], q[2];\n")),
               std::logic_error);
  EXPECT_THROW(mapper.map(qc::QuantumComputation::fromQASM(
                   "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[5];\ncx "
                   "q[0], q[1];\n")),
               std::logic_error);
}

TEST(NAMapper, QAOA10) {
  std::istringstream archIS(R"({
      "name": "Nature",
      "initialZones": [
          "storage"
      ],
      "zones": [
          {
              "name": "entangling",
              "xmin": -300,
              "xmax": 656,
              "ymin": -10,
              "ymax": 46,
              "fidelity": 0.9959
          },
          {
              "name": "storage",
              "xmin": -300,
              "xmax": 656,
              "ymin": 47,
              "ymax": 121,
              "fidelity": 1
          },
          {
              "name": "readout",
              "xmin": -300,
              "xmax": 656,
              "ymin": 122,
              "ymax": 156,
              "fidelity": 0.99
          }
      ],
      "operations": [
          {
              "name": "rz",
              "type": "local",
              "zones": [
                  "entangling",
                  "storage",
                  "readout"
              ],
              "time": 0.5,
              "fidelity": 0.999
          },
          {
              "name": "ry",
              "type": "global",
              "zones": [
                  "entangling",
                  "storage",
                  "readout"
              ],
              "time": 0.5,
              "fidelity": 0.999
          },
          {
              "name": "cz",
              "type": "global",
              "zones": [
                  "entangling"
              ],
              "time": 0.2,
              "fidelity": 0.9959
          },
          {
              "name": "measure",
              "type": "global",
              "zones": [
                  "readout"
              ],
              "time": 0.2,
              "fidelity": 0.95
          }
      ],
      "decoherence": {
          "t1": 100000000,
          "t2": 1500000
      },
      "interactionRadius": 2,
      "noInteractionRadius": 5,
      "minAtomDistance": 1,
      "shuttling": [
          {
              "rows": 5,
              "columns": 5,
              "xmin": -2.5,
              "xmax": 2.5,
              "ymin": -2.5,
              "ymax": 2.5,
              "move": {
                  "speed": 0.55,
                  "fidelity": 1
              },
              "load": {
                  "time": 20,
                  "fidelity": 1
              },
              "store": {
                  "time": 20,
                  "fidelity": 1
              }
          }
      ]
  })");
  std::stringstream  gridSS;
  gridSS << "x,y\n";
  // entangling zone (4 x 36 = 144 sites)
  for (std::size_t y = 0; y <= 36; y += 12) {
    for (std::size_t x = 3; x <= 353; x += 10) {
      gridSS << x << "," << y << "\n";
    }
  }
  // storage zone (12 x 72 = 864 sites)
  for (std::size_t y = 56; y <= 111; y += 5) {
    for (std::size_t x = 0; x <= 355; x += 5) {
      gridSS << x << "," << y << "\n";
    }
  }
  // readout zone (4 x 72 = 288 sites)
  for (std::size_t y = 131; y <= 146; y += 5) {
    for (std::size_t x = 0; x <= 355; x += 5) {
      gridSS << x << "," << y << "\n";
    }
  }
  // total: 1296 sites
  // For the test, we removed all rz gates because the mapping task remains the
  // same
  const std::string qasm = R"(OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
ry(-pi/4) q;
ry(pi/4) q;
cp(pi) q[0],q[1];
cp(pi) q[5],q[7];
ry(-pi/2) q;
ry(pi/2) q;
cp(pi) q[0],q[1];
cp(pi) q[5],q[7];
ry(-pi/4) q;
ry(pi/4) q;
cp(pi) q[1],q[3];
cp(pi) q[0],q[6];
cp(pi) q[9],q[5];
ry(-pi/2) q;
ry(pi/2) q;
cp(pi) q[1],q[3];
cp(pi) q[0],q[6];
cp(pi) q[9],q[5];
ry(-1.0312062) q;
ry(1.0312062) q;
cp(pi) q[3],q[4];
cp(pi) q[0],q[1];
cp(pi) q[2],q[6];
cp(pi) q[8],q[2];
cp(pi) q[9],q[7];
ry(-pi/2) q;
ry(pi/2) q;
cp(pi) q[3],q[4];
cp(pi) q[0],q[1];
cp(pi) q[2],q[6];
cp(pi) q[8],q[2];
cp(pi) q[9],q[7];
ry(-1.0312062) q;
ry(1.0312062) q;
cp(pi) q[1],q[3];
cp(pi) q[0],q[6];
cp(pi) q[8],q[4];
cp(pi) q[5],q[7];
cp(pi) q[9],q[5];
ry(-0.92609333) q;
ry(0.92609333) q;
cp(pi) q[1],q[3];
cp(pi) q[0],q[6];
cp(pi) q[2],q[6];
cp(pi) q[8],q[4];
cp(pi) q[5],q[7];
cp(pi) q[9],q[5];
ry(-pi/4) q;
ry(pi/4) q;
cp(pi) q[2],q[6];
cp(pi) q[3],q[4];
cp(pi) q[8],q[2];
cp(pi) q[9],q[7];
ry(-0.3223636) q;
ry(0.3223636) q;
cp(pi) q[3],q[4];
cp(pi) q[8],q[2];
cp(pi) q[9],q[7];
ry(-pi/4) q;
ry(pi/4) q;
cp(pi) q[8],q[4];
ry(-0.3223636) q;
ry(0.3223636) q;
cp(pi) q[8],q[4];
ry(-pi/4) q;
ry(pi/4) q;)";
  const auto&       circ = qc::QuantumComputation::fromQASM(qasm);
  const auto&       arch = na::Architecture(archIS, gridSS);
  // ---------------------------------------------------------------------
  na::NAMapper mapper(
      arch, na::Configuration(3, 3, na::NAMappingMethod::MaximizeParallelism));
  mapper.map(circ);
  std::ignore = mapper.getResult();
  std::ignore = mapper.getStats();
  // ---------------------------------------------------------------------
  na::NAMapper mapper2(arch,
                       na::Configuration(1, 1, na::NAMappingMethod::Naive));
  mapper2.map(circ);
  // ---------------------------------------------------------------------
}

TEST(NAMapper, QAOA16Narrow) {
  std::istringstream archIS(R"({
      "name": "Nature",
      "initialZones": [
          "storage"
      ],
      "zones": [
          {
              "name": "entangling",
              "xmin": -300,
              "xmax": 656,
              "ymin": -10,
              "ymax": 46,
              "fidelity": 0.9959
          },
          {
              "name": "storage",
              "xmin": -300,
              "xmax": 656,
              "ymin": 47,
              "ymax": 421,
              "fidelity": 1
          },
          {
              "name": "readout",
              "xmin": -300,
              "xmax": 656,
              "ymin": 422,
              "ymax": 456,
              "fidelity": 0.99
          }
      ],
      "operations": [
          {
              "name": "rz",
              "type": "local",
              "zones": [
                  "entangling",
                  "storage",
                  "readout"
              ],
              "time": 0.5,
              "fidelity": 0.999
          },
          {
              "name": "ry",
              "type": "global",
              "zones": [
                  "entangling",
                  "storage",
                  "readout"
              ],
              "time": 0.5,
              "fidelity": 0.999
          },
          {
              "name": "cz",
              "type": "global",
              "zones": [
                  "entangling"
              ],
              "time": 0.2,
              "fidelity": 0.9959
          },
          {
              "name": "measure",
              "type": "global",
              "zones": [
                  "readout"
              ],
              "time": 0.2,
              "fidelity": 0.95
          }
      ],
      "decoherence": {
          "t1": 100000000,
          "t2": 1500000
      },
      "interactionRadius": 2,
      "noInteractionRadius": 5,
      "minAtomDistance": 1,
      "shuttling": [
          {
              "rows": 5,
              "columns": 5,
              "xmin": -2.5,
              "xmax": 2.5,
              "ymin": -2.5,
              "ymax": 2.5,
              "move": {
                  "speed": 0.55,
                  "fidelity": 1
              },
              "load": {
                  "time": 20,
                  "fidelity": 1
              },
              "store": {
                  "time": 20,
                  "fidelity": 1
              }
          }
      ]
  }
  )");
  std::stringstream  gridSS;
  gridSS << "x,y\n";
  // entangling zone (4 x 36 = 144 sites)
  for (std::size_t y = 0; y <= 36; y += 12) {
    for (std::size_t x = 3; x <= 353; x += 10) {
      gridSS << x << "," << y << "\n";
    }
  }
  // storage zone (72 x 12 = 864 sites)
  for (std::size_t y = 56; y <= 411; y += 5) {
    for (std::size_t x = 150; x <= 205; x += 5) {
      gridSS << x << "," << y << "\n";
    }
  }
  // readout zone (4 x 12 = 48 sites)
  for (std::size_t y = 431; y <= 446; y += 5) {
    for (std::size_t x = 150; x <= 205; x += 5) {
      gridSS << x << "," << y << "\n";
    }
  }
  // total: 1056 sites
  // For the test, we removed all rz gates because the mapping task remains the
  // same
  const std::string qasm = R"(OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
ry(-pi/4) q;
ry(pi/4) q;
cp(pi) q[0],q[2];
cp(pi) q[1],q[7];
cp(pi) q[8],q[3];
cp(pi) q[12],q[6];
ry(-pi/2) q;
ry(pi/2) q;
cp(pi) q[0],q[2];
cp(pi) q[1],q[7];
cp(pi) q[8],q[3];
cp(pi) q[12],q[6];
ry(-pi/2) q;
ry(pi/2) q;
cp(pi) q[0],q[4];
cp(pi) q[8],q[9];
cp(pi) q[1],q[10];
cp(pi) q[13],q[6];
cp(pi) q[2],q[14];
cp(pi) q[3],q[15];
ry(-pi/2) q;
ry(pi/2) q;
cp(pi) q[0],q[4];
cp(pi) q[8],q[9];
cp(pi) q[1],q[10];
cp(pi) q[13],q[6];
cp(pi) q[2],q[14];
cp(pi) q[3],q[15];
ry(-pi/4) q;
ry(pi/4) q;
cp(pi) q[4],q[5];
cp(pi) q[12],q[13];
cp(pi) q[0],q[2];
cp(pi) q[14],q[7];
cp(pi) q[10],q[15];
cp(pi) q[8],q[3];
ry(-pi/2) q;
ry(pi/2) q;
cp(pi) q[4],q[5];
cp(pi) q[12],q[13];
cp(pi) q[0],q[2];
cp(pi) q[14],q[7];
cp(pi) q[10],q[15];
cp(pi) q[8],q[3];
ry(-pi/4) q;
ry(pi/4) q;
cp(pi) q[11],q[5];
cp(pi) q[12],q[6];
cp(pi) q[13],q[6];
cp(pi) q[0],q[4];
cp(pi) q[2],q[14];
cp(pi) q[1],q[7];
cp(pi) q[1],q[10];
cp(pi) q[3],q[15];
ry(-pi/2) q;
ry(pi/2) q;
cp(pi) q[11],q[5];
cp(pi) q[12],q[6];
cp(pi) q[13],q[6];
cp(pi) q[0],q[4];
cp(pi) q[2],q[14];
cp(pi) q[1],q[7];
cp(pi) q[1],q[10];
cp(pi) q[3],q[15];
ry(-pi/4) q;
ry(pi/4) q;
cp(pi) q[9],q[11];
cp(pi) q[12],q[13];
cp(pi) q[4],q[5];
cp(pi) q[14],q[7];
cp(pi) q[10],q[15];
ry(-0.64469806) q;
ry(0.64469806) q;
cp(pi) q[9],q[11];
cp(pi) q[12],q[13];
cp(pi) q[4],q[5];
cp(pi) q[14],q[7];
cp(pi) q[10],q[15];
ry(-2.2154814) q;
ry(2.2154814) q;
cp(pi) q[11],q[5];
cp(pi) q[8],q[9];
ry(-0.3223291) q;
ry(0.3223291) q;
cp(pi) q[11],q[5];
cp(pi) q[8],q[9];
ry(-pi/4) q;
ry(pi/4) q;
cp(pi) q[9],q[11];
ry(-0.3223291) q;
ry(0.3223291) q;
cp(pi) q[9],q[11];
ry(-2.2154814) q;
ry(2.2154814) q;)";
  const auto&       circ = qc::QuantumComputation::fromQASM(qasm);
  const auto&       arch = na::Architecture(archIS, gridSS);
  // ---------------------------------------------------------------------
  na::NAMapper mapper(
      arch, na::Configuration(3, 3, na::NAMappingMethod::MaximizeParallelism));
  mapper.map(circ);
  std::ignore = mapper.getStats();
  std::ignore = mapper.getResult();
}
