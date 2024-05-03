//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "Architecture.hpp"
#include "Configuration.hpp"
#include "NAMapper.hpp"

#include "gtest/gtest.h"

TEST(NeutralAtomMapper, QAOA5) {
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
  const std::string qasm = R"(OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
rz(5.7203885) q[0];
rz(3*pi/2) q[1];
ry(-0.86876537) q;
rz(3.9153065) q[0];
rz(pi) q[1];
ry(0.86876537) q;
rz(4.6666067) q[0];
rz(pi/2) q[1];
cp(pi) q[0],q[1];
rz(3*pi/2) q[0];
rz(pi) q[1];
ry(-pi/2) q;
rz(pi) q[0];
rz(5.5433691) q[1];
ry(pi/2) q;
rz(3*pi/2) q[0];
rz(3*pi/2) q[1];
cp(pi) q[0],q[1];
rz(0.51701458) q[0];
rz(pi) q[1];
rz(0) q[3];
rz(pi) q[4];
ry(-pi/2) q;
rz(pi) q[0];
rz(3*pi/2) q[1];
rz(3.2944225) q[3];
rz(3.408658) q[4];
ry(pi/2) q;
rz(3*pi/2) q[0];
rz(2.9748582) q[1];
rz(pi) q[3];
rz(pi) q[4];
cp(pi) q[0],q[2];
cp(pi) q[1],q[3];
rz(3*pi/2) q[1];
rz(pi) q[2];
rz(pi) q[3];
ry(-pi/2) q;
rz(pi) q[1];
rz(3.8814089) q[2];
rz(5.5433691) q[3];
ry(pi/2) q;
rz(3*pi/2) q[1];
rz(pi) q[2];
rz(3*pi/2) q[3];
cp(pi) q[0],q[2];
cp(pi) q[1],q[3];
rz(3*pi/2) q[0];
rz(5.3006653) q[1];
rz(2.1590727) q[2];
ry(-1.0641775) q;
rz(pi) q[0];
rz(4.3993567) q[1];
rz(4.3993567) q[2];
ry(1.0641775) q;
rz(3*pi/2) q[0];
rz(4.2874278) q[1];
rz(0.58827634) q[2];
cp(pi) q[0],q[1];
cp(pi) q[2],q[4];
rz(pi) q[1];
rz(pi) q[3];
rz(pi) q[4];
ry(-pi/2) q;
rz(6.0088424) q[1];
rz(pi) q[2];
rz(3*pi/2) q[3];
rz(5.5433691) q[4];
ry(pi/2) q;
rz(pi) q[1];
rz(3*pi/2) q[2];
rz(6.1303554) q[3];
cp(pi) q[0],q[1];
cp(pi) q[2],q[4];
rz(3.276332) q[4];
ry(-pi/4) q;
rz(pi) q[2];
rz(5.9043581) q[4];
ry(pi/4) q;
rz(2.1283551) q[2];
rz(3.276332) q[4];
cp(pi) q[0],q[2];
cp(pi) q[3],q[4];
rz(5.1105864) q[1];
rz(3.1951436) q[2];
rz(3*pi/2) q[4];
ry(-1.2008882) q;
rz(4.5614608) q[1];
rz(5.9887998) q[2];
rz(pi) q[4];
ry(1.2008882) q;
rz(5.1105864) q[1];
rz(3.1951436) q[2];
rz(3*pi/2) q[4];
cp(pi) q[0],q[2];
cp(pi) q[3],q[4];
rz(pi/2) q[0];
rz(3.9006217) q[2];
rz(2.3298254) q[3];
rz(2.1283551) q[4];
ry(-2.1984352) q;
rz(pi) q[0];
rz(4.1580788) q[2];
rz(4.1580788) q[3];
ry(2.1984352) q;
rz(pi/2) q[0];
rz(3.9006217) q[2];
rz(6.0289768) q[3];
cp(pi) q[1],q[3];
cp(pi) q[2],q[4];
rz(3*pi/2) q[3];
rz(3*pi/2) q[4];
ry(-0.13717147) q;
rz(pi) q[3];
rz(pi) q[4];
ry(0.13717147) q;
rz(3*pi/2) q[3];
rz(3*pi/2) q[4];
cp(pi) q[1],q[3];
cp(pi) q[2],q[4];
rz(pi/2) q[1];
rz(pi/2) q[2];
rz(3.9006217) q[3];
ry(-2.1984352) q;
rz(pi) q[1];
rz(pi) q[2];
rz(4.1580788) q[3];
ry(2.1984352) q;
rz(pi/2) q[1];
rz(pi/2) q[2];
rz(3.9006217) q[3];
cp(pi) q[3],q[4];
rz(3*pi/2) q[4];
ry(-0.13717147) q;
rz(pi) q[4];
ry(0.13717147) q;
rz(3*pi/2) q[4];
cp(pi) q[3],q[4];
rz(pi/2) q[3];
rz(5.7869366) q[4];
ry(-2.1984352) q;
rz(pi) q[3];
rz(4.1580788) q[4];
ry(2.1984352) q;
rz(pi/2) q[3];
rz(3.9006217) q[4];)";
  const auto&       circ = qc::QuantumComputation::fromQASM(qasm);
  const auto&       arch = na::Architecture(archIS, gridSS);
  // ---------------------------------------------------------------------
  na::NAMapper mapper(
      arch, na::Configuration(3, 3, na::NAMappingMethod::MaximizeParallelism));
  EXPECT_THROW(std::ignore = mapper.getResult(), std::logic_error);
  EXPECT_THROW(std::ignore = mapper.getStats(), std::logic_error);
  mapper.map(circ);
  std::ignore = mapper.getStats();
  std::ignore = mapper.getResult();
  // ---------------------------------------------------------------------
  EXPECT_THROW(
      mapper.map(qc::QuantumComputation::fromQASM(
          "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[5];\nx q[0];\n")),
      std::invalid_argument);
  // ---------------------------------------------------------------------
  EXPECT_THROW(mapper.map(qc::QuantumComputation::fromQASM(
                   "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg "
                   "q[5];\nry(pi/2) q[0];\n")),
               std::invalid_argument);
  // ---------------------------------------------------------------------
  na::NAMapper mapper2(arch,
                       na::Configuration(1, 1, na::NAMappingMethod::Naive));
  mapper2.map(circ);
  // ---------------------------------------------------------------------
}
