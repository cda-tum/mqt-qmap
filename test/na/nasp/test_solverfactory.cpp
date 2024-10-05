#include "ir/QuantumComputation.hpp"
#include "ir/operations/OpType.hpp"
#include "na/Architecture.hpp"
#include "na/nasp/SolverFactory.hpp"

#include <cstddef>
#include <cstdint>
#include <gtest/gtest.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>

TEST(SolverFactory, Create) {
  na::Architecture arch;
  // write content to a file
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
            "ymin": -20,
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
  std::stringstream gridSS;
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
  arch.fromFileStream(archIS, gridSS);
  // create solver
  auto solver = na::SolverFactory::create(arch);
  const std::string qasm = R"(OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
cz q[0],q[6];
cz q[1],q[3];
cz q[4],q[5];
cz q[0],q[4];
cz q[5],q[6];
cz q[1],q[2];
cz q[0],q[2];
cz q[3],q[5];
cz q[1],q[4];
h q[2];
h q[3];
h q[4];
h q[6];
)";
  const auto& circ = qc::QuantumComputation::fromQASM(qasm);
  // get operations for solver
  const auto& pairs =
      na::SolverFactory::getOpsForSolver(circ, {qc::Z, 1}, true);
  // solve
  const auto result = solver.solve(
      pairs, static_cast<std::uint16_t>(circ.getNqubits()), 5, false, true);
  EXPECT_TRUE(result.isSat());
}

TEST(SolverFactory, CreateExceptions) {
  na::Architecture arch;
  // write content to a file
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
            "ymin": -20,
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
  std::stringstream gridSS;
  gridSS << "x,y\n";
  // entangling zone (1 x 36 = 36 sites)
  for (std::size_t y = 0; y <= 0; y += 12) {
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
  arch.fromFileStream(archIS, gridSS);
  EXPECT_THROW(std::ignore = na::SolverFactory::create(arch),
               std::invalid_argument);

  auto circ = qc::QuantumComputation(3);
  circ.h(0);
  circ.cz(0, 1);
  circ.cecr(0, 1, 2);
  // get operations for solver
  EXPECT_THROW(std::ignore =
                   na::SolverFactory::getOpsForSolver(circ, {qc::Z, 1}, false),
               std::invalid_argument);
  EXPECT_THROW(std::ignore =
                   na::SolverFactory::getOpsForSolver(circ, {qc::ECR, 1}, true),
               std::invalid_argument);
}
