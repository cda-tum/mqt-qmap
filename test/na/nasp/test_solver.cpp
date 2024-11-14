#include "ir/QuantumComputation.hpp"
#include "ir/operations/OpType.hpp"
#include "na/nasp/Solver.hpp"
#include "na/nasp/SolverFactory.hpp"

#include <algorithm>
#include <cstdint>
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
#include <tuple>
#include <yaml-cpp/node/parse.h>

TEST(Solver, SteaneDoubleSidedStorage) {
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
  // create solver
  na::NASolver solver;
  solver.init(3, 7, 2, 3, 2, 2, 2, 2, 2, 4);
  // get operations for solver
  const auto& pairs =
      na::SolverFactory::getOpsForSolver(circ, {qc::Z, 1}, true);
  // solve
  const auto result = solver.solve(
      pairs, static_cast<std::uint16_t>(circ.getNqubits()), 4, false, true);
  EXPECT_TRUE(result.sat);
  EXPECT_EQ(result.stages.size(), 4);
}

TEST(Solver, SteaneBottomStorage) {
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
  // create solver
  na::NASolver solver;
  solver.init(3, 7, 2, 3, 2, 2, 2, 2, 0, 4);
  // get operations for solver
  const auto& pairs =
      na::SolverFactory::getOpsForSolver(circ, {qc::Z, 1}, true);
  // solve
  const auto resultUnsat = solver.solve(
      pairs, static_cast<std::uint16_t>(circ.getNqubits()), 4, false, true);
  EXPECT_FALSE(resultUnsat.sat);
  const auto resultSat = solver.solve(
      pairs, static_cast<std::uint16_t>(circ.getNqubits()), 5, false, true);
  EXPECT_TRUE(resultSat.sat);
  EXPECT_TRUE(resultSat.stages.front().rydberg);
  for (const auto& q : resultSat.stages.front().qubits) {
    EXPECT_GE(q.x, 0);
    EXPECT_LE(q.x, 3);
    EXPECT_GE(q.y, 0);
    EXPECT_LE(q.y, 7);
    EXPECT_GE(q.c, 0);
    EXPECT_LE(q.c, 2);
    EXPECT_GE(q.r, 0);
    EXPECT_LE(q.r, 3);
    EXPECT_GE(q.h, -2);
    EXPECT_LE(q.h, 2);
    EXPECT_GE(q.v, -2);
    EXPECT_LE(q.v, 2);
    EXPECT_NO_THROW(std::ignore = q.a);
  }
  for (const auto& g : resultSat.stages.front().gates) {
    EXPECT_TRUE(std::find(pairs.cbegin(), pairs.cend(), g.qubits) !=
                pairs.cend());
  }
}

TEST(Solver, NoShieldingFixedOrder) {
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
  // create solver
  na::NASolver solver;
  solver.init(3, 7, 2, 3, 2, 2, 2, 2, 0, 7);
  // get operations for solver
  const auto& pairs =
      na::SolverFactory::getOpsForSolver(circ, {qc::Z, 1}, true);
  // solve
  const auto result = solver.solve(
      pairs, static_cast<std::uint16_t>(circ.getNqubits()), 3, true, false);
  EXPECT_TRUE(result.sat);
}

TEST(Solver, FixedTransfer) {
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
  // create solver
  na::NASolver solver;
  solver.init(3, 7, 2, 3, 2, 2, 2, 2, 2, 4);
  // get operations for solver
  const auto& pairs =
      na::SolverFactory::getOpsForSolver(circ, {qc::Z, 1}, true);
  // solve
  const auto result = solver.solve(
      pairs, static_cast<std::uint16_t>(circ.getNqubits()), 5, 2, false, true);
  EXPECT_TRUE(result.sat);
}

TEST(Solver, Unsat) {
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
  // create solver
  na::NASolver solver;
  solver.init(3, 7, 2, 3, 2, 2, 2, 2, 2, 4);
  // get operations for solver
  const auto& pairs =
      na::SolverFactory::getOpsForSolver(circ, {qc::Z, 1}, true);
  // solve
  const auto result = solver.solve(
      pairs, static_cast<std::uint16_t>(circ.getNqubits()), 3, false, true);
  EXPECT_FALSE(result.sat);
}

TEST(Solver, Exceptions) {
  na::NASolver solver;
  EXPECT_THROW(solver.init(3, 7, 2, 3, 2, 2, 2, 2, 2, 7),
               std::invalid_argument);
  solver.init(3, 7, 2, 3, 2, 2, 2, 2, 0, 7);
  EXPECT_THROW(std::ignore = solver.solve({{0, 1}}, 3, 1, false, true),
               std::invalid_argument);
}

TEST(Solver, YAMLRoundTrip) {
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
  // create solver
  na::NASolver solver;
  solver.init(3, 7, 2, 3, 2, 2, 2, 2, 2, 4);
  // get operations for solver
  const auto& pairs =
      na::SolverFactory::getOpsForSolver(circ, {qc::Z, 1}, true);
  // solve
  const auto result = solver.solve(
      pairs, static_cast<std::uint16_t>(circ.getNqubits()), 4, false, true);
  const auto resultRT = na::NASolver::Result::fromYAML(
      YAML::Load(result.yaml())); // Round-Tripped result
  EXPECT_EQ(resultRT, result);
}
