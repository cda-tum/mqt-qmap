#include "QuantumComputation.hpp"
#include "Solver.hpp"
#include "SolverFactory.hpp"
#include "operations/OpType.hpp"

#include <algorithm>
#include <stdexcept>
#include <cstdint>
#include <gtest/gtest.h>
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
  const auto&       circ = qc::QuantumComputation::fromQASM(qasm);
  // create solver
  na::NASolver solver;
  solver.init(3, 7, 2, 3, 2, 2, 2, 2, 2, 4);
  // get operations for solver
  const auto& pairs =
      na::SolverFactory::getOpsForSolver(circ, {qc::Z, 1}, true);
  // solve
  const auto result = solver.solve(
      pairs, static_cast<std::uint16_t>(circ.getNqubits()), 4, false, true);
  EXPECT_TRUE(result.isSat());
  EXPECT_EQ(result.numStages(), 4);
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
  const auto&       circ = qc::QuantumComputation::fromQASM(qasm);
  // create solver
  na::NASolver solver;
  solver.init(3, 7, 2, 3, 2, 2, 2, 2, 0, 4);
  // get operations for solver
  const auto& pairs =
      na::SolverFactory::getOpsForSolver(circ, {qc::Z, 1}, true);
  // solve
  const auto resultUnsat = solver.solve(
      pairs, static_cast<std::uint16_t>(circ.getNqubits()), 4, false, true);
  EXPECT_FALSE(resultUnsat.isSat());
  const auto resultSat = solver.solve(
      pairs, static_cast<std::uint16_t>(circ.getNqubits()), 5, false, true);
  EXPECT_TRUE(resultSat.isSat());
  EXPECT_TRUE(resultSat.front().isRydberg());
  for (const auto& q : resultSat.front().getQubits()) {
    EXPECT_GE(q.getX(), 0);
    EXPECT_LE(q.getX(), 3);
    EXPECT_GE(q.getY(), 0);
    EXPECT_LE(q.getY(), 7);
    EXPECT_GE(q.getC(), 0);
    EXPECT_LE(q.getC(), 2);
    EXPECT_GE(q.getR(), 0);
    EXPECT_LE(q.getR(), 3);
    EXPECT_GE(q.getH(), -2);
    EXPECT_LE(q.getH(), 2);
    EXPECT_GE(q.getV(), -2);
    EXPECT_LE(q.getV(), 2);
    EXPECT_NO_THROW(std::ignore = q.isAOD());
  }
  for (const auto& g : resultSat.front().getGates()) {
    EXPECT_TRUE(std::find(pairs.cbegin(), pairs.cend(), g.getQubits()) !=
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
  const auto&       circ = qc::QuantumComputation::fromQASM(qasm);
  // create solver
  na::NASolver solver;
  solver.init(3, 7, 2, 3, 2, 2, 2, 2, 0, 7);
  // get operations for solver
  const auto& pairs =
      na::SolverFactory::getOpsForSolver(circ, {qc::Z, 1}, true);
  // solve
  const auto result = solver.solve(
      pairs, static_cast<std::uint16_t>(circ.getNqubits()), 3, true, false);
  EXPECT_TRUE(result.isSat());
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
  const auto&       circ = qc::QuantumComputation::fromQASM(qasm);
  // create solver
  na::NASolver solver;
  solver.init(3, 7, 2, 3, 2, 2, 2, 2, 2, 4);
  // get operations for solver
  const auto& pairs =
      na::SolverFactory::getOpsForSolver(circ, {qc::Z, 1}, true);
  // solve
  const auto result = solver.solve(
      pairs, static_cast<std::uint16_t>(circ.getNqubits()), 5, 2, false, true);
  EXPECT_TRUE(result.isSat());
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
  const auto&       circ = qc::QuantumComputation::fromQASM(qasm);
  // create solver
  na::NASolver solver;
  solver.init(3, 7, 2, 3, 2, 2, 2, 2, 2, 4);
  // get operations for solver
  const auto& pairs =
      na::SolverFactory::getOpsForSolver(circ, {qc::Z, 1}, true);
  // solve
  const auto result = solver.solve(
      pairs, static_cast<std::uint16_t>(circ.getNqubits()), 3, false, true);
  EXPECT_FALSE(result.isSat());
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
  const auto&       circ = qc::QuantumComputation::fromQASM(qasm);
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
