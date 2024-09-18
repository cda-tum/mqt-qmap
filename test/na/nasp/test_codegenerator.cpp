#include "na/nasp/CodeGenerator.hpp"
#include "ir/QuantumComputation.hpp"
#include "na/nasp/Solver.hpp"
#include "na/nasp/SolverFactory.hpp"
#include "na/NAComputation.hpp"
#include "ir/operations/OpType.hpp"

#include <cstdint>
#include <gtest/gtest.h>
#include <string>

TEST(CodeGenerator, Generate) {
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
  const auto& comp = na::CodeGenerator::generate(circ, result, 2, 2, 2, 4);
  EXPECT_TRUE(comp.validateAODConstraints());
}
