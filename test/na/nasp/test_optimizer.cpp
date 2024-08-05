#include "Architecture.hpp"
#include "Configuration.hpp"
#include "Definitions.hpp"
#include "NAMapper.hpp"
#include "QuantumComputation.hpp"
#include "datastructures/Layer.hpp"
#include "na/NAComputation.hpp"
#include "na/operations/NAGlobalOperation.hpp"
#include "na/operations/NALocalOperation.hpp"
#include "na/operations/NAShuttlingOperation.hpp"
#include "nasp/Optimizer.hpp"
#include "nasp/Solver.hpp"
#include "nasp/SolverFactory.hpp"
#include "operations/CompoundOperation.hpp"
#include "operations/OpType.hpp"
#include "operations/StandardOperation.hpp"

#include <algorithm>
#include <cstddef>
#include <gtest/gtest.h>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <vector>

TEST(NASPOptimizer, Steane) {
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
  const auto& pairs = na::SolverFactory::getOpsForSolver(circ, {qc::Z, 1}, true);
  // Setup optimizer
  na::Optimizer optimizer;
  optimizer.setTimeout(std::chrono::minutes(1));
  optimizer.setMaxNSubProcs(2);
  optimizer.setInitialValue(3);
  optimizer.setMaxValue(17);
  optimizer.setObjectiveFunction([&solver, &circ, &pairs](const auto i) {
    return solver.solve(pairs, circ.getNqubits(), i, false, true);
  });
  optimizer.minimize();
  const auto resultOpt = optimizer.getExtremumOpt();
  EXPECT_TRUE(resultOpt.has_value());
  EXPECT_TRUE(resultOpt->isSat());
  EXPECT_EQ(resultOpt->numStages(), 4);
}