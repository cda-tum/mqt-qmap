#pragma once

#include "Solver.hpp"
#include "ir/QuantumComputation.hpp"
#include "na/Architecture.hpp"
#include "na/NADefinitions.hpp"

#include <utility>
#include <vector>

namespace na {
class SolverFactory {
public:
  [[nodiscard]] static auto create(const Architecture& arch) -> NASolver;

  [[nodiscard]] static auto getOpsForSolver(
      const qc::QuantumComputation& circ, FullOpType opType,
      bool quiet = false) -> std::vector<std::pair<unsigned int, unsigned int>>;
};
} // namespace na
