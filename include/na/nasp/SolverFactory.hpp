#pragma once

#include "Architecture.hpp"
#include "QuantumComputation.hpp"
#include "Solver.hpp"
#include <cstddef>
#include <utility>
#include <vector>
#include "na/NADefinitions.hpp"

namespace na {
class SolverFactory {
public:
  [[nodiscard]] static auto create(const Architecture& arch) -> NASolver;

  [[nodiscard]] static auto getOpsForSolver(const qc::QuantumComputation& circ,
                                            FullOpType opType,
                                            bool quiet = false)
    -> std::vector<std::pair<unsigned int, unsigned int> >;
};
} // namespace na