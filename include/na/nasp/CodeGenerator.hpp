# pragma once

#include "QuantumComputation.hpp"
#include "Solver.hpp"
#include "na/NAComputation.hpp"
#include "na/NADefinitions.hpp"
#include <cstddef>
#include <cstdint>

namespace na {
using namespace qc;

class CodeGenerator {
private:
  static auto coordFromDiscrete(
      std::size_t x, std::size_t y, std::int32_t h, std::int32_t v,
      std::size_t maxHOffset, std::size_t maxVOffset,
      std::size_t minEntanglingY, std::size_t maxEntanglingY) -> Point;

public:
  [[nodiscard]] static auto generate(
      const QuantumComputation& input, const NASolver::Result& result,
      std::size_t maxHOffset, std::size_t maxVOffset,
      std::size_t minEntanglingY, std::size_t maxEntanglingY) ->
    NAComputation;
};
} // namespace na