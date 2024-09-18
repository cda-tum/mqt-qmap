#pragma once

#include "ir/QuantumComputation.hpp"
#include "Solver.hpp"
#include "na/NAComputation.hpp"
#include "na/NADefinitions.hpp"

#include <cstdint>

namespace na {
using namespace qc;

class CodeGenerator {
private:
  static auto coordFromDiscrete(std::int32_t x, std::int32_t y, std::int32_t h,
                                std::int32_t v, std::int32_t maxHOffset,
                                std::int32_t maxVOffset,
                                std::int32_t minEntanglingY,
                                std::int32_t maxEntanglingY) -> Point;

public:
  [[nodiscard]] static auto
  generate(const QuantumComputation& input, const NASolver::Result& result,
           std::uint16_t maxHOffset, std::uint16_t maxVOffset,
           std::uint16_t minEntanglingY,
           std::uint16_t maxEntanglingY) -> NAComputation;
};
} // namespace na
