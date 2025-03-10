#pragma once

#include "Definitions.hpp"
#include "na/azac/Architecture.hpp"

#include <nlohmann/json_fwd.hpp>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na {
/**
 * The class BMReuseAnalyzer implements the default reuse analysis for the
 * zoned neutral atom compiler that uses a bipartite maximum matching.
 */
class BMReuseAnalyzer {
public:
  /**
   * Create a new BMReuseAnalyzer.
   * @note Both parameters are unused. Hence, the constructor does nothing
   * and the function @ref analyzeReuse is a static function.
   */
  BMReuseAnalyzer(const Architecture& /* unused */,
                  const nlohmann::json& /* unused */) {}
  /// Analyze the reuse of qubits in the given two-qubit gate layers.
  [[nodiscard]] static auto
  analyzeReuse(const std::vector<std::vector<std::pair<qc::Qubit, qc::Qubit>>>&
                   twoQubitGateLayers)
      -> std::vector<std::unordered_set<qc::Qubit>>;
};
} // namespace na
