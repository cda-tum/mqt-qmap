#pragma once

#include "Definitions.hpp"
#include "na/azac/Architecture.hpp"

#include <cstddef>
#include <nlohmann/json_fwd.hpp>
#include <optional>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na::azac {
/**
 * The class VMReuseAnalyzer implements the default reuse analysis for the
 * zoned neutral atom compiler that uses a bipartite maximum matching.
 */
class VMReuseAnalyzer {
  friend class VMReuseAnalyzerMaximumBipartiteMatchingTest_Direct_Test;
  friend class VMReuseAnalyzerMaximumBipartiteMatchingTest_Inverse_Test;
  friend class VMReuseAnalyzerMaximumBipartiteMatchingInvertedTest_Direct_Test;

public:
  /**
   * Create a new VMReuseAnalyzer.
   * @note Both parameters are unused. Hence, the constructor does nothing
   * and the function @ref analyzeReuse is a static function.
   */
  VMReuseAnalyzer(const Architecture& /* unused */,
                  const nlohmann::json& config);
  /// Analyze the reuse of qubits in the given two-qubit gate layers.
  [[nodiscard]] static auto
  analyzeReuse(const std::vector<std::vector<std::array<qc::Qubit, 2>>>&
                   twoQubitGateLayers)
      -> std::vector<std::unordered_set<qc::Qubit>>;

private:
  /// Computes a maximum matching in a bipartite graph
  /// @note implemented pseudocode from
  /// https://epubs.siam.org/doi/pdf/10.1137/0202019?download=true
  [[nodiscard]] static auto maximumBipartiteMatching(
      const std::vector<std::vector<std::size_t>>& sparseMatrix,
      bool inverted = false) -> std::vector<std::optional<std::size_t>>;
};
} // namespace na
