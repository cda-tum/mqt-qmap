#pragma once

#include "Definitions.hpp"
#include "na/azac/Architecture.hpp"

#include <nlohmann/json_fwd.hpp>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na::azac {
/**
 * The class NoReuseAnalyzer implements a fake reuse analyser that basically
 * disables the reuse analysis.
 */
class NoReuseAnalyzer {
protected:
  /**
   * Create a new NoReuseAnalyzer.
   * @note Both parameters are unused. Hence, the constructor does nothing
   * and the function @ref analyzeReuse is a static function.
   */
  NoReuseAnalyzer(const Architecture& /* unused */,
                  const nlohmann::json& config);
  /// Analyze the reuse of qubits in the given two-qubit gate layers.
  [[nodiscard]] static auto
  analyzeReuse(const std::vector<std::vector<std::array<qc::Qubit, 2>>>&
                   twoQubitGateLayers)
      -> std::vector<std::unordered_set<qc::Qubit>>;
};
} // namespace na::azac
