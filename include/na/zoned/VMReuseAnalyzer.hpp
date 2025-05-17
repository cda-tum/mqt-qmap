/*
 * Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
 * Copyright (c) 2025 Munich Quantum Software Company GmbH
 * All rights reserved.
 *
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License
 */

#pragma once

#include "ir/Definitions.hpp"
#include "na/zoned/Architecture.hpp"
#include "na/zoned/Types.hpp"

#include <cstddef>
#include <nlohmann/json.hpp>
#include <optional>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na::zoned {
/**
 * The class VMReuseAnalyzer implements the default reuse analysis for the
 * zoned neutral atom compiler that uses a bipartite maximum matching.
 */
class VMReuseAnalyzer {
  friend class VMReuseAnalyzerMaximumBipartiteMatchingTest_Direct_Test;
  friend class VMReuseAnalyzerMaximumBipartiteMatchingTest_Inverse_Test;
  friend class VMReuseAnalyzerMaximumBipartiteMatchingInvertedTest_Direct_Test;

public:
  /// The configuration of the VMReuseAnalyzer
  /// @note VMReuseAnalyzer does not have any configuration parameters.
  struct Config {
    template <typename BasicJsonType>
    friend void to_json(BasicJsonType& /* unused */,
                        const Config& /* unused */) {}
    template <typename BasicJsonType>
    friend void from_json(const BasicJsonType& /* unused */,
                          Config& /* unused */) {}
  };
  /**
   * Create a new VMReuseAnalyzer.
   * @note Both parameters are unused. Hence, the constructor does nothing
   * and the function @ref analyzeReuse is a static function.
   */
  VMReuseAnalyzer(const Architecture& /* unused */,
                  const Config& /* unused */) {}
  /// Analyze the reuse of qubits in the given two-qubit gate layers.
  [[nodiscard]] static auto
  analyzeReuse(const std::vector<TwoQubitGateLayer>& twoQubitGateLayers)
      -> std::vector<std::unordered_set<qc::Qubit>>;

private:
  /// Computes a maximum matching in a bipartite graph
  /// @note implemented pseudocode from
  /// https://epubs.siam.org/doi/pdf/10.1137/0202019?download=true
  [[nodiscard]] static auto maximumBipartiteMatching(
      const std::vector<std::vector<std::size_t>>& sparseMatrix,
      bool inverted = false) -> std::vector<std::optional<std::size_t>>;
};
} // namespace na::zoned
