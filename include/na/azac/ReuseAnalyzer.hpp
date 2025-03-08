#pragma once

namespace na {
/**
 * The class ReuseAnalyzer implements the default reuse analysis for the
 * zoned neutral atom compiler.
 */
class ReuseAnalyzer {
protected:
  /**
   * Create a new ReuseAnalyzer.
   * @details Both parameters are unused. Hence, the constructor does nothing
   * and the function @ref analyzeReuse is a static function.
   */
  ReuseAnalyzer(const Architecture& /* unused */,
                const nlohmann::json& /* unused */) {}
  /// Analyze the reuse of qubits in the given two-qubit gate layers.
  [[nodiscard]] static auto
  analyzeReuse(const std::vector<std::vector<std::pair<qc::Qubit, qc::Qubit>>>&
                   twoQubitGateLayers)
      -> std::vector<std::unordered_set<qc::Qubit>> {
    std::vector<std::unordered_map<qc::Qubit, size_t>> usedQubitsInLayers;
    usedQubitsInLayers.reserve(twoQubitGateLayers.size());
    const auto& firstGateLayer = twoQubitGateLayers.front();
    auto& usedQubitsInFirstLayer = usedQubitsInLayers.emplace_back();
    for (size_t gateIdx = 0; gateIdx < firstGateLayer.size(); ++gateIdx) {
      const auto& gate = firstGateLayer[gateIdx];
      usedQubitsInFirstLayer[gate.first] = gateIdx;
      usedQubitsInFirstLayer[gate.second] = gateIdx;
    }
    std::vector<std::unordered_set<qc::Qubit>> reuseQubits;
    reuseQubits.reserve(twoQubitGateLayers.size());
    for (auto layer = twoQubitGateLayers.begin();
         layer != twoQubitGateLayers.end();) {
      const auto& twoQubitGatesInPreviousLayer = *layer;
      const auto& twoQubitGatesInCurrentLayer = *++layer;
      std::vector matrix(
          twoQubitGatesInCurrentLayer.size(),
          std::vector(twoQubitGatesInPreviousLayer.size(), false));
      const auto& usedQubitsInPreviousLayer = usedQubitsInLayers.back();
      auto& usedQubitsInCurrentLayer = usedQubitsInLayers.emplace_back();
      auto& reuseQubitsInCurrentLayer = reuseQubits.emplace_back();
      for (size_t gateIdx = 0; gateIdx < twoQubitGatesInCurrentLayer.size();
           ++gateIdx) {
        const auto& gate = twoQubitGatesInCurrentLayer[gateIdx];
        const auto& itFirst = usedQubitsInPreviousLayer.find(gate.first);
        const auto& itSecond = usedQubitsInPreviousLayer.find(gate.second);
        if (itFirst != usedQubitsInPreviousLayer.end()) {
          // If the both qubits of the gate are used in the previous layer also
          // by the identical gate, then both qubits can stay at their location
          // and be reused.
          if (itSecond != usedQubitsInPreviousLayer.end() &&
              itFirst->second == itSecond->second) {
            reuseQubitsInCurrentLayer.emplace(gate.first);
            reuseQubitsInCurrentLayer.emplace(gate.second);
          } else {
            matrix[gateIdx][itFirst->second] = true;
          }
        }
        if (itSecond != usedQubitsInPreviousLayer.end()) {
          matrix[gateIdx][itSecond->second] = true;
        }
        usedQubitsInCurrentLayer[gate.first] = gateIdx;
        usedQubitsInCurrentLayer[gate.second] = gateIdx;
      }
      std::vector sparseMatrix(matrix.size(), std::vector<std::size_t>{});
      for (std::size_t r = 0; r < matrix.size(); ++r) {
        for (std::size_t c = 0; c < matrix[r].size(); ++c) {
          if (matrix[r][c]) {
            sparseMatrix[r].emplace_back(c);
          }
        }
      }
      const auto& matching = maximumBipartiteMatching(sparseMatrix, true);
      for (std::size_t gateIdx = 0; gateIdx < matching.size(); ++gateIdx) {
        if (const auto& reuseGate = matching[gateIdx]; reuseGate) {
          const auto& gate = twoQubitGatesInCurrentLayer[gateIdx];
          if (usedQubitsInPreviousLayer.at(gate.first) == reuseGate) {
            reuseQubitsInCurrentLayer.emplace(gate.first);
          } else {
            assert(usedQubitsInPreviousLayer.at(gate.second) == reuseGate);
            reuseQubitsInCurrentLayer.emplace(gate.second);
          }
        }
      }
    }
    return reuseQubits;
  }
};
} // namespace na
