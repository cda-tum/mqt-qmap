#include "na/azac/NoReuseAnalyzer.hpp"

#include "Definitions.hpp"
#include "na/azac/Architecture.hpp"

#include <iostream>
#include <sstream>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na {
NoReuseAnalyzer::NoReuseAnalyzer(const Architecture&,
                                 const nlohmann::json& config) {
  if (const auto& configIt = config.find("no_reuse_analyzer");
      configIt != config.end() && configIt->is_object()) {
    for (const auto& [key, value] : configIt.value().items()) {
      std::ostringstream oss;
      oss << "[WARN] Configuration for NoReuseAnalyzer contains an unknown "
             "key: "
          << key << ". Ignoring.\n";
      std::cout << oss.str();
    }
  }
}
auto NoReuseAnalyzer::analyzeReuse(
    const std::vector<std::vector<std::pair<qc::Qubit, qc::Qubit>>>&
        twoQubitGateLayers) -> std::vector<std::unordered_set<qc::Qubit>> {
  return std::vector<std::unordered_set<qc::Qubit>>(twoQubitGateLayers.size() -
                                                    1);
}
}; // namespace na
