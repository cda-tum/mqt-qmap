#include "na/azac/CompilerBase.hpp"

#include "Definitions.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/StandardOperation.hpp"
#include "na/azac/Utils.hpp"

#include <unordered_map>

namespace na {
auto CompilerBase::loadSettings(std::istream&& is) -> void {
  nlohmann::json settingsJson{};
  std::move(is) >> settingsJson;
  if (settingsJson.contains("name")) {
    result.name = settingsJson["name"];
  }
  if (settingsJson.contains("dir")) {
    dir = std::filesystem::path(settingsJson["dir"]);
  }
  if (settingsJson.contains("dependency")) {
    hasDependency = settingsJson["dependency"];
  }
  if (settingsJson.contains("routing_strategy")) {
    routingStrategy = toRoutingStrategy(settingsJson["routing_strategy"]);
  }
  if (settingsJson.contains("trivial_placement")) {
    trivialPlacement = settingsJson["trivial_placement"];
  }
  if (settingsJson.contains("dynamic_placement")) {
    dynamicPlacement = settingsJson["dynamic_placement"];
  }
  if (settingsJson.contains("use_window")) {
    useWindow = settingsJson["use_window"];
  }
  if (settingsJson.contains("use_verifier")) {
    toVerify = settingsJson["use_verifier"];
  }
  if (settingsJson.contains("window_size")) {
    windowSize = settingsJson["window_size"];
  }
  if (settingsJson.contains("l2")) {
    l2 = settingsJson["l2"];
  }
  if (settingsJson.contains("reuse")) {
    reuse = settingsJson["reuse"];
  }
  if (settingsJson.contains("scheduling")) {
    schedulingStrategy = toSchedulingStrategy(settingsJson["scheduling"]);
  }
  if (settingsJson.contains("arch_spec")) {
    if (settingsJson["arch_spec"].is_string()) {
      const auto& path = std::filesystem::path(settingsJson["arch_spec"]);
      if (exists(path)) {
        result.architectureSpecPath = path;
        architecture = Architecture(result.architectureSpecPath);
      } else {
        throw std::filesystem::filesystem_error(
            "File with architecture specification not found", path,
            std::make_error_code(std::errc::no_such_file_or_directory));
      }
    } else if (settingsJson["arch_spec"].is_object()) {
      architecture = Architecture(settingsJson["arch_spec"]);
    } else {
      throw std::invalid_argument("Architecture specification is invalid");
    }
  } else {
    throw std::invalid_argument("Architecture specification is missing");
  }
}
std::string CompilerBase::toString() const {
  std::stringstream ss{};
  ss << "[INFO] AZAC: Settings\n";
  ss << "[INFO]           Result directory: " << dir << "\n";
  if (hasDependency) {
    ss << "[INFO]           Scheduling strategy: " << schedulingStrategy
       << "\n";
  } else {
    ss << "[INFO]           Scheduling strategy: edge coloring\n";
  }
  if (trivialPlacement) {
    ss << "[INFO]           Placement strategy: trivial placement\n";
  } else if (givenInitialMapping) {
    ss << "[INFO]           Initial placement strategy: user-defined "
          "placement\n";
  } else if (l2) {
    ss << "[INFO]           Initial placement strategy: SA-based placement "
          "with L2 distance model\n";
  } else {
    ss << "[INFO]           Initial placement strategy: SA-based placement "
          "with Euclidean distance model\n";
  }
  if (dynamicPlacement) {
    ss << "[INFO]           Intermediate placement strategy: minimal "
          "weighted matching\n";
  } else {
    ss << "[INFO]           Intermediate placement strategy: return to "
          "initial mapping\n";
  }
  if (reuse) {
    ss << "[INFO]                                            reuse aware\n";
  } else {
    ss << "[INFO]                                            no reuse\n";
  }
  ss << "[INFO]           Routing strategy: " << routingStrategy;
  if (useWindow) {
    ss << " with window size " << windowSize << "\n";
  } else {
    ss << " without window\n";
  }
  if (toVerify) {
    ss << "[INFO]           Verifier: enable\n";
  } else {
    ss << "[INFO]           Verifier: disable\n";
  }
  return ss.str();
}
auto CompilerBase::collectReuseQubit() -> void {
  reuseQubits.clear();
  std::vector qubitIsUsed(
      gateScheduling.size(),
      std::vector<std::optional<std::size_t>>(nQubits, std::nullopt));
  for (std::size_t gateIdx = 0; gateIdx < gateScheduling.front().size();
       ++gateIdx) {
    const auto& gate = gateScheduling.front()[gateIdx];
    qubitIsUsed[0][gate->first] = gateIdx;
    qubitIsUsed[0][gate->second] = gateIdx;
  }
  for (std::size_t i = 1; i < gateScheduling.size(); ++i) {
    reuseQubits.emplace_back();
    std::vector matrix(gateScheduling[i].size(),
                       std::vector(gateScheduling[i - 1].size(), false));
    for (std::size_t gateIdx = 0; gateIdx < gateScheduling[i].size();
         ++gateIdx) {
      const auto& gate = gateScheduling[i][gateIdx];
      if (qubitIsUsed[i - 1][gate->first] != -1 &&
          qubitIsUsed[i - 1][gate->first] == qubitIsUsed[i - 1][gate->second]) {
        reuseQubits.back().emplace(gate->first);
        reuseQubits.back().emplace(gate->second);
      } else {
        if (qubitIsUsed[i - 1][gate->first]) {
          matrix[gateIdx][*qubitIsUsed[i - 1][gate->first]] = true;
        }
        if (qubitIsUsed[i - 1][gate->second]) {
          matrix[gateIdx][*qubitIsUsed[i - 1][gate->second]] = true;
        }
      }
      qubitIsUsed[i][gate->first] = gateIdx;
      qubitIsUsed[i][gate->second] = gateIdx;
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
      if (const auto reuseGate = matching[gateIdx]; reuseGate) {
        const auto& gate = gateScheduling[i][gateIdx];
        if (qubitIsUsed[i - 1][gate->first] == reuseGate) {
          reuseQubits.back().emplace(gate->first);
        }
        if (qubitIsUsed[i - 1][gate->second] == reuseGate) {
          reuseQubits.back().emplace(gate->second);
        }
      }
    }
  }
  reuseQubits.emplace_back();
}
std::string CompilerBase::toString(const RoutingStrategy strategy) {
  switch (strategy) {
  case RoutingStrategy::MaximalIs:
    return "maximalis";
  case RoutingStrategy::MaximalIsSort:
  default:
    return "maximalis_sort";
  }
}
CompilerBase::RoutingStrategy
CompilerBase::toRoutingStrategy(const std::string& strategy) {
  static const std::unordered_map<std::string, RoutingStrategy>
      STRINGTostrategy{{"maximalis_sort", RoutingStrategy::MaximalIsSort},
                       {"maximalis", RoutingStrategy::MaximalIs}};
  const auto it = STRINGTostrategy.find(strategy);
  if (it == STRINGTostrategy.end()) {
    throw std::invalid_argument("Unknown routing strategy");
  }
  return it->second;
}
std::string CompilerBase::toString(const SchedulingStrategy strategy) {
  switch (strategy) {
  case SchedulingStrategy::ASAP:
    return "asap";
  case SchedulingStrategy::TRIVIAL:
  default:
    return "trivial";
  }
}
auto CompilerBase::toSchedulingStrategy(const std::string& strategy)
    -> CompilerBase::SchedulingStrategy {
  static const std::unordered_map<std::string, SchedulingStrategy>
      STRINGTostrategy{{"asap", SchedulingStrategy::ASAP},
                       {"trivial", SchedulingStrategy::TRIVIAL}};
  const auto it = STRINGTostrategy.find(strategy);
  if (it == STRINGTostrategy.end()) {
    throw std::invalid_argument("Unknown scheduling strategy");
  }
  return it->second;
}
auto CompilerBase::setProgram(const qc::QuantumComputation& qc) -> void {
  twoQubitGates.clear();
  dictG1QParent.emplace(nullptr, std::vector<qc::StandardOperation>{});
  std::cout << "[INFO] AZAC: Read in Circuit\n";
  nQubits = qc.getNqubits();
  /// array that stores the index of the last 2-qubit gate acting on each
  /// qubit
  std::vector<const std::pair<qc::Qubit, qc::Qubit>*> listQubitLast2qGate(
      nQubits, nullptr);
  std::size_t nSingleQubitGate = 0;
  for (const auto& op : qc) {
    if (op->isStandardOperation()) {
      const auto& stdop = dynamic_cast<qc::StandardOperation&>(*op);
      if (stdop.getNcontrols() == 1 && stdop.getNtargets() == 1 &&
          stdop.getType() == qc::Z) {
        const auto& usedQubits = stdop.getUsedQubits();
        const std::vector qubits(usedQubits.cbegin(), usedQubits.cend());
        listQubitLast2qGate[qubits[0]] = listQubitLast2qGate[qubits[1]] =
            twoQubitGates
                .emplace_back(std::make_unique<std::pair<qc::Qubit, qc::Qubit>>(
                    std::min(qubits[0], qubits[1]),
                    std::max(qubits[0], qubits[1])))
                .get();
      } else if (stdop.getNcontrols() == 0 && stdop.getNtargets() == 1) {
        const auto& qubit = *stdop.getTargets().cbegin();
        if (dictG1QParent.find(listQubitLast2qGate[qubit]) ==
            dictG1QParent.cend()) {
          dictG1QParent[listQubitLast2qGate[qubit]] =
              std::vector<qc::StandardOperation>{};
        }
        dictG1QParent[listQubitLast2qGate[qubit]].emplace_back(stdop);
        ++nSingleQubitGate;
      } else {
        std::stringstream ss{};
        ss << "Standard operation ";
        stdop.print(ss, stdop.getNqubits());
        ss << " is not supported";
        throw std::invalid_argument(ss.str());
      }
    } else {
      throw std::invalid_argument("Non-standard operation is not supported");
    }
  }
  nTwoQubitGates = twoQubitGates.size();
  std::cout << "[INFO]           Number of qubits: " << nQubits << "\n";
  std::cout << "[INFO]           Number of two-qubit gates: " << nTwoQubitGates
            << "\n";
  std::cout << "[INFO]           Number of single-qubit gates: "
            << nSingleQubitGate << "\n";
}
} // namespace na
