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
    routing_strategy = toRoutingStrategy(settingsJson["routing_strategy"]);
  }
  if (settingsJson.contains("trivial_placement")) {
    trivial_placement = settingsJson["trivial_placement"];
  }
  if (settingsJson.contains("dynamic_placement")) {
    dynamic_placement = settingsJson["dynamic_placement"];
  }
  if (settingsJson.contains("use_window")) {
    use_window = settingsJson["use_window"];
  }
  if (settingsJson.contains("use_verifier")) {
    to_verify = settingsJson["use_verifier"];
  }
  if (settingsJson.contains("window_size")) {
    window_size = settingsJson["window_size"];
  }
  if (settingsJson.contains("l2")) {
    l2 = settingsJson["l2"];
  }
  if (settingsJson.contains("reuse")) {
    reuse = settingsJson["reuse"];
  }
  if (settingsJson.contains("scheduling")) {
    scheduling_strategy = toSchedulingStrategy(settingsJson["scheduling"]);
  }
  if (settingsJson.contains("arch_spec")) {
    if (settingsJson["arch_spec"].is_string()) {
      const auto& path = std::filesystem::path(settingsJson["arch_spec"]);
      if (exists(path)) {
        result.architecture_spec_path = path;
        architecture = Architecture(result.architecture_spec_path);
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
  ss << "[INFO] ZAC: Setting\n";
  ss << "[INFO]           Result directory: " << dir << "\n";
  if (hasDependency) {
    ss << "[INFO]           Scheduling strategy: " << scheduling_strategy
       << "\n";
  } else {
    ss << "[INFO]           Scheduling strategy: edge coloring\n";
  }
  if (trivial_placement) {
    ss << "[INFO]           Placement strategy: trivial placement\n";
  } else if (given_initial_mapping) {
    ss << "[INFO]           Initial placement strategy: user-defined "
          "placement\n";
  } else if (l2) {
    ss << "[INFO]           Initial placement strategy: SA-based placement "
          "with L2 distance model\n";
  } else {
    ss << "[INFO]           Initial placement strategy: SA-based placement "
          "with Euclidean distance model\n";
  }
  if (dynamic_placement) {
    ss << "[INFO]           Intermediate placement strategy: minimal "
          "weighted matching\n";
  } else {
    ss << "[INFO]           Intermediate placement strategy: return to "
          "intial mapping\n";
  }
  if (reuse) {
    ss << "[INFO]                                            reuse aware\n";
  } else {
    ss << "[INFO]                                            no reuse\n";
  }
  ss << "[INFO]           Routing strategy: " << routing_strategy;
  if (use_window) {
    ss << " with window size " << window_size << "\n";
  } else {
    ss << " without window\n";
  }
  if (to_verify) {
    ss << "[INFO]           Verifier: enable\n";
  } else {
    ss << "[INFO]           Verifier: disable\n";
  }
  return ss.str();
}
auto CompilerBase::collect_reuse_qubit() -> void {
  reuse_qubits.clear();
  std::vector qubitIsUsed(
      gate_scheduling.size(),
      std::vector<std::optional<std::size_t>>(n_q, std::nullopt));
  for (std::size_t gateIdx = 0; gateIdx < gate_scheduling.front().size();
       ++gateIdx) {
    const auto& gate = gate_scheduling.front()[gateIdx];
    qubitIsUsed[0][gate->first] = gateIdx;
    qubitIsUsed[0][gate->second] = gateIdx;
  }
  for (std::size_t i = 1; i < gate_scheduling.size(); ++i) {
    reuse_qubits.emplace_back();
    std::vector matrix(gate_scheduling[i].size(),
                       std::vector(gate_scheduling[i - 1].size(), false));
    for (std::size_t gateIdx = 0; gateIdx < gate_scheduling[i].size();
         ++gateIdx) {
      const auto& gate = gate_scheduling[i][gateIdx];
      if (qubitIsUsed[i - 1][gate->first] != -1 &&
          qubitIsUsed[i - 1][gate->first] == qubitIsUsed[i - 1][gate->second]) {
        reuse_qubits.back().emplace(gate->first);
        reuse_qubits.back().emplace(gate->second);
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
      if (const auto reuseGate = matching[gateIdx]; reuseGate != -1) {
        const auto& gate = gate_scheduling[i][gateIdx];
        if (qubitIsUsed[i - 1][gate->first] == reuseGate) {
          reuse_qubits.back().emplace(gate->first);
        }
        if (qubitIsUsed[i - 1][gate->second] == reuseGate) {
          reuse_qubits.back().emplace(gate->second);
        }
      }
    }
  }
  reuse_qubits.emplace_back();
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
      STRING_TO_STRATEGY{{"maximalis_sort", RoutingStrategy::MaximalIsSort},
                       {"maximalis", RoutingStrategy::MaximalIs}};
  const auto it = STRING_TO_STRATEGY.find(strategy);
  if (it == STRING_TO_STRATEGY.end()) {
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
auto
CompilerBase::toSchedulingStrategy(const std::string& strategy) -> CompilerBase::SchedulingStrategy {
  static const std::unordered_map<std::string, SchedulingStrategy>
      STRING_TO_STRATEGY{{"asap", SchedulingStrategy::ASAP},
                       {"trivial", SchedulingStrategy::TRIVIAL}};
  const auto it = STRING_TO_STRATEGY.find(strategy);
  if (it == STRING_TO_STRATEGY.end()) {
    throw std::invalid_argument("Unknown scheduling strategy");
  }
  return it->second;
}
auto CompilerBase::setProgram(const qc::QuantumComputation& qc) {
  g_q.clear();
  n_q = qc.getNqubits();
  dict_g_1q_parent.emplace(nullptr, std::vector<qc::StandardOperation>{});
  /// array that stores the index of the last 2-qubit gate acting on each
  /// qubit
  std::vector<const std::pair<qc::Qubit, qc::Qubit>*> listQubitLast2qGate(
      n_q, nullptr);
  std::size_t nSingleQubitGate = 0;
  for (const auto& op : qc) {
    if (op->isStandardOperation()) {
      const auto& stdop = dynamic_cast<qc::StandardOperation&>(*op);
      if (stdop.getNcontrols() == 1 && stdop.getNtargets() == 1 &&
          stdop.getType() == qc::Z) {
        const auto& usedQubits = stdop.getUsedQubits();
        const std::vector qubits(usedQubits.cbegin(), usedQubits.cend());
        if (qubits[0] < qubits[1]) {
          g_q.emplace_back(qubits[0], qubits[1]);
        } else {
          g_q.emplace_back(qubits[1], qubits[0]);
        }
      } else if (stdop.getNcontrols() == 0 && stdop.getNtargets() == 1) {
        const auto& qubit = *stdop.getTargets().cbegin();
        if (dict_g_1q_parent.find(listQubitLast2qGate[qubit]) ==
            dict_g_1q_parent.cend()) {
          dict_g_1q_parent[listQubitLast2qGate[qubit]] =
              std::vector<qc::StandardOperation>{};
        }
        dict_g_1q_parent[listQubitLast2qGate[qubit]].emplace_back(stdop);
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
  n_g = g_q.size();
  std::cout << "[INFO]           number of qubits: " << n_q << "\n";
  std::cout << "[INFO]           number of two-qubit gates: " << n_g << "\n";
  std::cout << "[INFO]           number of single-qubit gates: "
            << nSingleQubitGate << "\n";
}
} // namespace na