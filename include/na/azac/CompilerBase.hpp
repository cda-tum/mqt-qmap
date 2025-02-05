#pragma once

#include "Definitions.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/StandardOperation.hpp"
#include "na/NAComputation.hpp"
#include "na/azac/Architecture.hpp"
#include "na/azac/Utils.hpp"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <nlohmann/json_fwd.hpp>
#include <optional>
#include <ostream>
#include <queue>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na {

template <typename T> class CompilerBase {
protected:
  std::filesystem::path dir = "./result/";
  std::size_t n_q = 0; ///< number of qubits
  std::size_t n_g = 0; ///< number of gates
  Architecture architecture;
  struct Result {
    std::string name = "Untitled";
    std::string architecture_spec_path = "inline";
    nlohmann::json instructions;
    double runtime = 0;
  };
  Result result{};
  struct RuntimeAnalysis {
    std::chrono::microseconds scheduling;
    std::chrono::microseconds initial_placement;
    std::chrono::microseconds intermediate_placement;
    std::chrono::microseconds routing;
    std::chrono::microseconds total;
  };
  RuntimeAnalysis runtime_analysis{};
  bool to_verify = true;
  /// trivial placement, i.e., place qubits in the order they appear in the
  /// circuit if false, a simulated annealing-based placement is chosen
  /// @see place_trivial
  bool trivial_placement = true;

public:
  enum class RoutingStrategy : std::uint8_t { MAXIMAL_IS, MAXIMAL_IS_SORT };
  static std::string toString(const RoutingStrategy strategy) {
    switch (strategy) {
    case RoutingStrategy::MAXIMAL_IS:
      return "maximal_is";
    case RoutingStrategy::MAXIMAL_IS_SORT:
    default:
      return "maximal_is_sort";
    }
  }
  static RoutingStrategy toRoutingStrategy(const std::string& strategy) {
    static const std::unordered_map<std::string, RoutingStrategy>
        stringToStrategy{{"maximal_is_sort", RoutingStrategy::MAXIMAL_IS_SORT},
                         {"maximal_is", RoutingStrategy::MAXIMAL_IS}};
    const auto it = stringToStrategy.find(strategy);
    if (it == stringToStrategy.end()) {
      throw std::invalid_argument("Unknown routing strategy");
    }
    return it->second;
  }
  friend std::ostream& operator<<(std::ostream& os, const RoutingStrategy rs) {
    os << toString(rs);
    return os;
  }

protected:
  RoutingStrategy routing_strategy = RoutingStrategy::MAXIMAL_IS_SORT;

public:
  enum class SchedulingStrategy : std::uint8_t {
    ASAP, ///< Schedule gates in groups of commutative gates as soon as possible
    TRIVIAL ///< Schedule gates in the order they appear in the circuit, esp.,
            ///< every group contains only one gate
  };
  static std::string toString(const SchedulingStrategy strategy) {
    switch (strategy) {
    case SchedulingStrategy::ASAP:
      return "asap";
    case SchedulingStrategy::TRIVIAL:
    default:
      return "trivial";
    }
  }
  static SchedulingStrategy toSchedulingStrategy(const std::string& strategy) {
    static const std::unordered_map<std::string, SchedulingStrategy>
        stringToStrategy{{"asap", SchedulingStrategy::ASAP},
                         {"trivial", SchedulingStrategy::TRIVIAL}};
    const auto it = stringToStrategy.find(strategy);
    if (it == stringToStrategy.end()) {
      throw std::invalid_argument("Unknown scheduling strategy");
    }
    return it->second;
  }
  friend std::ostream& operator<<(std::ostream& os,
                                  const SchedulingStrategy ss) {
    os << toString(ss);
    return os;
  }

protected:
  SchedulingStrategy scheduling_strategy = SchedulingStrategy::ASAP;
  bool dynamic_placement = true;
  /// initial mapping of qubits to SLM sites, if this is not given either a
  /// trivial placement is chosen, see @ref trivial_placement, or a simulated
  /// annealing-based placement is chosen
  std::optional<std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>
      given_initial_mapping = std::nullopt;
  /// Mind the dependency between gates, i.e., do not allow changing there order
  /// if they are commutative
  bool hasDependency = true;
  bool l2 = false;
  bool use_window = true;
  std::size_t window_size = 0;
  bool reuse = true;
  std::size_t common1q = 0;
  /// list of 2-qubit CZ gates as a list of pairs of qubits
  std::vector<std::pair<qc::Qubit, qc::Qubit>> g_q{};
  /// map that stores the 1-qubit gates that act on a qubit after the
  /// respective 2-qubit gate
  std::unordered_map<const std::pair<qc::Qubit, qc::Qubit>*,
                     std::vector<qc::StandardOperation>>
      dict_g_1q_parent{};
  /// list of qubit placements for all layers
  std::vector<std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>
      qubit_mapping{};
  /// list of qubit lists that are reused in each layer
  std::vector<std::vector<qc::Qubit>> reuse_qubits;
  /// list of 2-qubit gates that are executed in each layer
  /// @note Computed by the scheduler
  std::vector<std::vector<const std::pair<qc::Qubit, qc::Qubit>*>>
      gate_scheduling{};
  /// list of 1-qubit gates that are executed in each layer
  /// @note Computed by the scheduler
  std::vector<std::vector<const qc::StandardOperation*>> gate_1q_scheduling{};
  std::vector<std::vector<std::size_t>> reuse_qubit{};

public:
  CompilerBase() = default;
  explicit CompilerBase(const std::string& filename)
      : CompilerBase(std::filesystem::path(filename)) {}
  explicit CompilerBase(const std::filesystem::path& path)
      : CompilerBase(std::ifstream(path)) {}
  explicit CompilerBase(std::istream& is)
      : CompilerBase(std::move(is)) {}
  explicit CompilerBase(std::ifstream&& is) {
    loadSettings(is);
  }
  auto loadSettings(const std::string& filename) -> void {
    loadSettings(std::filesystem::path(filename));
  }
  auto loadSettings(const std::filesystem::path& filepath) -> void {
    loadSettings(std::ifstream(filepath));
  }
  auto loadSettings(std::istream& is) -> void {
    loadSettings(std::move(is));
  }
  auto loadSettings(std::istream&& is) -> void {
    nlohmann::json settings_json{};
    is >> settings_json;
    if (settings_json.contains("name")) {
      result.name = settings_json["name"];
    }
    if (settings_json.contains("dir")) {
      dir = std::filesystem::path(settings_json["dir"]);
    }
    if (settings_json.contains("dependency")) {
      hasDependency = settings_json["dependency"];
    }
    if (settings_json.contains("routing_strategy")) {
      routing_strategy = toRoutingStrategy(settings_json["routing_strategy"]);
    }
    if (settings_json.contains("trivial_placement")) {
      trivial_placement = settings_json["trivial_placement"];
    }
    if (settings_json.contains("dynamic_placement")) {
      dynamic_placement = settings_json["dynamic_placement"];
    }
    if (settings_json.contains("use_window")) {
      use_window = settings_json["use_window"];
    }
    if (settings_json.contains("use_verifier")) {
      to_verify = settings_json["use_verifier"];
    }
    if (settings_json.contains("window_size")) {
      window_size = settings_json["window_size"];
    }
    if (settings_json.contains("l2")) {
      l2 = settings_json["l2"];
    }
    if (settings_json.contains("reuse")) {
      reuse = settings_json["reuse"];
    }
    if (settings_json.contains("scheduling")) {
      scheduling_strategy = toSchedulingStrategy(settings_json["scheduling"]);
    }
    if (settings_json.contains("arch_spec")) {
      if (settings_json["arch_spec"].is_string()) {
        if (std::filesystem::exists(settings_json["arch_spec"])) {
          result.architecture_spec_path =
              std::filesystem::path(settings_json["arch_spec"]);
          architecture = Architecture(result.architecture_spec_path);
        } else {
          throw std::invalid_argument("Architecture specification is missing");
        }
      } else if (settings_json["arch_spec"].is_object()) {
        architecture = Architecture(settings_json["arch_spec"]);
      } else {
        throw std::invalid_argument("Architecture specification is invalid");
      }
    } else {
      throw std::invalid_argument("Architecture specification is missing");
    }
  }
  [[nodiscard]] std::string toString() const {
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
      ss << "[INFO]                                         : reuse aware\n";
    } else {
      ss << "[INFO]                                         : no reuse\n";
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

  friend std::ostream& operator<<(std::ostream& os, const CompilerBase& cb) {
    os << cb.toString();
    return os;
  }

  auto setProgram(const qc::QuantumComputation& qc) {
    g_q.clear();
    n_q = qc.getNqubits();
    dict_g_1q_parent.emplace(nullptr, std::vector<qc::StandardOperation>{});
    /// array that stores the index of the last 2-qubit gate acting on each
    /// qubit
    std::vector<const std::pair<qc::Qubit, qc::Qubit>*> list_qubit_last_2q_gate(
        n_q, nullptr);
    std::size_t n_single_qubit_gate = 0;
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
          if (dict_g_1q_parent.find(list_qubit_last_2q_gate[qubit]) ==
              dict_g_1q_parent.cend()) {
            dict_g_1q_parent[list_qubit_last_2q_gate[qubit]] =
                std::vector<qc::StandardOperation>{};
          }
          dict_g_1q_parent[list_qubit_last_2q_gate[qubit]].emplace_back(stdop);
          ++n_single_qubit_gate;
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
              << n_single_qubit_gate << "\n";
  }

  auto solve(bool save_file = true) -> Result {
    // member to hold intermediate results
    gate_scheduling.clear();
    gate_1q_scheduling.clear();
    qubit_mapping.clear();
    reuse_qubit.clear();
    qubit_mapping.clear();

    std::cout << "[INFO] ZAC: A compiler for neutral atom-based compute-store "
                 "architecture\n";
    std::cout << *this;
    // todo: check if the program input is valid, i.e., #q < #p
    const auto t_s = std::chrono::system_clock::now();
    // gate scheduling with graph coloring
    std::cout << "[INFO] ZAC: Run scheduling\n";
    static_cast<T*>(this)->scheduling();

    if (reuse) {
      collect_reuse_qubit();
    } else {
      reuse_qubit.reserve(gate_scheduling.size());
      for (std::size_t i = 0; i < gate_scheduling.size(); ++i) {
        reuse_qubit.emplace_back(std::vector<std::size_t>{});
      }
    }

    static_cast<T*>(this)->place_qubit_initial();
    std::cout << "[INFO]               Time for initial placement: "
              << runtime_analysis.initial_placement << "s\n";
    static_cast<T*>(this)->place_qubit_intermedeiate();
    std::cout << "[INFO]               Time for intermediate placement: "
              << runtime_analysis.intermediate_placement << "s\n";
    static_cast<T*>(this)->route_qubit();
    runtime_analysis.total = std::chrono::system_clock::now() - t_s;
    std::cout << "[INFO]               Time for routing: "
              << runtime_analysis.routing << "s\n";
    std::cout << "[INFO] ZAC: Toal Time: " << runtime_analysis.total << "s\n";
    if (save_file) {
      if (dir.empty()) {
        dir = "./result/";
      }
      const auto code_filename = dir / "code" / (result.name + "_code.na");
      std::ofstream code_ofs(code_filename);
      code_ofs << result.instructions;

      const auto result_json_filename =
          dir / "time" / (result.name + "_time.json");
      std::ofstream result_json_ofs(result_json_filename);
      nlohmann::json result_json{};
      result_json["scheduling"] = runtime_analysis.scheduling;
      result_json["initial_placement"] = runtime_analysis.initial_placement;
      result_json["intermediate_placement"] =
          runtime_analysis.intermediate_placement;
      result_json["routing"] = runtime_analysis.routing;
      result_json["total"] = runtime_analysis.total;
      result_json_ofs << result_json;
    }

    if (to_verify) {
      std::cout << "[INFO] ZAC: Start Verification\n";
      throw std::invalid_argument("Verification is not implemented yet");
    }

    return result;
  }

  /// collect qubits that will remain in Rydberg zone between two Rydberg stages
  auto collect_reuse_qubit() -> void {
    reuse_qubit.clear();
    std::vector qubit_is_used(gate_scheduling.size(),
                              std::vector<std::size_t>(n_q, -1));
    for (std::size_t gate_idx = 0; gate_idx < gate_scheduling.front().size();
         ++gate_idx) {
      const auto& gate = gate_scheduling.front()[gate_idx];
      qubit_is_used[0][gate->first] = gate_idx;
      qubit_is_used[0][gate->second] = gate_idx;
    }
    for (std::size_t i = 1; i < gate_scheduling.size(); ++i) {
      reuse_qubit.emplace_back();
      std::vector matrix(gate_scheduling[i].size(),
                         std::vector(gate_scheduling[i - 1].size(), false));
      for (std::size_t gate_idx = 0; gate_idx < gate_scheduling[i].size();
           ++gate_idx) {
        const auto& gate = gate_scheduling[i][gate_idx];
        if (qubit_is_used[i - 1][gate->first] != -1 &&
            qubit_is_used[i - 1][gate->first] ==
                qubit_is_used[i - 1][gate->second]) {
          reuse_qubit.back().emplace_back(gate->first);
          reuse_qubit.back().emplace_back(gate->second);
        } else {
          if (qubit_is_used[i - 1][gate->first] > -1) {
            matrix[gate_idx][qubit_is_used[i - 1][gate->first]] = true;
          }
          if (qubit_is_used[i - 1][gate->second] > -1) {
            matrix[gate_idx][qubit_is_used[i - 1][gate->second]] = true;
          }
        }
        qubit_is_used[i][gate->first] = gate_idx;
        qubit_is_used[i][gate->second] = gate_idx;
      }
      std::vector sparse_matrix(matrix.size(), std::vector<std::size_t>{});
      for (std::size_t r = 0; r < matrix.size(); ++r) {
        for (std::size_t c = 0; c < matrix[r].size(); ++c) {
          if (matrix[r][c]) {
            sparse_matrix[r].emplace_back(c);
          }
        }
      }
      const auto& matching = maximumBipartiteMatching(sparse_matrix, true);
      for (std::size_t gate_idx = 0; gate_idx < matching.size(); ++gate_idx) {
        const auto reuse_gate = matching[gate_idx];
        if (reuse_gate != -1) {
          const auto& gate = gate_scheduling[i][gate_idx];
          if (qubit_is_used[i - 1][gate->first] == reuse_gate) {
            reuse_qubit.back().emplace_back(gate->first);
          }
          if (qubit_is_used[i - 1][gate->second] == reuse_gate) {
            reuse_qubit.back().emplace_back(gate->second);
          }
        }
      }
    }
    reuse_qubit.emplace_back();
  }
};

} // namespace na
