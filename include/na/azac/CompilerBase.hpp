#pragma once

#include "Definitions.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/StandardOperation.hpp"
#include "na/azac/Architecture.hpp"

#include <chrono>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <list>
#include <nlohmann/json.hpp>
#include <nlohmann/json_fwd.hpp>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na {

class CompilerBase {
protected:
  std::filesystem::path dir = "./result/";
  std::size_t n_q = 0; ///< number of qubits
  std::size_t n_g = 0; ///< number of gates
  Architecture architecture;
  struct Result {
    std::string name{};
    std::filesystem::path architecture_spec_path{};
    std::vector<std::string> instructions{};
    double runtime = 0;
  };
  Result result{};
  struct RuntimeAnalysis {
    std::chrono::microseconds scheduling;
    std::chrono::microseconds initial_placement;
    std::chrono::microseconds intermediate_placement;
  };
  RuntimeAnalysis runtime_analysis{};
  bool to_verify = true;
  /// trivial placement, i.e., place qubits in the order they appear in the
  /// circuit if false, a simulated annealing-based placement is chosen
  /// @see place_trivial
  bool trivial_placement = true;

public:
  enum class RoutingStrategy : std::uint8_t { MAXIMAL_IS_SORT };
  static std::string toString(const RoutingStrategy strategy) {
    static const std::map<RoutingStrategy, std::string> strategyToString{
        {RoutingStrategy::MAXIMAL_IS_SORT, "maximal_is_sort"}};
    const auto it = strategyToString.find(strategy);
    if (it == strategyToString.end()) {
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
    static const std::map<SchedulingStrategy, std::string> strategyToString{
        {SchedulingStrategy::ASAP, "asap"},
        {SchedulingStrategy::TRIVIAL, "trivial"}};
    const auto it = strategyToString.find(strategy);
    if (it == strategyToString.end()) {
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
                     std::unordered_set<qc::StandardOperation>>
      dict_g_1q_parent{};
  /// list of qubit placements for all layers
  std::vector<std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>
      qubit_mapping{};
  /// list of qubit lists that are reused in each layer
  std::vector<std::vector<qc::Qubit>> reuse_qubits;

public:
  explicit CompilerBase(std::ifstream& settingsIfs) {
    nlohmann::json settings_json{};
    settingsIfs >> settings_json;
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
      routing_strategy = settings_json["routing_strategy"];
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
      scheduling_strategy = settings_json["scheduling"];
    }
    if (settings_json.contains("arch_spec")) {
      result.architecture_spec_path =
          std::filesystem::path(settings_json["arch_spec"]);
    }
    architecture = Architecture(result.architecture_spec_path);
  }
  explicit CompilerBase(std::ifstream&& settingsIfs)
      : CompilerBase(std::move(settingsIfs)) {}
  explicit CompilerBase(const std::filesystem::path& settingsPath)
      : CompilerBase(std::ifstream(settingsPath)) {}
  explicit CompilerBase(const std::string& settingsFilename)
      : CompilerBase(std::filesystem::path(settingsFilename)) {}
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
    n_q = qc.getNqubits();
    /// array that stores the index of the last 2-qubit gate acting on each
    /// qubit
    std::vector<const std::pair<qc::Qubit, qc::Qubit>*> list_qubit_last_2q_gate(
        n_q, nullptr);

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
                std::unordered_set<qc::StandardOperation>{};
          }
          dict_g_1q_parent.at(list_qubit_last_2q_gate[qubit]).insert(stdop);
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

        selfg_q = []
        dict_g_1q_parent = {-1: []}
            n_single_qubit_gate = 0

            list_qubit_last_2q_gate = [-1 for i in range(0, self.n_q)]
            instruction = cz_circuit.data
            for ins in instruction:
                if ins.operation.num_qubits == 2:
                    offset = 0
                    if ins.qubits[0]._register != None:
                        offset = register_idx[ins.qubits[0]._register]
                    q0 = offset + ins.qubits[0]._index
                    offset = 0
                    if ins.qubits[1]._register != None:
                        offset = register_idx[ins.qubits[1]._register]
                    q1 = offset + ins.qubits[1]._index
                    list_qubit_last_2q_gate[q0] = len(self.g_q)
                    list_qubit_last_2q_gate[q1] = len(self.g_q)
                    if q0 < q1:
                        self.g_q.append([q0, q1])
                    else:
                        self.g_q.append([q1, q0])
                elif ins.operation.name != "measure" and ins.operation.name != "barrier":
                    offset = 0
                    if ins.qubits[0]._register != None:
                        offset = register_idx[ins.qubits[0]._register]
                    q0 = offset + ins.qubits[0]._index
                    if list_qubit_last_2q_gate[q0] not in self.dict_g_1q_parent:
                        self.dict_g_1q_parent[list_qubit_last_2q_gate[q0]] = []
                    self.dict_g_1q_parent[list_qubit_last_2q_gate[q0]].append((ins.operation.name, q0))
                    n_single_qubit_gate += 1

        self.n_g = len(self.g_q)
        self.g_s = tuple(['CRZ' for _ in range(self.n_g)])
  };

} // namespace na
