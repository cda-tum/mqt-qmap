#pragma once

#include "Definitions.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/StandardOperation.hpp"
#include "na/azac/Architecture.hpp"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <nlohmann/json_fwd.hpp>
#include <optional>
#include <ostream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace na {

class CompilerBase {
public:
  enum class RoutingStrategy : std::uint8_t { MaximalIs, MaximalIsSort };
  enum class SchedulingStrategy : std::uint8_t {
    ASAP, ///< Schedule gates in groups of commutative gates as soon as possible
    TRIVIAL ///< Schedule gates in the order they appear in the circuit, esp.,
    ///< every group contains only one gate
  };
  struct Result {
    std::string name = "Untitled";
    std::string architecture_spec_path = "inline";
    nlohmann::json instructions;
    double runtime = 0;
  };
  struct RuntimeAnalysis {
    std::chrono::microseconds scheduling;
    std::chrono::microseconds initial_placement;
    std::chrono::microseconds intermediate_placement;
    std::chrono::microseconds routing;
    std::chrono::microseconds total;
  };

protected:
  std::filesystem::path dir = "./result/";
  std::size_t n_q = 0; ///< number of qubits
  std::size_t n_g = 0; ///< number of gates
  Architecture architecture;
  Result result{};
  RuntimeAnalysis runtime_analysis{};
  bool to_verify = true;
  /// trivial placement, i.e., place qubits in the order they appear in the
  /// circuit if false, a simulated annealing-based placement is chosen
  /// @see place_trivial
  bool trivial_placement = true;
  RoutingStrategy routing_strategy = RoutingStrategy::MaximalIsSort;
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
  std::vector<std::unordered_set<std::size_t>> reuse_qubits{};
  /// list of 2-qubit gates that are executed in each layer
  /// @note Computed by the scheduler
  std::vector<std::vector<const std::pair<qc::Qubit, qc::Qubit>*>>
      gate_scheduling{};
  std::vector<std::vector<std::size_t>>
      gate_scheduling_idx{};
  /// list of 1-qubit gates that are executed in each layer
  /// @note Computed by the scheduler
  std::vector<std::vector<const qc::StandardOperation*>> gate_1q_scheduling{};

public:
  CompilerBase() = default;
  explicit CompilerBase(const std::string& filename)
      : CompilerBase(std::filesystem::path(filename)) {}
  explicit CompilerBase(const std::filesystem::path& path)
      : CompilerBase(std::ifstream(path)) {}
  explicit CompilerBase(std::istream& is)
      : CompilerBase(std::move(is)) {}
  explicit CompilerBase(std::istream&& is) {
    loadSettings(std::move(is));
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
  auto loadSettings(std::istream&& is) -> void;
  auto setProgram(const qc::QuantumComputation& qc);
  [[nodiscard]] auto toString() const -> std::string;
  friend auto operator<<(std::ostream& os, const CompilerBase& cb) -> std::ostream& {
    os << cb.toString();
    return os;
  }

protected:
  /// collect qubits that will remain in Rydberg zone between two Rydberg stages
  auto collect_reuse_qubit() -> void;
  static auto toString(RoutingStrategy strategy) -> std::string;
  static auto toRoutingStrategy(const std::string& strategy) -> RoutingStrategy;
  friend auto operator<<(std::ostream& os, const RoutingStrategy rs) -> std::ostream& {
    os << toString(rs);
    return os;
  }
  static auto toString(SchedulingStrategy strategy) -> std::string;
  static auto toSchedulingStrategy(const std::string& strategy) -> SchedulingStrategy;
  friend auto operator<<(std::ostream& os,
                                  const SchedulingStrategy ss) -> std::ostream& {
    os << toString(ss);
    return os;
  }
public:
  [[nodiscard]] auto get_dir() const -> const std::filesystem::path& { return dir; }
  [[nodiscard]] auto get_n_q() const -> std::size_t { return n_q; }
  [[nodiscard]] auto get_n_g() const -> std::size_t { return n_g; }
  [[nodiscard]] auto get_architecture() const -> const Architecture& {
    return architecture;
  }
  [[nodiscard]] auto get_result() -> Result& { return result; }
  [[nodiscard]] auto get_runtime_analysis() -> RuntimeAnalysis& {
    return runtime_analysis;
  }
  [[nodiscard]] auto is_to_verify() const -> bool { return to_verify; }
  [[nodiscard]] auto is_trivial_placement() const -> bool { return trivial_placement; }
  [[nodiscard]] auto get_routing_strategy() const -> RoutingStrategy {
    return routing_strategy;
  }
  [[nodiscard]] auto get_scheduling_strategy() const -> SchedulingStrategy {
    return scheduling_strategy;
  }
  [[nodiscard]] auto is_dynamic_placement() const -> bool { return dynamic_placement; }
  [[nodiscard]] auto get_given_initial_mapping() const -> const std::optional<
      std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>& {
    return given_initial_mapping;
  }
  [[nodiscard]] auto is_has_dependency() const -> bool { return hasDependency; }
  [[nodiscard]] auto is_l2() const -> bool { return l2; }
  [[nodiscard]] auto is_use_window() const -> bool { return use_window; }
  [[nodiscard]] auto get_window_size() const -> std::size_t { return window_size; }
  [[nodiscard]] auto is_reuse() const -> bool { return reuse; }
  [[nodiscard]] auto get_common1_q() const -> std::size_t { return common1q; }
  [[nodiscard]] auto get_g_q() const
      -> const std::vector<std::pair<qc::Qubit, qc::Qubit>>& {
    return g_q;
  }
  [[nodiscard]] auto get_dict_g_1q_parent() const
      -> const std::unordered_map<const std::pair<qc::Qubit, qc::Qubit>*,
                                  std::vector<qc::StandardOperation>>& {
    return dict_g_1q_parent;
  }
  [[nodiscard]] auto get_qubit_mapping() -> std::vector<
      std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>& {
    return qubit_mapping;
  }
  [[nodiscard]] auto get_reuse_qubits() const
      -> const std::vector<std::unordered_set<std::size_t>>& {
    return reuse_qubits;
  }
  [[nodiscard]] auto get_gate_scheduling()
      -> std::vector<std::vector<const std::pair<qc::Qubit, qc::Qubit>*>>& {
    return gate_scheduling;
  }
  [[nodiscard]] auto get_gate_scheduling_idx() const
      -> const std::vector<std::vector<std::size_t>>& {
    return gate_scheduling_idx;
  }
  [[nodiscard]] auto get_gate_1q_scheduling()
      -> std::vector<std::vector<const qc::StandardOperation*>>& {
    return gate_1q_scheduling;
  }
  auto set_dict_g_1q_parent(
      const std::unordered_map<const std::pair<qc::Qubit, qc::Qubit>*,
                               std::vector<qc::StandardOperation>>&
          new_dict_g_1q_parent) -> void {
    dict_g_1q_parent = new_dict_g_1q_parent;
  }
  auto set_qubit_mapping(
      const std::vector<
          std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>&
          new_qubit_mapping) -> void {
    qubit_mapping = new_qubit_mapping;
  }
  auto set_gate_scheduling_idx(
      const std::vector<std::vector<std::size_t>>&
          new_gate_scheduling_idx) -> void {
    gate_scheduling_idx = new_gate_scheduling_idx;
  }
  auto set_gate_1_q_scheduling(
      const std::vector<std::vector<const qc::StandardOperation*>>&
          new_gate_1q_scheduling) -> void {
    gate_1q_scheduling = new_gate_1q_scheduling;
  }
};

} // namespace na
