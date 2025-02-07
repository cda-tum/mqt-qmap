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
    std::string architectureSpecPath = "inline";
    nlohmann::json instructions;
    double runtime = 0;
  };
  struct RuntimeAnalysis {
    std::chrono::microseconds scheduling;
    std::chrono::microseconds initialPlacement;
    std::chrono::microseconds intermediatePlacement;
    std::chrono::microseconds routing;
    std::chrono::microseconds total;
  };

protected:
  std::filesystem::path dir = "./result/";
  std::size_t nQubits = 0; ///< number of qubits
  std::size_t nTwoQubitGates = 0; ///< number of 2-qubit gates
  Architecture architecture;
  Result result{};
  RuntimeAnalysis runtimeAnalysis{};
  bool toVerify = true;
  /// trivial placement, i.e., place qubits in the order they appear in the
  /// circuit if false, a simulated annealing-based placement is chosen
  /// @see placeTrivial
  bool trivialPlacement = true;
  RoutingStrategy routingStrategy = RoutingStrategy::MaximalIsSort;
  SchedulingStrategy schedulingStrategy = SchedulingStrategy::ASAP;
  bool dynamicPlacement = true;
  /// initial mapping of qubits to SLM sites, if this is not given either a
  /// trivial placement is chosen, see @ref trivialPlacement, or a simulated
  /// annealing-based placement is chosen
  std::optional<std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>
      givenInitialMapping = std::nullopt;
  /// Mind the dependency between gates, i.e., do not allow changing there order
  /// if they are commutative
  bool hasDependency = true;
  bool l2 = false;
  bool useWindow = true;
  std::size_t windowSize = 0;
  bool reuse = true;
  std::size_t common1Q = 0;
  /// list of 2-qubit CZ gates as a list of pairs of qubits
  std::vector<std::pair<qc::Qubit, qc::Qubit>> twoQubitGates{};
  /// map that stores the 1-qubit gates that act on a qubit after the
  /// respective 2-qubit gate
  std::unordered_map<const std::pair<qc::Qubit, qc::Qubit>*,
                     std::vector<qc::StandardOperation>>
      dictG1QParent{};
  /// list of qubit placements for all layers
  std::vector<std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>
      qubitMapping{};
  /// list of qubit lists that are reused in each layer
  std::vector<std::unordered_set<std::size_t>> reuseQubits{};
  /// list of 2-qubit gates that are executed in each layer
  /// @note Computed by the scheduler
  std::vector<std::vector<const std::pair<qc::Qubit, qc::Qubit>*>>
      gateScheduling{};
  std::vector<std::vector<std::size_t>>
      gateSchedulingIdx{};
  /// list of 1-qubit gates that are executed in each layer
  /// @note Computed by the scheduler
  std::vector<std::vector<const qc::StandardOperation*>> gate1QScheduling{};

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
  auto loadSettings(const std::filesystem::path& path) -> void {
    loadSettings(std::ifstream(path));
  }
  auto loadSettings(std::istream& is) -> void {
    loadSettings(std::move(is));
  }
  auto loadSettings(std::istream&& is) -> void;
  auto setProgram(const qc::QuantumComputation& qc) -> void;
  [[nodiscard]] auto toString() const -> std::string;
  friend auto operator<<(std::ostream& os, const CompilerBase& cb) -> std::ostream& {
    os << cb.toString();
    return os;
  }

protected:
  /// collect qubits that will remain in Rydberg zone between two Rydberg stages
  auto collectReuseQubit() -> void;
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
  [[nodiscard]] auto getDir() const -> const std::filesystem::path& { return dir; }
  [[nodiscard]] auto getNQubits() const -> std::size_t { return nQubits; }
  [[nodiscard]] auto getNTwoQubitGates() const -> std::size_t { return nTwoQubitGates; }
  [[nodiscard]] auto getArchitecture() const -> const Architecture& {
    return architecture;
  }
  [[nodiscard]] auto getResult() -> Result& { return result; }
  [[nodiscard]] auto getRuntimeAnalysis() -> RuntimeAnalysis& {
    return runtimeAnalysis;
  }
  [[nodiscard]] auto isToVerify() const -> bool { return toVerify; }
  [[nodiscard]] auto isTrivialPlacement() const -> bool { return trivialPlacement; }
  [[nodiscard]] auto getRoutingStrategy() const -> RoutingStrategy {
    return routingStrategy;
  }
  [[nodiscard]] auto getSchedulingStrategy() const -> SchedulingStrategy {
    return schedulingStrategy;
  }
  [[nodiscard]] auto isDynamicPlacement() const -> bool { return dynamicPlacement; }
  [[nodiscard]] auto getGivenInitialMapping() const -> const std::optional<
      std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>& {
    return givenInitialMapping;
  }
  [[nodiscard]] auto isHasDependency() const -> bool { return hasDependency; }
  [[nodiscard]] auto isL2() const -> bool { return l2; }
  [[nodiscard]] auto isUseWindow() const -> bool { return useWindow; }
  [[nodiscard]] auto getWindowSize() const -> std::size_t { return windowSize; }
  [[nodiscard]] auto isReuse() const -> bool { return reuse; }
  [[nodiscard]] auto getCommon1Q() const -> std::size_t { return common1Q; }
  [[nodiscard]] auto getTwoQubitGates() const
      -> const std::vector<std::pair<qc::Qubit, qc::Qubit>>& {
    return twoQubitGates;
  }
  [[nodiscard]] auto getDictG1QParent() const
      -> const std::unordered_map<const std::pair<qc::Qubit, qc::Qubit>*,
                                  std::vector<qc::StandardOperation>>& {
    return dictG1QParent;
  }
  [[nodiscard]] auto getQubitMapping() -> std::vector<
      std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>& {
    return qubitMapping;
  }
  [[nodiscard]] auto getReuseQubits() const
      -> const std::vector<std::unordered_set<std::size_t>>& {
    return reuseQubits;
  }
  [[nodiscard]] auto getGateScheduling()
      -> std::vector<std::vector<const std::pair<qc::Qubit, qc::Qubit>*>>& {
    return gateScheduling;
  }
  [[nodiscard]] auto getGateSchedulingIdx() const
      -> const std::vector<std::vector<std::size_t>>& {
    return gateSchedulingIdx;
  }
  [[nodiscard]] auto getGate1QScheduling()
      -> std::vector<std::vector<const qc::StandardOperation*>>& {
    return gate1QScheduling;
  }
  auto setDictg1QParent(
      const std::unordered_map<const std::pair<qc::Qubit, qc::Qubit>*,
                               std::vector<qc::StandardOperation>>&
          newDictg1QParent) -> void {
    dictG1QParent = newDictg1QParent;
  }
  auto setQubitMapping(
      const std::vector<
          std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>&
          newQubitMapping) -> void {
    qubitMapping = newQubitMapping;
  }
  auto setGateSchedulingIdx(
      const std::vector<std::vector<std::size_t>>&
          newGateSchedulingIdx) -> void {
    gateSchedulingIdx = newGateSchedulingIdx;
  }
  auto setGate1QScheduling(
      const std::vector<std::vector<const qc::StandardOperation*>>&
          newGate1QScheduling) -> void {
    gate1QScheduling = newGate1QScheduling;
  }
};

} // namespace na
