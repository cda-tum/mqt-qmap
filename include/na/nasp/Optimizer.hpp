#pragma once

#include "Solver.hpp"

#include <chrono>
#include <cstddef>
#include <functional>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace na {

class Optimizer {
public:
  using ObjectiveFunction = std::function<NASolver::Result(std::size_t)>;

private:
  struct OptimizerProcess {
    std::size_t arg        = 0;
    int         readPipeFd = 0;
  };

  std::chrono::minutes            timeout      = std::chrono::minutes::zero();
  std::size_t                     maxNSubProcs = 1;
  std::size_t                     initialValue = 0;
  std::size_t                     maxValue     = 0;
  ObjectiveFunction               objective    = nullptr;
  std::optional<NASolver::Result> extremum     = std::nullopt;
  bool                            quiet        = true;

  // auxiliary variables for minimize
  std::optional<std::uint16_t>                maxUnsat = std::nullopt;
  std::optional<std::uint16_t>                minSat   = std::nullopt;
  std::unordered_map<pid_t, OptimizerProcess> processData;

  [[nodiscard]] auto getNSubProcsRunning() const -> std::size_t {
    return processData.size();
  }

  [[nodiscard]] auto isSubProcRunning() const -> bool {
    return !processData.empty();
  }

  auto forkChildProcess(std::uint16_t        arg,
                        std::chrono::seconds childTimeout) -> void;

  auto waitForChildProcess() -> void;

  auto killAllChildProcesses() -> void;

public:
  [[nodiscard]] auto getExtremum() const -> NASolver::Result {
    if (!extremum.has_value()) {
      throw std::runtime_error("No extremum found.");
    }
    return extremum.value();
  }

  [[nodiscard]] auto getExtremumOpt() const -> std::optional<NASolver::Result> {
    return extremum;
  }

  auto setTimeout(const std::chrono::minutes newTimeout) -> void {
    timeout = newTimeout;
  }

  auto setMaxNSubProcs(const std::uint16_t newMaxNSubProcs) -> void {
    maxNSubProcs = newMaxNSubProcs;
  }

  auto setInitialValue(const std::uint16_t newInitialValue) -> void {
    initialValue = newInitialValue;
  }

  auto setMaxValue(const std::uint16_t newMaxValue) -> void {
    maxValue = newMaxValue;
  }

  auto setQuiet(const bool newQuiet) -> void { quiet = newQuiet; }

  auto setObjectiveFunction(
      const std::function<NASolver::Result(std::uint16_t)>& newObjective)
      -> void {
    objective = newObjective;
  }

  auto minimize() -> void;
};
} // namespace na
