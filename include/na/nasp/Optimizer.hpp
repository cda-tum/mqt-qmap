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
  std::optional<std::size_t>                  maxUnsat = std::nullopt;
  std::optional<std::size_t>                  minSat   = std::nullopt;
  std::unordered_map<pid_t, OptimizerProcess> processData;

  [[nodiscard]] auto getNSubProcsRunning() const -> std::size_t {
    return processData.size();
  }

  [[nodiscard]] auto isSubProcRunning() const -> bool {
    return !processData.empty();
  }

  auto forkChildProcess(std::size_t arg, std::chrono::seconds timeout) -> void;

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

  auto setTimeout(const std::chrono::minutes timeout) -> void {
    this->timeout = timeout;
  }

  auto setMaxNSubProcs(const std::size_t maxNSubProcs) -> void {
    this->maxNSubProcs = maxNSubProcs;
  }

  auto setInitialValue(const std::size_t initialValue) -> void {
    this->initialValue = initialValue;
  }

  auto setMaxValue(const std::size_t maxValue) -> void {
    this->maxValue = maxValue;
  }

  auto setQuiet(const bool quiet) -> void { this->quiet = quiet; }

  auto setObjectiveFunction(
      const std::function<NASolver::Result(std::size_t)>& objective) -> void {
    this->objective = objective;
  }

  auto minimize() -> void;
};

/**
 * @brief This auxiliary function maps every integer uniquely to a pair of
 * psoitive integers @code std::pair{x, y}@endcode, such that always
 * @code x ≤ y@endcode holds.
 * @details The function is defined as follows:
 * @code i = ceil((x + y - 2) * (x + y) / 4) + x@endcode.
 *
 * @details Example: @code
 * x \ y |  1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 | ...
 *     1 |  1 |  2 |  3 |  5 |  7 | 10 | 13 | 17 | ...
 *     2 |  - |  4 |  6 |  8 | 11 | 14 | 18 | 22 | ...
 *     3 |  - |  - |  9 | 12 | 15 | 19 | 23 | 28 | ...
 *     4 |  - |  - |  - | 16 | 20 | 24 | 29 | 34 | ...
 *     5 |  - |  - |  - |  - | 25 | 30 | 35 | 41 | ...
 *     6 |  - |  - |  - |  - |  - | 36 | 42 | 48 | ...
 *     7 |  - |  - |  - |  - |  - |  - | 49 | 56 | ...
 *     8 |  - |  - |  - |  - |  - |  - |  - | 64 | ...
 * @endcode
 * @param i a positive integer
 * @return a pair of integers
 */
static auto reverseUpperPairingFunction(const std::size_t i)
    -> std::pair<std::size_t, std::size_t> {
  const double      w = floor(1.0 + sqrt(1.0 + 4.0 * static_cast<double>(i)));
  const double      t = ceil((w - 2.0) * w / 4.0);
  const std::size_t x = i - static_cast<std::size_t>(t);
  const std::size_t y = static_cast<std::size_t>(w) - x - 2;
  return std::pair{x, y};
}

/**
 * @brief This auxiliary function maps a pair of positive integers
 * @code std::pair{x, y}@endcode, where always @code x ≤ y@endcode, to a
 * unique positive integer i.
 * @details The function is defined as follows:
 * @code i = ceil((x + y) * (x + y + 2) / 4) + x@endcode.
 *
 * @details Example: @code
 * x \ y |  0 |  1 |  2 |  3 |  4 |  5 |  6 |  7 | ...
 * 0     |  0 |  1 |  3 |  4 |  6 |  9 | 12 | 16 | ...
 * 1     |  - |  3 |  5 |  7 | 10 | 13 | 17 | 21 | ...
 * 2     |  - |  - |  8 | 11 | 14 | 18 | 22 | 27 | ...
 * 3     |  - |  - |  - | 15 | 19 | 23 | 28 | 33 | ...
 * 4     |  - |  - |  - |  - | 24 | 29 | 34 | 40 | ...
 * 5     |  - |  - |  - |  - |  - | 35 | 41 | 47 | ...
 * 6     |  - |  - |  - |  - |  - |  - | 48 | 55 | ...
 * 7     |  - |  - |  - |  - |  - |  - |  - | 63 | ...
 * @endcode
 * @param x a positive integer
 * @param y a positive integer
 * @return a positive integer
 */
static auto upperPairingFunction(const std::size_t x,
                                 const std::size_t y) -> std::size_t {
  if (x > y) {
    throw std::invalid_argument("x must be smaller or equal y.");
  }
  const auto w = static_cast<double>(x + y);
  const auto i = static_cast<std::size_t>(ceil((w + 2.0) * w / 4.0)) + x;
  return i;
}
} // namespace na
