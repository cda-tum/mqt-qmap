#pragma once

#include "Definitions.hpp"
#include "yaml-cpp/node/node.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <z3++.h>

namespace na {
using namespace z3;

class NASolver {
protected:
  static context ctx;

private:
  std::int32_t maxX = 0;
  std::int32_t maxY = 0;
  std::int32_t minEntanglingY = 0;
  std::int32_t maxEntanglingY = 0;
  std::int32_t maxC = 0;
  std::int32_t maxR = 0;
  std::int32_t maxHOffset = 0;
  std::int32_t maxVOffset = 0;
  std::int32_t maxHDist = 0;
  std::int32_t maxVDist = 0;

  enum class Storage : std::uint8_t { None, Bottom, TwoSided };

  Storage storage = Storage::None;

  static auto minBitsToRepresentUInt(std::int32_t num) -> std::uint32_t;
  static auto minBitsToRepresentInt(std::int32_t num) -> std::uint32_t;

  /// A class to collect all variables associated with one qubit
  class Qubit {
  private:
    std::uint16_t id;
    expr x;
    expr y;
    expr a;
    expr c;
    expr r;
    expr h;
    expr v;

  public:
    [[nodiscard]] Qubit(const std::uint16_t id, const std::uint16_t t,
                        const std::uint16_t maxX, const std::uint16_t maxY,
                        const std::uint16_t maxC, const std::uint16_t maxR,
                        const std::uint16_t maxHOffset,
                        const std::uint16_t maxVOffset)
        : id(id),
          x(ctx.bv_const(
              ("x" + std::to_string(t) + "^" + std::to_string(id)).c_str(),
              minBitsToRepresentUInt(maxX))),
          y(ctx.bv_const(
              ("y" + std::to_string(t) + "^" + std::to_string(id)).c_str(),
              minBitsToRepresentUInt(maxY))),
          a(ctx.bool_const(
              ("a" + std::to_string(t) + "^" + std::to_string(id)).c_str())),
          c(ctx.bv_const(
              ("c" + std::to_string(t) + "^" + std::to_string(id)).c_str(),
              minBitsToRepresentUInt(maxC))),
          r(ctx.bv_const(
              ("r" + std::to_string(t) + "^" + std::to_string(id)).c_str(),
              minBitsToRepresentUInt(maxR))),
          h(ctx.bv_const(
              ("h" + std::to_string(t) + "^" + std::to_string(id)).c_str(),
              minBitsToRepresentInt(maxHOffset))),
          v(ctx.bv_const(
              ("v" + std::to_string(t) + "^" + std::to_string(id)).c_str(),
              minBitsToRepresentInt(maxVOffset))) {}

    /// unique identifier of the qubit
    [[nodiscard]] std::uint16_t getId() const { return id; }

    /// x-coordinate of the site the atom is loaded in
    [[nodiscard]] const expr& getX() const { return x; }

    /// y-coordinate of the site the atom is loaded in
    [[nodiscard]] const expr& getY() const { return y; }

    /// boolean variable to indicate whether the atom is loaded in an AOD, SLM
    /// otherwise
    [[nodiscard]] const expr& getA() const { return a; }

    /// if the atom is loaded in an AOD, this is the index of the AOD column,
    /// otherwise it has no meaning
    [[nodiscard]] const expr& getC() const { return c; }

    /// if the atom is loaded in an AOD, this is the index of the AOD row,
    /// otherwise it has no meaning
    [[nodiscard]] const expr& getR() const { return r; }

    /// denotes the horizontal offset from the SLM trap
    [[nodiscard]] const expr& getH() const { return h; }

    /// denotes the vertical offset from the SLM trap
    [[nodiscard]] const expr& getV() const { return v; }
  };

  class Stage {
  private:
    std::uint16_t t;
    std::vector<Qubit> qubits;
    std::vector<expr> loadCols;
    std::vector<expr> loadRows;
    std::vector<expr> storeCols;
    std::vector<expr> storeRows;

  public:
    [[nodiscard]] explicit Stage(
        const std::uint16_t t, const std::uint16_t numQubits,
        const std::uint16_t maxX, const std::uint16_t maxY,
        const std::uint16_t maxC, const std::uint16_t maxR,
        const std::uint16_t maxHOffset, const std::uint16_t maxVOffset)
        : t(t) {
      qubits.reserve(numQubits);
      for (std::uint16_t id = 0; id < numQubits; ++id) {
        qubits.emplace_back(id, t, maxX, maxY, maxC, maxR, maxHOffset,
                            maxVOffset);
      }
      loadCols.reserve(maxC);
      storeCols.reserve(maxC);
      for (std::uint16_t c = 0; c <= maxC; ++c) {
        std::stringstream suffixStream;
        suffixStream << "_" << t << "^c" << c;
        const auto& suffix = suffixStream.str();
        loadCols.emplace_back(ctx.bool_const(("load" + suffix).c_str()));
        storeCols.emplace_back(ctx.bool_const(("store" + suffix).c_str()));
      }
      loadRows.reserve(maxR);
      storeRows.reserve(maxR);
      for (std::uint16_t r = 0; r <= maxR; ++r) {
        std::stringstream suffixStream;
        suffixStream << "_" << t << "^r" << r;
        const auto& suffix = suffixStream.str();
        loadRows.emplace_back(ctx.bool_const(("load" + suffix).c_str()));
        storeRows.emplace_back(ctx.bool_const(("store" + suffix).c_str()));
      }
    }

    [[nodiscard]] std::uint16_t getT() const { return t; }

    [[nodiscard]] auto getQubit(const size_t i) const -> const Qubit& {
      return qubits[i];
    }

    [[nodiscard]] auto numQubits() const {
      return static_cast<std::uint16_t>(qubits.size());
    }

    [[nodiscard]] auto getLoadCol(const size_t i) const -> const expr& {
      return loadCols[i];
    }

    [[nodiscard]] auto getLoadRow(const size_t i) const -> const expr& {
      return loadRows[i];
    }

    [[nodiscard]] auto getStoreCol(const size_t i) const -> const expr& {
      return storeCols[i];
    }

    [[nodiscard]] auto getStoreRow(const size_t i) const -> const expr& {
      return storeRows[i];
    }
  };

  std::uint16_t numQubits = 0;
  std::uint16_t numStages = 0;
  std::optional<std::uint16_t> numTransfers = std::nullopt;
  std::vector<Stage> stages;
  std::vector<expr> transfers;
  std::vector<expr> gates;

  /// Initializes the variables for all stages and all qubits
  auto initVariables() -> void;

  /// Return constraints ensuring that exactly @code numTransfers@endcode
  /// transfers trake place
  [[nodiscard]] auto
  getExactNumTransfersConstraints() const -> std::vector<expr>;

  /// Returns the constraint @code (x_t^(q0) = x_t^(q1)) ∧ (y_t^(q0) = y_t^(q1))
  /// @endcode
  [[nodiscard]] auto
  getHaveSamePositionConstraint(std::uint16_t q0, std::uint16_t q1,
                                std::uint16_t t) const -> expr;

  /// Returns the constraint @code (x_t^(q0) ≠ x_t^(q1)) ∨ (y_t^(q0) ≠ y_t^(q1))
  /// @endcode
  [[nodiscard]] auto
  getHaveDifferentPositionConstraint(std::uint16_t q0, std::uint16_t q1,
                                     std::uint16_t t) const -> expr;

  /// Return constraints ensuring that the qubits is in the entangling zone at
  /// stage t
  [[nodiscard]] auto
  getAffectedByRydbergBeamConstraint(std::uint16_t q,
                                     std::uint16_t t) const -> expr;

  /// Return constraints ensuring that the qubits is in the entangling zone at
  /// stage t
  [[nodiscard]] auto
  getShieldedFromRydbergBeamConstraint(std::uint16_t q,
                                       std::uint16_t t) const -> expr;

  /// Returns a vector of constraints ensuring that transition from a Rydberg
  /// stage to the next stage is valid
  [[nodiscard]] auto getValidRydbergTransitionConstraints(std::uint16_t t) const
      -> std::vector<expr>;

  /// Returns a vector of constraints ensuring that transition from a Transfer
  /// stage to the next stage is valid
  [[nodiscard]] auto getValidTransferTransitionConstraints(
      std::uint16_t t) const -> std::vector<expr>;

  /**
   * @brief Returns the constraints extracted from the quantum circuit to ensure
   * execution of each gate and not more.
   * @details Creates the variables @code gate_i@endcode for every gate between
   * the qubits @code q0, q1@endcode and returns the following constraints:
   * @code
   * (0 ≤ gate_i) ∧ (gate_i < numStages) for all i
   * (gate_i = t) ⟷ rydbergStage(t) ∧ haveSamePosition(q0, q1, t) for all i, t
   * and (q0, q1) is a gate_i rydbergStage(t) ⟶ haveDifferentPosition(q, q', t)
   * for all q, q', t where (q, q') is not a gate
   * @endcode
   * @return a vector of the constraints described above
   */
  [[nodiscard]] auto getCircuitExecutionConstraints(
      const std::vector<std::pair<qc::Qubit, qc::Qubit>>& ops,
      bool mindOpsOrder, bool shieldIdleAtoms) -> std::vector<expr>;

  /// Returns a constraint expressing that this stage is a Rydberg stage, that
  /// is if
  /// @code numTransfers_{t-1} = numTransfers_t @endcode
  [[nodiscard]] auto getRydbergStageConstraint(std::uint16_t t) const -> expr;

  /// Returns a constraint expressing that this stage is a Transfer stage, that
  /// is if
  /// @code numTransfers_{t-1} + 1 = numTransfers_t @endcode
  [[nodiscard]] auto getTransferStageConstraint(std::uint16_t t) const -> expr;

  /// Returns constraints esnuring that the state at the given stage is valid
  [[nodiscard]] auto
  getValidStageConstraints(std::uint16_t t) const -> std::vector<expr>;

public:
  [[nodiscard]] NASolver() = default;

  auto init(std::uint16_t newMaxX, std::uint16_t newMaxY, std::uint16_t newMaxC,
            std::uint16_t newMaxR, std::uint16_t newMaxHOffset,
            std::uint16_t newMaxVOffset, std::uint16_t newMaxHDist,
            std::uint16_t newMaxVDist, std::uint16_t newMinEntanglingY,
            std::uint16_t newMaxEntanglingY) -> void;

  class Result {
  public:
    class Qubit {
    private:
      std::int32_t x;
      std::int32_t y;
      bool a;
      std::int32_t c;
      std::int32_t r;
      std::int32_t h;
      std::int32_t v;

    public:
      [[nodiscard]] Qubit() = default;

      [[nodiscard]] Qubit(const std::uint16_t x, const std::uint16_t y,
                          const bool a, const std::uint16_t c,
                          const std::uint16_t r, const std::int32_t h,
                          const std::int32_t v)
          : x(x), y(y), a(a), c(c), r(r), h(h), v(v) {}

      [[nodiscard]] static auto fromYAML(const YAML::Node& yaml) -> Qubit;

      [[nodiscard]] auto getX() const -> std::int32_t { return x; }
      [[nodiscard]] auto getY() const -> std::int32_t { return y; }
      [[nodiscard]] auto isAOD() const -> bool { return a; }
      [[nodiscard]] auto getC() const -> std::int32_t { return c; }
      [[nodiscard]] auto getR() const -> std::int32_t { return r; }
      [[nodiscard]] auto getH() const -> std::int32_t { return h; }
      [[nodiscard]] auto getV() const -> std::int32_t { return v; }

      [[nodiscard]] auto yaml(std::size_t indent, bool item = true,
                              bool compact = true) const -> std::string;
      [[nodiscard]] auto operator==(const Qubit& other) const -> bool;
    };

    class Gate {
    private:
      std::uint16_t stage = 0;
      std::pair<qc::Qubit, qc::Qubit> qubits;

    public:
      [[nodiscard]] Gate() = default;

      [[nodiscard]] Gate(const std::uint16_t stage,
                         const std::pair<qc::Qubit, qc::Qubit>& qubits)
          : stage(stage), qubits(qubits) {}

      [[nodiscard]]
      Gate(const std::uint16_t stage,
           std::pair<qc::Qubit, qc::Qubit>&& qubits) noexcept
          : stage(stage), qubits(std::move(qubits)) {}

      [[nodiscard]] static auto fromYAML(const YAML::Node& yaml) -> Gate;

      [[nodiscard]] auto getStage() const -> std::uint16_t { return stage; }

      [[nodiscard]] auto
      getQubits() const -> const std::pair<qc::Qubit, qc::Qubit>& {
        return qubits;
      }

      [[nodiscard]] auto yaml(std::size_t indent, bool item = true,
                              bool compact = true) const -> std::string;

      [[nodiscard]] auto operator==(const Gate& other) const -> bool;
    };

    class Stage {
    private:
      bool rydberg = true;
      std::vector<Qubit> qubits;
      std::vector<Gate> gates;

    public:
      [[nodiscard]] Stage() = default;

      [[nodiscard]] Stage(const bool rydberg, const std::vector<Qubit>& qubits,
                          const std::vector<Gate>& gates)
          : rydberg(rydberg), qubits(qubits), gates(gates) {}

      [[nodiscard]] Stage(const bool rydberg, std::vector<Qubit>&& qubits,
                          std::vector<Gate>&& gates) noexcept
          : rydberg(rydberg), qubits(std::move(qubits)),
            gates(std::move(gates)) {}

      [[nodiscard]] static auto fromYAML(const YAML::Node& yaml) -> Stage;

      [[nodiscard]] auto isRydberg() const -> bool { return rydberg; }

      [[nodiscard]] auto getQubit(const std::uint16_t i) const -> const Qubit& {
        return qubits[i];
      }

      [[nodiscard]] auto numQubits() const {
        return static_cast<std::uint16_t>(qubits.size());
      }

      [[nodiscard]] auto getQubits() const -> const std::vector<Qubit>& {
        return qubits;
      }

      [[nodiscard]] auto getGate(const std::uint16_t i) const -> const Gate& {
        return gates[i];
      }

      [[nodiscard]] auto numGates() const {
        return static_cast<std::uint16_t>(gates.size());
      }

      [[nodiscard]] auto getGates() const -> const std::vector<Gate>& {
        return gates;
      }

      [[nodiscard]] auto yaml(std::size_t indent, bool item = true,
                              bool compact = true) const -> std::string;

      [[nodiscard]] auto operator==(const Stage& other) const -> bool;
    };

  private:
    bool sat = false;
    std::vector<Stage> stages;

  public:
    [[nodiscard]] Result() = default;

    [[nodiscard]] explicit Result(const bool sat) : sat(sat) {}

    [[nodiscard]] Result(const bool sat, const std::vector<Stage>& stages)
        : sat(sat), stages(stages) {}

    [[nodiscard]] Result(const bool sat, std::vector<Stage>&& stages) noexcept
        : sat(sat), stages(std::move(stages)) {}

    [[nodiscard]] static auto fromYAML(const YAML::Node& yaml) -> Result;

    [[nodiscard]] auto getStage(const std::uint16_t i) const -> const Stage& {
      return stages[i];
    }

    [[nodiscard]] auto numStages() const {
      return static_cast<std::uint16_t>(stages.size());
    }

    [[nodiscard]] auto isSat() const -> bool { return sat; }

    [[nodiscard]] auto front() const -> const Stage& { return stages.front(); }

    [[nodiscard]] auto begin() const -> decltype(stages.begin()) {
      return stages.begin();
    }

    [[nodiscard]] auto end() const -> decltype(stages.end()) {
      return stages.end();
    }

    [[nodiscard]] auto yaml(std::size_t indent = 0,
                            bool compact = true) const -> std::string;

    [[nodiscard]] auto operator==(const Result& other) const -> bool;
  };

  [[nodiscard]] auto
  solve(const std::vector<std::pair<qc::Qubit, qc::Qubit>>& ops,
        std::uint16_t newNumQubits, std::uint16_t newNumStages,
        std::uint16_t newNumTransfers = 0, bool mindOpsOrder = false,
        bool shieldIdleQubits = true) -> Result;

  [[nodiscard]] auto
  solve(const std::vector<std::pair<qc::Qubit, qc::Qubit>>& ops,
        std::uint16_t newNumQubits, std::uint16_t newNumStages,
        bool mindOpsOrder = false, bool shieldIdleQubits = true) -> Result;
};

struct ExprHash {
  std::uint32_t operator()(const expr& e) const { return e.hash(); }
};
} // namespace na
