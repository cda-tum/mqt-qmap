#pragma once

#include "yaml-cpp/node/node.h"
#include <cstddef>
#include <cstdint>
#include <optional>
#include <z3++.h>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace na {
using namespace z3;

class NASolver {
protected:
  static context ctx;

private:
  unsigned int maxX = 0;
  unsigned int maxY = 0;
  unsigned int minEntanglingY = 0;
  unsigned int maxEntanglingY = 0;
  unsigned int maxC = 0;
  unsigned int maxR = 0;
  unsigned int maxHOffset = 0;
  unsigned int maxVOffset = 0;
  unsigned int maxHDist = 0;
  unsigned int maxVDist = 0;

  enum class Storage : std::uint8_t { None, Bottom, TwoSided };

  Storage storage = Storage::None;

  static auto minBitsToRepresentUInt(unsigned int num) -> unsigned int;
  static auto minBitsToRepresentInt(unsigned int num) -> unsigned int;

  /// A class to collect all varaibles associated with one qubit
  class Qubit {
  private:
    unsigned int id;
    expr x;
    expr y;
    expr a;
    expr c;
    expr r;
    expr h;
    expr v;

  public:
    [[nodiscard]] Qubit(const unsigned int id, const unsigned int t,
                        const unsigned int maxX,
                        const unsigned int maxY, const unsigned int maxC,
                        const unsigned int maxR, const unsigned int maxHOffset,
                        const unsigned int maxVOffset)
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
            minBitsToRepresentInt(maxVOffset))) {
    }

    /// unique identifier of the qubit
    [[nodiscard]] unsigned int getId() const {
      return id;
    }

    /// x-coordinate of the site the atom is loaded in
    [[nodiscard]] const expr& getX() const {
      return x;
    }

    /// y-coordinate of the site the atom is loaded in
    [[nodiscard]] const expr& getY() const {
      return y;
    }

    /// boolean variable to indicate whether the atom is loaded in an AOD, SLM otherwise
    [[nodiscard]] const expr& getA() const {
      return a;
    }

    /// if the atom is loaded in an AOD, this is the index of the AOD column, otherwise it has no meaning
    [[nodiscard]] const expr& getC() const {
      return c;
    }

    /// if the atom is loaded in an AOD, this is the index of the AOD row, otherwise it has no meaning
    [[nodiscard]] const expr& getR() const {
      return r;
    }

    /// denotes the horizontal offset from the SLM trap
    [[nodiscard]] const expr& getH() const {
      return h;
    }

    /// denotes the vertical offset from the SLM trap
    [[nodiscard]] const expr& getV() const {
      return v;
    }
  };

  class Stage {
  private:
    unsigned int t;
    std::vector<Qubit> qubits;
    std::vector<expr> loadCols;
    std::vector<expr> loadRows;
    std::vector<expr> storeCols;
    std::vector<expr> storeRows;

  public:
    [[nodiscard]] explicit
    Stage(const size_t t, const unsigned int numQubits, const unsigned int maxX,
          const unsigned int maxY, const unsigned int maxC,
          const unsigned int maxR, const unsigned int maxHOffset,
          const unsigned int maxVOffset
        ) : t(t) {
      qubits.reserve(numQubits);
      for (unsigned int id = 0; id < numQubits; ++id) {
        qubits.emplace_back(id, t, maxX, maxY, maxC, maxR, maxHOffset,
                            maxVOffset);
      }
      loadCols.reserve(maxC);
      storeCols.reserve(maxC);
      for (unsigned int c = 0; c <= maxC; ++c) {
        std::stringstream suffixStream;
        suffixStream << "_" << t << "^c" << c;
        const auto& suffix = suffixStream.str();
        loadCols.emplace_back(ctx.bool_const(("load" + suffix).c_str()));
        storeCols.emplace_back(ctx.bool_const(("store" + suffix).c_str()));
      }
      loadRows.reserve(maxR);
      storeRows.reserve(maxR);
      for (unsigned int r = 0; r <= maxR; ++r) {
        std::stringstream suffixStream;
        suffixStream << "_" << t << "^r" << r;
        const auto& suffix = suffixStream.str();
        loadRows.emplace_back(ctx.bool_const(("load" + suffix).c_str()));
        storeRows.emplace_back(ctx.bool_const(("store" + suffix).c_str()));
      }
    }

    [[nodiscard]] unsigned int getT() const {
      return t;
    }

    [[nodiscard]] auto getQubit(const size_t i) const -> const Qubit& {
      return qubits[i];
    }

    [[nodiscard]] auto numQubits() const -> size_t {
      return qubits.size();
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

  unsigned int numQubits = 0;
  unsigned int numStages = 0;
  std::optional<unsigned int> numTransfers = std::nullopt;
  std::vector<Stage> stages;
  std::vector<expr> transfers;
  std::vector<expr> gates;

  /// Initializes the variables for all stages and all qubits
  auto initVariables() -> void;

  /// Return constraints ensuring that exactly @code numTransfers@endcode transfers trake place
  [[nodiscard]] auto
  getExactNumTransfersConstraints() const -> std::vector<expr>;

  /// Returns the constraint @code (x_t^(q0) = x_t^(q1)) ∧ (y_t^(q0) = y_t^(q1)) @endcode
  [[nodiscard]] auto getHaveSamePositionConstraint(
      unsigned int q0, unsigned int q1,
      unsigned int t) const -> expr;

  /// Returns the constraint @code (x_t^(q0) ≠ x_t^(q1)) ∨ (y_t^(q0) ≠ y_t^(q1)) @endcode
  [[nodiscard]] auto getHaveDifferentPositionConstraint(
      unsigned int q0, unsigned int q1,
      unsigned int t) const -> expr;

  /// Return constraints ensuring that the qubits is in the entangling zone at stage t
  [[nodiscard]] auto getAffectedByRydbergBeamConstraint(
      unsigned int q, unsigned int t) const -> expr;

  /// Return constraints ensuring that the qubits is in the entangling zone at stage t
  [[nodiscard]] auto getShieldedFromRydbergBeamConstraint(
      unsigned int q, unsigned int t) const -> expr;

  /// Retruns a vecotr of constraints ensuring that transition from a Rydberg stage to the next stage is valid
  [[nodiscard]] auto getValidRydbergTransitionConstraints(
      unsigned int t) const -> std::vector<expr>;

  /// Retruns a vecotr of constraints ensuring that transition from a Transfer stage to the next stage is valid
  [[nodiscard]] auto getValidTransferTransitionConstraints(
      unsigned int t) const -> std::vector<expr>;

  /**
   * @brief Returns the constraints extracted from the quantum circuit to ensure execution of each gate and not more.
   * @details Creates the variables @code gate_i@endcode for every gate between the qubits @code q0, q1@endcode
   * and returns the following constraints:
   * @code
   * (0 ≤ gate_i) ∧ (gate_i < numStages) for all i
   * (gate_i = t) ⟷ rydbergStage(t) ∧ haveSamePosition(q0, q1, t) for all i, t and (q0, q1) is a gate_i
   * rydbergStage(t) ⟶ haveDifferentPosition(q, q', t) for all q, q', t where (q, q') is not a gate
   * @endcode
   * @return a vector of the constraints described above
   */
  [[nodiscard]] auto getCircuitExecutionConstraints(
      const std::vector<std::pair<unsigned int, unsigned int> >& ops,
      bool mindOpsOrder, bool shieldIdleAtoms) ->
    std::vector<expr>;

  /// Returns a constraint expressing that this stage is a Rydberg stage, that is if
  /// @code numTransfers_{t-1} = numTransfers_t @endcode
  [[nodiscard]] auto getRydbergStageConstraint(unsigned int t) const -> expr;

  /// Returns a constraint expressing that this stage is a Transfer stage, that is if
  /// @code numTransfers_{t-1} + 1 = numTransfers_t @endcode
  [[nodiscard]] auto getTransferStageConstraint(unsigned int t) const -> expr;

  /// Returns constraints esnuring that the state at the given stage is valid
  [[nodiscard]] auto getValidStageConstraints(
      unsigned int t) const -> std::vector<expr>;

public:
  [[nodiscard]] NASolver() = default;

  auto init(unsigned int maxX, unsigned int maxY, unsigned int maxC,
            unsigned int maxR, unsigned int maxHOffset,
            unsigned int maxVOffset, unsigned int maxHDist,
            unsigned int maxVDist, unsigned int minEntanglingY,
            unsigned int maxEntanglingY) -> void;

  class Result {
  public:
    class Qubit {
    private:
      unsigned int x;
      unsigned int y;
      bool a;
      unsigned int c;
      unsigned int r;
      int h;
      int v;

    public:
      [[nodiscard]] Qubit() = default;

      [[nodiscard]] Qubit(const unsigned int x, const unsigned int y,
                          const bool a, const unsigned int c,
                          const unsigned int r, const int h,
                          const int v)
        : x(x), y(y), a(a), c(c), r(r), h(h), v(v) {
      }

      [[nodiscard]] static auto fromYAML(const YAML::Node& yaml) -> Qubit;

      [[nodiscard]] auto getX() const -> unsigned int { return x; }
      [[nodiscard]] auto getY() const -> unsigned int { return y; }
      [[nodiscard]] auto isAOD() const -> bool { return a; }
      [[nodiscard]] auto getC() const -> unsigned int { return c; }
      [[nodiscard]] auto getR() const -> unsigned int { return r; }
      [[nodiscard]] auto getH() const -> int { return h; }
      [[nodiscard]] auto getV() const -> int { return v; }

      [[nodiscard]] auto yaml(unsigned int indent,
                              bool item = true,
                              bool compact = true) const -> std::string;
    };

    class Gate {
    private:
      unsigned int stage = 0;
      std::pair<unsigned int, unsigned int> qubits;

    public:
      [[nodiscard]] Gate() = default;

      [[nodiscard]] Gate(const unsigned int stage,
                         const std::pair<unsigned int, unsigned int>& qubits)
        : stage(stage), qubits(qubits) {
      }

      [[nodiscard]] Gate(const unsigned int stage,
                         std::pair<unsigned int, unsigned int>&& qubits)
        noexcept
        : stage(stage), qubits(std::move(qubits)) {
      }

      [[nodiscard]] static auto fromYAML(const YAML::Node& yaml) -> Gate;

      [[nodiscard]] auto getStage() const -> unsigned int {
        return stage;
      }

      [[nodiscard]] auto
      getQubits() const -> const std::pair<unsigned int, unsigned int>& {
        return qubits;
      }

      [[nodiscard]] auto yaml(unsigned int indent,
                              bool item = true,
                              bool compact = true) const -> std::string;
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
        : rydberg(rydberg), qubits(qubits), gates(gates) {
      }

      [[nodiscard]] Stage(const bool rydberg,
                          std::vector<Qubit>&& qubits,
                          std::vector<Gate>&& gates) noexcept
        : rydberg(rydberg), qubits(std::move(qubits)), gates(std::move(gates)) {
      }

      [[nodiscard]] static auto fromYAML(const YAML::Node& yaml) -> Stage;

      [[nodiscard]] auto isRydberg() const -> bool { return rydberg; }

      [[nodiscard]] auto getQubit(const unsigned int i) const -> const Qubit& {
        return qubits[i];
      }

      [[nodiscard]] auto numQubits() const -> unsigned int {
        return qubits.size();
      }

      [[nodiscard]] auto getQubits() const -> const std::vector<Qubit>& {
        return qubits;
      }

      [[nodiscard]] auto getGate(const unsigned int i) const -> const Gate& {
        return gates[i];
      }

      [[nodiscard]] auto numGates() const -> unsigned int {
        return gates.size();
      }

      [[nodiscard]] auto getGates() const -> const std::vector<Gate>& {
        return gates;
      }

      [[nodiscard]] auto yaml(unsigned int indent,
                              bool item = true,
                              bool compact = true) const -> std::string;
    };

  private:
    bool sat = false;
    std::vector<Stage> stages;

  public:
    [[nodiscard]] Result() = default;

    [[nodiscard]] explicit Result(const bool sat)
      : sat(sat) {
    }

    [[nodiscard]] Result(const bool sat,
                         const std::vector<Stage>& stages)
      : sat(sat), stages(stages) {
    }

    [[nodiscard]] Result(const bool sat,
                         std::vector<Stage>&& stages) noexcept
      : sat(sat), stages(std::move(stages)) {
    }

    [[nodiscard]] static auto fromYAML(const YAML::Node& yaml) -> Result;

    [[nodiscard]] auto getStage(const unsigned int i) const -> const Stage& {
      return stages[i];
    }

    [[nodiscard]] auto numStages() const -> unsigned int {
      return stages.size();
    }

    [[nodiscard]] auto isSat() const -> bool {
      return sat;
    }

    [[nodiscard]] auto front() const -> const Stage& {
      return stages.front();
    }

    [[nodiscard]] auto begin() const -> decltype(stages.begin()) {
      return stages.begin();
    }

    [[nodiscard]] auto end() const -> decltype(stages.end()) {
      return stages.end();
    }

    [[nodiscard]] auto yaml(unsigned int indent = 0,
                            bool compact = true) const -> std::string;
  };

  [[nodiscard]] auto solve(
      const std::vector<std::pair<unsigned int, unsigned int> >& ops,
      unsigned int numQubits,
      unsigned int numStages,
      unsigned int numTransfers = 0,
      bool mindOpsOrder = false,
      bool shieldIdleQubits = true) -> Result;

  [[nodiscard]] auto solve(
      const std::vector<std::pair<unsigned int, unsigned int> >& ops,
      unsigned int numQubits,
      unsigned int numStages,
      bool mindOpsOrder = false,
      bool shieldIdleQubits = true) -> Result;
};

struct ExprHash {
  unsigned int operator()(const expr& e) const {
    return e.hash();
  }
};
} // namespace na