//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/mqt-qmap for more
// information.
//

#pragma once

#include "Architecture.hpp"
#include "Definitions.hpp"
#include "Graph.hpp"
#include "QuantumComputation.hpp"
#include "na/GlobalOperation.hpp"
#include "operations/OpType.hpp"
#include "operations/Operation.hpp"

#include <cassert>
#include <iostream>
#include <iterator>
#include <memory>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na {
static constexpr std::array<qc::OpType, 10> DIAGONAL_GATES = {
    qc::Barrier, qc::I,   qc::Z, qc::S,  qc::Sdg,
    qc::T,       qc::Tdg, qc::P, qc::RZ, qc::RZZ};
/**
 * @brief Class to manage the creation of layers when traversing a quantum
 * circuit.
 * @details The class uses the DAG of the circuit to create layers of gates that
 * can be executed at the same time. It can be used to create the front or look
 * ahead layer.
 */
class Layer {
public:
  class DAGVertex : public std::enable_shared_from_this<DAGVertex> {
  protected:
    // if the executableCounter becomes equal to the executableThreshold the
    // vertex becomes executable
    std::size_t                             executableThreshold = 0;
    std::size_t                             executableCounter   = 0;
    std::vector<std::shared_ptr<DAGVertex>> enabledSuccessors   = {};
    std::vector<std::shared_ptr<DAGVertex>> disabledSuccessors  = {};
    bool                                    executed            = false;
    const std::unique_ptr<qc::Operation>*   operation;
    std::unique_ptr<std::unordered_set<std::shared_ptr<DAGVertex>>>*
        executableSet;

    /**
     * @brief Construct a new DAGVertex
     * @details The constructor initializes the vertex with the given operation
     * and the given executable set, which is shared among all vertices.
     *
     * @param operation
     * @param executableSet
     */
    DAGVertex(const std::unique_ptr<qc::Operation>* operation,
              std::unique_ptr<std::unordered_set<std::shared_ptr<DAGVertex>>>*
                  executableSet)
        : operation(operation), executableSet(executableSet) {}

  public:
    [[nodiscard]] static auto
    create(const std::unique_ptr<qc::Operation>* operation,
           std::unique_ptr<std::unordered_set<std::shared_ptr<DAGVertex>>>*
               executableSet) -> std::shared_ptr<DAGVertex> {
      std::shared_ptr<DAGVertex> v(new DAGVertex(operation, executableSet));
      v->updateExecutableSet();
      return v;
    }
    [[nodiscard]] auto isExecutable() const {
      assert(executableCounter <= executableThreshold);
      return (!executed) && executableCounter == executableThreshold;
    }
    [[nodiscard]] auto isExecuted() const { return executed; }
    [[nodiscard]] auto getOperation() const
        -> const std::unique_ptr<qc::Operation>* {
      return operation;
    }

  protected:
    auto incExecutableCounter() {
      executableCounter++;
      updateExecutableSet();
    }
    auto decExecutableCounter() {
      executableCounter--;
      updateExecutableSet();
    }
    /**
     * @brief Inserts or removes the vertex from the executable set.
     * @warning May not be called from the constructor.
     */
    auto updateExecutableSet() -> void {
      if (isExecutable()) {
        if (const auto& it = (*executableSet)->find(shared_from_this());
            it == (*executableSet)->end()) {
          (*executableSet)->insert(shared_from_this());
        }
      } else {
        if (const auto& it = (*executableSet)->find(shared_from_this());
            it != (*executableSet)->end()) {
          (*executableSet)->erase(shared_from_this());
        }
      }
    }

  public:
    auto execute() -> void {
      if (isExecutable()) {
        executed = true;
        for (const auto& successor : disabledSuccessors) {
          successor->decExecutableCounter();
        }
        for (const auto& successor : enabledSuccessors) {
          successor->incExecutableCounter();
        }
        updateExecutableSet();
      } else {
        throw std::logic_error(
            "The vertex is not executable and cannot be executed.");
      }
    }
    auto addEnabledSuccessor(std::shared_ptr<DAGVertex> successor) {
      enabledSuccessors.emplace_back(successor);
      ++successor->executableThreshold;
      successor->updateExecutableSet();
    }
    auto addDisabledSuccessor(std::shared_ptr<DAGVertex> successor) {
      disabledSuccessors.emplace_back(successor);
      --successor->executableThreshold;
      successor->updateExecutableSet();
    }
  };

protected:
  std::unique_ptr<std::unordered_set<std::shared_ptr<DAGVertex>>>
      executableSet =
          std::make_unique<std::unordered_set<std::shared_ptr<DAGVertex>>>();

  /// Checks whether two operations commute on a given qubit.
  [[nodiscard]] static auto
  commutesAtQubit(const std::unique_ptr<qc::Operation>* op1,
                  const std::unique_ptr<qc::Operation>* op2,
                  const qc::Qubit&                      qubit) -> bool;
  /// Checks whether two consecutive operations cancel each other out.
  [[nodiscard]] static auto isInverse(const std::unique_ptr<qc::Operation>* op1,
                                      const std::unique_ptr<qc::Operation>* op2)
      -> bool;
  auto constructDAG(const qc::QuantumComputation& qc) -> void;

public:
  Layer() = default;
  Layer(const Layer& other) {
    executableSet =
        std::make_unique<std::unordered_set<std::shared_ptr<DAGVertex>>>(
            other.executableSet->cbegin(), other.executableSet->cend());
  }
  virtual ~Layer() = default;
  explicit Layer(const qc::QuantumComputation& qc) { constructDAG(qc); }
  Layer& operator=(const Layer& other) {
    executableSet =
        std::make_unique<std::unordered_set<std::shared_ptr<DAGVertex>>>(
            other.executableSet->cbegin(), other.executableSet->cend());
    return *this;
  }
  [[nodiscard]] auto getExecutableSet() const -> const
      std::unique_ptr<std::unordered_set<std::shared_ptr<DAGVertex>>>* {
    return &executableSet;
  }
  auto setCircuit(const qc::QuantumComputation& qc) -> void {
    executableSet->clear();
    constructDAG(qc);
  }
  static auto execute(const std::shared_ptr<DAGVertex>& vertex) -> void {
    vertex->execute();
  }
  [[nodiscard]] auto constructInteractionGraph(qc::OpType opType, Number nctrl) const
      -> Graph<std::shared_ptr<DAGVertex>>;
};

} // namespace na
