#pragma once

#include "Definitions.hpp"
#include "Solver.hpp"
#include "ir/QuantumComputation.hpp"
#include "na/nalac/datastructures/Architecture.hpp"

#include <utility>
#include <vector>

namespace na {
class SolverFactory {
public:
  /**
   * @brief Create a NASolver from a given architecture.
   *
   * @details The solver considers an abstraction of the real architecture. This
   * function derives the parameters for the abstraction from the given
   * architecture. These parameters are used to initialize the solver, which
   * is then returned.
   *
   * @param arch The architecture to create the solver for.
   * @return The created NASolver.
   */
  [[nodiscard]] static auto create(const nalac::Architecture& arch) -> NASolver;

  /**
   * @brief Get the list of entangling operations that the solver takes as
   * input.
   *
   * @details The solver only considers the entangling operations of a circuit.
   * For that it receives a list of qubit pairs that represent each one
   * entangling operation. This function generates this list from a given
   * QuantumComputation and a FullOpType that specifies the entangling
   * operation.
   *
   * @warning This function expects a QuantumComputation that was used as input
   * for the NASolver. Additionally, this function assumes the quantum circuit
   * represented by the QuantumComputation to be of the following form:
   * First, all qubits are initialized in the |+> state by applying a Hadamard
   * gate to each qubit. Then, a set of entangling gates (CZ) is applied to the
   * qubits. Finally, hadamard gates are applied to some qubits. Unfortunately,
   * the function cannot deal with arbitrary quantum circuits as the NASolver
   * cannot do either.
   *
   * @param circ
   * @param opType
   * @param ctrls
   * @param quiet
   * @return
   */
  [[nodiscard]] static auto
  getOpsForSolver(const qc::QuantumComputation& circ, qc::OpType opType,
                  std::size_t ctrls, bool quiet = false)
      -> std::vector<std::pair<qc::Qubit, qc::Qubit>>;
};
} // namespace na
