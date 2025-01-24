//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "Definitions.hpp"
#include "hybridmap/NeutralAtomDefinitions.hpp"

#include <cstdint>
#include <set>
#include <utility>
#include <vector>

namespace na {
/**
 * @brief Class to manage the creation of layers when traversing a quantum
 * circuit.
 * @details The class uses the qc::DAG of the circuit to create layers of gates
 * that can be executed at the same time. It can be used to create the front or
 * look ahead layer.
 */

class NeutralAtomLayer {
protected:
  qc::DAG dag;
  qc::DAGIterators iterators;
  qc::DAGIterators ends;
  GateList gates;
  GateList newGates;
  GateList mappedSingleQubitGates;
  GateLists candidates;

  /**
   * @brief Updates the gates for the given qubits
   * @details The function iterates over the qc::DAG and updates the gates for
   * the given qubits as far es possible.
   * @param qubitsToUpdate The qubits that have been updated
   * @param commuteWith Gates the new gates should commute with
   */
  void updateByQubits(const std::set<qc::Qubit>& qubitsToUpdate);

  /**
   * @brief Updates the candidates for the given qubits
   */
  void updateCandidatesByQubits(const std::set<qc::Qubit>& qubitsToUpdate);
  /**
   * @brief Checks the candidates and add them to the gates if possible
   * @param qubitsToUpdate The qubits that have been updated
   */
  void candidatesToGates(const std::set<qc::Qubit>& qubitsToUpdate);

public:
  // Constructor
  explicit NeutralAtomLayer(qc::DAG graph) : dag(std::move(graph)) {
    iterators.reserve(dag.size());
    candidates.reserve(dag.size());
    for (auto& i : dag) {
      auto it = i.begin();
      iterators.emplace_back(it);
      ends.emplace_back(i.end());
      candidates.emplace_back();
    }
  }

  /**
   * @brief Returns the current layer of gates
   * @return The current layer of gates
   */
  GateList getGates() const { return gates; }
  GateList getNewGates() const { return newGates; }
  /**
   * @brief Returns a vector of the iterator indices for debugging
   * @return A copy of the current iterator indices
   */
  std::vector<uint32_t> getIteratorOffset();
  /**
   * @brief Initializes the layer by updating all qubits starting
   */
  void initAllQubits();
  /**
   * @brief Removes the provided gates from the current layer and update the
   * the layer depending on the qubits of the gates.
   * @param gatesToRemove Gates to remove from the current layer
   * @param commuteWith Gates the new gates should commute with
   */
  void removeGatesAndUpdate(const GateList& gatesToRemove);

  /**
   * @brief Returns the mapped single qubit gates
   * @return The mapped single qubit gates
   */
  GateList getMappedSingleQubitGates() { return mappedSingleQubitGates; }
};

// Commutation checks
bool commutesWithAtQubit(const GateList& layer, const qc::Operation* opPointer,
                         const qc::Qubit& qubit);
bool commuteAtQubit(const qc::Operation* opPointer1,
                    const qc::Operation* opPointer2, const qc::Qubit& qubit);
} // namespace na
