//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "Definitions.hpp"
#include "hybridmap/NeutralAtomDefinitions.hpp"
#include "vector"

#include <cstdint>
#include <set>
#include <utility>

namespace qc {
/**
 * @brief Class to manage the creation of layers when traversing a quantum
 * circuit.
 * @details The class uses the DAG of the circuit to create layers of gates that
 * can be executed at the same time. It can be used to create the front or look
 * ahead layer.
 */

class NeutralAtomLayer {
protected:
  DAG                   dag;
  DAGIterators          iterators;
  GateList              gates;
  GateList              mappedSingleQubitGates;
  std::vector<GateList> candidates;

  /**
   * @brief Updates the gates for the given qubits
   * @details The function iterates over the DAG and updates the gates for the
   * given qubits as far es possible.
   * @param qubitsToUpdate The qubits that have been updated
   * @param commuteWith Gates the new gates should commute with
   */
  void updateByQubits(const std::set<Qubit>& qubitsToUpdate);
  /**
   * @brief Updates the candidates for the given qubits
   */
  void updateCandidatesByQubits(const std::set<Qubit>& qubitsToUpdate);
  /**
   * @brief Checks the candidates and add them to the gates if possible
   * @param qubitsToUpdate The qubits that have been updated
   */
  void candidatesToGates(const std::set<Qubit>& qubitsToUpdate);

  // Commutation checks
  static bool commutesWithAtQubit(const GateList&  layer,
                                  const Operation* opPointer,
                                  const Qubit&     qubit);
  static bool commuteAtQubit(const Operation* opPointer1,
                             const Operation* opPointer2, const Qubit& qubit);

public:
  // Constructor
  explicit NeutralAtomLayer(DAG dag) : dag(std::move(dag)) {
    for (auto& i : this->dag) {
      auto it = i.begin();
      this->iterators.emplace_back(it);
      this->candidates.emplace_back();
    }
  }

  /**
   * @brief Returns the current layer of gates
   * @return The current layer of gates
   */
  GateList getGates() { return gates; }
  /**
   * @brief Returns a vector of the iterator indices
   * @return A copy of the current iterator indices
   */
  std::vector<uint32_t> getIteratorOffset();
  /**
   * @brief Initializes the layer by updating all qubits starting from the
   * iterators
   * @param The iterator offset to start from
   */
  void initLayerOffset(const std::vector<uint32_t>& iteratorOffset = {});
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

} // namespace qc
