//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/mqt-qmap for more
// information.
//

#include "Layer.hpp"

#include "Definitions.hpp"
#include "Graph.hpp"
#include "operations/OpType.hpp"

#include <algorithm>
#include <iostream>
#include <memory>
#include <set>

namespace na {

/// Checks whether the two operations commute on the given qubit.
[[nodiscard]] auto
Layer::commutesAtQubit(const std::unique_ptr<qc::Operation>* op1,
                       const std::unique_ptr<qc::Operation>* op2,
                       const qc::Qubit&                      qubit) -> bool {
  // FIXME: Cleanup this conditional branch, but how?
  // check whether both operations act on the given qubit
  if (!(*op1)->actsOn(qubit) or !(*op2)->actsOn(qubit)) {
    throw std::invalid_argument(
        "The operations do not act on the given qubit.");
  }
  if (const auto& controls1 = (*op1)->getControls();
      controls1.find(qubit) != controls1.end()) {
    // if op1 is controlled on the given qubit
    if (const auto& controls2 = (*op2)->getControls();
        controls2.find(qubit) != controls2.end()) {
      // if op2 is controlled on the given qubit
      // q: ──■────■──
      //      |    |
      return true;
    }
    // here: qubit is a target of op2
    return isDiagonal((*op2)->getType());
    // true, iff qubit is a target and op2 is a diagonal gate
    //         ┌────┐
    // q: ──■──┤ RZ ├
    //      |  └────┘
  }
  // here: qubit is a target of op1
  if (const auto& controls2 = (*op2)->getControls();
      controls2.find(qubit) != controls2.end()) {
    return isDiagonal((*op1)->getType());
    // true, iff qubit is a target and op1 is a diagonal gate and op2 is
    // controlled
    //    ┌────┐
    // q: ┤ RZ ├──■──
    //    └────┘  |
  }
  // here: qubit is a target of both operations
  if (isDiagonal((*op1)->getType()) and isDiagonal((*op2)->getType())) {
    // if both operations are diagonal gates
    //    ┌────┐┌────┐
    // q: ┤ RZ ├┤ RZ ├
    //    └────┘└────┘
    return true;
  }
  return (*op1)->getType() == (*op2)->getType() and
         (*op1)->getType() < qc::OpType::Compound;
  // FIXME: Beautify the hack with '< qc::OpType::Compound' to identify
  // FIXME: standard gates
  // true, iff both operations are of the same type
  //    ┌───┐┌───┐
  // q: ┤ A ├┤ A ├
  //    └───┘└───┘
}
[[nodiscard]] auto Layer::isInverse(const std::unique_ptr<qc::Operation>* op1,
                                    const std::unique_ptr<qc::Operation>* op2)
    -> bool {
  // TODO: Add check for remaining operations
  if ((*op1)->getControls() == (*op2)->getControls() and
      (*op1)->getTargets() == (*op2)->getTargets()) {
    auto result =
        ((*op1)->getType() == qc::OpType::I and
         (*op2)->getType() == qc::OpType::I) or
        ((*op1)->getType() == qc::OpType::X and
         (*op2)->getType() == qc::OpType::X) or
        ((*op1)->getType() == qc::OpType::Y and
         (*op2)->getType() == qc::OpType::Y) or
        ((*op1)->getType() == qc::OpType::Z and
         (*op2)->getType() == qc::OpType::Z) or
        ((*op1)->getType() == qc::OpType::S and
         (*op2)->getType() == qc::OpType::Sdg) or
        ((*op1)->getType() == qc::OpType::Sdg and
         (*op2)->getType() == qc::OpType::S) or
        ((*op1)->getType() == qc::OpType::SX and
         (*op2)->getType() == qc::OpType::SXdg) or
        ((*op1)->getType() == qc::OpType::SXdg and
         (*op2)->getType() == qc::OpType::SX) or
        ((*op1)->getType() == qc::OpType::T and
         (*op2)->getType() == qc::OpType::Tdg) or
        ((*op1)->getType() == qc::OpType::Tdg and
         (*op2)->getType() == qc::OpType::T) or
        ((*op1)->getType() == qc::OpType::H and
         (*op2)->getType() == qc::OpType::H) or
        ((*op1)->getType() == qc::OpType::P and
         (*op2)->getType() == qc::OpType::P and
         std::abs((*op1)->getParameter()[0] + (*op2)->getParameter()[0]) <
             qc::PARAMETER_TOLERANCE) or
        ((*op1)->getType() == qc::OpType::RX and
         (*op2)->getType() == qc::OpType::RX and
         std::abs((*op1)->getParameter()[0] + (*op2)->getParameter()[0]) <
             qc::PARAMETER_TOLERANCE) or
        ((*op1)->getType() == qc::OpType::RY and
         (*op2)->getType() == qc::OpType::RY and
         std::abs((*op1)->getParameter()[0] + (*op2)->getParameter()[0]) <
             qc::PARAMETER_TOLERANCE) or
        ((*op1)->getType() == qc::OpType::RZ and
         (*op2)->getType() == qc::OpType::RZ and
         std::abs((*op1)->getParameter()[0] + (*op2)->getParameter()[0]) <
             qc::PARAMETER_TOLERANCE);
    return result;
  }
  return false;
}
auto Layer::constructDAG(const qc::QuantumComputation& qc) -> void {
  const auto nQubits = qc.getNqubits();
  // For a pair of self-canceling operations like two consecutive X operations
  // or RY rotations with opposite angles the first operations is a
  // destructive operation that disables operations until the consecutive
  // constructive operation enables them again
  // ---
  // those that add a (+) edge to the current group members
  std::vector<std::vector<std::shared_ptr<DAGVertex>>> constructive(nQubits);
  // those that add a (-) edge to the current group members
  std::vector<std::vector<std::shared_ptr<DAGVertex>>> destructive(nQubits);
  // those that are already in the current group where all gates commute on
  // this qubit
  std::vector<std::vector<std::shared_ptr<DAGVertex>>> currentGroup(nQubits);
  // lookahead of 1 serves as a buffer for the next operation on each qubit
  std::vector<std::shared_ptr<DAGVertex>> lookahead(nQubits);
  // the predecessor of the current group members
  std::vector<std::vector<std::shared_ptr<DAGVertex>>> predecessorGroup(
      nQubits);
  // all operations acting on a qubit (processed so far) excluding
  // constructive and destructive operations
  std::vector<std::vector<std::shared_ptr<DAGVertex>>> qubitOperations(nQubits);
  // iterate over all operations in the quantum circuit
  for (const auto& op : qc) {
    // create a vertex for the current operation
    std::shared_ptr<DAGVertex> const vertex =
        DAGVertex::create(&op, &executableSet);
    // iterate over all qubits the operation acts on
    for (const auto& qubit : op->getUsedQubits()) {
      // check whether the lookahead is empty
      if (lookahead[qubit] == nullptr) {
        // here: the lookahead is empty
        // add the current operation to the lookahead
        lookahead[qubit] = vertex;
      } else {
        // here: the lookahead is not empty
        // get the current vertex from the lookahead and store the new
        // vertex in the lookahead
        std::shared_ptr<DAGVertex> const current = lookahead[qubit];
        lookahead[qubit]                         = vertex;
        // check whether the current operation is the inverse of the
        // lookahead
        if (isInverse(current->getOperation(),
                      lookahead[qubit]->getOperation())) {
          // here: the current operation is the inverse of the lookahead
          // add an enabling edge from the lookahead to all operations on this
          // qubit including the destructive ones
          for (const auto& qubitOperation : qubitOperations[qubit]) {
            lookahead[qubit]->addEnabledSuccessor(qubitOperation);
          }
          for (const auto& qubitOperation : destructive[qubit]) {
            lookahead[qubit]->addEnabledSuccessor(qubitOperation);
          }
          // add the lookahead to the constructive group
          constructive[qubit].emplace_back(lookahead[qubit]);
          // add a disabling edge to all operations on this qubit including
          // the destructive ones
          for (const auto& qubitOperation : qubitOperations[qubit]) {
            current->addDisabledSuccessor(qubitOperation);
          }
          for (const auto& qubitOperation : destructive[qubit]) {
            current->addDisabledSuccessor(qubitOperation);
          }
          // add an enabling edge to the lookahead
          current->addEnabledSuccessor(lookahead[qubit]);
          // add the current vertex to the destructive group
          destructive[qubit].emplace_back(current);
          // clear the lookahead
          lookahead[qubit] = nullptr;
        } else {
          // add an enabling edge from each constructive operation
          for (const auto& constructiveOp : constructive[qubit]) {
            constructiveOp->addEnabledSuccessor(current);
          }
          // add a disabling edge from each destructive operation
          for (const auto& destructiveOp : destructive[qubit]) {
            destructiveOp->addDisabledSuccessor(current);
          }
          // check whether the current operation commutes with the current
          // group members
          if (!currentGroup[qubit].empty() and
              !commutesAtQubit(currentGroup[qubit][0]->getOperation(),
                               current->getOperation(), qubit)) {
            // here: the current operation does not commute with the current
            // group members and is not the inverse of the lookahead
            // --> start a new group
            predecessorGroup[qubit].clear();
            predecessorGroup[qubit] = currentGroup[qubit];
            currentGroup[qubit].clear();
          }
          // add an enabling edge from each predecessor
          for (const auto& predecessor : predecessorGroup[qubit]) {
            predecessor->addEnabledSuccessor(current);
          }
          // add the current vertex to the current group
          currentGroup[qubit].emplace_back(current);
          qubitOperations[qubit].emplace_back(current);
        }
      }
    }
  }
  // process the remaining lookahead for every qubit
  for (qc::Qubit qubit = 0; qubit < nQubits; ++qubit) {
    if (lookahead[qubit] != nullptr) {
      auto const current = lookahead[qubit];
      lookahead[qubit]   = nullptr;
      // add an enabling edge from each constructive operation
      for (const auto& constructiveOp : constructive[qubit]) {
        constructiveOp->addEnabledSuccessor(current);
      }
      // add a disabling edge from each destructive operation
      for (const auto& destructiveOp : destructive[qubit]) {
        destructiveOp->addDisabledSuccessor(current);
      }
      // check whether the current operation commutes with the current
      // group members
      if (!currentGroup[qubit].empty() and
          !commutesAtQubit(currentGroup[qubit][0]->getOperation(),
                           current->getOperation(), qubit)) {
        // here: the current operation does not commute with the current
        // group members and is not the inverse of the lookahead
        // --> start a new group
        predecessorGroup[qubit].clear();
        predecessorGroup[qubit] = currentGroup[qubit];
        currentGroup[qubit].clear();
      }
      // add an enabling edge from each predecessor
      for (const auto& predecessor : predecessorGroup[qubit]) {
        predecessor->addEnabledSuccessor(current);
      }
      // add the current vertex to the current group
      currentGroup[qubit].emplace_back(current);
      qubitOperations[qubit].emplace_back(current);
    }
  }
}
[[nodiscard]] auto Layer::constructInteractionGraph(OpType opType) const
    -> Graph<std::shared_ptr<DAGVertex>> {
  switch (opType.type) {
  case qc::X:
  case qc::Y:
  case qc::Z:
  case qc::RX:
  case qc::RY:
  case qc::RZ:
    if (opType.nctrl == 1) {
      break;
    }
    [[fallthrough]];
  default:
    std::stringstream ss;
    ss << "The operation type ";
    for (std::size_t i = 0; i < opType.nctrl; ++i) {
      ss << "c";
    }
    ss << opType << " is not supported for constructing an interaction graph.";
    throw std::invalid_argument(ss.str());
  }
  Graph<std::shared_ptr<DAGVertex>> graph;
  for (const auto& vertex : *executableSet) {
    const auto& gate = *vertex->getOperation();
    if (gate->getType() == opType.type and
        gate->getNcontrols() == opType.nctrl) {
      const auto& usedQubits = gate->getUsedQubits();
      if (usedQubits.size() != 2) {
        throw std::invalid_argument(
            "The interaction graph can only be constructed for two-qubit "
            "gates.");
      }
      std::vector<qc::Qubit> q(usedQubits.cbegin(), usedQubits.cend());
      graph.addEdge(q[0], q[1], vertex);
    }
  }
  return graph;
}
[[nodiscard]] auto Layer::getExecutablesOfType(OpType opType) const
    -> std::vector<std::shared_ptr<DAGVertex>> {
  std::vector<std::shared_ptr<DAGVertex>> executables;
  for (const auto& vertex : *executableSet) {
    if ((*vertex->getOperation())->getType() == opType.type and
        (*vertex->getOperation())->getNcontrols() == opType.nctrl) {
      executables.emplace_back(vertex);
    }
  }
  return executables;
}
} // namespace na
