//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/mqt-qmap for more
// information.
//

#include "Layer.hpp"
#include "Graph.hpp"

#include "Definitions.hpp"

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
  if (!(*op1)->actsOn(qubit) || !(*op2)->actsOn(qubit)) {
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
    return std::find(DIAGONAL_GATES.begin(), DIAGONAL_GATES.end(),
                     (*op2)->getType()) != DIAGONAL_GATES.end() ||
           ((*op2)->getType() == qc::OpType::Global &&
            std::find(
                DIAGONAL_GATES.begin(), DIAGONAL_GATES.end(),
                (dynamic_cast<GlobalOperation*>(op2->get())->getOpsType())) !=
                DIAGONAL_GATES.end());
    // true, iff qubit is a target and op2 is a diagonal gate
    //         ┌────┐
    // q: ──■──┤ RZ ├
    //      |  └────┘
  }
  // here: qubit is a target of op1
  if (const auto& controls2 = (*op2)->getControls();
      controls2.find(qubit) != controls2.end()) {
    return std::find(DIAGONAL_GATES.begin(), DIAGONAL_GATES.end(),
                     (*op1)->getType()) != DIAGONAL_GATES.end() ||
           ((*op1)->getType() == qc::OpType::Global &&
            std::find(
                DIAGONAL_GATES.begin(), DIAGONAL_GATES.end(),
                (dynamic_cast<GlobalOperation*>(op1->get())->getOpsType())) !=
                DIAGONAL_GATES.end());
    // true, iff qubit is a target and op1 is a diagonal gate and op2 is
    // controlled
    //    ┌────┐
    // q: ┤ RZ ├──■──
    //    └────┘  |
  }
  // here: qubit is a target of both operations
  if ((std::find(DIAGONAL_GATES.begin(), DIAGONAL_GATES.end(),
                 (*op1)->getType()) != DIAGONAL_GATES.end() ||
       ((*op1)->getType() == qc::OpType::Global &&
        std::find(DIAGONAL_GATES.begin(), DIAGONAL_GATES.end(),
                  (dynamic_cast<GlobalOperation*>(op1->get())->getOpsType())) !=
            DIAGONAL_GATES.end())) &&
      (std::find(DIAGONAL_GATES.begin(), DIAGONAL_GATES.end(),
                 (*op2)->getType()) != DIAGONAL_GATES.end() ||
       ((*op2)->getType() == qc::OpType::Global &&
        std::find(DIAGONAL_GATES.begin(), DIAGONAL_GATES.end(),
                  (dynamic_cast<GlobalOperation*>(op2->get())->getOpsType())) !=
            DIAGONAL_GATES.end()))) {
    // if both operations are diagonal gates
    //    ┌────┐┌────┐
    // q: ┤ RZ ├┤ RZ ├
    //    └────┘└────┘
    return true;
  }
  return ((*op1)->getType() == (*op2)->getType() &&
          (*op1)->getType() < qc::OpType::Compound) ||
         // FIXME: Beautify the hack with '< qc::OpType::Compound' to identify
         // FIXME: standard gates
         ((*op1)->getType() == qc::OpType::Global &&
          (*op2)->getType() == qc::OpType::Global &&
          (dynamic_cast<GlobalOperation*>(op1->get())->getOpsType() ==
           dynamic_cast<GlobalOperation*>(op2->get())->getOpsType()));
  // true, iff both operations are of the same type
  //    ┌───┐┌───┐
  // q: ┤ A ├┤ A ├
  //    └───┘└───┘
}
[[nodiscard]] auto Layer::isInverse(const std::unique_ptr<qc::Operation>* op1,
                                    const std::unique_ptr<qc::Operation>* op2)
    -> bool {
  // TODO: Add check for remaining operations
  if ((*op1)->getControls() == (*op2)->getControls() &&
      (*op1)->getTargets() == (*op2)->getTargets()) {
    if ((*op1)->getType() == qc::OpType::Global &&
        (*op2)->getType() == qc::OpType::Global) {
      const auto& gop1 = dynamic_cast<GlobalOperation*>(op1->get());
      const auto& gop2 = dynamic_cast<GlobalOperation*>(op2->get());
      return (gop1->getOpsType() == qc::OpType::I &&
              gop2->getOpsType() == qc::OpType::I) ||
             (gop1->getOpsType() == qc::OpType::X &&
              gop2->getOpsType() == qc::OpType::X) ||
             (gop1->getOpsType() == qc::OpType::Y &&
              gop2->getOpsType() == qc::OpType::Y) ||
             (gop1->getOpsType() == qc::OpType::Z &&
              gop2->getOpsType() == qc::OpType::Z) ||
             (gop1->getOpsType() == qc::OpType::S &&
              gop2->getOpsType() == qc::OpType::Sdg) ||
             (gop1->getOpsType() == qc::OpType::Sdg &&
              gop2->getOpsType() == qc::OpType::S) ||
             (gop1->getOpsType() == qc::OpType::SX &&
              gop2->getOpsType() == qc::OpType::SXdg) ||
             (gop1->getOpsType() == qc::OpType::SXdg &&
              gop2->getOpsType() == qc::OpType::SX) ||
             (gop1->getOpsType() == qc::OpType::T &&
              gop2->getOpsType() == qc::OpType::Tdg) ||
             (gop1->getOpsType() == qc::OpType::Tdg &&
              gop2->getOpsType() == qc::OpType::T) ||
             (gop1->getOpsType() == qc::OpType::H &&
              gop2->getOpsType() == qc::OpType::H) ||
             (gop1->getOpsType() == qc::OpType::P &&
              gop2->getOpsType() == qc::OpType::P &&
              std::abs(gop1->getParameter()[0] + gop2->getParameter()[0]) <
                  qc::PARAMETER_TOLERANCE) ||
             (gop1->getOpsType() == qc::OpType::RX &&
              gop2->getOpsType() == qc::OpType::RX &&
              std::abs(gop1->getParameter()[0] + gop2->getParameter()[0]) <
                  qc::PARAMETER_TOLERANCE) ||
             (gop1->getOpsType() == qc::OpType::RY &&
              gop2->getOpsType() == qc::OpType::RY &&
              std::abs(gop1->getParameter()[0] + gop2->getParameter()[0]) <
                  qc::PARAMETER_TOLERANCE) ||
             (gop1->getOpsType() == qc::OpType::RZ &&
              gop2->getOpsType() == qc::OpType::RZ &&
              std::abs(gop1->getParameter()[0] + gop2->getParameter()[0]) <
                  qc::PARAMETER_TOLERANCE);
    }
    return ((*op1)->getType() == qc::OpType::I &&
            (*op2)->getType() == qc::OpType::I) ||
           ((*op1)->getType() == qc::OpType::X &&
            (*op2)->getType() == qc::OpType::X) ||
           ((*op1)->getType() == qc::OpType::Y &&
            (*op2)->getType() == qc::OpType::Y) ||
           ((*op1)->getType() == qc::OpType::Z &&
            (*op2)->getType() == qc::OpType::Z) ||
           ((*op1)->getType() == qc::OpType::S &&
            (*op2)->getType() == qc::OpType::Sdg) ||
           ((*op1)->getType() == qc::OpType::Sdg &&
            (*op2)->getType() == qc::OpType::S) ||
           ((*op1)->getType() == qc::OpType::SX &&
            (*op2)->getType() == qc::OpType::SXdg) ||
           ((*op1)->getType() == qc::OpType::SXdg &&
            (*op2)->getType() == qc::OpType::SX) ||
           ((*op1)->getType() == qc::OpType::T &&
            (*op2)->getType() == qc::OpType::Tdg) ||
           ((*op1)->getType() == qc::OpType::Tdg &&
            (*op2)->getType() == qc::OpType::T) ||
           ((*op1)->getType() == qc::OpType::H &&
            (*op2)->getType() == qc::OpType::H) ||
           ((*op1)->getType() == qc::OpType::P &&
            (*op2)->getType() == qc::OpType::P &&
            std::abs((*op1)->getParameter()[0] + (*op2)->getParameter()[0]) <
                qc::PARAMETER_TOLERANCE) ||
           ((*op1)->getType() == qc::OpType::RX &&
            (*op2)->getType() == qc::OpType::RX &&
            std::abs((*op1)->getParameter()[0] + (*op2)->getParameter()[0]) <
                qc::PARAMETER_TOLERANCE) ||
           ((*op1)->getType() == qc::OpType::RY &&
            (*op2)->getType() == qc::OpType::RY &&
            std::abs((*op1)->getParameter()[0] + (*op2)->getParameter()[0]) <
                qc::PARAMETER_TOLERANCE) ||
           ((*op1)->getType() == qc::OpType::RZ &&
            (*op2)->getType() == qc::OpType::RZ &&
            std::abs((*op1)->getParameter()[0] + (*op2)->getParameter()[0]) <
                qc::PARAMETER_TOLERANCE);
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
          if (!currentGroup[qubit].empty() &&
              !commutesAtQubit(currentGroup[qubit][0]->getOperation(),
                               current->getOperation(), qubit)) {
            // here: the current operation does not commute with the current
            // group members and is not the inverse of the lookahead
            // --> start a new group
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
      if (!currentGroup[qubit].empty() &&
          !commutesAtQubit(currentGroup[qubit][0]->getOperation(),
                           current->getOperation(), qubit)) {
        // here: the current operation does not commute with the current
        // group members and is not the inverse of the lookahead
        // --> start a new group
        // TODO: Can this be beautified with the copy-assign operator?
        predecessorGroup[qubit].clear();
        for (const auto& v : currentGroup[qubit]) {
          predecessorGroup[qubit].emplace_back(v);
        }
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
[[nodiscard]] auto Layer::constructInteractionGraph(qc::OpType opType, Number nctrl) const
    -> Graph<std::shared_ptr<DAGVertex>> {
  switch (opType) {
  case qc::X:
  case qc::Y:
  case qc::Z:
  case qc::RX:
  case qc::RY:
  case qc::RZ:
    if (nctrl == 1) {
      break;
    }
    [[fallthrough]];
  default:
    std::stringstream ss;
    ss << "The operation type ";
    for (std::size_t i = 0; i < nctrl; ++i) {
      ss << "c";
    }
    ss << qc::toString(opType)
       << " is not supported for constructing an interaction graph.";
    throw std::invalid_argument(ss.str());
  }
  Graph<std::shared_ptr<DAGVertex>> graph;
  for (const auto& vertex : *executableSet) {
    const auto& gate = *vertex->getOperation();
    if (gate->getType() == opType && gate->getNcontrols() == nctrl) {
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
} // namespace na
