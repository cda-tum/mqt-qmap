//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "hybridmap/AodScheduler.hpp"

#include "Definitions.hpp"
#include "QuantumComputation.hpp"
#include "hybridmap/NeutralAtomDefinitions.hpp"
#include "hybridmap/NeutralAtomUtils.hpp"
#include "operations/AodOperation.hpp"
#include "operations/OpType.hpp"
#include "utils.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <set>
#include <utility>
#include <vector>

namespace qc {

QuantumComputation AodScheduler::schedule(QuantumComputation& qc) {
  initMoveGroups(qc);
  if (moveGroups.empty()) {
    return qc;
  }
  processMoveGroups();

  // create new quantum circuit and insert AOD operations at the correct
  // indices
  auto     groupIt = moveGroups.begin();
  uint32_t idx     = 0;
  for (const auto& op : qc) {
    if (groupIt != moveGroups.end() && idx == groupIt->getFirstIdx()) {
      // add move group
      for (auto& aodOp : groupIt->processedOpsInit) {
        qcScheduled.emplace_back(std::make_unique<AodOperation>(aodOp));
      }
      qcScheduled.emplace_back(
          std::make_unique<AodOperation>(groupIt->processedOpShuttle));
      for (auto& aodOp : groupIt->processedOpsFinal) {
        qcScheduled.emplace_back(std::make_unique<AodOperation>(aodOp));
      }
      groupIt++;
    } else if (op->getType() != OpType::Move) {
      qcScheduled.emplace_back(op->clone());
    }
    idx++;
  }

  return qcScheduled;
}

void AodScheduler::initMoveGroups(QuantumComputation& qc) {
  MoveGroup       currentMoveGroup{arch};
  MoveGroup const lastMoveGroup{arch};
  uint32_t        idx = 0;
  for (auto& op : qc) {
    if (op->getType() == OpType::Move) {
      AtomMove const move{op->getTargets()[0], op->getTargets()[1]};
      if (currentMoveGroup.canAdd(move)) {
        currentMoveGroup.add(move, idx);
      } else {
        moveGroups.push_back(std::move(currentMoveGroup));
        currentMoveGroup = MoveGroup{arch};
        currentMoveGroup.add(move, idx);
      }
    } else if (op->getNqubits() > 1 && !currentMoveGroup.moves.empty()) {
      for (const auto& qubit : op->getUsedQubits()) {
        if (std::find(currentMoveGroup.qubitsUsedByGates.begin(),
                      currentMoveGroup.qubitsUsedByGates.end(),
                      qubit) == currentMoveGroup.qubitsUsedByGates.end()) {
          currentMoveGroup.qubitsUsedByGates.push_back(qubit);
        }
      }
    }
    idx++;
  }
  if (!currentMoveGroup.moves.empty()) {
    moveGroups.push_back(std::move(currentMoveGroup));
  }
}

bool AodScheduler::MoveGroup::canAdd(const AtomMove& move) {
  // if move would move a qubit that is used by a gate in this move group
  // return false
  if (std::find(qubitsUsedByGates.begin(), qubitsUsedByGates.end(),
                move.first) != qubitsUsedByGates.end()) {
    return false;
  }
  // checks if the op can be executed in parallel
  auto moveVector = arch.getVector(move.first, move.second);
  return std::all_of(
      moves.begin(), moves.end(),
      [&moveVector, this](const std::pair<AtomMove, uint32_t> opPair) {
        auto moveGroup = opPair.first;
        auto opVector  = arch.getVector(moveGroup.first, moveGroup.second);
        return parallelCheck(moveVector, opVector);
      });
}

bool AodScheduler::MoveGroup::parallelCheck(const MoveVector& v1,
                                            const MoveVector& v2) {
  if (!v1.overlap(v2)) {
    return true;
  }
  // overlap -> check if same direction
  if (v1.direction != v2.direction) {
    return false;
  }
  // same direction -> check if include
  if (v1.include(v2) || v2.include(v1)) {
    return false;
  }
  return true;
}

void AodScheduler::MoveGroup::add(const AtomMove& move, const uint32_t idx) {
  moves.emplace_back(move, idx);
  qubitsUsedByGates.push_back(move.second);
}

void AodScheduler::AodActivationHelper::addActivation(
    std::pair<ActivationMergeType, ActivationMergeType> merge,
    const Coordinate& origin, const AtomMove& move, MoveVector v) {
  const auto x         = origin.getX();
  const auto y         = origin.getY();
  const auto signX     = v.direction.getSignX();
  const auto signY     = v.direction.getSignY();
  const auto deltaX    = v.xEnd - v.xStart;
  const auto deltaY    = v.yEnd - v.yStart;
  auto       aodMovesX = getAodMovesFromInit(Dimension::X, x);
  auto       aodMovesY = getAodMovesFromInit(Dimension::Y, y);

  switch (merge.first) {
  case ActivationMergeType::Trivial:
    switch (merge.second) {
    case ActivationMergeType::Trivial:
      allActivations.emplace_back(
          AodActivation{{x, deltaX, signX}, {y, deltaY, signY}, move});
      break;
    case ActivationMergeType::Merge:
      mergeActivationDim(Dimension::Y,
                         AodActivation{Dimension::Y, {y, deltaY, signY}, move},
                         AodActivation{Dimension::X, {x, deltaX, signX}, move});
      break;
    case ActivationMergeType::Append:
      allActivations.emplace_back(
          AodActivation{{x, deltaX, signX}, {y, deltaY, signY}, move});
      aodMovesY = getAodMovesFromInit(Dimension::Y, y);
      reAssignOffsets(aodMovesY, signY);
      break;
    default:
      break;
    }
    break;
  case ActivationMergeType::Merge:
    switch (merge.second) {
    case ActivationMergeType::Trivial:
      mergeActivationDim(Dimension::X,
                         AodActivation{Dimension::X, {x, deltaX, signX}, move},
                         AodActivation{Dimension::Y, {y, deltaY, signY}, move});
      break;
    case ActivationMergeType::Merge:
      throw QMAPException("Merge in both dimensions should never happen.");
    case ActivationMergeType::Append:
      mergeActivationDim(Dimension::X,
                         AodActivation{Dimension::X, {x, deltaX, signX}, move},
                         AodActivation{Dimension::Y, {y, deltaY, signY}, move});
      //      allActivations.emplace_back(
      //          AodActivation{Dimension::Y, {y, deltaY, signY}, move});
      aodMovesY = getAodMovesFromInit(Dimension::Y, y);
      reAssignOffsets(aodMovesY, signY);
      break;
    default:
      break;
    }
    break;
  case ActivationMergeType::Append:
    switch (merge.second) {
    case ActivationMergeType::Trivial:
      allActivations.emplace_back(
          AodActivation{{x, deltaX, signX}, {y, deltaY, signY}, move});
      aodMovesX = getAodMovesFromInit(Dimension::X, x);
      reAssignOffsets(aodMovesX, signX);
      break;
    case ActivationMergeType::Merge:
      mergeActivationDim(Dimension::Y,
                         AodActivation{Dimension::Y, {y, deltaY, signY}, move},
                         AodActivation{Dimension::X, {x, deltaX, signX}, move});
      aodMovesX = getAodMovesFromInit(Dimension::X, x);
      reAssignOffsets(aodMovesX, signX);
      break;
    case ActivationMergeType::Append:
      allActivations.emplace_back(
          AodActivation{{x, deltaX, signX}, {y, deltaY, signY}, move});
      aodMovesX = getAodMovesFromInit(Dimension::X, x);
      reAssignOffsets(aodMovesX, signX);
      aodMovesY = getAodMovesFromInit(Dimension::Y, y);
      reAssignOffsets(aodMovesY, signY);
      break;
    default:
      break;
    }
    break;
  default:
    break;
  }
}

std::pair<ActivationMergeType, ActivationMergeType>
AodScheduler::AodActivationHelper::canAddActivation(
    const qc::Coordinate& origin, qc::MoveVector v) const {
  auto aodMovesX = getAodMovesFromInit(Dimension::X, origin.getX());
  auto aodMovesY = getAodMovesFromInit(Dimension::Y, origin.getY());

  auto const canX = canAddActivationDim(Dimension::X, origin, v);
  auto const canY = canAddActivationDim(Dimension::Y, origin, v);
  return std::make_pair(canX, canY);
}

ActivationMergeType AodScheduler::AodActivationHelper::canAddActivationDim(
    Dimension dim, const Coordinate& origin, MoveVector v) const {
  auto x = dim == Dimension::X ? origin.getX() : origin.getY();
  auto sign =
      dim == Dimension::X ? v.direction.getSignX() : v.direction.getSignY();
  auto delta    = dim == Dimension::X ? v.xEnd - v.xStart : v.yEnd - v.yStart;
  auto aodMoves = getAodMovesFromInit(dim, x);
  if (aodMoves.empty()) {
    // return false as no merge is required
    return ActivationMergeType::Trivial;
  }
  // check if it can be combined with existing activations
  for (auto& aodMove : aodMoves) {
    if (aodMove->init == x && std::abs(aodMove->delta - delta) < 0.0001 &&
        aodMove->offset == sign) {
      // combine activations
      return ActivationMergeType::Merge;
    }
  }

  // check if increase in x direction possible
  if (checkIntermediateSpaceAtInit(dim, x, sign)) {
    return ActivationMergeType::Append;
  }
  return ActivationMergeType::Impossible;
}

void AodScheduler::AodActivationHelper::reAssignOffsets(
    std::vector<std::shared_ptr<AodMove>>& aodMoves, int32_t sign) {
  std::sort(
      aodMoves.begin(), aodMoves.end(),
      [](const std::shared_ptr<AodMove>& a, const std::shared_ptr<AodMove> b) {
        return std::abs(a->delta) < std::abs(b->delta);
      });
  int32_t offset = sign;
  for (auto& aodMove : aodMoves) {
    // same sign
    if (aodMove->delta * sign >= 0) {
      aodMove->offset = offset;
      offset += sign;
    }
  }
}

void AodScheduler::processMoveGroups() {
  // convert the moves from MoveGroup to AodOperations
  for (auto groupIt = moveGroups.begin(); groupIt != moveGroups.end();
       ++groupIt) {
    AodActivationHelper   aodActivationHelper{arch, OpType::AodActivate};
    AodActivationHelper   aodDeactivationHelper{arch, OpType::AodDeactivate};
    MoveGroup             possibleNewMoveGroup{arch};
    std::vector<AtomMove> movesToRemove;
    for (auto& movePair : groupIt->moves) {
      auto& move              = movePair.first;
      auto  idx               = movePair.second;
      auto  origin            = arch.getCoordinate(move.first);
      auto  target            = arch.getCoordinate(move.second);
      auto  v                 = arch.getVector(move.first, move.second);
      auto  vReverse          = arch.getVector(move.second, move.first);
      auto activationCanAddXY = aodActivationHelper.canAddActivation(origin, v);
      auto deactivationCanAddXY =
          aodDeactivationHelper.canAddActivation(target, vReverse);
      if (activationCanAddXY.first == ActivationMergeType::Impossible ||
          activationCanAddXY.second == ActivationMergeType::Impossible ||
          deactivationCanAddXY.first == ActivationMergeType::Impossible ||
          deactivationCanAddXY.second == ActivationMergeType::Impossible) {
        // move could not be added as not sufficient intermediate levels
        // add new move group and add move to it
        possibleNewMoveGroup.add(move, idx);
        movesToRemove.push_back(move);
      } else {
        aodActivationHelper.addActivation(activationCanAddXY, origin, move, v);
        aodDeactivationHelper.addActivation(deactivationCanAddXY, target, move,
                                            vReverse);
      }
    }
    // remove from current move group
    for (const auto& moveToRemove : movesToRemove) {
      groupIt->moves.erase(
          std::remove_if(groupIt->moves.begin(), groupIt->moves.end(),
                         [&moveToRemove](const auto& movePair) {
                           return movePair.first == moveToRemove;
                         }),
          groupIt->moves.end());
    }
    if (!possibleNewMoveGroup.moves.empty()) {
      groupIt = moveGroups.insert(groupIt + 1, std::move(possibleNewMoveGroup));
      possibleNewMoveGroup = MoveGroup{arch};
      groupIt--;
    }
    groupIt->processedOpsInit   = aodActivationHelper.getAodOperations();
    groupIt->processedOpsFinal  = aodDeactivationHelper.getAodOperations();
    groupIt->processedOpShuttle = MoveGroup::connectAodOperations(
        groupIt->processedOpsInit, groupIt->processedOpsFinal);
  }
}

AodOperation AodScheduler::MoveGroup::connectAodOperations(
    const std::vector<AodOperation>& opsInit,
    const std::vector<AodOperation>& opsFinal) {
  // for each init operation find the corresponding final operation
  // and connect with an aod move operations
  // all can be done in parallel in a single move
  std::vector<SingleOperation> aodOperations;
  std::set<CoordIndex>         targetQubits;

  for (const auto& opInit : opsInit) {
    if (opInit.getType() == OpType::AodMove) {
      for (const auto& opFinal : opsFinal) {
        if (opFinal.getType() == OpType::AodMove) {
          if (opInit.getTargets().size() <= 1 ||
              opFinal.getTargets().size() <= 1) {
            throw QFRException("AodScheduler::MoveGroup::connectAodOperations: "
                               "AodMove operation with less than 2 targets");
          }
          if (opInit.getTargets() == opFinal.getTargets()) {
            targetQubits.insert(opInit.getTargets().begin(),
                                opInit.getTargets().end());
            // found corresponding final operation
            // connect with aod move
            const auto startXs = opInit.getEnds(Dimension::X);
            const auto endXs   = opFinal.getStarts(Dimension::X);
            const auto startYs = opInit.getEnds(Dimension::Y);
            const auto endYs   = opFinal.getStarts(Dimension::Y);
            if (!startXs.empty() && !endXs.empty()) {
              for (size_t i = 0; i < startXs.size(); i++) {
                const auto startX = startXs[i];
                const auto endX   = endXs[i];
                if (std::abs(startX - endX) > 0.0001) {
                  aodOperations.emplace_back(Dimension::X, startX, endX);
                }
              }
            }
            if (!startYs.empty() && !endYs.empty()) {
              for (size_t i = 0; i < startYs.size(); i++) {
                const auto startY = startYs[i];
                const auto endY   = endYs[i];
                if (std::abs(startY - endY) > 0.0001) {
                  aodOperations.emplace_back(Dimension::Y, startY, endY);
                }
              }
            }
          }
        }
      }
    }
  }
  std::vector<CoordIndex> targetQubitsVec;
  targetQubitsVec.reserve(targetQubits.size());
  for (const auto& qubit : targetQubits) {
    targetQubitsVec.push_back(qubit);
  }
  return {OpType::AodMove, targetQubitsVec, aodOperations};
}

std::vector<std::shared_ptr<AodScheduler::AodActivationHelper::AodMove>>
AodScheduler::AodActivationHelper::getAodMovesFromInit(Dimension dim,
                                                       uint32_t  init) const {
  std::vector<std::shared_ptr<AodMove>> aodMoves;
  for (const auto& activation : allActivations) {
    for (auto& aodMove : activation.getActivates(dim)) {
      if (aodMove->init == init) {
        aodMoves.push_back(aodMove);
      }
    }
  }
  return aodMoves;
}

uint32_t AodScheduler::AodActivationHelper::getMaxOffsetAtInit(
    Dimension dim, uint32_t init, int32_t sign) const {
  auto aodMoves = getAodMovesFromInit(dim, init);
  if (aodMoves.empty()) {
    return 0;
  }
  uint32_t maxOffset = 0;
  for (const auto& aodMove : aodMoves) {
    auto offset = aodMove->offset;
    if (offset * sign >= 0) {
      maxOffset = std::max(maxOffset, static_cast<uint32_t>(std::abs(offset)));
    }
  }
  return maxOffset;
}

bool AodScheduler::AodActivationHelper::checkIntermediateSpaceAtInit(
    Dimension dim, uint32_t init, int32_t sign) const {
  uint32_t neighborX = init;
  if (sign > 0) {
    neighborX += 1;
  } else {
    neighborX -= 1;
  }
  auto aodMoves         = getAodMovesFromInit(dim, init);
  auto aodMovesNeighbor = getAodMovesFromInit(dim, neighborX);
  if (aodMoves.empty() && aodMovesNeighbor.empty()) {
    return true;
  }
  if (aodMoves.empty()) {
    return getMaxOffsetAtInit(dim, neighborX, sign) <
           arch.getNAodIntermediateLevels();
  }
  if (aodMovesNeighbor.empty()) {
    return getMaxOffsetAtInit(dim, init, sign) <
           arch.getNAodIntermediateLevels();
  }
  return getMaxOffsetAtInit(dim, init, sign) +
             getMaxOffsetAtInit(dim, neighborX, sign) <
         arch.getNAodIntermediateLevels();
}

void AodScheduler::AodActivationHelper::mergeActivationDim(
    Dimension dim, const AodActivation& activationDim,
    const AodActivation& activationOtherDim) {
  // merge activations
  for (auto& activationCurrent : allActivations) {
    auto activates = activationDim.getActivates(dim);
    for (auto& aodMove : activates) {
      if (aodMove->init == activationDim.getActivates(dim)[0]->init &&
          aodMove->delta == activationDim.getActivates(dim)[0]->delta &&
          aodMove->offset == activationDim.getActivates(dim)[0]->offset) {
        // append move
        activationCurrent.moves.push_back(activationDim.moves[0]);
        // add activation in the other dimension
        if (dim == Dimension::X) {
          activationCurrent.activateYs.push_back(
              activationOtherDim.activateYs[0]);
        } else {
          activationCurrent.activateXs.push_back(
              activationOtherDim.activateXs[0]);
        }
        return;
      }
    }
  }
}

std::pair<AodOperation, AodOperation>
AodScheduler::AodActivationHelper::getAodOperation(
    const AodActivationHelper::AodActivation& activation) const {
  std::vector<CoordIndex> qubitsActivation;
  qubitsActivation.reserve(activation.moves.size());
  for (const auto& move : activation.moves) {
    if (type == OpType::AodActivate) {
      qubitsActivation.push_back(move.first);
    } else {
      qubitsActivation.push_back(move.second);
    }
  }
  std::vector<CoordIndex> qubitsMove;
  qubitsMove.reserve(activation.moves.size() * 2);
  for (const auto& move : activation.moves) {
    if (std::find(qubitsMove.begin(), qubitsMove.end(), move.first) ==
        qubitsMove.end()) {
      qubitsMove.push_back(move.first);
    }
    if (std::find(qubitsMove.begin(), qubitsMove.end(), move.second) ==
        qubitsMove.end()) {
      qubitsMove.push_back(move.second);
    }
  }

  std::vector<SingleOperation> initOperations;
  std::vector<SingleOperation> offsetOperations;

  auto d      = this->arch.getInterQubitDistance();
  auto interD = this->arch.getInterQubitDistance() /
                this->arch.getNAodIntermediateLevels();

  for (const auto& aodMove : activation.activateXs) {
    initOperations.emplace_back(Dimension::X,
                                static_cast<fp>(aodMove->init) * d,
                                static_cast<fp>(aodMove->init) * d);
    if (type == OpType::AodActivate) {
      offsetOperations.emplace_back(
          Dimension::X, static_cast<fp>(aodMove->init) * d,
          static_cast<fp>(aodMove->init) * d +
              static_cast<fp>(aodMove->offset) * interD);
    } else {
      offsetOperations.emplace_back(Dimension::X,
                                    static_cast<fp>(aodMove->init) * d +
                                        static_cast<fp>(aodMove->offset) *
                                            interD,
                                    static_cast<fp>(aodMove->init) * d);
    }
  }
  for (const auto& aodMove : activation.activateYs) {
    initOperations.emplace_back(Dimension::Y,
                                static_cast<fp>(aodMove->init) * d,
                                static_cast<fp>(aodMove->init) * d);
    if (type == OpType::AodActivate) {
      offsetOperations.emplace_back(
          Dimension::Y, static_cast<fp>(aodMove->init) * d,
          static_cast<fp>(aodMove->init) * d +
              static_cast<fp>(aodMove->offset) * interD);
    } else {
      offsetOperations.emplace_back(Dimension::Y,
                                    static_cast<fp>(aodMove->init) * d +
                                        static_cast<fp>(aodMove->offset) *
                                            interD,
                                    static_cast<fp>(aodMove->init) * d);
    }
  }

  auto initOp   = AodOperation(type, qubitsActivation, initOperations);
  auto offsetOp = AodOperation(OpType::AodMove, qubitsMove, offsetOperations);
  if (this->type == OpType::AodActivate) {
    return std::make_pair(initOp, offsetOp);
  }
  return std::make_pair(offsetOp, initOp);
}

std::vector<AodOperation>
AodScheduler::AodActivationHelper::getAodOperations() const {
  std::vector<AodOperation> aodOperations;
  for (const auto& activation : allActivations) {
    auto operations = getAodOperation(activation);
    aodOperations.emplace_back(std::move(operations.first));
    aodOperations.emplace_back(std::move(operations.second));
  }
  return aodOperations;
}

} // namespace qc
