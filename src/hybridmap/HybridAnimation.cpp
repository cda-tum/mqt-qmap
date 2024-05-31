//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "hybridmap/HybridAnimation.hpp"

#include "Definitions.hpp"
#include "hybridmap/MoveToAodConverter.hpp"
#include "hybridmap/NeutralAtomArchitecture.hpp"
#include "operations/OpType.hpp"
#include "utils.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>

qc::AnimationAtoms::AnimationAtoms(const std::map<HwQubit, HwQubit>&  initHwPos,
                                   const qc::NeutralAtomArchitecture& arch) {
  auto nCols = arch.getNcolumns();

  for (const auto& [id, coord] : initHwPos) {
    coordIdxToId[coord] = id;
    auto column         = coord % nCols;
    auto row            = coord / nCols;
    idToCoord[id]       = {column * arch.getInterQubitDistance(),
                           row * arch.getInterQubitDistance()};
  }
}

std::string qc::AnimationAtoms::getInitString() {
  std::string initString;
  initString +=
      "time;id;x;y;size;fill;color;axes;axesId;margin;marginId;marginSize\n";
  for (const auto& [id, coord] : idToCoord) {
    initString += "0.000;" + std::to_string(id) + ";" +
                  std::to_string(coord.first) + ";" +
                  std::to_string(coord.second) + ";1;0;0;0;0;0;0;0\n";
  }
  return initString;
}

std::string qc::AnimationAtoms::getEndString(fp endTime) {
  std::string initString;
  for (const auto& [id, coord] : idToCoord) {
    initString += std::to_string(endTime) + ";" + std::to_string(id) + ";" +
                  std::to_string(coord.first) + ";" +
                  std::to_string(coord.second) + ";1;0;0;0;0;0;0;0\n";
  }
  return initString;
}

qc::AnimationAtoms::axesId qc::AnimationAtoms::addAxis(qc::HwQubit id) {
  if (axesIds.find(id) == axesIds.end()) {
    axesIdCounter++;
    axesIds[id] = axesIdCounter;
  } else {
    throw QMAPException("Tried to add axis but axis already exists for qubit " +
                        std::to_string(id));
  }
  return axesIds[id];
}
qc::AnimationAtoms::marginId qc::AnimationAtoms::addMargin(qc::HwQubit id) {
  if (marginIds.find(id) == marginIds.end()) {
    marginIdCounter++;
    marginIds[id] = marginIdCounter;
  } else {
    throw QMAPException(
        "Tried to add margin but margin already exists for qubit " +
        std::to_string(id));
  }
  return marginIds[id];
}

std::string
qc::AnimationAtoms::createCsvOp(const std::unique_ptr<qc::Operation>& op,
                                fp startTime, fp endTime,
                                const qc::NeutralAtomArchitecture& arch) {
  std::string csvLine;

  for (const auto& coordIdx : op->getUsedQubits()) {
    // if coordIdx unmapped -> continue except it is an AodDeactivate
    if (OpType::AodDeactivate != op->getType() &&
        coordIdxToId.find(coordIdx) == coordIdxToId.end()) {
      continue;
    }
    if (op->getType() == OpType::AodDeactivate) {
      // check if there is a qubit at coordIdx
      // if yes -> update coordIdxToId with new coordIdx
      // if not -> throw exception
      for (const auto& idAndCoord : idToCoord) {
        auto id    = idAndCoord.first;
        auto coord = idAndCoord.second;
        auto col   = coordIdx % arch.getNcolumns();
        auto row   = coordIdx / arch.getNcolumns();
        if (std::abs(coord.first - col * arch.getInterQubitDistance()) <
                0.0001 &&
            std::abs(coord.second - row * arch.getInterQubitDistance()) <
                0.0001) {
          // remove old coordIdx with same id
          for (const auto& [oldCoordIdx, oldId] : coordIdxToId) {
            if (oldId == id) {
              coordIdxToId.erase(oldCoordIdx);
              break;
            }
          }
          // add new coordIdx with id
          coordIdxToId[coordIdx] = id;
          break;
        }
      }
    }
    if (coordIdxToId.find(coordIdx) == coordIdxToId.end() ||
        idToCoord.find(coordIdxToId.at(coordIdx)) == idToCoord.end()) {
      throw QMAPException("Tried to create csv line for qubit at coordIdx " +
                          std::to_string(coordIdx) +
                          " but there is no qubit at this coordIdx");
    }
    auto id    = coordIdxToId.at(coordIdx);
    auto coord = idToCoord.at(id);
    if (op->getType() == OpType::AodActivate) {
      addAxis(id);
      if (axesIds.find(id) == axesIds.end()) {
        throw QMAPException(
            "Tried to activate qubit at coordIdx " + std::to_string(coordIdx) +
            " but there is no axis for qubit " + std::to_string(id));
      }
      csvLine += createCsvLine(startTime, id, coord.first, coord.second, 1,
                               colorSlm, true, axesIds.at(id));
      csvLine += createCsvLine(endTime, id, coord.first, coord.second, 1,
                               colorAod, true, axesIds.at(id));
    } else if (op->getType() == OpType::AodDeactivate) {
      if (axesIds.find(id) == axesIds.end()) {
        throw QMAPException("Tried to deactivate qubit at coordIdx " +
                            std::to_string(coordIdx) +
                            " but there is no axis for qubit " +
                            std::to_string(id));
      }
      csvLine += createCsvLine(startTime, id, coord.first, coord.second, 1,
                               colorAod, true, axesIds.at(id));
      csvLine += createCsvLine(endTime, id, coord.first, coord.second, 1,
                               colorSlm, true, axesIds.at(id));
      removeAxis(id);

    } else if (op->getType() == OpType::AodMove) {
      if (axesIds.find(id) == axesIds.end()) {
        throw QMAPException(
            "Tried to move qubit at coordIdx " + std::to_string(coordIdx) +
            " but there is no axis for qubit " + std::to_string(id));
      }
      csvLine += createCsvLine(startTime, id, coord.first, coord.second, 1,
                               colorAod, true, axesIds.at(id));

      // update atom coordinates
      auto startsX =
          dynamic_cast<AodOperation*>(op.get())->getStarts(Dimension::X);
      auto endsX = dynamic_cast<AodOperation*>(op.get())->getEnds(Dimension::X);
      auto startsY =
          dynamic_cast<AodOperation*>(op.get())->getStarts(Dimension::Y);
      auto endsY = dynamic_cast<AodOperation*>(op.get())->getEnds(Dimension::Y);
      for (size_t i = 0; i < startsX.size(); i++) {
        if (std::abs(startsX[i] - coord.first) < 0.0001) {
          coord.first = endsX[i];
        }
      }
      for (size_t i = 0; i < startsY.size(); i++) {
        if (std::abs(startsY[i] - coord.second) < 0.0001) {
          coord.second = endsY[i];
        }
      }
      // save new coordinates
      idToCoord[id] = coord;
      csvLine += createCsvLine(endTime, id, coord.first, coord.second, 1,
                               colorAod, true, axesIds.at(id));
    } else if (op->getUsedQubits().size() > 1) { // multi qubit gates
      addMargin(id);
      csvLine += createCsvLine(startTime, id, coord.first, coord.second, 1,
                               colorSlm, false, 0, false, marginIds.at(id));
      auto midTime = (startTime + endTime) / 2;
      csvLine +=
          createCsvLine(midTime, id, coord.first, coord.second, 1, colorCz,
                        false, 0, true, marginIds.at(id),
                        arch.getBlockingFactor() * arch.getInteractionRadius() *
                            arch.getInterQubitDistance());
      csvLine += createCsvLine(endTime, id, coord.first, coord.second, 1,
                               colorSlm, false, 0, false, marginIds.at(id));
      removeMargin(id);

    } else { // single qubit gates
      csvLine +=
          createCsvLine(startTime, id, coord.first, coord.second, 1, colorSlm);
      auto midTime = (startTime + endTime) / 2;
      csvLine +=
          createCsvLine(midTime, id, coord.first, coord.second, 1, colorLocal);
      csvLine +=
          createCsvLine(endTime, id, coord.first, coord.second, 1, colorSlm);
    }
  }
  return csvLine;
}
std::string qc::AnimationAtoms::createCsvLine(
    qc::fp startTime, qc::HwQubit id, qc::fp x, qc::fp y, uint32_t size,
    uint32_t color, bool axes, qc::AnimationAtoms::axesId axId, bool margin,
    qc::AnimationAtoms::marginId marginId, fp marginSize) {
  std::string csvLine;
  csvLine +=
      std::to_string(startTime) + ";" + std::to_string(id) + ";" +
      std::to_string(x) + ";" + std::to_string(y) + ";" + std::to_string(size) +
      ";" + std::to_string(color) + ";" + std::to_string(color) + ";" +
      std::to_string(static_cast<int>(axes)) + ";" + std::to_string(axId) +
      ";" + std::to_string(static_cast<int>(margin)) + ";" +
      std::to_string(marginId) + ";" + std::to_string(marginSize) + "\n";
  return csvLine;
}
