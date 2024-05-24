//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "hybridmap/NeutralAtomScheduler.hpp"

#include "CircuitOptimizer.hpp"

#include <string>

qc::SchedulerResults
qc::NeutralAtomScheduler::schedule(const qc::QuantumComputation&     qc,
                                   const std::map<HwQubit, HwQubit>& initHwPos,
                                   bool verbose, bool createAnimationCsv,
                                   fp shuttlingSpeedFactor) {
  if (verbose) {
    std::cout << "\n* schedule start!\n";
  }

  std::vector<fp> totalExecutionTimes(arch.getNpositions(), 0);
  // saves for each coord the time slots that are blocked by a multi qubit gate
  std::vector<std::deque<std::pair<fp, fp>>> rydbergBlockedQubitsTimes(
      arch.getNpositions(), std::deque<std::pair<fp, fp>>());
  fp aodLastBlockedTime  = 0;
  fp totalGateTime       = 0;
  fp totalGateFidelities = 1;

  AnimationAtoms animationAtoms(initHwPos, arch);
  if (createAnimationCsv) {
    animationCsv += animationAtoms.getInitString();
    animationArchitectureCsv = arch.getAnimationCsv();
  }

  int      index        = 0;
  int      nAodActivate = 0;
  uint32_t nCZs         = 0;
  for (const auto& op : qc) {
    index++;
    if (verbose) {
      std::cout << "\n" << index << "\n";
    }
    if (op->getType() == qc::AodActivate) {
      nAodActivate++;
    } else if (op->getType() == qc::OpType::Z && op->getNcontrols() == 1) {
      nCZs++;
    }

    auto qubits = op->getUsedQubits();
    auto opTime = arch.getOpTime(op.get());
    if (op->getType() == qc::AodMove || op->getType() == qc::AodActivate ||
        op->getType() == qc::AodDeactivate) {
      opTime *= shuttlingSpeedFactor;
    }
    auto opFidelity = arch.getOpFidelity(op.get());

    // DEBUG info
    if (verbose) {
      std::cout << op->getName() << "  ";
      for (const auto& qubit : qubits) {
        std::cout << "q" << qubit << " ";
      }
      std::cout << "-> time: " << opTime << ", fidelity: " << opFidelity
                << "\n";
    }

    fp maxTime = 0;
    if (op->getType() == qc::AodMove || op->getType() == qc::AodActivate ||
        op->getType() == qc::AodDeactivate) {
      // AodBlocking
      maxTime = aodLastBlockedTime;
      for (const auto& qubit : qubits) {
        maxTime = std::max(maxTime, totalExecutionTimes[qubit]);
      }
      aodLastBlockedTime = maxTime + opTime;
    } else if (qubits.size() > 1) {
      // multi qubit gates -> take into consideration blocking
      auto rydbergBlockedQubits = arch.getBlockedCoordIndices(op.get());
      // get max execution time over all blocked qubits
      bool rydbergBlocked = true;
      while (rydbergBlocked) {
        // get regular max execution time
        for (const auto& qubit : qubits) {
          maxTime = std::max(maxTime, totalExecutionTimes[qubit]);
        }
        // check if all blocked qubits are free at maxTime
        rydbergBlocked = false;
        for (const auto& qubit : rydbergBlockedQubits) {
          // check if qubit is blocked at maxTime
          for (const auto& startEnd : rydbergBlockedQubitsTimes[qubit]) {
            auto start = startEnd.first;
            auto end   = startEnd.second;
            if ((start <= maxTime && end > maxTime) ||
                (start <= maxTime + opTime && end > maxTime + opTime)) {
              rydbergBlocked = true;
              // update maxTime to the end of the blocking
              maxTime = end;
              // remove the blocking
              break;
            }
            // skip rest of the times
            if (end > maxTime) {
              break;
            }
          }
        }
      }

      for (const auto& qubit : rydbergBlockedQubits) {
        rydbergBlockedQubitsTimes[qubit].emplace_back(maxTime,
                                                      maxTime + opTime);
      }

    } else {
      // other operations -> no blocking
      // get max execution time over all qubits
      for (const auto& qubit : qubits) {
        maxTime = std::max(maxTime, totalExecutionTimes[qubit]);
        // remove all blocked times that are smaller than maxTime
        while (!rydbergBlockedQubitsTimes[qubit].empty() &&
               rydbergBlockedQubitsTimes[qubit].front().second < maxTime) {
          rydbergBlockedQubitsTimes[qubit].pop_front();
        }
      }
    }
    // update total execution times
    for (const auto& qubit : qubits) {
      totalExecutionTimes[qubit] = maxTime + opTime;
    }

    totalGateFidelities *= opFidelity;
    totalGateTime += opTime;
    if (verbose) {
      std::cout << "\n";
      printTotalExecutionTimes(totalExecutionTimes, rydbergBlockedQubitsTimes);
    }

    // update animation
    if (createAnimationCsv) {
      animationCsv +=
          animationAtoms.createCsvOp(op, maxTime, maxTime + opTime, arch);
    }
  }
  if (verbose) {
    std::cout << "\n* schedule end!\n";
    std::cout << "nAodActivate: " << nAodActivate << "\n";
  }

  const auto maxExecutionTime =
      *std::max_element(totalExecutionTimes.begin(), totalExecutionTimes.end());
  const auto totalIdleTime =
      maxExecutionTime * arch.getNqubits() - totalGateTime;
  const auto totalFidelities =
      totalGateFidelities *
      std::exp(-totalIdleTime / arch.getDecoherenceTime());

  if (createAnimationCsv) {
    animationCsv += animationAtoms.getEndString(maxExecutionTime);
  }
  if (verbose) {
    printSchedulerResults(totalExecutionTimes, totalIdleTime,
                          totalGateFidelities, totalFidelities, nCZs);
  }
  return {maxExecutionTime, totalIdleTime, totalGateFidelities, totalFidelities,
          nCZs};
}

void qc::NeutralAtomScheduler::printSchedulerResults(
    std::vector<fp>& totalExecutionTimes, fp totalIdleTime,
    fp totalGateFidelities, fp totalFidelities, uint32_t nCZs) {
  auto totalExecutionTime =
      *std::max_element(totalExecutionTimes.begin(), totalExecutionTimes.end());
  std::cout << "\ntotalExecutionTimes: " << totalExecutionTime << "\n";
  std::cout << "totalIdleTime: " << totalIdleTime << "\n";
  std::cout << "totalGateFidelities: " << totalGateFidelities << "\n";
  std::cout << "totalFidelities: " << totalFidelities << "\n";
  std::cout << "totalNumCZs: " << nCZs << "\n";
}

void qc::NeutralAtomScheduler::printTotalExecutionTimes(
    std::vector<fp>&                            totalExecutionTimes,
    std::vector<std::deque<std::pair<fp, fp>>>& blockedQubitsTimes) {
  std::cout << "ExecutionTime: "
            << "\n";
  for (size_t qubit = 0; qubit < totalExecutionTimes.size(); qubit++) {
    std::cout << "[" << qubit << "] " << totalExecutionTimes[qubit] << " \t";
    for (const auto& blockedTime : blockedQubitsTimes[qubit]) {
      std::cout << blockedTime.first << "-" << blockedTime.second << " \t";
    }
    std::cout << "\n";
  }
}
qc::AnimationAtoms::AnimationAtoms(const std::map<HwQubit, HwQubit>&  initHwPos,
                                   const qc::NeutralAtomArchitecture& arch) {
  auto nCols = arch.getNcolumns();

  for (const auto& [id, coord] : initHwPos) {
    coordIdxToId[coord] = id;
    idToCoord[id]       = {coord % nCols * arch.getInterQubitDistance(),
                           static_cast<fp>(coord) / nCols *
                               arch.getInterQubitDistance()};
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
    throw std::runtime_error(
        "Tried to add axis but axis already exists for qubit " +
        std::to_string(id));
  }
  return axesIds[id];
}
qc::AnimationAtoms::marginId qc::AnimationAtoms::addMargin(qc::HwQubit id) {
  if (marginIds.find(id) == marginIds.end()) {
    marginIdCounter++;
    marginIds[id] = marginIdCounter;
  } else {
    throw std::runtime_error(
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
    if (coordIdxToId.find(coordIdx) == coordIdxToId.end()) {
      // check if qubit was moved to coordIdx, if yes, update coordIdxToId
      if (op->getType() == OpType::AodDeactivate) {
        for (const auto& [id, coord] : idToCoord) {
          if (std::abs(coord.first - coordIdx % arch.getNcolumns() *
                                         arch.getInterQubitDistance()) <
                  0.0001 &&
              std::abs(coord.second -
                       static_cast<fp>(coordIdx) / arch.getNcolumns() *
                           arch.getInterQubitDistance()) < 0.0001) {
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
      } else {
        continue;
      }
    }
    auto id    = coordIdxToId.at(coordIdx);
    auto coord = idToCoord.at(id);
    if (op->getType() == OpType::AodActivate) {
      addAxis(id);
      csvLine += createCsvLine(startTime, id, coord.first, coord.second, 1,
                               colorSlm, true, axesIds.at(id));
      csvLine += createCsvLine(endTime, id, coord.first, coord.second, 1,
                               colorAod, true, axesIds.at(id));
    } else if (op->getType() == OpType::AodDeactivate) {
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
