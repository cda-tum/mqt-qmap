//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "Architecture.hpp"
#include "Mapper.hpp"
#include "MappingResults.hpp"
#include "QuantumComputation.hpp"

#include <fstream>
#include <string>

class DataLogger {
public:
  DataLogger(std::string path, Architecture& arch, qc::QuantumComputation qc)
      : dataLoggingPath(std::move(path)), architecture(&arch),
        nqubits(arch.getNqubits()), inputCircuit(std::move(qc)) {
    initLog();
    logArchitecture();
    logInputCircuit(inputCircuit);
    for (std::size_t i = 0; i < inputCircuit.getNqubits(); ++i) {
      qregs.emplace_back("q", "q[" + std::to_string(i) + "]");
    }
    for (std::size_t i = 0; i < inputCircuit.getNcbits(); ++i) {
      cregs.emplace_back("c", "c[" + std::to_string(i) + "]");
    }
  }

  void initLog();
  void clearLog();
  void logArchitecture();
  void logSearchNode(std::size_t layer, std::size_t nodeId,
                     std::size_t parentId, double costFixed, double costHeur,
                     double lookaheadPenalty,
                     const std::array<std::int16_t, MAX_DEVICE_QUBITS>& qubits,
                     bool validMapping, const std::vector<Exchange>& swaps,
                     std::size_t depth);
  void logFinalizeLayer(
      std::size_t layer, const qc::CompoundOperation& ops,
      const std::vector<std::uint16_t>& singleQubitMultiplicity,
      const std::map<std::pair<std::uint16_t, std::uint16_t>,
                     std::pair<std::uint16_t, std::uint16_t>>&
          twoQubitMultiplicity,
      const std::array<std::int16_t, MAX_DEVICE_QUBITS>& initialLayout,
      std::size_t finalNodeId, double finalCostFixed, double finalCostHeur,
      double finalLookaheadPenalty,
      const std::array<std::int16_t, MAX_DEVICE_QUBITS>& finalLayout,
      const std::vector<Exchange>& finalSwaps, std::size_t finalSearchDepth);
  void splitLayer();
  void logMappingResult(MappingResults& result);
  void logInputCircuit(qc::QuantumComputation& qc) {
    if (deactivated) {
      return;
    }
    qc.dump(dataLoggingPath + "/input.qasm", qc::Format::OpenQASM3);
  };
  void logOutputCircuit(qc::QuantumComputation& qc) {
    if (deactivated) {
      return;
    }
    qc.dump(dataLoggingPath + "/output.qasm", qc::Format::OpenQASM3);
  }
  void close();

protected:
  std::string dataLoggingPath;
  Architecture* architecture;
  std::uint16_t nqubits;
  qc::QuantumComputation inputCircuit;
  qc::RegisterNames qregs{};
  qc::RegisterNames cregs{};
  std::vector<std::ofstream> searchNodesLogFiles; // 1 per layer
  bool deactivated = false;

  void openNewLayer(std::size_t layer);
};
