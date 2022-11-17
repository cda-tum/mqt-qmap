//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "Architecture.hpp"
#include "configuration/Configuration.hpp"

#include <iostream>
#include <sstream>
#include <string>

#pragma once

struct MappingResults {
  struct CircuitInfo {
    // general info
    std::string    name{};
    unsigned short qubits           = 0;
    unsigned long  gates            = 0;
    unsigned long  singleQubitGates = 0;
    unsigned long  cnots            = 0;
    unsigned long  layers           = 0;

    // info in output circuit
    unsigned long swaps            = 0;
    unsigned long directionReverse = 0;
    unsigned long teleportations   = 0;
  };

  CircuitInfo input{};

  std::string   architecture{};
  Configuration config{};

  double time    = 0.0;
  bool   timeout = true;

  CircuitInfo output{};
  std::string mappedCircuit{};

  std::string wcnf{};

  MappingResults()          = default;
  virtual ~MappingResults() = default;

  virtual void copyInput(const MappingResults& mappingResults) {
    input        = mappingResults.input;
    architecture = mappingResults.architecture;
    config       = mappingResults.config;
    output       = mappingResults.output;
    wcnf         = mappingResults.wcnf;
  }

  [[nodiscard]] std::string toString() const { return json().dump(2); }

  [[nodiscard]] virtual nlohmann::json json() const {
    nlohmann::json resultJSON{};
    auto&          circuit        = resultJSON["circuit"];
    circuit["name"]               = input.name;
    circuit["qubits"]             = input.qubits;
    circuit["gates"]              = input.gates;
    circuit["single_qubit_gates"] = input.singleQubitGates;
    circuit["cnots"]              = input.cnots;

    auto& mapped_circuit                 = resultJSON["mapped_circuit"];
    mapped_circuit["name"]               = output.name;
    mapped_circuit["qubits"]             = output.qubits;
    mapped_circuit["gates"]              = output.gates;
    mapped_circuit["single_qubit_gates"] = output.singleQubitGates;
    mapped_circuit["cnots"]              = output.cnots;
    if (!mappedCircuit.empty()) {
      mapped_circuit["qasm"] = mappedCircuit;
    }

    resultJSON["config"] = config.json();

    auto& stats           = resultJSON["statistics"];
    stats["timeout"]      = timeout;
    stats["mapping_time"] = time;
    stats["arch"]         = architecture;
    stats["layers"]       = input.layers;
    stats["swaps"]        = output.swaps;
    if (config.method == Method::Exact) {
      stats["direction_reverse"] = output.directionReverse;
      if (config.includeWCNF && !wcnf.empty()) {
        stats["WCNF"] = wcnf;
      }
    } else if (config.method == Method::Heuristic) {
      stats["teleportations"] = output.teleportations;
    }
    stats["additional_gates"] = static_cast<long>(output.gates) - input.gates;

    return resultJSON;
  }

  virtual std::string csv() {
    std::stringstream ss{};
    ss << input.name << ";" << input.qubits << ";" << input.gates << ";"
       << input.singleQubitGates << ";" << input.cnots << ";" << architecture
       << ";" << output.name << ";" << output.qubits << ";" << output.gates
       << ";" << output.singleQubitGates << ";" << output.cnots << ";"
       << output.swaps << ";" << output.directionReverse << ";"
       << output.teleportations << ";";
    if (timeout) {
      ss << "TO";
    } else {
      ss << time;
    }
    ss << ";";
    return ss.str();
  }
};
