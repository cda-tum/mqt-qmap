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
    std::string   name{};
    std::uint16_t qubits           = 0;
    std::size_t   gates            = 0;
    std::size_t   singleQubitGates = 0;
    std::size_t   cnots            = 0;
    std::size_t   layers           = 0;
    double        totalFidelity    = 1.;
    double        totalLogFidelity = 0.;

    // info in output circuit
    std::size_t swaps            = 0;
    std::size_t directionReverse = 0;
    std::size_t teleportations   = 0;
  };

  struct HeuristicBenchmarkInfo {
    std::size_t expandedNodes            = 0;
    std::size_t generatedNodes           = 0;
    std::size_t solutionDepth            = 0;
    double      timePerNode              = 0.;
    double      averageBranchingFactor   = 0.;
    double      effectiveBranchingFactor = 0.;
  };

  CircuitInfo input{};

  std::string   architecture{};
  Configuration config{};

  double time    = 0.0;
  bool   timeout = true;

  CircuitInfo output{};
  std::string mappedCircuit{};

  std::string wcnf{};

  HeuristicBenchmarkInfo              heuristicBenchmark{};
  std::vector<HeuristicBenchmarkInfo> layerHeuristicBenchmark{};

  MappingResults()          = default;
  virtual ~MappingResults() = default;

  virtual void copyInput(const MappingResults& mappingResults) {
    input                   = mappingResults.input;
    architecture            = mappingResults.architecture;
    config                  = mappingResults.config;
    output                  = mappingResults.output;
    wcnf                    = mappingResults.wcnf;
    heuristicBenchmark      = mappingResults.heuristicBenchmark;
    layerHeuristicBenchmark = mappingResults.layerHeuristicBenchmark;
  }

  [[nodiscard]] std::string toString() const { return json().dump(2); }

  [[nodiscard]] virtual nlohmann::json json() const {
    nlohmann::json resultJSON{};

    auto& circuit                 = resultJSON["circuit"];
    circuit["name"]               = input.name;
    circuit["qubits"]             = input.qubits;
    circuit["gates"]              = input.gates;
    circuit["single_qubit_gates"] = input.singleQubitGates;
    circuit["cnots"]              = input.cnots;

    auto& mappedCirc                 = resultJSON["mapped_circuit"];
    mappedCirc["name"]               = output.name;
    mappedCirc["qubits"]             = output.qubits;
    mappedCirc["gates"]              = output.gates;
    mappedCirc["single_qubit_gates"] = output.singleQubitGates;
    mappedCirc["cnots"]              = output.cnots;
    if (!mappedCircuit.empty()) {
      mappedCirc["qasm"] = mappedCircuit;
    }

    resultJSON["config"] = config.json();

    auto& stats             = resultJSON["statistics"];
    stats["timeout"]        = timeout;
    stats["mapping_time"]   = time;
    stats["arch"]           = architecture;
    stats["layers"]         = input.layers;
    stats["swaps"]          = output.swaps;
    stats["total_fidelity"] = output.totalFidelity;
    if (config.method == Method::Exact) {
      stats["direction_reverse"] = output.directionReverse;
      if (config.includeWCNF && !wcnf.empty()) {
        stats["WCNF"] = wcnf;
      }
    } else if (config.method == Method::Heuristic) {
      stats["teleportations"]      = output.teleportations;
      auto& benchmark              = stats["benchmark"];
      benchmark["expanded_nodes"]  = heuristicBenchmark.expandedNodes;
      benchmark["generated_nodes"] = heuristicBenchmark.generatedNodes;
      benchmark["time_per_node"]   = heuristicBenchmark.timePerNode;
      benchmark["average_branching_factor"] =
          heuristicBenchmark.averageBranchingFactor;
      benchmark["effective_branching_factor"] =
          heuristicBenchmark.effectiveBranchingFactor;
    }
    stats["additional_gates"] =
        static_cast<std::make_signed_t<decltype(output.gates)>>(output.gates) -
        static_cast<std::make_signed_t<decltype(input.gates)>>(input.gates);

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
