/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#include "configuration/Configuration.hpp"

#include <iostream>
#include <sstream>
#include <string>

#ifndef QMAP_MAPPINGRESULTS_HPP
    #define QMAP_MAPPINGRESULTS_HPP

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

    MappingResults()          = default;
    virtual ~MappingResults() = default;

    virtual void copyInput(const MappingResults& mappingResults) {
        input        = mappingResults.input;
        architecture = mappingResults.architecture;
        config       = mappingResults.config;
        output       = mappingResults.output;
    }

    [[nodiscard]] std::string toString() const {
        return json().dump(2);
    }

    [[nodiscard]] virtual nlohmann::json json() const {
        nlohmann::json resultJSON{};
        auto&          circuit        = resultJSON["circuit"];
        circuit["name"]               = input.name;
        circuit["qubits"]             = input.qubits;
        circuit["gates"]              = input.gates;
        circuit["single_qubit_gates"] = input.singleQubitGates;
        circuit["cnots"]              = input.cnots;
        circuit["layers"]             = input.layers;

        auto& mapped_circuit                 = resultJSON["mapped_circuit"];
        mapped_circuit["name"]               = output.name;
        mapped_circuit["qubits"]             = output.qubits;
        mapped_circuit["gates"]              = output.gates;
        mapped_circuit["single_qubit_gates"] = output.singleQubitGates;
        mapped_circuit["cnots"]              = output.cnots;
        mapped_circuit["swaps"]              = output.swaps;
        if (!mappedCircuit.empty()) {
            mapped_circuit["qasm"] = mappedCircuit;
        }
        if (config.method == Method::Exact) {
            mapped_circuit["direction_reverse"] = output.directionReverse;
        } else if (config.method == Method::Heuristic) {
            mapped_circuit["teleportations"] = output.teleportations;
        }

        resultJSON["config"] = config.json();

        auto& stats               = resultJSON["statistics"];
        stats["timeout"]          = timeout;
        stats["mapping_time"]     = time;
        stats["additional_gates"] = output.gates - input.gates;
        stats["arch"]             = architecture;

        return resultJSON;
    }

    virtual std::string csv() {
        std::stringstream ss{};
        ss << input.name << ";"
           << input.qubits << ";"
           << input.gates << ";"
           << input.singleQubitGates << ";"
           << input.cnots << ";"
           << architecture << ";"
           << output.name << ";"
           << output.qubits << ";"
           << output.gates << ";"
           << output.singleQubitGates << ";"
           << output.cnots << ";"
           << output.swaps << ";"
           << output.directionReverse << ";"
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

#endif //QMAP_MAPPINGRESULTS_HPP
