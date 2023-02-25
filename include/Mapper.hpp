//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "Architecture.hpp"
#include "MappingResults.hpp"
#include "QuantumComputation.hpp"
#include "configuration/Configuration.hpp"

#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_set>

constexpr std::int16_t  DEFAULT_POSITION  = -1;
constexpr double        INITIAL_FIDELITY  = 1.0;
constexpr std::uint16_t MAX_DEVICE_QUBITS = 128;

class Mapper {
protected:
  // internal structures

  /**
   * @brief Structure to store an operation on 1 or 2 logical qubits.
   *
   * For a single-qubit operation `control` is set to `-1`
   */
  struct Gate {
    std::int16_t  control = -1;
    std::uint16_t target  = 0;

    qc::Operation* op = nullptr;

    Gate(const std::int16_t c, const std::uint16_t t) : control(c), target(t){};
    Gate(const std::int16_t c, const std::uint16_t t, qc::Operation* operation)
        : control(c), target(t), op(operation){};

    [[nodiscard]] bool singleQubit() const { return control == -1; }
  };

  /**
   * @brief The quantum circuit to be mapped
   */
  qc::QuantumComputation qc;
  /**
   * @brief The quantum architecture on which to map the circuit
   */
  Architecture& architecture;

  /**
   * @brief The resulting quantum circuit after mapping
   */
  qc::QuantumComputation qcMapped;
  /**
   * @brief The gates of the circuit split into layers
   *
   * Each entry in the outer vector corresponds to 1 layer, containing all its
   * gates in an inner vector
   */
  std::vector<std::vector<Gate>> layers{};

  /**
   * @brief containing the logical qubit currently mapped to each physical
   * qubit. `qubits[physical_qubit] = logical_qubit`
   *
   * The inverse of `locations`
   */
  std::array<std::int16_t, MAX_DEVICE_QUBITS> qubits{};
  /**
   * @brief containing the logical qubit currently mapped to each physical
   * qubit. `locations[logical_qubit] = physical_qubit`
   *
   * The inverse of `qubits`
   */
  std::array<std::int16_t, MAX_DEVICE_QUBITS> locations{};
  std::array<double, MAX_DEVICE_QUBITS>       fidelities{};

  std::unordered_set<std::uint16_t> usedDeviceQubits{};

  MappingResults results{};

  /**
   * @brief Initialize the results structure with circuit names, registers in
   * the output circuit, gate counts, etc.
   */
  virtual void initResults();

  /**
   * @brief Splits the circuit into layers according to the method set in
   * `config.layering` and saves the result in `layers`
   *
   * methods of layering described in
   * https://iic.jku.at/files/eda/2019_dac_mapping_quantum_circuits_ibm_architectures_using_minimal_number_swap_h_gates.pdf
   *
   * Layering::IndividualGates/Layering::None -> each gate on separate layer
   * Layering::DisjointQubits -> each layer contains gates only acting on a
   * disjoint set of qubits 
   * Layering::OddGates -> always 2 gates per layer
   * (assigned by order of original gate index in the circuit)
   * Layering::QubitTriangle -> intended for architectures which contain
   * triangles of physical qubits, each layer only contains gates acting on 3
   * distinct qubits
   * Layering::Disjoint2qBlocks -> each layer contains 2Q-Blocks only acting on a
   * disjoint set of qubits
   */
  virtual void createLayers();

  /**
   * gates are put in the last layer (from the back of the circuit) in which
   * all of its qubits are not yet used by another gate in a circuit diagram
   * this can be thought of shifting all gates as far left as possible and
   * defining each column of gates as one layer.
   *
   * @param lastLayer the array storing the last layer each qubit is used in
   * @param control the (potential) control qubit of the gate
   * @param target the target qubit of the gate
   * @param gate the gate to be added to the layerh
   */
  void processDisjointQubitLayer(
      std::array<std::optional<std::size_t>, MAX_DEVICE_QUBITS>& lastLayer,
      const std::optional<std::uint16_t>& control, std::uint16_t target,
      qc::Operation* gate);
  
  /**
   * gates are put in the last layer (from the back of the circuit) in which
   * all of its qubits are not yet used by another gate in a circuit diagram,
   * except for gates acting on exactly matching qubit sets, which are 
   * collected in the same layer. Additionally also single qubit gates can be 
   * added to a layer already containing gates acting on that qubit.
   * This approach is equivalent to collecting all gates in 2Q-Blocks and then
   * applying DisjointQubitLayering.
   *
   * @param lastLayer the array storing the last layer each qubit is used in
   * @param control the (potential) control qubit of the gate
   * @param target the target qubit of the gate
   * @param gate the gate to be added to the layerh
   */
  void processDisjoint2qBlocksLayer(
      std::array<std::optional<std::size_t>, MAX_DEVICE_QUBITS>& lastLayer,
      const std::optional<std::uint16_t>& control, const std::uint16_t target,
      qc::Operation* gate);

  /**
   * @brief Get the index of the next layer after the given index containing a
   * gate acting on more than one qubit
   */
  [[gnu::pure]] virtual std::size_t getNextLayer(std::size_t idx);

  /**
   * @brief adding additional qubits to the result circuit if architecture has
   * more physical qubits than the original circuit has logical qubits
   */
  virtual void placeRemainingArchitectureQubits();

  /**
   * @brief finalizes the circuit after mapping
   * (e.g. adding unused qubits if architecture has more physical qubits than
   * mapped circuit has logical qubits)
   */
  virtual void finalizeMappedCircuit();

  /**
   * @brief count number of elementary gates and cnots in circuit and save the
   * results in `info.gates` and `info.cnots`
   */
  virtual void countGates(const qc::QuantumComputation& circuit,
                          MappingResults::CircuitInfo&  info) {
    countGates(circuit.cbegin(), circuit.cend(), info);
  }
  /**
   * @brief count number of elementary gates and cnots in circuit and save the
   * results in `info.gates` and `info.cnots`
   */
  virtual void countGates(decltype(qcMapped.cbegin())      it,
                          const decltype(qcMapped.cend())& end,
                          MappingResults::CircuitInfo&     info);

  /**
   * @brief performs optimizations on the circuit before mapping
   *
   * @param config contains settings of the current mapping run (e.g.
   * `config.preMappingOptimizations` controls if pre-mapping optimizations are
   * performed)
   */
  virtual void preMappingOptimizations(const Configuration& config);

  /**
   * @brief performs optimizations on the circuit before mapping
   *
   * @param config contains settings of the current mapping run (e.g.
   * `config.postMappingOptimizations` controls if post-mapping optimizations
   * are performed)
   */
  virtual void postMappingOptimizations(const Configuration& config);

public:
  Mapper(const qc::QuantumComputation& quantumComputation,
         Architecture&                 architecture);
  virtual ~Mapper() = default;

  /**
   * @brief map the circuit passed at initialization to the architecture
   *
   * @param config the settings for this mapping run (controls e.g. layering
   * methods, pre- and post-optimizations, etc.)
   */
  virtual void map(const Configuration& config) = 0;

  virtual void dumpResult(const std::string& outputFilename) {
    if (qcMapped.empty()) {
      std::cerr << "Mapped circuit is empty." << std::endl;
      return;
    }

    const std::size_t dot       = outputFilename.find_last_of('.');
    std::string       extension = outputFilename.substr(dot + 1);
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](unsigned char c) { return ::tolower(c); });
    if (extension == "real") {
      dumpResult(outputFilename, qc::Format::Real);
    } else if (extension == "qasm") {
      dumpResult(outputFilename, qc::Format::OpenQASM);
    } else {
      throw QMAPException("[dump] Extension " + extension +
                          " not recognized/supported for dumping.");
    }
  }

  virtual void dumpResult(const std::string& outputFilename,
                          qc::Format         format) {
    const std::size_t slash = outputFilename.find_last_of('/');
    const std::size_t dot   = outputFilename.find_last_of('.');
    results.output.name     = outputFilename.substr(slash + 1, dot - slash - 1);
    qcMapped.dump(outputFilename, format);
  }

  virtual void dumpResult(std::ostream& os, qc::Format format) {
    qcMapped.dump(os, format);
  }

  virtual std::ostream& printResult(std::ostream& out) {
    out << results.toString();
    return out;
  }

  virtual MappingResults& getResults() { return results; }

  virtual nlohmann::json json() { return results.json(); }

  virtual std::string csv() { return results.csv(); }

  std::ostream& printLayering(std::ostream& out) {
    out << "---------------- Layering -------------------" << std::endl;
    for (auto& layer : layers) {
      for (auto& gate : layer) {
        if (gate.singleQubit()) {
          out << "(" << gate.target << ") ";
        } else {
          out << "(" << gate.control << " " << gate.target << ") ";
        }
      }
      out << std::endl;
    }
    out << "---------------------------------------------" << std::endl;
    return out;
  }

  std::ostream& printLocations(std::ostream& out) {
    out << "---------------- Locations -------------------" << std::endl;
    for (std::size_t i = 0; i < qc.getNqubits(); ++i) {
      out << locations.at(i) << " ";
    }
    out << std::endl;
    out << "---------------------------------------------" << std::endl;
    return out;
  }
  std::ostream& printQubits(std::ostream& out) {
    out << "---------------- Qubits -------------------" << std::endl;
    for (std::size_t i = 0; i < architecture.getNqubits(); ++i) {
      out << qubits.at(i) << " ";
    }
    out << std::endl;
    out << "---------------------------------------------" << std::endl;
    return out;
  }

  virtual void reset() {
    architecture.reset();
    qc.reset();
    layers.clear();
    qubits.fill(DEFAULT_POSITION);
    locations.fill(DEFAULT_POSITION);
    usedDeviceQubits.clear();

    results = MappingResults();
  }
};
