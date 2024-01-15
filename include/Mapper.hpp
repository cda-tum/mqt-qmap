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
#include <set>
#include <map>

/**
 * number of two-qubit gates acting on pairs of logical qubits in some layer
 * where the keys correspond to logical qubit pairs ({q1, q2}, with q1<=q2)
 * and the values to the number of gates acting on a pair in each direction
 * (the first number with control=q1, target=q2 and the second the reverse).
 *
 * e.g., with multiplicity {{0,1},{2,3}} there are 2 gates with logical
 * qubit 0 as control and qubit 1 as target, and 3 gates with 1 as control
 * and 0 as target.
 */
using TwoQubitMultiplicity =
    std::map<Edge, std::pair<std::uint16_t, std::uint16_t>>;

/**
 * number of single-qubit gates acting on each logical qubit in some
 * layer.
 *
 * e.g. with multiplicity {1,0,2} there is 1 1Q-gate acting on q0, no 1Q-gates
 * acting on q1, and 2 1Q-gates acting on q2
 */
using SingleQubitMultiplicity = std::vector<std::uint16_t>;

constexpr std::int16_t DEFAULT_POSITION = -1;

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
  Architecture* architecture;

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
   * @brief The number of 1Q-gates acting on each logical qubit in each layer
   */
  std::vector<SingleQubitMultiplicity> singleQubitMultiplicities{};

  /**
   * @brief The number of 2Q-gates acting on each pair of logical qubits in each
   * layer
   */
  std::vector<TwoQubitMultiplicity> twoQubitMultiplicities{};

  /**
   * @brief For each layer the set of all logical qubits, which are acted on by
   * a gate in the layer
   */
  std::vector<std::set<std::uint16_t>> activeQubits{};

  /**
   * @brief For each layer the set of all logical qubits, which are acted on by
   * a 1Q-gate in the layer
   */
  std::vector<std::set<std::uint16_t>> activeQubits1QGates{};

  /**
   * @brief For each layer the set of all logical qubits, which are acted on by
   * a 2Q-gate in the layer
   */
  std::vector<std::set<std::uint16_t>> activeQubits2QGates{};

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
   * Layering::IndividualGates -> each gate on separate layer
   * Layering::DisjointQubits -> each layer contains gates only acting on a
   * disjoint set of qubits
   * Layering::OddGates -> always 2 gates per layer
   * (assigned by order of original gate index in the circuit)
   * Layering::QubitTriangle -> intended for architectures which contain
   * triangles of physical qubits, each layer only contains gates acting on 3
   * distinct qubits
   * Layering::Disjoint2qBlocks -> each layer contains 2Q-Blocks only acting on
   * a disjoint set of qubits
   */
  virtual void createLayers();

  /**
   * @brief Returns true if the layer at the given index can be split into two
   * without resulting in an empty layer (assuming the original layer only has
   * disjoint 2Q-gate-blocks)
   *
   * @param index the index of the layer to be split
   * @return true if the layer is splittable
   * @return false if splitting the layer will result in an empty layer
   */
  virtual bool isLayerSplittable(std::size_t index);

  /**
   * @brief Splits the layer at the given index into two layers with half as
   * many qubits acted on by gates in each layer
   *
   * @param index the index of the layer to be split
   * @param arch architecture on which the circuit is mapped
   */
  virtual void splitLayer(std::size_t index, Architecture& arch);

  /**
   * gates are put in the last layer (from the back of the circuit) in which
   * all of its qubits are not yet used by another gate in a circuit diagram
   * this can be thought of shifting all gates as far left as possible and
   * defining each column of gates as one layer.
   *
   * @param lastLayer the array storing the last layer each qubit is used in
   * @param control the (potential) control qubit of the gate
   * @param target the target qubit of the gate
   * @param gate the gate to be added to the layer
   * @param collect2qBlocks if true, gates are collected in 2Q-blocks, and
   * layering is performed on these blocks
   */
  void processDisjointQubitLayer(
      std::array<std::optional<std::size_t>, MAX_DEVICE_QUBITS>& lastLayer,
      const std::optional<std::uint16_t>& control, std::uint16_t target,
      qc::Operation* gate);

  /**
   * Similar to processDisjointQubitLayer, but instead of treating each gate
   * individually, gates are collected in 2Q-blocks, which are layered
   * disjointly to each other
   *
   * @param lastLayer the array storing the last layer each qubit is used in
   * @param control the (potential) control qubit of the gate
   * @param target the target qubit of the gate
   * @param gate the gate to be added to the layer
   */
  void processDisjoint2qBlockLayer(
      std::array<std::optional<std::size_t>, MAX_DEVICE_QUBITS>& lastLayer,
      const std::optional<std::uint16_t>& control, std::uint16_t target,
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
  Mapper(qc::QuantumComputation quantumComputation, Architecture& architecture);
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
    out << "---------------- Layering -------------------\n";
    for (auto& layer : layers) {
      for (auto& gate : layer) {
        if (gate.singleQubit()) {
          out << "(" << gate.target << ") ";
        } else {
          out << "(" << gate.control << " " << gate.target << ") ";
        }
      }
      out << "\n";
    }
    out << "---------------------------------------------\n";
    return out;
  }

  std::ostream& printLocations(std::ostream& out) {
    out << "---------------- Locations -------------------\n";
    for (std::size_t i = 0; i < qc.getNqubits(); ++i) {
      out << locations.at(i) << " ";
    }
    out << "\n---------------------------------------------\n";
    return out;
  }
  std::ostream& printQubits(std::ostream& out) {
    out << "---------------- Qubits -------------------\n";
    for (std::size_t i = 0; i < architecture->getNqubits(); ++i) {
      out << qubits.at(i) << " ";
    }
    out << "\n---------------------------------------------\n";
    return out;
  }

  virtual void reset() {
    architecture->reset();
    qc.reset();
    layers.clear();
    qubits.fill(DEFAULT_POSITION);
    locations.fill(DEFAULT_POSITION);

    results = MappingResults();
  }
};
