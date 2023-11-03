//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "Mapper.hpp"

#include "CircuitOptimizer.hpp"

void Mapper::initResults() {
  countGates(qc, results.input);
  results.input.name    = qc.getName();
  results.input.qubits  = static_cast<std::uint16_t>(qc.getNqubits());
  results.architecture  = architecture.getName();
  results.output.name   = qc.getName() + "_mapped";
  results.output.qubits = architecture.getNqubits();
  results.output.gates  = std::numeric_limits<std::size_t>::max();
  qcMapped.addQubitRegister(architecture.getNqubits());
}

Mapper::Mapper(const qc::QuantumComputation& quantumComputation,
               Architecture&                 arch)
    : qc(quantumComputation), architecture(arch) {
  qubits.fill(DEFAULT_POSITION);
  locations.fill(DEFAULT_POSITION);
  fidelities.fill(INITIAL_FIDELITY);

  // strip away qubits that are not used in the circuit
  qc.stripIdleQubits(true, true);
  // strip away final measurement gates
  qc::CircuitOptimizer::removeFinalMeasurements(qc);
}

void Mapper::processDisjointQubitLayer(
    std::array<std::optional<std::size_t>, MAX_DEVICE_QUBITS>& lastLayer,
    const std::optional<std::uint16_t>& control, const std::uint16_t target,
    qc::Operation* gate, bool collect2qBlocks) {
  std::size_t layer = 0;
  if (!control.has_value()) {
    if (lastLayer.at(target).has_value()) {
      layer = *lastLayer.at(target);
      // single qubit gates can always be added to the last 2Q block
      if (!collect2qBlocks) {
        layer++;
      }
    }
    if (!collect2qBlocks) {
      lastLayer.at(target) = layer;
    }
    lastLayer.at(target) = layer;
  } else {
    if (!lastLayer.at(*control).has_value() &&
        !lastLayer.at(target).has_value()) {
      layer = 0;
    } else if (!lastLayer.at(*control).has_value()) {
      layer = *lastLayer.at(target) + 1;
    } else if (!lastLayer.at(target).has_value()) {
      layer = *lastLayer.at(*control) + 1;
    } else {
      layer = std::max(*lastLayer.at(*control), *lastLayer.at(target)) + 1;

      if (collect2qBlocks &&
          (*lastLayer.at(*control) == *lastLayer.at(target))) {
        for (auto& g : layers.at(layer - 1)) {
          if ((g.control == *control && g.target == target) ||
              (g.control == target && g.target == *control)) {
            // if last layer contained gate with equivalent qubit set, use that
            // layer
            layer--;
            break;
          }
        }
      }
    }
    lastLayer.at(*control) = layer;
    lastLayer.at(target)   = layer;
  }

  if (layers.size() <= layer) {
    layers.emplace_back();
  }
  if (control.has_value()) {
    layers.at(layer).emplace_back(*control, target, gate);
  } else {
    layers.at(layer).emplace_back(-1, target, gate);
  }
}

void Mapper::createLayers() {
  const auto& config = results.config;
  std::array<std::optional<std::size_t>, MAX_DEVICE_QUBITS> lastLayer{};

  auto qubitsInLayer = std::set<std::uint16_t>{};

  bool even = true;
  for (auto& gate : qc) {
    // skip over barrier instructions
    if (gate->getType() == qc::Barrier || gate->getType() == qc::Measure) {
      continue;
    }

    if (!gate->isUnitary()) {
      throw QMAPException(
          "Mapping not possible: circuit contains non-unitary operation: " +
          std::string(gate->getName()));
    }

    if (gate->getControls().size() > 1 || gate->getTargets().size() > 1) {
      throw QMAPException("Circuit contains gates with more than one control. "
                          "Please make sure that the input circuit's gates are "
                          "decomposed to the appropriate gate set!");
    }

    const bool                   singleQubit = !gate->isControlled();
    std::optional<std::uint16_t> control     = std::nullopt;
    if (!singleQubit) {
      control = static_cast<std::uint16_t>(
          qc.initialLayout.at((*gate->getControls().begin()).qubit));
    }
    const auto target = static_cast<std::uint16_t>(
        qc.initialLayout.at(gate->getTargets().at(0)));

    // methods of layering described in
    // https://iic.jku.at/files/eda/2019_dac_mapping_quantum_circuits_ibm_architectures_using_minimal_number_swap_h_gates.pdf
    switch (config.layering) {
    case Layering::IndividualGates:
    case Layering::None:
      // each gate is put in a new layer
      layers.emplace_back();
      if (control.has_value()) {
        layers.back().emplace_back(*control, target, gate.get());
      } else {
        layers.back().emplace_back(-1, target, gate.get());
      }
      break;
    case Layering::DisjointQubits:
      processDisjointQubitLayer(lastLayer, control, target, gate.get(), false);
      break;
    case Layering::Disjoint2qBlocks:
      processDisjointQubitLayer(lastLayer, control, target, gate.get(), true);
      break;
    case Layering::OddGates:
      // every other gate is put in a new layer
      if (even) {
        layers.emplace_back();
        if (control.has_value()) {
          layers.back().emplace_back(*control, target, gate.get());
        } else {
          layers.back().emplace_back(-1, target, gate.get());
        }
      } else {
        if (control.has_value()) {
          layers.back().emplace_back(*control, target, gate.get());
        } else {
          layers.back().emplace_back(-1, target, gate.get());
        }
      }
      even = !even;
      break;
    case Layering::QubitTriangle:
      if (layers.empty()) {
        layers.emplace_back();
      }

      if (singleQubit) {
        // single qubit gates can be added in any layer
        layers.back().emplace_back(-1, target, gate.get());
      } else {
        qubitsInLayer.insert(*control);
        qubitsInLayer.insert(target);

        if (qubitsInLayer.size() <= 3) {
          layers.back().emplace_back(*control, target, gate.get());
        } else {
          layers.emplace_back();
          layers.back().emplace_back(*control, target, gate.get());
          qubitsInLayer.clear();
          qubitsInLayer.insert(*control);
          qubitsInLayer.insert(target);
        }
      }
      break;
    }
  }
  results.input.layers = layers.size();
}

std::size_t Mapper::getNextLayer(std::size_t idx) {
  auto next = idx + 1;
  while (next < layers.size()) {
    for (const auto& gate : layers.at(next)) {
      if (!gate.singleQubit()) {
        return next;
      }
    }
    next++;
  }
  return std::numeric_limits<std::size_t>::max();
}

void Mapper::finalizeMappedCircuit() {
  // add additional qubits if the architecture contains more qubits than the
  // circuit
  if (architecture.getNqubits() > qcMapped.getNqubits()) {
    for (auto logicalQubit = qcMapped.getNqubits();
         logicalQubit < architecture.getNqubits(); ++logicalQubit) {
      std::optional<qc::Qubit> physicalQubit = std::nullopt;

      // check if the corresponding physical qubit is already in use
      if (qcMapped.initialLayout.find(static_cast<qc::Qubit>(logicalQubit)) !=
          qcMapped.initialLayout.end()) {
        // get the next unused physical qubit
        for (physicalQubit = 0; *physicalQubit < architecture.getNqubits();
             ++(*physicalQubit)) {
          if (qcMapped.initialLayout.find(*physicalQubit) ==
              qcMapped.initialLayout.end()) {
            break;
          }
        }
      } else {
        physicalQubit = static_cast<qc::Qubit>(logicalQubit);
      }

      assert(physicalQubit.has_value());

      // the added logical qubits are not used in the circuit itself, so they
      // are regarded garbage
      qcMapped.addAncillaryQubit(*physicalQubit, std::nullopt);
    }
  }
  // unify quantum registers
  qcMapped.unifyQuantumRegisters();

  // append measurements according to output permutation
  if (results.config.addMeasurementsToMappedCircuit) {
    qcMapped.appendMeasurementsAccordingToOutputPermutation();
  }
}

void Mapper::placeRemainingArchitectureQubits() {
  if (qc.getNqubits() < architecture.getNqubits()) {
    for (auto logical = qc.getNqubits(); logical < architecture.getNqubits();
         ++logical) {
      std::optional<qc::Qubit> physical = std::nullopt;

      // check if the corresponding physical qubit is already in use
      if (qcMapped.initialLayout.find(static_cast<qc::Qubit>(logical)) !=
          qcMapped.initialLayout.end()) {
        // get the next unused physical qubit
        for (physical = 0; *physical < architecture.getNqubits();
             ++(*physical)) {
          if (qcMapped.initialLayout.find(*physical) ==
              qcMapped.initialLayout.end()) {
            break;
          }
        }
      } else {
        physical = static_cast<qc::Qubit>(logical);
      }

      assert(physical.has_value());
      qubits.at(*physical) = static_cast<std::int16_t>(logical);

      // mark architecture qubit as ancillary and garbage
      qcMapped.initialLayout[*physical] = static_cast<qc::Qubit>(logical);
      qcMapped.setLogicalQubitAncillary(static_cast<qc::Qubit>(logical));
      qcMapped.setLogicalQubitGarbage(static_cast<qc::Qubit>(logical));
    }
  }
}

void Mapper::preMappingOptimizations(const Configuration& config
                                     [[maybe_unused]]) {
  if (!config.preMappingOptimizations) {
    return;
  }

  // at the moment there are no pre-mapping optimizations
}

void Mapper::postMappingOptimizations(const Configuration& config) {
  if (!config.postMappingOptimizations) {
    return;
  }

  // try to cancel adjacent CNOT gates
  qc::CircuitOptimizer::cancelCNOTs(qcMapped);
}

void Mapper::countGates(decltype(qcMapped.cbegin())      it,
                        const decltype(qcMapped.cend())& end,
                        MappingResults::CircuitInfo&     info) {
  for (; it != end; ++it) {
    const auto& g = *it;
    if (g->getType() == qc::Teleportation) {
      info.gates += GATES_OF_TELEPORTATION;
      continue;
    }

    if (g->isStandardOperation()) {
      if (g->getType() == qc::SWAP) {
        if (architecture.bidirectional()) {
          info.gates += GATES_OF_BIDIRECTIONAL_SWAP;
          info.cnots += GATES_OF_BIDIRECTIONAL_SWAP;
        } else {
          info.gates += GATES_OF_UNIDIRECTIONAL_SWAP;
          info.cnots += GATES_OF_BIDIRECTIONAL_SWAP;
          info.singleQubitGates += GATES_OF_DIRECTION_REVERSE;
        }
      } else if (g->getControls().empty()) {
        ++info.singleQubitGates;
        ++info.gates;
      } else {
        assert(g->getType() == qc::X);
        ++info.cnots;
        ++info.gates;
      }
      continue;
    }

    if (g->isCompoundOperation()) {
      const auto& cg = dynamic_cast<const qc::CompoundOperation*>(g.get());
      countGates(cg->cbegin(), cg->cend(), info);
    }
  }
}
