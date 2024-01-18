//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "Mapper.hpp"

#include "CircuitOptimizer.hpp"

#include <utility>

void Mapper::initResults() {
  countGates(qc, results.input);
  results.input.name    = qc.getName();
  results.input.qubits  = static_cast<std::uint16_t>(qc.getNqubits());
  results.architecture  = architecture->getName();
  results.output.name   = qc.getName() + "_mapped";
  results.output.qubits = architecture->getNqubits();
  results.output.gates  = std::numeric_limits<std::size_t>::max();
  qcMapped.addQubitRegister(architecture->getNqubits());
}

Mapper::Mapper(qc::QuantumComputation quantumComputation, Architecture& arch)
    : qc(std::move(quantumComputation)), architecture(&arch) {
  qubits.fill(DEFAULT_POSITION);
  locations.fill(DEFAULT_POSITION);

  // strip away qubits that are not used in the circuit
  qc.stripIdleQubits(true, true);
  // strip away final measurement gates
  qc::CircuitOptimizer::removeFinalMeasurements(qc);
}

std::size_t Mapper::processDisjointQubitLayer(
    std::array<std::optional<std::size_t>, MAX_DEVICE_QUBITS>& lastLayer,
    const std::optional<std::uint16_t>& control, const std::uint16_t target) {
  std::size_t layer = 0;
  if (!control.has_value()) {
    if (lastLayer.at(target).has_value()) {
      layer = *lastLayer.at(target) + 1;
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
    }
    lastLayer.at(*control) = layer;
    lastLayer.at(target)   = layer;
  }
  return layer;
}

std::size_t Mapper::processDisjoint2qBlockLayer(
    std::array<std::optional<std::size_t>, MAX_DEVICE_QUBITS>& lastLayer,
    const std::optional<std::uint16_t>& control, const std::uint16_t target) {
  std::size_t layer = 0;
  if (!control.has_value()) {
    // single qubit gates can always be added to the last 2Q block and should
    // not affect placings of future 2Q blocks
    if (lastLayer.at(target).has_value()) {
      layer = *lastLayer.at(target);
    }
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

      if (*lastLayer.at(*control) == *lastLayer.at(target)) {
        for (auto& g : layers.at(layer - 1)) {
          if ((g.control == *control && g.target == target) ||
              (g.control == target && g.target == *control)) {
            // if last layer contained gate with equivalent qubit set, use that
            // layer
            --layer;
            break;
          }
        }
      }
    }
    lastLayer.at(*control) = layer;
    lastLayer.at(target)   = layer;
  }
  return layer;
}

std::size_t Mapper::processDisjointSameOpTypeBlockLayer(
      std::array<std::optional<std::size_t>, MAX_DEVICE_QUBITS>& lastLayer,
      std::array<qc::OpType, MAX_DEVICE_QUBITS>& lastLayerOpType,
      const std::optional<std::uint16_t>& control, std::uint16_t target,
      qc::Operation* gate) {
  std::size_t layer = 0;
  if (!control.has_value()) {
    if (lastLayer.at(target).has_value()) {
      layer = *lastLayer.at(target);
      if (lastLayerOpType.at(target) != gate->getType()) {
        ++layer;
      }
      lastLayer.at(target) = layer;
    }
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

      if (*lastLayer.at(*control) == *lastLayer.at(target) && 
          lastLayerOpType.at(*control) == gate->getType() && 
          lastLayerOpType.at(target) == gate->getType()) {
        for (auto& g : layers.at(layer - 1)) {
          if ((g.control == *control && g.target == target) ||
              (g.control == target && g.target == *control)) {
            // if last layer contained gate with equivalent qubit set, use that
            // layer
            --layer;
            break;
          }
        }
      }
    }
    lastLayer.at(*control) = layer;
    lastLayer.at(target)   = layer;
  }
  return layer;
}

void Mapper::createLayers() {
  const auto& config = results.config;
  std::array<std::optional<std::size_t>, MAX_DEVICE_QUBITS> lastLayer{};
  std::array<qc::OpType, MAX_DEVICE_QUBITS> lastLayerOpType{};
  lastLayerOpType.fill(qc::OpType::None);
  
  singleQubitMultiplicities.clear();
  twoQubitMultiplicities.clear();
  activeQubits.clear();
  activeQubits1QGates.clear();
  activeQubits2QGates.clear();
  layers.clear();
  
  maximumActiveQubits = config.maximumActiveQubits;
  maximumActiveQubits1QGates = config.maximumActiveQubits1QGates;
  maximumActiveQubits2QGates = config.maximumActiveQubits2QGates;
  // TODO: find maximum matching in architecture to prevent overfilling
  // maximumActiveQubits2QGates = std::max(config.maximumActiveQubits2QGates, 2*architecture->getMaxMatchingSize());
  
  if (maximumActiveQubits == 1) {
    throw QMAPException("maximumActiveQubits cannot be 1, since this prevents placing any 2q gates!");
  }
  if (maximumActiveQubits2QGates == 1) {
    throw QMAPException("maximumActiveQubits2QGates cannot be 1, since this prevents placing any 2q gates!");
  }

  auto qubitsInLayer = std::set<std::uint16_t>{};
  std::size_t layer = 0;

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
      // each gate is put in a new layer
      layer = layers.size();
      break;
    case Layering::DisjointQubits:
      layer = processDisjointQubitLayer(lastLayer, control, target);
      break;
    case Layering::Disjoint2qBlocks:
      layer = processDisjoint2qBlockLayer(lastLayer, control, target);
      break;
    case Layering::DisjointSameOpTypeBlocks:
      layer = processDisjointSameOpTypeBlockLayer(lastLayer, lastLayerOpType, control, target, gate.get());
      break;
    case Layering::OddGates:
      // every other gate is put in a new layer
      if (even) {
        layer = layers.size();
      } else {
        assert(!layers.empty());
        layer = layers.size()-1;
      }
      even = !even;
      break;
    case Layering::QubitTriangle:
      if (singleQubit) {
        // single qubit gates can be added in any layer
        layer = layers.empty() ? 0 : layers.size()-1;
      } else {
        qubitsInLayer.insert(*control);
        qubitsInLayer.insert(target);
        if (qubitsInLayer.size() <= 3) {
          layer = layers.size()-1;
        } else {
          layer = layers.size();
          qubitsInLayer.clear();
          qubitsInLayer.insert(*control);
          qubitsInLayer.insert(target);
        }
      }
      break;
    }
    
    // check if layer is full, if so step to the next layer
    bool layerChanged = false;
    while(layer < layers.size()) {
      if (maximumActiveQubits > 0 &&
        maximumActiveQubits < activeQubits.at(layer).size() + 
        static_cast<std::size_t>(
          activeQubits.at(layer).find(target) == 
          activeQubits.at(layer).end()) + 
        (singleQubit ? 0 : 
          static_cast<std::size_t>(
            activeQubits.at(layer).find(*control) == 
            activeQubits.at(layer).end()
          )
        )
      ) {
        // if activeQubits.size() would grow above maximumActiveQubits
        ++layer;
        layerChanged = true;
        continue;
      }
      if (control.has_value()) {
        // 2q gate
        if (maximumActiveQubits2QGates > 0 &&
          maximumActiveQubits2QGates < activeQubits2QGates.at(layer).size() + 
          static_cast<std::size_t>(
            activeQubits2QGates.at(layer).find(target) == 
            activeQubits2QGates.at(layer).end()) + 
          static_cast<std::size_t>(
            activeQubits2QGates.at(layer).find(*control) == 
            activeQubits2QGates.at(layer).end()
          )
        ) {
          // if activeQubits2QGates.size() would grow above 
          // maximumActiveQubits2QGates
          ++layer;
          layerChanged = true;
          continue;
        }
      } else {
        // 1q gate
        if (maximumActiveQubits1QGates > 0 &&
          maximumActiveQubits1QGates < activeQubits1QGates.at(layer).size() + 
          static_cast<std::size_t>(
            activeQubits1QGates.at(layer).find(target) == 
            activeQubits1QGates.at(layer).end())
        ) {
          // if activeQubits1QGates.size() would grow above 
          // maximumActiveQubits1QGates
          ++layer;
          layerChanged = true;
          continue;
        }
      }
      break;
    }
    
    if (layerChanged) {
      even = true;
      lastLayer.at(target)   = layer;
      lastLayerOpType.at(target) = gate->getType();
      if (!singleQubit) {
        lastLayer.at(*control) = layer;
        lastLayerOpType.at(*control) = gate->getType();
        qubitsInLayer.clear();
        qubitsInLayer.insert(*control);
        qubitsInLayer.insert(target);
      }
    }
    
    while(layers.size() <= layer) {
      singleQubitMultiplicities.emplace_back(architecture->getNqubits(), 0);
      twoQubitMultiplicities.emplace_back();
      activeQubits.emplace_back();
      activeQubits1QGates.emplace_back();
      activeQubits2QGates.emplace_back();
      layers.emplace_back();
    }
    
    if (singleQubit) {
      layers.at(layer).emplace_back(-1, target, gate.get());
      activeQubits.at(layer).emplace(target);
      activeQubits1QGates.at(layer).emplace(target);
      ++singleQubitMultiplicities.at(layer).at(target);
    } else {
      layers.at(layer).emplace_back(*control, target, gate.get());
      activeQubits.at(layer).emplace(*control);
      activeQubits.at(layer).emplace(target);
      activeQubits2QGates.at(layer).emplace(*control);
      activeQubits2QGates.at(layer).emplace(target);
      if (*control >= target) {
        auto insertResult = twoQubitMultiplicities.at(layer).insert({{target, *control}, {0, 1}});
        if (insertResult.second) {
          // if there is already an entry for this edge, the multiplicity in 
          // the backward direction is increased
          insertResult.first->second.second++;
        }
      } else {
        auto insertResult = twoQubitMultiplicities.at(layer).insert({{*control, target}, {1, 0}});
        if (insertResult.second) {
          // if there is already an entry for this edge, the multiplicity in
          // the forward direction is increased
          insertResult.first->second.first++;
        }
      }
    }
  }
  results.input.layers = layers.size();
}

bool Mapper::isLayerSplittable(std::size_t index) {
  if (twoQubitMultiplicities.at(index).size() > 1) {
    return true;
  }
  if (activeQubits1QGates.at(index).size() > 2) {
    return true;
  }
  if (twoQubitMultiplicities.at(index).empty()) {
    return false;
  }
  // check if there is a 1Q gate on a qubit that is not part of the 2Q gate
  return std::any_of(activeQubits1QGates.at(index).begin(),
                     activeQubits1QGates.at(index).end(),
                     [this, index](auto q) {
                       return activeQubits2QGates.at(index).find(q) ==
                              activeQubits2QGates.at(index).end();
                     });
}

void Mapper::splitLayer(std::size_t index, Architecture& arch) {
  const SingleQubitMultiplicity& singleQubitMultiplicity =
      singleQubitMultiplicities.at(index);
  const TwoQubitMultiplicity& twoQubitMultiplicity =
      twoQubitMultiplicities.at(index);
  std::vector<Gate>       layer0{};
  std::vector<Gate>       layer1{};
  SingleQubitMultiplicity singleQubitMultiplicity0(arch.getNqubits(), 0);
  SingleQubitMultiplicity singleQubitMultiplicity1(arch.getNqubits(), 0);
  TwoQubitMultiplicity    twoQubitMultiplicity0{};
  TwoQubitMultiplicity    twoQubitMultiplicity1{};
  std::set<std::uint16_t> activeQubits0{};
  std::set<std::uint16_t> activeQubits1QGates0{};
  std::set<std::uint16_t> activeQubits2QGates0{};
  std::set<std::uint16_t> activeQubits1{};
  std::set<std::uint16_t> activeQubits1QGates1{};
  std::set<std::uint16_t> activeQubits2QGates1{};

  // 2Q-gates
  bool even = false;
  for (const auto& edge : twoQubitMultiplicity) {
    if (even) {
      twoQubitMultiplicity0.insert(edge);
      activeQubits0.emplace(edge.first.first);
      activeQubits0.emplace(edge.first.second);
      activeQubits2QGates0.emplace(edge.first.first);
      activeQubits2QGates0.emplace(edge.first.second);
    } else {
      twoQubitMultiplicity1.insert(edge);
      activeQubits1.emplace(edge.first.first);
      activeQubits1.emplace(edge.first.second);
      activeQubits2QGates1.emplace(edge.first.first);
      activeQubits2QGates1.emplace(edge.first.second);
    }
    even = !even;
  }

  // 1Q-gates
  even = true;
  for (std::size_t q = 0; q < singleQubitMultiplicity.size(); ++q) {
    if (singleQubitMultiplicity[q] == 0) {
      continue;
    }
    // if a qubit is also acted on by a 2Q-gate, put it on the same layer as
    // the 2Q-gate
    if (activeQubits2QGates0.find(static_cast<std::uint16_t>(q)) !=
        activeQubits2QGates0.end()) {
      singleQubitMultiplicity0[q] = singleQubitMultiplicity[q];
      activeQubits0.emplace(static_cast<std::uint16_t>(q));
      activeQubits1QGates0.emplace(static_cast<std::uint16_t>(q));
      continue;
    }
    if (activeQubits2QGates1.find(static_cast<std::uint16_t>(q)) !=
        activeQubits2QGates1.end()) {
      singleQubitMultiplicity1[q] = singleQubitMultiplicity[q];
      activeQubits1.emplace(static_cast<std::uint16_t>(q));
      activeQubits1QGates1.emplace(static_cast<std::uint16_t>(q));
      continue;
    }

    if (even) {
      singleQubitMultiplicity0[q] = singleQubitMultiplicity[q];
      activeQubits0.emplace(q);
      activeQubits1QGates0.emplace(q);
    } else {
      singleQubitMultiplicity1[q] = singleQubitMultiplicity[q];
      activeQubits1.emplace(q);
      activeQubits1QGates1.emplace(q);
    }
    even = !even;
  }

  for (auto& gate : layers[index]) {
    if (gate.singleQubit()) {
      if (singleQubitMultiplicity0[gate.target] > 0) {
        layer0.emplace_back(gate);
      } else {
        layer1.emplace_back(gate);
      }
    } else {
      if (activeQubits2QGates0.find(gate.target) !=
          activeQubits2QGates0.end()) {
        layer0.emplace_back(gate);
      } else {
        layer1.emplace_back(gate);
      }
    }
  }

  // insert new layers
  layers[index] = layer0;
  layers.insert(
      layers.begin() +
          static_cast<std::vector<std::vector<Mapper::Gate>>::difference_type>(
              index) +
          1,
      layer1);
  singleQubitMultiplicities[index] = singleQubitMultiplicity0;
  singleQubitMultiplicities.insert(
      singleQubitMultiplicities.begin() +
          static_cast<std::vector<SingleQubitMultiplicity>::difference_type>(
              index) +
          1,
      singleQubitMultiplicity1);
  twoQubitMultiplicities[index] = twoQubitMultiplicity0;
  twoQubitMultiplicities.insert(
      twoQubitMultiplicities.begin() +
          static_cast<std::vector<TwoQubitMultiplicity>::difference_type>(
              index) +
          1,
      twoQubitMultiplicity1);
  activeQubits[index] = activeQubits0;
  activeQubits.insert(
      activeQubits.begin() +
          static_cast<std::vector<std::set<std::uint16_t>>::difference_type>(
              index) +
          1,
      activeQubits1);
  activeQubits1QGates[index] = activeQubits1QGates0;
  activeQubits1QGates.insert(
      activeQubits1QGates.begin() +
          static_cast<std::vector<std::set<std::uint16_t>>::difference_type>(
              index) +
          1,
      activeQubits1QGates1);
  activeQubits2QGates[index] = activeQubits2QGates0;
  activeQubits2QGates.insert(
      activeQubits2QGates.begin() +
          static_cast<std::vector<std::set<std::uint16_t>>::difference_type>(
              index) +
          1,
      activeQubits2QGates1);
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
  if (architecture->getNqubits() > qcMapped.getNqubits()) {
    for (auto logicalQubit = qcMapped.getNqubits();
         logicalQubit < architecture->getNqubits(); ++logicalQubit) {
      std::optional<qc::Qubit> physicalQubit = std::nullopt;

      // check if the corresponding physical qubit is already in use
      if (qcMapped.initialLayout.find(static_cast<qc::Qubit>(logicalQubit)) !=
          qcMapped.initialLayout.end()) {
        // get the next unused physical qubit
        for (physicalQubit = 0; *physicalQubit < architecture->getNqubits();
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
  if (qc.getNqubits() < architecture->getNqubits()) {
    for (auto logical = qc.getNqubits(); logical < architecture->getNqubits();
         ++logical) {
      std::optional<qc::Qubit> physical = std::nullopt;

      // check if the corresponding physical qubit is already in use
      if (qcMapped.initialLayout.find(static_cast<qc::Qubit>(logical)) !=
          qcMapped.initialLayout.end()) {
        // get the next unused physical qubit
        for (physical = 0; *physical < architecture->getNqubits();
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
        auto q1 = static_cast<std::uint16_t>(g->getTargets()[0]);
        auto q2 = static_cast<std::uint16_t>(g->getTargets()[1]);
        if (architecture->isFidelityAvailable()) {
          info.totalLogFidelity += architecture->getSwapFidelityCost(q1, q2);
        }
        if (architecture->getCouplingMap().find({q1, q2}) !=
                architecture->getCouplingMap().end() &&
            architecture->getCouplingMap().find({q2, q1}) !=
                architecture->getCouplingMap().end()) {
          // bidirectional edge
          info.gates += GATES_OF_BIDIRECTIONAL_SWAP;
          info.cnots += GATES_OF_BIDIRECTIONAL_SWAP;
        } else {
          // unidirectional edge
          info.gates += GATES_OF_UNIDIRECTIONAL_SWAP;
          info.cnots += GATES_OF_BIDIRECTIONAL_SWAP;
          info.singleQubitGates += GATES_OF_DIRECTION_REVERSE;
        }
      } else if (g->getControls().empty()) {
        ++info.singleQubitGates;
        ++info.gates;
        auto q1 = static_cast<std::uint16_t>(g->getTargets()[0]);
        if (architecture->isFidelityAvailable()) {
          info.totalLogFidelity += architecture->getSingleQubitFidelityCost(q1);
        }
      } else {
        assert(g->getType() == qc::X);
        ++info.cnots;
        ++info.gates;
        auto q1 =
            static_cast<std::uint16_t>((*(g->getControls().begin())).qubit);
        auto q2 = static_cast<std::uint16_t>(g->getTargets()[0]);
        if (architecture->isFidelityAvailable()) {
          info.totalLogFidelity +=
              architecture->getTwoQubitFidelityCost(q1, q2);
        }
      }
      continue;
    }

    if (g->isCompoundOperation()) {
      const auto& cg = dynamic_cast<const qc::CompoundOperation*>(g.get());
      countGates(cg->cbegin(), cg->cend(), info);
    }
  }
  info.totalFidelity = std::pow(2, -info.totalLogFidelity);
}
