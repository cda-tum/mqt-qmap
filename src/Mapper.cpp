/*
 * This file is part of the MQT QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#include "Mapper.hpp"

#include "CircuitOptimizer.hpp"

void Mapper::initResults() {
    countGates(qc, results.input);
    results.input.name    = qc.getName();
    results.input.qubits  = qc.getNqubits();
    results.architecture  = architecture.getArchitectureName();
    results.output.name   = qc.getName() + "_mapped";
    results.output.qubits = architecture.getNqubits();
    results.output.gates  = std::numeric_limits<unsigned long>::max();
    qcMapped.addQubitRegister(architecture.getNqubits());
}

Mapper::Mapper(const qc::QuantumComputation& quantumComputation, Architecture& arch):
    qc(quantumComputation.clone()), architecture(arch) {
    qubits.fill(DEFAULT_POSITION);
    locations.fill(DEFAULT_POSITION);
    fidelities.fill(INITIAL_FIDELITY);

    // strip away qubits that are not used in the circuit
    qc.stripIdleQubits(true, true);
    // strip away final measurement gates
    qc::CircuitOptimizer::removeFinalMeasurements(qc);
}

void Mapper::createLayers() {
    const auto&                          config = results.config;
    std::array<short, MAX_DEVICE_QUBITS> lastLayer{};
    lastLayer.fill(DEFAULT_POSITION);

    auto qubitsInLayer = std::set<unsigned short>{};

    bool even = true;
    for (auto& gate: qc) {
        // skip over barrier instructions
        if (gate->getType() == qc::Barrier || gate->getType() == qc::Measure) {
            continue;
        }

        if (!gate->isUnitary()) {
            throw QMAPException("Mapping not possible: circuit contains non-unitary operation: " + std::string(gate->getName()));
        }

        if (gate->getControls().size() > 1 || gate->getTargets().size() > 1) {
            throw QMAPException("Circuit contains gates with more than one control. Please make sure that the input circuit's gates are decomposed to the appropriate gate set!");
        }

        bool  singleQubit = gate->getControls().empty();
        short control     = -1;
        if (!singleQubit) {
            control = static_cast<short>(qc.initialLayout.at((*gate->getControls().begin()).qubit));
        }
        unsigned short target = qc.initialLayout.at(gate->getTargets().at(0));
        size_t         layer  = 0;

        switch (config.layering) {
            case Layering::IndividualGates:
            case Layering::None:
                layers.emplace_back();
                layers.back().emplace_back(control, target, gate.get());
                break;
            case Layering::DisjointQubits:
                if (singleQubit) {
                    layer                = lastLayer.at(target) + 1;
                    lastLayer.at(target) = layer;
                } else {
                    layer                 = std::max(lastLayer.at(control), lastLayer.at(target)) + 1;
                    lastLayer.at(control) = lastLayer.at(target) = layer;
                }

                if (layers.size() <= layer) {
                    layers.emplace_back();
                }
                layers.at(layer).emplace_back(control, target, gate.get());
                break;
            case Layering::OddGates:
                if (even) {
                    layers.emplace_back();
                    layers.back().emplace_back(control, target, gate.get());
                } else {
                    layers.back().emplace_back(control, target, gate.get());
                }
                even = !even;
                break;
            case Layering::QubitTriangle:
                if (layers.empty()) {
                    layers.emplace_back();
                }

                if (singleQubit) {
                    // single qubit gates can be added in any layer
                    layers.back().emplace_back(control, target, gate.get());
                } else {
                    qubitsInLayer.insert(control);
                    qubitsInLayer.insert(target);

                    if (qubitsInLayer.size() <= 3) {
                        layers.back().emplace_back(control, target, gate.get());
                    } else {
                        layers.emplace_back();
                        layers.back().emplace_back(control, target, gate.get());
                        qubitsInLayer.clear();
                        qubitsInLayer.insert(control);
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
        for (const auto& gate: layers.at(next)) {
            if (!gate.singleQubit()) {
                return next;
            }
        }
        next++;
    }
    return -1;
}

void Mapper::finalizeMappedCircuit() {
    // add additional qubits if the architecture contains more qubits than the circuit
    if (architecture.getNqubits() > qcMapped.getNqubits()) {
        for (auto logicalQubit = qcMapped.getNqubits(); logicalQubit < static_cast<dd::QubitCount>(architecture.getNqubits()); ++logicalQubit) {
            dd::Qubit physicalQubit = -1;

            // check if the corresponding physical qubit is already in use
            if (qcMapped.initialLayout.find(static_cast<dd::Qubit>(logicalQubit)) != qcMapped.initialLayout.end()) {
                // get the next unused physical qubit
                for (physicalQubit = 0; physicalQubit < static_cast<dd::Qubit>(architecture.getNqubits()); ++physicalQubit) {
                    if (qcMapped.initialLayout.find(physicalQubit) == qcMapped.initialLayout.end()) {
                        break;
                    }
                }
            } else {
                physicalQubit = static_cast<dd::Qubit>(logicalQubit);
            }

            assert(physicalQubit != -1);

            // the added logical qubits are not used in the circuit itself, so they are regarded garbage
            qcMapped.addAncillaryQubit(physicalQubit, -1);
        }
    }
    // unify quantum registers
    qcMapped.unifyQuantumRegisters();

    // append measurements according to output permutation
    qcMapped.appendMeasurementsAccordingToOutputPermutation();
}

void Mapper::placeRemainingArchitectureQubits() {
    if (qc.getNqubits() < architecture.getNqubits()) {
        for (auto logical = qc.getNqubits(); logical < static_cast<decltype(logical)>(architecture.getNqubits()); ++logical) {
            dd::Qubit physical = -1;

            // check if the corresponding physical qubit is already in use
            if (qcMapped.initialLayout.find(static_cast<dd::Qubit>(logical)) != qcMapped.initialLayout.end()) {
                // get the next unused physical qubit
                for (physical = 0; physical < static_cast<dd::Qubit>(architecture.getNqubits()); ++physical) {
                    if (qcMapped.initialLayout.find(physical) == qcMapped.initialLayout.end()) {
                        break;
                    }
                }
            } else {
                physical = static_cast<dd::Qubit>(logical);
            }

            qubits.at(physical) = logical;

            // mark architecture qubit as ancillary and garbage
            qcMapped.initialLayout[physical] = static_cast<dd::Qubit>(logical);
            qcMapped.setLogicalQubitAncillary(logical);
            qcMapped.setLogicalQubitGarbage(logical);
        }
    }
}

void Mapper::preMappingOptimizations(const Configuration& config [[maybe_unused]]) {
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

void Mapper::countGates(decltype(qcMapped.cbegin()) it, const decltype(qcMapped.cend())& end, MappingResults::CircuitInfo& info) {
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
