//
// Created by Sarah on 21.09.2022.
//

#include "Synthesizer.hpp"

#include <utility>
void Synthesizer::initResults() {
    results                  = SynthesisResults();
    results.architectureName = architecture.getName();
    resultCircuit.addQubitRegister(architecture.getNqubits());
}

Synthesizer::Synthesizer(Architecture architecture):
    architecture(std::move(architecture)) {}

void Synthesizer::initCouplingMap(std::uint32_t nqubits) {
    if (architecture.isArchitectureAvailable()){
        auto& cm = highestFidelityCouplingMap.emplace_back();
        architecture.getHighestFidelityCouplingMap(nqubits, cm);
    }
}
