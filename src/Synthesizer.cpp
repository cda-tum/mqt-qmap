//
// Created by Sarah on 21.09.2022.
//

#include "Synthesizer.hpp"
void Synthesizer::initResults() {
    results                  = SynthesisResults();
    results.architectureName = architecture.getName();
    resultCircuit.addQubitRegister(architecture.getNqubits());
}

Synthesizer::Synthesizer(Architecture architecture):
    architecture(architecture) {}
