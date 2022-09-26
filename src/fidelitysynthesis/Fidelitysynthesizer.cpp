#include "fidelitysynthesis/Fidelitysynthesizer.hpp"
void Fidelitysynthesizer::makeSynthesis(const CliffordSynthesizer::SynthesisData& data) {
    if (!architecture.isArchitectureAvailable()) {
        util::fatal("No fidelity architecture specified in coupling map.");
    }
    makeMultipleGateConstraints(data);
    // COST
    if (strategy == SynthesisStrategy::UseMinimizer ||
        strategy == SynthesisStrategy::SplitIter) {
        logicbase::LogicTerm cost = logicbase::LogicTerm(0);
        // For each edge in the coupling map, get the fidelity cost
        for (const auto& edge: data.reducedCM) {
            logicbase::LogicTerm fidelity = logicbase::LogicTerm(
                    (1 - std::log(architecture.getFidelityTable()[edge.first][edge.second])) * 1000);
            auto a = std::find(data.qubitChoice.begin(), data.qubitChoice.end(), edge.first);
            auto b = std::find(data.qubitChoice.begin(), data.qubitChoice.end(), edge.second);
            if (a == data.qubitChoice.end() || b == data.qubitChoice.end()) {
                util::fatal("Coupling map contains invalid qubit.");
            }
            // at each time t if there is a gate on the edge, add the cost
            for (unsigned int gateStep = 0; gateStep < data.timesteps; ++gateStep) {
                cost = cost + (data.gC[gateStep][std::distance(data.qubitChoice.begin(), a)][std::distance(data.qubitChoice.begin(), b)] * fidelity);
            }
        }
        // For each qubit, get the fidelity cost
        for (unsigned int a = 0; a < data.nqubits; ++a) {
            logicbase::LogicTerm fidelity =
                    logicbase::LogicTerm((1 - std::log(architecture.getSingleQubitFidelities()[a])) * 1000);
            // at each time t if there is a gate on a, add the cost
            for (int gateStep = 0; gateStep < data.timesteps; ++gateStep) {
                for (auto gate: Gates::SINGLE_QUBIT_WITHOUT_NOP) {
                    cost = cost + (data.gS[gateStep][Gates::toIndex(gate)][a] * fidelity);
                }
            }
        }
        dynamic_cast<logicbase::LogicBlockOptimizer*>(data.lb.get())->minimize(cost);
        cost = logicbase::LogicTerm(0);
        for (unsigned int gateStep = 0; gateStep < data.timesteps; ++gateStep) {
            for (unsigned int a = 0; a < data.nqubits; ++a) {
                cost = cost + data.gS[gateStep][1][a] + data.gS[gateStep][2][a];
                for (unsigned int b = 0; b < data.nqubits; ++b) {
                    cost = cost + data.gC[gateStep][a][b];
                }
            }
        }
        dynamic_cast<logicbase::LogicBlockOptimizer*>(data.lb.get())->maximize(cost);
    }
}
