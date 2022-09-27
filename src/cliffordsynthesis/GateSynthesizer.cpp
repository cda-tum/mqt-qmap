#include "cliffordsynthesis/GateSynthesizer.hpp"
void GateSynthesizer::makeSpecificEncoding(const CliffordSynthesizer::SynthesisData& data, const SynthesisConfiguration& configuration) {
    logicbase::LogicTerm changes = logicbase::LogicTerm(true);
    makeSingleGateConstraints(data);
    // COST
    if (configuration.strategy == SynthesisStrategy::UseMinimizer) {
        logicbase::LogicTerm cost = logicbase::LogicTerm(0);
        for (int gateStep = 1; gateStep < data.timesteps + 1; ++gateStep) {
            for (int a = 0; a < configuration.nqubits; ++a) {
                if (configuration.target != SynthesisTarget::GATES_ONLY_CNOT) {
                    for (auto gate: Gates::SINGLE_QUBIT_WITHOUT_NOP) {
                        cost = cost + data.gS[gateStep][Gates::toIndex(gate)][a];
                    }
                }
                for (int b = 0; b <= a; ++b) {
                    if (a == b) {
                        continue;
                    }
                    cost = cost + data.gC[gateStep][a][b] + data.gC[gateStep][b][a];
                }
            }
        }
        dynamic_cast<logicbase::LogicBlockOptimizer*>(data.lb.get())->minimize(cost);
    }
}
void GateSynthesizer::updateResults(SynthesisResults& results) {
    if (results.gateCount < results.gateCount ||
        results.gateCount == 0) {
        results = results;
    }
}
