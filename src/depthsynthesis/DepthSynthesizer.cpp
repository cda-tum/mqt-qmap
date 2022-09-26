#include "depthsynthesis/DepthSynthesizer.hpp"
void DepthSynthesizer::makeSpecificEncoding(const CliffordSynthesizer::SynthesisData& data) {
    using namespace logicbase;
    makeMultipleGateConstraints(data);
    // COST
    if (strategy == SynthesisStrategy::UseMinimizer ||
        strategy == SynthesisStrategy::SplitIter) {
        LogicTerm cost = LogicTerm(0);
        for (int gateStep = 1; gateStep < data.timesteps + 1; ++gateStep) {
            LogicTerm anyGate = LogicTerm(true);
            for (int a = 0; a < nqubits; ++a) {
                for (auto gate: Gates::SINGLE_QUBIT_WITHOUT_NOP) {
                    anyGate = anyGate && !data.gS[gateStep][Gates::toIndex(gate)][a];
                }
                for (int b = 0; b <= a; ++b) {
                    if (a == b) {
                        continue;
                    }
                    anyGate = anyGate && !data.gC[gateStep][a][b] && !data.gC[gateStep][b][a];
                }
            }
            cost = cost + LogicTerm::ite(anyGate, LogicTerm(5), LogicTerm(0));
        }
        dynamic_cast<LogicBlockOptimizer*>(data.lb.get())->maximize(cost);
    }
}
void DepthSynthesizer::updateResults(SynthesisResults& results) {
    if (results.depth < optimalResults.depth || optimalResults.depth == 0) {
        optimalResults = results;
    }
}
