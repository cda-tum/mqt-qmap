/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
*/

#include "cliffordsynthesis/TargetMetricHandler.hpp"

#include "cliffordsynthesis/GateEncoding.hpp"
namespace cs {

    void TargetMetricHandler::makeTargetMetric(const SynthesisData& data, const Configuration& configuration) {
        bool useMaxSat = configuration.strategy == OptimizationStrategy::UseMinimizer || configuration.strategy == OptimizationStrategy::SplitIter;
        bool onlyCnot  = configuration.target == TargetMetric::TWO_QUBIT_GATES;

        switch (configuration.target) {
            case TargetMetric::GATES:
            case TargetMetric::TWO_QUBIT_GATES:
                makeGateMetric(data, useMaxSat, onlyCnot);
                break;
            case TargetMetric::DEPTH:
                makeDepthMetric(data, useMaxSat);
                break;
            case TargetMetric::FIDELITY:
                makeFidelityMetric(data, useMaxSat, configuration.architecture);
                break;
        }
    }
    void TargetMetricHandler::makeGateMetric(const SynthesisData& data, bool useMaxSat, bool onlyCNOT) {
        logicbase::LogicTerm changes = logicbase::LogicTerm(true);
        // COST
        if (useMaxSat) {
            logicbase::LogicTerm cost = logicbase::LogicTerm(0);
            for (unsigned int gateStep = 1; gateStep < data.timesteps + 1; ++gateStep) {
                for (unsigned int a = 0; a < data.nqubits; ++a) {
                    if (!onlyCNOT) {
                        for (auto gate: Gates::SINGLE_QUBIT_WITHOUT_NOP) {
                            cost = cost + data.gS[gateStep][Gates::toIndex(gate)][a];
                        }
                    }
                    for (unsigned int b = 0; b <= a; ++b) {
                        if (a == b) {
                            continue;
                        }
                        cost = cost + data.gTwoQubit[gateStep][a][b] + data.gTwoQubit[gateStep][b][a];
                    }
                }
            }
            dynamic_cast<logicbase::LogicBlockOptimizer*>(data.lb.get())->minimize(cost);
        }
    }

    void TargetMetricHandler::makeDepthMetric(const SynthesisData& data, bool useMaxSat) {
        using namespace logicbase;
        // COST
        if (useMaxSat) {
            LogicTerm cost = LogicTerm(0);
            for (unsigned int gateStep = 1; gateStep < data.timesteps + 1; ++gateStep) {
                LogicTerm anyGate = LogicTerm(true);
                for (unsigned int a = 0; a < data.nqubits; ++a) {
                    for (auto gate: Gates::SINGLE_QUBIT_WITHOUT_NOP) {
                        anyGate = anyGate && !data.gS[gateStep][Gates::toIndex(gate)][a];
                    }
                    for (unsigned int b = 0; b <= a; ++b) {
                        if (a == b) {
                            continue;
                        }
                        anyGate = anyGate && !data.gTwoQubit[gateStep][a][b] && !data.gTwoQubit[gateStep][b][a];
                    }
                }
                cost = cost + LogicTerm::ite(anyGate, LogicTerm(5), LogicTerm(0));
            }
            dynamic_cast<LogicBlockOptimizer*>(data.lb.get())->maximize(cost);
        }
    }

    void TargetMetricHandler::makeFidelityMetric(const SynthesisData& data, bool useMaxSat, const Architecture& architecture) {
        if (!architecture.isArchitectureAvailable() || !architecture.isCalibrationDataAvailable()) {
            util::fatal("No fidelity architecture specified in coupling map.");
        }
        // COST
        if (useMaxSat) {
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
                    cost = cost + (data.gTwoQubit[gateStep][std::distance(data.qubitChoice.begin(), a)][std::distance(data.qubitChoice.begin(), b)] * fidelity);
                }
            }
            // For each qubit, get the fidelity cost
            for (unsigned int a = 0; a < data.nqubits; ++a) {
                logicbase::LogicTerm fidelity =
                        logicbase::LogicTerm((1 - std::log(architecture.getSingleQubitFidelities()[a])) * 1000);
                // at each time t if there is a gate on a, add the cost
                for (unsigned int gateStep = 0; gateStep < data.timesteps; ++gateStep) {
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
                        cost = cost + data.gTwoQubit[gateStep][a][b];
                    }
                }
            }
            dynamic_cast<logicbase::LogicBlockOptimizer*>(data.lb.get())->maximize(cost);
        }
    }
    void TargetMetricHandler::updateResults(const Configuration& configuration, Results& results, Results& currentResults) {
        switch (configuration.target) {
            case TargetMetric::GATES:
            case TargetMetric::TWO_QUBIT_GATES:
                if ((results.sat && results.gateCount < currentResults.gateCount) || currentResults.gateCount == 0) {
                    currentResults = results;
                }
                break;
            case TargetMetric::DEPTH:
                if ((results.sat && results.depth < currentResults.depth) || currentResults.depth == 0) {
                    currentResults = results;
                }
                break;
            case TargetMetric::FIDELITY:
                if ((results.sat && results.fidelity < currentResults.fidelity) || currentResults.fidelity == 0) {
                    currentResults = results;
                }
                break;
        }
    }

} // namespace cs
