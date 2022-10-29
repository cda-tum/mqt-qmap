/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
*/

#include "cliffordsynthesis/HeuristicStrategy.hpp"

#include <cstddef>
namespace cs {
    void cs::HeuristicStrategy::runHeuristicStrategy(const CouplingMap& reducedCM, const QubitSubset& qubitChoice, const Configuration& configuration, CliffordSynthesizer& synthesizer) {
        if (configuration.strategy == OptimizationStrategy::SplitIter) {
            runSplitIter(reducedCM, qubitChoice, configuration, synthesizer);
        } else {
            throw std::runtime_error("Unknown optimization strategy");
        }
    }

    void HeuristicStrategy::runSplitIter(const CouplingMap& reducedCM, const QubitSubset& qubitChoice, const Configuration& configuration, CliffordSynthesizer& synthesizer) {
        if (configuration.targetCircuit->size() < 2) {
            return;
        }
        DEBUG() << "Running split iter" << std::endl;
        Tableau                   fullTableau  = *configuration.targetTableau;
        auto                      circuitSplit = static_cast<std::size_t>(std::log(configuration.targetCircuit->getNindividualOps()));
        std::uint16_t             split        = std::min(5, configuration.nqubits / 2);
        std::vector<std::shared_ptr<std::thread>> threads;
        std::vector<std::shared_ptr<Results>>     results;
        std::uint32_t             nThreads = configuration.nThreads;
        qc::QuantumComputation    circuit  = configuration.targetCircuit->clone();
        bool                      stopping = false;
        while (!stopping) {
            results.clear();
            DEBUG() << "Current split size: " << split << std::endl;
            DEBUG() << "Current circuit split size: " << circuitSplit << std::endl;
            auto    start = std::chrono::high_resolution_clock::now();
            Results totalResult;
            totalResult.result = logicbase::Result::SAT;
            qc::QuantumComputation localResultCircuit;
            localResultCircuit.addQubitRegister(configuration.nqubits);
            for (size_t i = 0; i * circuitSplit < circuit.getNindividualOps();
                 i += nThreads) {
                threads.clear();
                DEBUG() << "Currently at " << i * circuitSplit << " of "
                        << circuit.getNindividualOps() << std::endl;
                for (std::uint32_t j = 0; j < nThreads; j++) {
                    std::shared_ptr<Results> r = std::make_shared<Results>();
                    std::shared_ptr<std::thread> t = std::make_unique<std::thread>(runSplinter, i, circuitSplit,
                                              split, std::ref(reducedCM), std::ref(qubitChoice),
                                              std::ref(circuit), r, &synthesizer, configuration);
                    threads.emplace_back(t);
                    results.push_back(r);
                }
                for (const auto& t: threads) {
                    t->join();
                }
                for (const auto& r: results) {
                    if (r->result == logicbase::Result::UNSAT) {
                        totalResult.result = logicbase::Result::UNSAT;
                        break;
                    }
                }
                if (totalResult.result == logicbase::Result::UNSAT) {
                    DEBUG() << "UNSAT, increasing split size." << std::endl;
                    split += std::max(1.0, split * configuration.circuitSplittingIncrease);
                    circuitSplit += std::max(1.0, circuitSplit * configuration.circuitSplittingIncrease);
                    if (circuitSplit > configuration.targetCircuit->getNindividualOps()) {
                        stopping = true;
                        break;
                    }
                    break;
                }
            }
            for (const auto& r: results) {
                if (r->resultStringCircuit.empty()) {
                    continue;
                }
                qc::QuantumComputation splitResult;
                std::istringstream     iss(r->resultStringCircuit);
                splitResult.import(iss, qc::OpenQASM);
                for (const auto& gate: splitResult) {
                    localResultCircuit.insert(localResultCircuit.end(),
                                              gate->clone());
                }
            }
            if (totalResult.result == logicbase::Result::SAT) {
                Tableau resultingTableau{localResultCircuit};
                DEBUG() << "Equality (Results): "
                        << ((fullTableau == resultingTableau) ? "True" : "False")
                        << std::endl;
                DEBUG() << "Original Circuit size: " << circuit.getNindividualOps()
                        << std::endl;
                DEBUG() << "Optimized Circuit size: "
                        << localResultCircuit.getNindividualOps() << std::endl;
                TRACE() << "Resulting Circuit: " << std::endl;
                std::ostringstream ss;
                localResultCircuit.dump(ss, qc::Format::OpenQASM);
                TRACE() << ss.str() << std::endl;
                auto                          end  = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> diff = end - start;
                INFO() << "Time for complete run: " << diff.count() << std::endl;
                if (circuit.getNindividualOps() ==
                    localResultCircuit.getNindividualOps()) {
                    break;
                }
                circuit = localResultCircuit.clone();
            }
        }
        std::stringstream ss;
        circuit.dump(ss, qc::Format::OpenQASM);
        synthesizer.optimalResults.resultStringCircuit = ss.str();
        synthesizer.optimalResults.resultTableaus.emplace_back(*configuration.targetTableau);
        synthesizer.optimalResults.gateCount = circuit.getNindividualOps();
        synthesizer.optimalResults.result    = logicbase::Result::SAT;
    }
    void HeuristicStrategy::runSplinter(int i, std::size_t circSplit, std::size_t split, const CouplingMap& reducedCM, const QubitSubset& qubitChoice, qc::QuantumComputation& circuit, std::shared_ptr<Results> r, CliffordSynthesizer* opt, const Configuration& configuration) {
        const Tableau targetTableau{circuit, 0, static_cast<std::size_t>((i + 1U)) * circSplit};
        const Tableau initTableau{circuit, 0, static_cast<std::size_t>(i) * circSplit};
        (*r) = opt->mainOptimization(split, reducedCM, qubitChoice, targetTableau,
                                     initTableau, configuration);
    }
} // namespace cs
