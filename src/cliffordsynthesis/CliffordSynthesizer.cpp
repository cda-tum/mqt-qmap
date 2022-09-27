/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
*/

#include "cliffordsynthesis/CliffordSynthesizer.hpp"
/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
*/

#include "LogicBlock/LogicBlock.hpp"
#include "LogicTerm/LogicTerm.hpp"
#include "LogicUtil/util_logicblock.hpp"
#include "cliffordsynthesis/GateEncoding.hpp"
#include "cliffordsynthesis/OptimizationStrategy.hpp"
#include "operations/OpType.hpp"
#include "operations/StandardOperation.hpp"
#include "utils.hpp"
#include "utils/logging.hpp"

namespace cs {
    void CliffordSynthesizer::synthesize(const Configuration& configuration) {
        TRACE() << "OptimizationStrategy: " << toString(configuration.strategy) << std::endl;
        TRACE() << "Target: " << TargetMetric::toString(configuration.target) << std::endl;
        TRACE() << "ReasoningEngine: " << toString(configuration.method) << std::endl;

        initResults();

        initCouplingMap(configuration.nqubits);

        auto                     totalStart = std::chrono::high_resolution_clock::now();
        std::vector<CouplingMap> reducedMaps;
        architecture.getReducedCouplingMaps(configuration.nqubits, reducedMaps);
        auto subsets =
                (configuration.chooseBest ? highestFidelityCouplingMap : reducedMaps);
        for (const auto& subset: subsets) {
            std::vector<std::uint16_t> qubitMap{Architecture::getQubitList(subset)};

            DEBUG() << "Reduced Coupling Map" << (configuration.chooseBest ? " (best)" : "") << ": ";
            std::stringstream strings;
            Architecture::printCouplingMap(subset, strings);
            DEBUG() << strings.str();
            DEBUG() << "Qubit Map: " << qubitMap;
            DEBUG() << "Coupling Map Fidelity: "
                    << Architecture::getAverageArchitectureFidelity(architecture.getCouplingMap(),
                                                                    std::set<std::uint16_t>(qubitMap.begin(), qubitMap.end()),
                                                                    architecture.getProperties());
            int timesteps =
                    configuration.initialTimesteps == 0 ? configuration.nqubits * configuration.nqubits : configuration.initialTimesteps;
            if (configuration.strategy == OptimizationStrategy::UseMinimizer) {
                runMinimizer(timesteps, subset, qubitMap, configuration);
            }
            if (configuration.strategy == OptimizationStrategy::StartLow) {
                runStartLow(timesteps, subset, qubitMap, configuration);
            }
            if (configuration.strategy == OptimizationStrategy::StartHigh) {
                runStartHigh(timesteps, subset, qubitMap, configuration);
            }
            if (configuration.strategy == OptimizationStrategy::MinMax) {
                runMinMax(timesteps, subset, qubitMap, configuration);
            }
            if (configuration.strategy == OptimizationStrategy::SplitIter) {
                runSplitIter(subset, qubitMap, configuration);
            }
            if (configuration.chooseBest && optimalResults.sat) {
                break;
            }
        }

        auto                          totalEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff     = totalEnd - totalStart;
        INFO() << "Time: " << diff.count() << std::endl;
        optimalResults.totalSeconds = diff.count();
    }

    Results CliffordSynthesizer::mainOptimization(
            std::uint32_t                                            timesteps,
            const std::set<std::pair<std::uint16_t, std::uint16_t>>& reducedCM,
            const std::vector<std::uint16_t>&                        qubitChoice,
            const Tableau& targetTableau, const Tableau& initialTableau,
            const Configuration& configuration) {
        std::unique_ptr<logicbase::LogicBlock> lb;
        using namespace logicbase;
        bool success = false;
        if (configuration.method == ReasoningEngine::Z3) {
            logicbase::LogicTerm::termType = TermType::BASE;
            if (configuration.strategy == OptimizationStrategy::UseMinimizer || configuration.strategy == OptimizationStrategy::SplitIter) {
                logicutil::Params params;
                params.addParam("pb.compile_equality", true);
                params.addParam("maxres.hill_climb", true);
                params.addParam("maxres.pivot_on_correction_set", false);
                lb = logicutil::getZ3LogicOptimizer(success, true, params);
            } else {
                logicutil::Params params;
                params.addParam("threads", static_cast<unsigned>(configuration.nThreads / 2));
                lb = logicutil::getZ3LogicBlock(success, true, params);
            }
        } else {
            return Results{};
        }
        if (!success) {
            throw QMAPException("Could not initialize Z3 logic block optimizer");
        }
        DEBUG() << "lb 1: " << lb.get() << std::endl;
        logicbase::LogicMatrix   x{};
        logicbase::LogicMatrix   z{};
        logicbase::LogicVector   r{};
        logicbase::LogicMatrix3D gS{};
        logicbase::LogicMatrix3D gC{};

        logicbase::LogicTerm changes = logicbase::LogicTerm(true);

        auto start = std::chrono::high_resolution_clock::now();
        /*
   * Tableau Variables x/z
   * k before gate k is applied
   * i column
   */
        std::stringstream xName{};
        std::stringstream zName{};
        std::stringstream rName{};
        for (int k = 0; k < timesteps + 1; ++k) {
            x.emplace_back();
            z.emplace_back();
            for (int i = 0; i < configuration.nqubits; ++i) {
                xName.str("");
                zName.str("");
                xName << "x_" << k << "_" << i;
                zName << "z_" << k << "_" << i;
                x.back().push_back(
                        lb->makeVariable(xName.str(), CType::BITVECTOR, configuration.nqubits));
                z.back().push_back(
                        lb->makeVariable(zName.str(), CType::BITVECTOR, configuration.nqubits));
            }
            rName.str("");
            rName << "r_" << k;
            r.push_back(lb->makeVariable(rName.str(), CType::BITVECTOR, configuration.nqubits));
        }

        /*
   * Gate Variables
   * k before gates k are applied
   * i qubit 1
   * j qubit
   */
        std::stringstream gName{};
        for (int gateStep = 0; gateStep < timesteps + 1; ++gateStep) {
            gS.emplace_back();
            for (auto gate: Gates::SINGLE_QUBIT) {
                gS.back().emplace_back();
                for (int j = 0; j < configuration.nqubits; ++j) {
                    gName.str("");
                    gName << "g_" << gateStep << "_" << Gates::gateName(gate) << "_" << j;
                    gS.back().back().push_back(lb->makeVariable(gName.str()));
                }
            }
        }
        for (int gateStep = 0; gateStep < timesteps + 1; ++gateStep) {
            gC.emplace_back();
            for (int j = 0; j < configuration.nqubits; ++j) {
                gC.back().emplace_back();
                for (int l = 0; l < configuration.nqubits; ++l) {
                    gName.str("");
                    gName << "g_" << gateStep << "_CNOT_" << j << "_" << l;
                    gC.back().back().push_back(lb->makeVariable(gName.str()));
                }
            }
        }

        assertTableau(SynthesisData{configuration.nqubits, timesteps, reducedCM, qubitChoice, lb, x, z, r, gS, gC}, initialTableau, 0);
        assertTableau(SynthesisData{configuration.nqubits, timesteps, reducedCM, qubitChoice, lb, x, z, r, gS, gC}, targetTableau, timesteps);

        makeSpecificEncoding(SynthesisData{configuration.nqubits, timesteps, reducedCM, qubitChoice, lb, x, z, r, gS, gC}, configuration);

        auto                          formulation = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff        = formulation - start;
        INFO() << "Time to prduce Formulation: " << diff.count() << std::endl;

        lb->produceInstance();

        auto modGen = std::chrono::high_resolution_clock::now();
        diff        = modGen - formulation;
        INFO() << "Time to generate Model: " << diff.count() << std::endl;

        TRACE() << "Clauses: " << TermImpl::getNextId(lb.get()) << std::endl;
        TRACE() << "None Terms: " << TermImpl::getNextId() << std::endl;

        Result result = lb->solve();

        auto end = std::chrono::high_resolution_clock::now();
        diff     = end - modGen;
        INFO() << "Time to solve Model: " << diff.count() << std::endl;
        std::chrono::duration<double, std::milli> elapsedMilliseconds = end - start;
        Results                                   results{};
        results.verbose          = configuration.verbosity;
        results.chooseBest       = configuration.chooseBest;
        results.nqubits          = configuration.nqubits;
        results.initialTimesteps = timesteps;
        results.strategy         = configuration.strategy;
        results.target           = configuration.target;
        results.totalSeconds     = elapsedMilliseconds.count();
        results.sat              = result == Result::SAT;
        results.doubleFidelity   = architecture.getFidelityTable();
        results.singleFidelity   = architecture.getSingleQubitFidelities();
        results.resultCM         = architecture.getCouplingMap();
        results.resultTableaus.clear();

        if (result == Result::SAT) {
            results.result               = logicbase::Result::SAT;
            Model*                 model = lb->getModel();
            qc::QuantumComputation resultCircuit;
            resultCircuit.addQubitRegister(configuration.nqubits);
            results.gateCount = 0;
            results.depth     = 0;
            results.fidelity  = 1;
            for (int gateStep = 0; gateStep < timesteps + 1; ++gateStep) {
                int oldGateCount = results.gateCount;
                TRACE() << "Gate Step: " << gateStep << std::endl
                        << " Actual gate count: " << results.gateCount << std::endl
                        << " Depth: " << results.depth << std::endl
                        << " Fidelity: " << results.fidelity << std::endl;
                if (gateStep > 0) {
                    for (int a = 0; a < configuration.nqubits; ++a) {
                        for (auto gate: Gates::SINGLE_QUBIT_WITHOUT_NOP) {
                            if (model->getBoolValue(gS[gateStep][Gates::toIndex(gate)][a], lb.get())) {
                                resultCircuit.emplace_back<qc::StandardOperation>(configuration.nqubits, a, Gates::toOpType(gate));
                                if (architecture.isCalibrationDataAvailable()) {
                                    results.fidelity *= (architecture.getSingleQubitFidelities()[a]);
                                }
                                TRACE() << Gates::gateName(gate) << "(" << a << ")" << std::endl;
                                if (architecture.isCalibrationDataAvailable()) {
                                    TRACE() << " Fidelity: " << architecture.getSingleQubitFidelities()[a]
                                            << std::endl;
                                }
                                ++results.gateCount;
                            }
                        }
                        for (int b = 0; b < configuration.nqubits; ++b) {
                            if (model->getBoolValue(gC[gateStep][a][b], lb.get())) {
                                results.gateCount++;
                                resultCircuit.emplace_back<qc::StandardOperation>(
                                        configuration.nqubits, dd::Control{static_cast<dd::Qubit>(a)}, b, qc::X);
                                if (architecture.isCalibrationDataAvailable()) {
                                    results.fidelity *=
                                            (1 - std::log(architecture.getFidelityTable()[qubitChoice.at(a)]
                                                                                         [qubitChoice.at(b)]));
                                }
                                TRACE() << "X(" << a << "," << b << ")" << std::endl;
                                if (architecture.isCalibrationDataAvailable()) {
                                    TRACE() << "Fidelity: "
                                            << (1 - std::log(architecture.getFidelityTable()[qubitChoice.at(a)]
                                                                                            [qubitChoice.at(b)]))
                                            << std::endl;
                                }
                            }
                        }
                    }
                }
                if (oldGateCount < results.gateCount) {
                    results.depth++;
                }
                auto tableau = results.resultTableaus.emplace_back();
                Tableau::generateTableau(tableau, resultCircuit);
                results.resultTableaus.back() = tableau;
                Tableau::initTableau(modelTableau, configuration.nqubits);
                for (int i = 0; i < configuration.nqubits; ++i) {
                    modelTableau.populateTableauFrom(model->getBitvectorValue(x[gateStep][i], lb.get()),
                                                     configuration.nqubits, i);
                    modelTableau.populateTableauFrom(model->getBitvectorValue(z[gateStep][i], lb.get()),
                                                     configuration.nqubits, i + configuration.nqubits);
                }
                modelTableau.populateTableauFrom(model->getBitvectorValue(r[gateStep], lb.get()),
                                                 configuration.nqubits, 2 * configuration.nqubits);
                if (configuration.verbosity >= 5) {
                    TRACE() << modelTableau;
                }
            }
            results.resultCircuit = resultCircuit.clone();
        }
        lb->reset();
        if (result == Result::SAT) {
            DEBUG() << "SAT" << std::endl;
            return results;
        }
        {
            results.result = logicbase::Result::UNSAT;
            DEBUG() << "UNSAT" << std::endl;
            return results;
        }
        return results;
    }

    void CliffordSynthesizer::runMinimizer(
            int timesteps, const CouplingMap& reducedCM,
            const std::vector<std::uint16_t>& qubitChoice, const Configuration& configuration) {
        DEBUG() << "Running minimizer" << std::endl;
        Results r = mainOptimization(timesteps, reducedCM, qubitChoice, configuration.targetTableau, configuration.initialTableau,
                                     configuration);
        updateResults(r);
    }
    void CliffordSynthesizer::runStartLow(
            int timesteps, const CouplingMap& reducedCM,
            const std::vector<std::uint16_t>& qubitChoice, const Configuration& configuration) {
        DEBUG() << "Running start low" << std::endl;
        Results r;
        while (r.result != logicbase::Result::SAT || r.result == logicbase::Result::NDEF) {
            DEBUG() << "Current t=" << timesteps << std::endl;
            r = mainOptimization(timesteps, reducedCM, qubitChoice, configuration.targetTableau, configuration.initialTableau, configuration);
            updateResults(r);
            if (r.result == logicbase::Result::UNSAT) {
                timesteps *= 1.5;
            }
        }
    }
    void CliffordSynthesizer::runStartHigh(
            int timesteps, const CouplingMap& reducedCM,
            const std::vector<std::uint16_t>& qubitChoice, const Configuration& configuration) {
        DEBUG() << "Running start high" << std::endl;
        Results r;
        int     oldTimesteps = timesteps;
        while (r.result == logicbase::Result::SAT || r.result == logicbase::Result::NDEF) {
            DEBUG() << "Current t=" << timesteps << std::endl;
            r = mainOptimization(timesteps, reducedCM, qubitChoice, configuration.targetTableau, configuration.initialTableau, configuration);
            updateResults(r);
            if (r.result == logicbase::Result::SAT) {
                oldTimesteps = timesteps;
                timesteps *= 0.5;
            } else {
                timesteps = oldTimesteps;
            }
        }
    }
    void CliffordSynthesizer::runMinMax(
            int timesteps, const CouplingMap& reducedCM,
            const std::vector<std::uint16_t>& qubitChoice, const Configuration& configuration) {
        DEBUG() << "Running minmax" << std::endl;
        Results r;
        int     t     = timesteps;
        int     upper = timesteps;
        int     lower = 0;
        while (std::abs(upper - lower) > 1) {
            DEBUG() << "Current t=" << t << std::endl;
            r = mainOptimization(t, reducedCM, qubitChoice, configuration.targetTableau, configuration.initialTableau, configuration);
            updateResults(r);
            if (r.result == logicbase::Result::SAT) {
                upper = t;
            } else if (r.result == logicbase::Result::UNSAT) {
                lower = t;
            } else {
                break;
            }
            if (upper - lower < 1 && r.result == logicbase::Result::UNSAT) {
                upper *= 1.5;
            }
            t = lower + std::abs(upper - lower) / 2;
        }
    }

    void CliffordSynthesizer::runSplinter(
            int i, unsigned int circuitSplit, unsigned int split,
            const CouplingMap& reducedCM, const std::vector<std::uint16_t>& qubitChoice,
            qc::QuantumComputation& circuit, Results* r,
            CliffordSynthesizer* opt, const Configuration& configuration) {
        Tableau targetTableau{};
        Tableau::generateTableau(targetTableau, circuit, 0, (i + 1U) * circuitSplit);
        Tableau initTableau{};
        Tableau::generateTableau(initTableau, circuit, 0, i * circuitSplit);
        (*r) = opt->mainOptimization(split, reducedCM, qubitChoice, targetTableau,
                                     initTableau, configuration);
    };

    void CliffordSynthesizer::runSplitIter(
            const CouplingMap&                reducedCM,
            const std::vector<std::uint16_t>& qubitChoice, const Configuration& configuration) {
        if (configuration.targetCircuit.size() < 2) {
            return;
        }
        DEBUG() << "Running split iter" << std::endl;
        Tableau                   fullTableau  = configuration.targetTableau;
        auto                      circuitSplit = static_cast<unsigned int>(std::log(configuration.targetCircuit.getNindividualOps()));
        int                       split        = std::min(5, configuration.nqubits / 2);
        std::vector<std::thread*> threads;
        std::vector<Results*>     results;
        int                       nThreads = 4;
        qc::QuantumComputation    circuit  = configuration.targetCircuit.clone();
        while (true) {
            results.clear();
            DEBUG() << "Current split size: " << split << std::endl;
            DEBUG() << "Current circuit split size: " << circuitSplit << std::endl;
            auto    start = std::chrono::high_resolution_clock::now();
            Results totalResult;
            totalResult.result = logicbase::Result::SAT;
            totalResult.resultCircuit.addQubitRegister(configuration.nqubits);
            for (size_t i = 0; i * circuitSplit < configuration.targetCircuit.getNindividualOps();
                 i += nThreads) {
                threads.clear();
                DEBUG() << "Currently at " << i * circuitSplit << " of "
                        << configuration.targetCircuit.getNindividualOps() << std::endl;
                for (int j = 0; j < nThreads; j++) {
                    auto* r = new Results();
                    auto* t = new std::thread(CliffordSynthesizer::runSplinter, i, circuitSplit,
                                              split, std::ref(reducedCM), std::ref(qubitChoice),
                                              std::ref(circuit), r, this, configuration);
                    threads.push_back(t);
                    results.push_back(r);
                }
                for (auto* t: threads) {
                    t->join();
                }
                for (auto* t: threads) {
                    delete t;
                }
                for (auto* r: results) {
                    if (r->result == logicbase::Result::UNSAT) {
                        totalResult.result = logicbase::Result::UNSAT;
                        break;
                    }
                }
                if (totalResult.result == logicbase::Result::UNSAT) {
                    DEBUG() << "UNSAT, increasing split size." << std::endl;
                    split += std::max(1.0, split * 0.2);
                    break;
                }
            }
            for (auto* r: results) {
                for (const auto& gate: r->resultCircuit) {
                    totalResult.resultCircuit.insert(totalResult.resultCircuit.end(),
                                                     gate->clone());
                }
                delete r;
            }
            if (totalResult.result == logicbase::Result::SAT) {
                Tableau resultingTableau{};
                Tableau::generateTableau(resultingTableau, totalResult.resultCircuit);
                DEBUG() << "Equality (Results): "
                        << ((fullTableau == resultingTableau) ? "True" : "False")
                        << std::endl;
                DEBUG() << "Original Circuit size: " << configuration.targetCircuit.getNindividualOps()
                        << std::endl;
                DEBUG() << "Optimized Circuit size: "
                        << totalResult.resultCircuit.getNindividualOps() << std::endl;
                TRACE() << "Resulting Circuit: " << std::endl;
                std::ostringstream ss;
                totalResult.resultCircuit.dump(ss, qc::Format::OpenQASM);
                TRACE() << ss.str() << std::endl;
                auto                          end  = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> diff = end - start;
                INFO() << "Time for complete run: " << diff.count() << std::endl;
                if (configuration.targetCircuit.getNindividualOps() ==
                    totalResult.resultCircuit.getNindividualOps()) {
                    split *= 1.2;
                    break;
                }
                circuit = totalResult.resultCircuit.clone();
            }
        }
        optimalResults.resultCircuit = configuration.targetCircuit.clone();
        optimalResults.resultTableaus.emplace_back(configuration.targetTableau);
        optimalResults.gateCount = configuration.targetCircuit.getNindividualOps();
        optimalResults.result    = logicbase::Result::SAT;
    }

    void CliffordSynthesizer::assertTableau(const SynthesisData& data, const Tableau& tableau, std::uint32_t position) {
        for (unsigned int a = 0; a < data.nqubits; ++a) {
            data.lb->assertFormula(data.x[position][a] ==
                                   logicbase::LogicTerm(tableau.getBVFrom(a), data.nqubits));
            data.lb->assertFormula(
                    data.z[position][a] ==
                    logicbase::LogicTerm(tableau.getBVFrom(a + data.nqubits), data.nqubits));
        }
        data.lb->assertFormula(
                data.r[position] ==
                logicbase::LogicTerm(tableau.getBVFrom(2 * data.nqubits), data.nqubits));
    }

    void cs::CliffordSynthesizer::setArchitecture(Architecture lArchitecture) {
        architecture = std::move(lArchitecture);
    }
    void CliffordSynthesizer::initCouplingMap(std::uint32_t nqubits) {
        if (architecture.isArchitectureAvailable()) {
            auto& cm = highestFidelityCouplingMap.emplace_back();
            architecture.getHighestFidelityCouplingMap(nqubits, cm);
        }
    }
    void CliffordSynthesizer::initResults() {
        optimalResults                  = Results();
        optimalResults.architectureName = architecture.getName();
        resultCircuit.addQubitRegister(architecture.getNqubits());
    }
    void CliffordSynthesizer::makeGateEncoding(const SynthesisData& data, const Configuration& configuration) {
        GateEncoding::makeGateEncoding(data, configuration);
    }
}