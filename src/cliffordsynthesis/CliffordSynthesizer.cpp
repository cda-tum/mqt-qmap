/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
*/

#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#include "LogicBlock/LogicBlock.hpp"
#include "LogicTerm/LogicTerm.hpp"
#include "LogicUtil/util_logicblock.hpp"
#include "operations/OpType.hpp"
#include "operations/StandardOperation.hpp"
#include "utils.hpp"
#include "utils/logging.hpp"

void CliffordSynthesizer::optimize() {
    TRACE() << "Strategy: " << toString(strategy) << std::endl;
    TRACE() << "Target: " << toString(target) << std::endl;
    TRACE() << "Method: " << toString(method) << std::endl;

    auto                     totalStart = std::chrono::high_resolution_clock::now();
    std::vector<CouplingMap> reducedMaps;
    architecture.getReducedCouplingMaps(nqubits, reducedMaps);
    auto subsets =
            (chooseBest ? highestFidelityMap : reducedMaps);
    for (const auto& subset: subsets) {
        std::vector<std::uint16_t> qubitMap{Architecture::getQubitList(subset)};

        DEBUG() << "Reduced Coupling Map" << (chooseBest ? " (best)" : "") << ": ";
        std::stringstream strings;
        Architecture::printCouplingMap(subset, strings);
        DEBUG() << strings.str();
        DEBUG() << "Qubit Map: " << qubitMap;
        DEBUG() << "Coupling Map Fidelity: "
                << Architecture::getAverageArchitectureFidelity(architecture.getCouplingMap(),
                                                                std::set<std::uint16_t>(qubitMap.begin(), qubitMap.end()),
                                                                architecture.getProperties());
        int timesteps =
                initialTimesteps == 0 ? nqubits * nqubits : initialTimesteps;
        if (strategy == SynthesisStrategy::UseMinimizer) {
            runMinimizer(timesteps, subset, qubitMap);
        }
        if (strategy == SynthesisStrategy::StartLow) {
            runStartLow(timesteps, subset, qubitMap);
        }
        if (strategy == SynthesisStrategy::StartHigh) {
            runStartHigh(timesteps, subset, qubitMap);
        }
        if (strategy == SynthesisStrategy::MinMax) {
            runMinMax(timesteps, subset, qubitMap);
        }
        if (strategy == SynthesisStrategy::SplitIter) {
            runSplitIter(subset, qubitMap);
        }
        if (chooseBest && optimalResults.sat) {
            break;
        }
    }

    auto                          totalEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff     = totalEnd - totalStart;
    INFO() << "Time: " << diff.count() << std::endl;
    optimalResults.totalSeconds = diff.count();
}

CliffordOptimizationResults CliffordSynthesizer::mainOptimization(
        int                                            timesteps,
        const std::set<std::pair<std::uint16_t, std::uint16_t>>& reducedCM,
        const std::vector<std::uint16_t>& qubitChoice, Tableau& initialTab,
        Tableau& targetTab) {
    std::unique_ptr<LogicBlock> lb;
    using namespace logicbase;
    bool success = false;
    if (method == SynthesisMethod::Z3) {
        LogicTerm::termType = TermType::BASE;
        if (strategy == SynthesisStrategy::UseMinimizer || strategy == SynthesisStrategy::SplitIter) {
            logicutil::Params params;
            params.addParam("pb.compile_equality", true);
            params.addParam("maxres.hill_climb", true);
            params.addParam("maxres.pivot_on_correction_set", false);
            lb = logicutil::getZ3LogicOptimizer(success, true, params);
        } else {
            logicutil::Params params;
            params.addParam("threads", static_cast<unsigned>(nthreads / 2));
            lb = logicutil::getZ3LogicBlock(success, true, params);
        }
    } else {
        return CliffordOptimizationResults{};
    }
    if (!success) {
        throw QMAPException("Could not initialize Z3 logic block optimizer");
    }
    DEBUG() << "lb 1: " << lb.get() << std::endl;
    LogicMatrix   x{};
    LogicMatrix   z{};
    LogicVector   r{};
    LogicMatrix3D gS{};
    LogicMatrix3D gC{};

    LogicTerm changes = LogicTerm(true);

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
        for (int i = 0; i < nqubits; ++i) {
            xName.str("");
            zName.str("");
            xName << "x_" << k << "_" << i;
            zName << "z_" << k << "_" << i;
            x.back().push_back(
                    lb->makeVariable(xName.str(), CType::BITVECTOR, nqubits));
            z.back().push_back(
                    lb->makeVariable(zName.str(), CType::BITVECTOR, nqubits));
        }
        rName.str("");
        rName << "r_" << k;
        r.push_back(lb->makeVariable(rName.str(), CType::BITVECTOR, nqubits));
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
            for (int j = 0; j < nqubits; ++j) {
                gName.str("");
                gName << "g_" << gateStep << "_" << Gates::gateName(gate) << "_" << j;
                gS.back().back().push_back(lb->makeVariable(gName.str()));
            }
        }
    }
    for (int gateStep = 0; gateStep < timesteps + 1; ++gateStep) {
        gC.emplace_back();
        for (int j = 0; j < nqubits; ++j) {
            gC.back().emplace_back();
            for (int l = 0; l < nqubits; ++l) {
                gName.str("");
                gName << "g_" << gateStep << "_CNOT_" << j << "_" << l;
                gC.back().back().push_back(lb->makeVariable(gName.str()));
            }
        }
    }

    assertTableau(initialTab, lb, x, z, r, nqubits, 0);
    assertTableau(targetTab, lb, x, z, r, nqubits, timesteps);

    if (target == SynthesisTarget::DEPTH) {
        makeDepthOptimizer(timesteps, reducedCM, qubitChoice, lb, x, z, r, gS,
                           gC);
    } else if (target == SynthesisTarget::GATES || target == SynthesisTarget::GATES_ONLY_CNOT) {
        makeGateOptimizer(timesteps, reducedCM, qubitChoice, lb, x, z, r, gS,
                          gC);
    } else if (target == SynthesisTarget::FIDELITY) {
        makeFidelityOptimizer(timesteps, reducedCM, qubitChoice, lb, x, z, r, gS,
                              gC);
    } else {
        ERROR() << "Unknown target" << std::endl;
        return CliffordOptimizationResults{};
    }
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
    CliffordOptimizationResults               results{};
    results.verbose          = verbose;
    results.chooseBest       = chooseBest;
    results.nqubits          = nqubits;
    results.initialTimesteps = timesteps;
    results.strategy         = strategy;
    results.target           = target;
    results.totalSeconds     = elapsedMilliseconds.count();
    results.sat              = result == Result::SAT;
    results.doubleFidelity   = architecture.getFidelityTable();
    results.singleFidelity   = architecture.getSingleQubitFidelities();
    results.resultCM         = architecture.getCouplingMap();
    results.resultTableaus.clear();

    if (result == Result::SAT) {
        results.result               = SynthesisResult::SAT;
        Model*                 model = lb->getModel();
        qc::QuantumComputation resultCircuit;
        resultCircuit.addQubitRegister(nqubits);
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
                for (int a = 0; a < nqubits; ++a) {
                    for (auto gate: Gates::SINGLE_QUBIT_WITHOUT_NOP) {
                        if (model->getBoolValue(gS[gateStep][Gates::toIndex(gate)][a], lb.get())) {
                            resultCircuit.emplace_back<qc::StandardOperation>(nqubits, a, Gates::toOpType(gate));
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
                    for (int b = 0; b < nqubits; ++b) {
                        if (model->getBoolValue(gC[gateStep][a][b], lb.get())) {
                            results.gateCount++;
                            resultCircuit.emplace_back<qc::StandardOperation>(
                                    nqubits, dd::Control{static_cast<dd::Qubit>(a)}, b, qc::X);
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
            Tableau::initTableau(modelTableau, nqubits);
            for (int i = 0; i < nqubits; ++i) {
                modelTableau.populateTableauFrom(model->getBitvectorValue(x[gateStep][i], lb.get()),
                                                 nqubits, i);
                modelTableau.populateTableauFrom(model->getBitvectorValue(z[gateStep][i], lb.get()),
                                                 nqubits, i + nqubits);
            }
            modelTableau.populateTableauFrom(model->getBitvectorValue(r[gateStep], lb.get()),
                                             nqubits, 2 * nqubits);
            if (verbose >= 5) {
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
        results.result = SynthesisResult::UNSAT;
        DEBUG() << "UNSAT" << std::endl;
        return results;
    }
    return results;
}

void CliffordSynthesizer::makeDepthOptimizer(int timesteps, const std::set<std::pair<std::uint16_t, std::uint16_t>>& reducedCM, const std::vector<std::uint16_t>& qubitChoice, std::unique_ptr<LogicBlock>& lb, const logicbase::LogicMatrix& x, const logicbase::LogicMatrix& z, const logicbase::LogicVector& r, const logicbase::LogicMatrix3D& gS, const logicbase::LogicMatrix3D& gC) const {
    makeMultipleGateConstraints(lb, x, z, r, nqubits, timesteps, reducedCM,
                                qubitChoice, gS, gC);
    // COST
    if (strategy == SynthesisStrategy::UseMinimizer ||
        strategy == SynthesisStrategy::SplitIter) {
        LogicTerm cost = LogicTerm(0);
        for (int gateStep = 1; gateStep < timesteps + 1; ++gateStep) {
            LogicTerm anyGate = LogicTerm(true);
            for (int a = 0; a < nqubits; ++a) {
                for (auto gate: Gates::SINGLE_QUBIT_WITHOUT_NOP) {
                    anyGate = anyGate && !gS[gateStep][Gates::toIndex(gate)][a];
                }
                for (int b = 0; b <= a; ++b) {
                    if (a == b) {
                        continue;
                    }
                    anyGate = anyGate && !gC[gateStep][a][b] && !gC[gateStep][b][a];
                }
            }
            cost = cost + LogicTerm::ite(anyGate, LogicTerm(5), LogicTerm(0));
        }
        dynamic_cast<LogicBlockOptimizer*>(lb.get())->maximize(cost);
    }
}

void CliffordSynthesizer::makeFidelityOptimizer(int timesteps, const std::set<std::pair<std::uint16_t, std::uint16_t>>& reducedCM, const std::vector<std::uint16_t>& qubitChoice, std::unique_ptr<LogicBlock>& lb, const logicbase::LogicMatrix& x, const logicbase::LogicMatrix& z, const logicbase::LogicVector& r, const logicbase::LogicMatrix3D& gS, const logicbase::LogicMatrix3D& gC) const {
    if (!architecture.isArchitectureAvailable()) {
        util::fatal("No fidelity architecture specified in coupling map.");
    }
    makeMultipleGateConstraints(lb, x, z, r, nqubits, timesteps, reducedCM,
                                qubitChoice, gS, gC);
    // COST
    if (strategy == SynthesisStrategy::UseMinimizer ||
        strategy == SynthesisStrategy::SplitIter) {
        LogicTerm cost = LogicTerm(0);
        // For each edge in the coupling map, get the fidelity cost
        for (const auto& edge: reducedCM) {
            LogicTerm fidelity = LogicTerm(
                    (1 - std::log(architecture.getFidelityTable()[edge.first][edge.second])) * 1000);
            auto a = std::find(qubitChoice.begin(), qubitChoice.end(), edge.first);
            auto b = std::find(qubitChoice.begin(), qubitChoice.end(), edge.second);
            if (a == qubitChoice.end() || b == qubitChoice.end()) {
                util::fatal("Coupling map contains invalid qubit.");
            }
            // at each time t if there is a gate on the edge, add the cost
            for (int gateStep = 0; gateStep < timesteps; ++gateStep) {
                cost = cost + (gC[gateStep][std::distance(qubitChoice.begin(), a)][std::distance(qubitChoice.begin(), b)] * fidelity);
            }
        }
        // For each qubit, get the fidelity cost
        for (int a = 0; a < nqubits; ++a) {
            LogicTerm fidelity =
                    LogicTerm((1 - std::log(architecture.getSingleQubitFidelities()[a])) * 1000);
            // at each time t if there is a gate on a, add the cost
            for (int gateStep = 0; gateStep < timesteps; ++gateStep) {
                for (auto gate: Gates::SINGLE_QUBIT_WITHOUT_NOP) {
                    cost = cost + (gS[gateStep][Gates::toIndex(gate)][a] * fidelity);
                }
            }
        }
        dynamic_cast<LogicBlockOptimizer*>(lb.get())->minimize(cost);
        cost = LogicTerm(0);
        for (int gateStep = 0; gateStep < timesteps; ++gateStep) {
            for (int a = 0; a < nqubits; ++a) {
                cost = cost + gS[gateStep][1][a] + gS[gateStep][2][a];
                for (int b = 0; b < nqubits; ++b) {
                    cost = cost + gC[gateStep][a][b];
                }
            }
        }
        dynamic_cast<LogicBlockOptimizer*>(lb.get())->maximize(cost);
    }
}

void CliffordSynthesizer::makeGateOptimizer(int timesteps, const std::set<std::pair<std::uint16_t, std::uint16_t>>& reducedCM, const std::vector<std::uint16_t>& qubitChoice, std::unique_ptr<LogicBlock>& lb, const logicbase::LogicMatrix& x, const logicbase::LogicMatrix& z, const logicbase::LogicVector& r, const logicbase::LogicMatrix3D& gS, const logicbase::LogicMatrix3D& gC) const {
    LogicTerm changes = LogicTerm(true);
    makeSingleGateConstraints(lb, x, z, r, nqubits, timesteps, reducedCM,
                              qubitChoice, gS, gC);
    // COST
    if (strategy == SynthesisStrategy::UseMinimizer) {
        LogicTerm cost = LogicTerm(0);
        for (int gateStep = 1; gateStep < timesteps + 1; ++gateStep) {
            for (int a = 0; a < nqubits; ++a) {
                if (target != SynthesisTarget::GATES_ONLY_CNOT) {
                    for (auto gate: Gates::SINGLE_QUBIT_WITHOUT_NOP) {
                        cost = cost + gS[gateStep][Gates::toIndex(gate)][a];
                    }
                }
                for (int b = 0; b <= a; ++b) {
                    if (a == b) {
                        continue;
                    }
                    cost = cost + gC[gateStep][a][b] + gC[gateStep][b][a];
                }
            }
        }
        dynamic_cast<LogicBlockOptimizer*>(lb.get())->minimize(cost);
    }
}

void CliffordSynthesizer::runMinimizer(
        int timesteps, const CouplingMap& reducedCM,
        const std::vector<std::uint16_t>& qubitChoice) {
    DEBUG() << "Running minimizer" << std::endl;
    CliffordOptimizationResults r = mainOptimization(timesteps, reducedCM, qubitChoice,
                                                     initialTableau, targetTableau);
    updateResults(r);
}
void CliffordSynthesizer::runStartLow(
        int timesteps, const CouplingMap& reducedCM,
        const std::vector<std::uint16_t>& qubitChoice) {
    DEBUG() << "Running start low" << std::endl;
    CliffordOptimizationResults r;
    while (r.result != SynthesisResult::SAT || r.result == SynthesisResult::UNDEF) {
        DEBUG() << "Current t=" << timesteps << std::endl;
        r = mainOptimization(timesteps, reducedCM, qubitChoice, initialTableau,
                             targetTableau);
        updateResults(r);
        if (r.result == SynthesisResult::UNSAT) {
            timesteps *= 1.5;
        }
    }
}
void CliffordSynthesizer::runStartHigh(
        int timesteps, const CouplingMap& reducedCM,
        const std::vector<std::uint16_t>& qubitChoice) {
    DEBUG() << "Running start high" << std::endl;
    CliffordOptimizationResults r;
    int                         oldTimesteps = timesteps;
    while (r.result == SynthesisResult::SAT || r.result == SynthesisResult::UNDEF) {
        DEBUG() << "Current t=" << timesteps << std::endl;
        r = mainOptimization(timesteps, reducedCM, qubitChoice, initialTableau,
                             targetTableau);
        updateResults(r);
        if (r.result == SynthesisResult::SAT) {
            oldTimesteps = timesteps;
            timesteps *= 0.5;
        } else {
            timesteps = oldTimesteps;
        }
    }
}
void CliffordSynthesizer::runMinMax(
        int timesteps, const CouplingMap& reducedCM,
        const std::vector<std::uint16_t>& qubitChoice) {
    DEBUG() << "Running minmax" << std::endl;
    CliffordOptimizationResults r;
    int                         t     = timesteps;
    int                         upper = timesteps;
    int                         lower = 0;
    while (std::abs(upper - lower) > 1) {
        DEBUG() << "Current t=" << t << std::endl;
        r = mainOptimization(t, reducedCM, qubitChoice, initialTableau,
                             targetTableau);
        updateResults(r);
        if (r.result == SynthesisResult::SAT) {
            upper = t;
        } else if (r.result == SynthesisResult::UNSAT) {
            lower = t;
        } else {
            break;
        }
        if (upper - lower < 1 && r.result == SynthesisResult::UNSAT) {
            upper *= 1.5;
        }
        t = lower + std::abs(upper - lower) / 2;
    }
}

void CliffordSynthesizer::runSplinter(
        int i, unsigned int circuitSplit, unsigned int split,
        const CouplingMap& reducedCM, const std::vector<std::uint16_t>& qubitChoice,
        qc::QuantumComputation& circuit, CliffordOptimizationResults* r,
        CliffordSynthesizer* opt) {
    Tableau targetTableau{};
    Tableau::generateTableau(targetTableau, circuit, 0, (i + 1U) * circuitSplit);
    Tableau initTableau{};
    Tableau::generateTableau(initTableau, circuit, 0, i * circuitSplit);
    (*r) = opt->mainOptimization(split, reducedCM, qubitChoice, targetTableau,
                                 initTableau);
};

void CliffordSynthesizer::runSplitIter(
        const CouplingMap&           reducedCM,
        const std::vector<std::uint16_t>& qubitChoice) {
    if (circuit.size() < 2) {
        return;
    }
    DEBUG() << "Running split iter" << std::endl;
    Tableau                                   fullTableau  = targetTableau;
    auto                                      circuitSplit = static_cast<unsigned int>(std::log(circuit.getNindividualOps()));
    int                                       split        = std::min(5, nqubits / 2);
    std::vector<std::thread*>                 threads;
    std::vector<CliffordOptimizationResults*> results;
    int                                       nThreads = nthreads;
    while (true) {
        results.clear();
        DEBUG() << "Current split size: " << split << std::endl;
        DEBUG() << "Current circuit split size: " << circuitSplit << std::endl;
        auto                        start = std::chrono::high_resolution_clock::now();
        CliffordOptimizationResults totalResult;
        totalResult.result = SynthesisResult::SAT;
        totalResult.resultCircuit.addQubitRegister(nqubits);
        for (size_t i = 0; i * circuitSplit < circuit.getNindividualOps();
             i += nThreads) {
            threads.clear();
            DEBUG() << "Currently at " << i * circuitSplit << " of "
                    << circuit.getNindividualOps() << std::endl;
            for (int j = 0; j < nThreads; j++) {
                auto* r = new CliffordOptimizationResults();
                auto* t = new std::thread(CliffordSynthesizer::runSplinter, i, circuitSplit,
                                          split, std::ref(reducedCM), std::ref(qubitChoice),
                                          std::ref(circuit), r, this);
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
                if (r->result == SynthesisResult::UNSAT) {
                    totalResult.result = SynthesisResult::UNSAT;
                    break;
                }
            }
            if (totalResult.result == SynthesisResult::UNSAT) {
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
        if (totalResult.result == SynthesisResult::SAT) {
            Tableau resultingTableau{};
            Tableau::generateTableau(resultingTableau, totalResult.resultCircuit);
            DEBUG() << "Equality (Results): "
                    << ((fullTableau == resultingTableau) ? "True" : "False")
                    << std::endl;
            DEBUG() << "Original Circuit size: " << circuit.getNindividualOps()
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
            if (circuit.getNindividualOps() ==
                totalResult.resultCircuit.getNindividualOps()) {
                split *= 1.2;
                break;
            }
            circuit = totalResult.resultCircuit.clone();
        }
    }
    optimalResults.resultCircuit = circuit.clone();
    optimalResults.resultTableaus.emplace_back(targetTableau);
    optimalResults.gateCount = circuit.getNindividualOps();
    optimalResults.result    = SynthesisResult::SAT;
}

void CliffordSynthesizer::updateResults(CliffordOptimizationResults& results) {
    if (!results.sat) {
        return;
    }
    switch (target) {
        case SynthesisTarget::GATES:
        case SynthesisTarget::GATES_ONLY_CNOT:
            if (results.gateCount < optimalResults.gateCount ||
                optimalResults.gateCount == 0) {
                optimalResults = results;
            }
            break;
        case SynthesisTarget::DEPTH:
            if (results.depth < optimalResults.depth || optimalResults.depth == 0) {
                optimalResults = results;
            }
            break;
        case SynthesisTarget::FIDELITY:
            if (results.fidelity >= optimalResults.fidelity ||
                optimalResults.fidelity == 0) {
                optimalResults = results;
            }
            break;
        default:
            break;
    }
}

void CliffordSynthesizer::assertTableau(const Tableau& tableau, std::unique_ptr<LogicBlock>& lb,
                                      const LogicMatrix& x,
                                      const LogicMatrix& z,
                                      const LogicVector& r, int nqubits,
                                      int position) {
    for (int a = 0; a < nqubits; ++a) {
        lb->assertFormula(x[position][a] ==
                          LogicTerm(tableau.getBVFrom(a), nqubits));
        lb->assertFormula(
                z[position][a] ==
                LogicTerm(tableau.getBVFrom(a + nqubits), nqubits));
    }
    lb->assertFormula(
            r[position] ==
            LogicTerm(tableau.getBVFrom(2 * nqubits), nqubits));
}

void CliffordSynthesizer::makeSingleGateConstraints(
        std::unique_ptr<LogicBlock>& lb, const LogicMatrix& x, const LogicMatrix& z,
        const LogicVector& r, int nqubits, int timesteps,
        const std::set<std::pair<std::uint16_t, std::uint16_t>>& reducedCM,
        const std::vector<std::uint16_t>& qubitChoice, const LogicMatrix3D& gS,
        const LogicMatrix3D& gC) {
    LogicTerm changes = LogicTerm(true);
    // CONSISTENCY
    // One gate per qubit, per step
    for (int gateStep = 1; gateStep < timesteps + 1; ++gateStep) {
        std::vector<LogicTerm> vars{};
        for (int a = 0; a < nqubits; ++a) {
            for (auto gate: Gates::SINGLE_QUBIT) {
                vars.emplace_back(gS[gateStep][Gates::toIndex(gate)][a]);
            }
            for (int b = 0; b < nqubits; ++b) {
                if (a == b || reducedCM.find({qubitChoice.at(a), qubitChoice.at(b)}) ==
                                      reducedCM.end()) {
                    continue;
                }
                vars.emplace_back(gC[gateStep][a][b]);
            }
        }
        lb->assertFormula(ExactlyOneCMDR(
                groupVars(vars, static_cast<std::size_t>(vars.size() / 2U)),
                LogicTerm::noneTerm(), lb.get()));
    }

    // GATE CONSTRAINTS
    for (int gateStep = 1; gateStep < timesteps + 1; ++gateStep) {
        for (int a = 0; a < nqubits; ++a) {
            // NO GATE
            changes = (x[gateStep][a] == x[gateStep - 1][a]);
            changes = changes && (z[gateStep][a] == z[gateStep - 1][a]);

            for (int b = 0; b < nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (x[gateStep][b] == x[gateStep - 1][b]);
                changes = changes && (z[gateStep][b] == z[gateStep - 1][b]);
            }

            changes = changes && (r[gateStep] == r[gateStep - 1]);
            changes = LogicTerm::implies(gS[gateStep][0][a], changes);
            lb->assertFormula(changes);

            // H
            changes = (z[gateStep][a] == x[gateStep - 1][a]);
            changes = changes && (x[gateStep][a] == z[gateStep - 1][a]);

            for (int b = 0; b < nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (x[gateStep][b] == x[gateStep - 1][b]);
                changes = changes && (z[gateStep][b] == z[gateStep - 1][b]);
            }

            changes = changes &&
                      (r[gateStep] == (r[gateStep - 1] ^
                                       (x[gateStep - 1][a] & z[gateStep - 1][a])));
            changes = LogicTerm::implies(gS[gateStep][1][a], changes);

            lb->assertFormula(changes);

            // S
            changes =
                    (z[gateStep][a] == (z[gateStep - 1][a] ^ x[gateStep - 1][a]));
            changes = changes && (x[gateStep][a] == x[gateStep - 1][a]);

            for (int b = 0; b < nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (x[gateStep][b] == x[gateStep - 1][b]);
                changes = changes && (z[gateStep][b] == z[gateStep - 1][b]);
            }

            changes = changes &&
                      (r[gateStep] == (r[gateStep - 1] ^
                                       (x[gateStep - 1][a] & z[gateStep - 1][a])));
            changes = LogicTerm::implies(gS[gateStep][2][a], changes);
            lb->assertFormula(changes);

            // Sdag
            changes =
                    (z[gateStep][a] == (z[gateStep - 1][a] ^ x[gateStep - 1][a]));
            changes = changes && (x[gateStep][a] == x[gateStep - 1][a]);

            for (int b = 0; b < nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (x[gateStep][b] == x[gateStep - 1][b]);
                changes = changes && (z[gateStep][b] == z[gateStep - 1][b]);
            }

            changes = changes &&
                      (r[gateStep] == (r[gateStep - 1] ^
                                       (x[gateStep - 1][a] & (x[gateStep - 1][a] ^ z[gateStep - 1][a]))));
            changes = LogicTerm::implies(gS[gateStep][Gates::toIndex(Gates::GATES::Sdag)][a], changes);
            lb->assertFormula(changes);

            // Z
            changes =
                    (z[gateStep][a] == z[gateStep - 1][a]);
            changes = changes && (x[gateStep][a] == x[gateStep - 1][a]);

            for (int b = 0; b < nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (x[gateStep][b] == x[gateStep - 1][b]);
                changes = changes && (z[gateStep][b] == z[gateStep - 1][b]);
            }

            changes = changes &&
                      (r[gateStep] == (r[gateStep - 1] ^ x[gateStep - 1][a]));
            changes = LogicTerm::implies(gS[gateStep][Gates::toIndex(Gates::GATES::Z)][a], changes);
            lb->assertFormula(changes);

            // X
            changes =
                    (z[gateStep][a] == z[gateStep - 1][a]);
            changes = changes && (x[gateStep][a] == x[gateStep - 1][a]);

            for (int b = 0; b < nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (x[gateStep][b] == x[gateStep - 1][b]);
                changes = changes && (z[gateStep][b] == z[gateStep - 1][b]);
            }

            changes = changes &&
                      (r[gateStep] == (r[gateStep - 1] ^ z[gateStep - 1][a]));
            changes = LogicTerm::implies(gS[gateStep][Gates::toIndex(Gates::GATES::X)][a], changes);
            lb->assertFormula(changes);

            // Y
            changes =
                    (z[gateStep][a] == z[gateStep - 1][a]);
            changes = changes && (x[gateStep][a] == x[gateStep - 1][a]);

            for (int b = 0; b < nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (x[gateStep][b] == x[gateStep - 1][b]);
                changes = changes && (z[gateStep][b] == z[gateStep - 1][b]);
            }

            changes = changes &&
                      (r[gateStep] == (r[gateStep - 1] ^ (z[gateStep - 1][a]) ^ x[gateStep - 1][a]));
            changes = LogicTerm::implies(gS[gateStep][Gates::toIndex(Gates::GATES::Y)][a], changes);
            lb->assertFormula(changes);

            // CNOT
            for (int b = 0; b < nqubits; ++b) {
                if (reducedCM.find({qubitChoice.at(a), qubitChoice.at(b)}) ==
                    reducedCM.end()) {
                    lb->assertFormula(!gC[gateStep][a][b]);
                } else {
                    changes =
                            (r[gateStep] == (r[gateStep - 1] ^
                                             ((x[gateStep - 1][a] & z[gateStep - 1][b]) &
                                              ((x[gateStep - 1][b] ^ z[gateStep - 1][a]) ^
                                               LogicTerm((1 << nqubits) - 1, nqubits)))));
                    changes = changes && (x[gateStep][b] ==
                                          (x[gateStep - 1][b] ^ x[gateStep - 1][a]));
                    changes = changes && (z[gateStep][a] ==
                                          (z[gateStep - 1][a] ^ z[gateStep - 1][b]));
                    changes = changes && (x[gateStep][a] == x[gateStep - 1][a]);
                    changes = changes && (z[gateStep][b] == z[gateStep - 1][b]);

                    for (int c = 0; c < nqubits; ++c) { // All other entries do not change
                        if (a == c || b == c) {
                            continue;
                        }
                        changes = changes && (x[gateStep][c] == x[gateStep - 1][c]);
                        changes = changes && (z[gateStep][c] == z[gateStep - 1][c]);
                    }

                    changes = LogicTerm::implies(gC[gateStep][a][b], changes);
                    lb->assertFormula(changes);
                }
            }
        }
    }
}
void CliffordSynthesizer::makeMultipleGateConstraints(
        std::unique_ptr<LogicBlock>& lb, const LogicMatrix& x, const LogicMatrix& z,
        const LogicVector& r, int nqubits, int timesteps,
        const std::set<std::pair<std::uint16_t, std::uint16_t>>& reducedCM,
        const std::vector<std::uint16_t>& qubitChoice, const LogicMatrix3D& gS,
        const LogicMatrix3D& gC) {
    LogicTerm changes = LogicTerm(true);
    // CONSISTENCY
    // One gate per qubit, per step
    for (int gateStep = 1; gateStep < timesteps + 1; ++gateStep) {
        for (int a = 0; a < nqubits; ++a) {
            std::vector<LogicTerm> vars{};
            for (auto gate: Gates::SINGLE_QUBIT) {
                vars.emplace_back(gS[gateStep][Gates::toIndex(gate)][a]);
            }
            for (int b = 0; b < nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                vars.emplace_back(gC[gateStep][a][b]);
                vars.emplace_back(gC[gateStep][b][a]);
            }
            lb->assertFormula(ExactlyOneCMDR(
                    groupVars(vars, static_cast<std::size_t>(vars.size() / 2)),
                    LogicTerm::noneTerm(), lb.get()));
        }
    }
    // Maximum any combination of 1 and 2 qubit gates adding up to n
    for (int gateStep = 1; gateStep < timesteps + 1; ++gateStep) {
        changes = LogicTerm(0);
        for (int a = 0; a < nqubits; ++a) {
            for (auto gate: Gates::SINGLE_QUBIT) {
                changes = changes + gS[gateStep][Gates::toIndex(gate)][a];
            }
            for (int b = 0; b < nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes + gC[gateStep][a][b] + gC[gateStep][a][b];
            }
        }
        changes = changes < LogicTerm(static_cast<int>(nqubits + 1));
        lb->assertFormula(changes);
    }

    // GATE CONSTRAINTS
    for (int gateStep = 1; gateStep < timesteps + 1; ++gateStep) {
        LogicTerm rChanges = r[gateStep - 1];
        for (int a = 0; a < nqubits; ++a) {
            // NO GATE
            changes = LogicTerm(true);
            changes = changes && (x[gateStep][a] == x[gateStep - 1][a]);
            changes = changes && (z[gateStep][a] == z[gateStep - 1][a]);
            changes = LogicTerm::implies(gS[gateStep][0][a], changes);
            lb->assertFormula(changes);

            // H
            changes = LogicTerm(true);
            changes = changes && (z[gateStep][a] == x[gateStep - 1][a]);
            changes = changes && (x[gateStep][a] == z[gateStep - 1][a]);

            rChanges = LogicTerm::ite(
                    gS[gateStep][1][a],
                    rChanges ^ (x[gateStep - 1][a] & z[gateStep - 1][a]), rChanges);
            changes = LogicTerm::implies(gS[gateStep][1][a], changes);

            lb->assertFormula(changes);

            // S
            changes = LogicTerm(true);
            changes = changes && (z[gateStep][a] ==
                                  (z[gateStep - 1][a] ^ x[gateStep - 1][a]));
            changes = changes && (x[gateStep][a] == x[gateStep - 1][a]);

            rChanges = LogicTerm::ite(
                    gS[gateStep][2][a],
                    rChanges ^ (x[gateStep - 1][a] & z[gateStep - 1][a]), rChanges);
            changes = LogicTerm::implies(gS[gateStep][2][a], changes);
            lb->assertFormula(changes);

            // Sdag
            changes =
                    (z[gateStep][a] == (z[gateStep - 1][a] ^ x[gateStep - 1][a]));
            changes = changes && (x[gateStep][a] == x[gateStep - 1][a]);

            rChanges = LogicTerm::ite(
                    gS[gateStep][Gates::toIndex(Gates::GATES::Sdag)][a],
                    rChanges ^ (x[gateStep - 1][a] & (x[gateStep - 1][a] ^ z[gateStep - 1][a])), rChanges);
            changes = LogicTerm::implies(gS[gateStep][Gates::toIndex(Gates::GATES::Sdag)][a], changes);
            lb->assertFormula(changes);

            // Z
            changes =
                    (z[gateStep][a] == z[gateStep - 1][a]);
            changes = changes && (x[gateStep][a] == x[gateStep - 1][a]);

            rChanges = LogicTerm::ite(
                    gS[gateStep][Gates::toIndex(Gates::GATES::Z)][a],
                    rChanges ^ (x[gateStep - 1][a]), rChanges);
            changes = LogicTerm::implies(gS[gateStep][Gates::toIndex(Gates::GATES::Z)][a], changes);
            lb->assertFormula(changes);

            // X
            changes =
                    (z[gateStep][a] == z[gateStep - 1][a]);
            changes = changes && (x[gateStep][a] == x[gateStep - 1][a]);

            rChanges = LogicTerm::ite(
                    gS[gateStep][Gates::toIndex(Gates::GATES::X)][a],
                    rChanges ^ (z[gateStep - 1][a]), rChanges);
            changes = LogicTerm::implies(gS[gateStep][Gates::toIndex(Gates::GATES::X)][a], changes);
            lb->assertFormula(changes);

            // Y
            changes =
                    (z[gateStep][a] == z[gateStep - 1][a]);
            changes = changes && (x[gateStep][a] == x[gateStep - 1][a]);

            rChanges = LogicTerm::ite(
                    gS[gateStep][Gates::toIndex(Gates::GATES::Y)][a],
                    rChanges ^ (z[gateStep - 1][a] ^ x[gateStep - 1][a]), rChanges);
            changes = LogicTerm::implies(gS[gateStep][Gates::toIndex(Gates::GATES::Y)][a], changes);
            lb->assertFormula(changes);

            // CNOT
            for (int b = 0; b < nqubits; ++b) {
                if (reducedCM.find({qubitChoice.at(a), qubitChoice.at(b)}) ==
                    reducedCM.end()) {
                    lb->assertFormula(!gC[gateStep][a][b]);
                } else {
                    changes  = LogicTerm(true);
                    rChanges = LogicTerm::ite(
                            gC[gateStep][a][b],
                            (rChanges ^ ((x[gateStep - 1][a] & z[gateStep - 1][b]) &
                                         ((x[gateStep - 1][b] ^ z[gateStep - 1][a]) ^
                                          LogicTerm((1 << nqubits) - 1, nqubits)))),
                            rChanges);
                    changes = changes && (x[gateStep][b] ==
                                          (x[gateStep - 1][b] ^ x[gateStep - 1][a]));
                    changes = changes && (z[gateStep][a] ==
                                          (z[gateStep - 1][a] ^ z[gateStep - 1][b]));
                    changes = changes && (x[gateStep][a] == x[gateStep - 1][a]);
                    changes = changes && (z[gateStep][b] == z[gateStep - 1][b]);
                    changes = LogicTerm::implies(gC[gateStep][a][b], changes);
                    lb->assertFormula(changes);
                }
            }
        }
        lb->assertFormula(r[gateStep] == rChanges);
    }
}


void CliffordSynthesizer::setCircuit(const qc::QuantumComputation& qc) {
    this->circuit = qc.clone();
    Tableau::initTableau(this->initialTableau, this->nqubits);
    Tableau::generateTableau(this->targetTableau, this->circuit);
}

void CliffordSynthesizer::setTargetTableau(Tableau& targetTabl) {
    Tableau::initTableau(this->initialTableau, this->nqubits);
    this->targetTableau = targetTabl;
}
void CliffordSynthesizer::setInitialTableau(Tableau& initialTabl) {
    this->initialTableau = initialTabl;
}


