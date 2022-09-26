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

void CliffordSynthesizer::synthesize() {
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

SynthesisResults CliffordSynthesizer::mainOptimization(
        std::uint32_t                                            timesteps,
        const std::set<std::pair<std::uint16_t, std::uint16_t>>& reducedCM,
        const std::vector<std::uint16_t>& qubitChoice, Tableau& initialTab,
        Tableau& targetTab) {
    std::unique_ptr<logicbase::LogicBlock> lb;
    using namespace logicbase;
    bool success = false;
    if (method == SynthesisMethod::Z3) {
        logicbase::LogicTerm::termType = TermType::BASE;
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
        return SynthesisResults{};
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

    assertTableau(SynthesisData{nqubits, timesteps, reducedCM, qubitChoice, lb, x, z, r, gS, gC}, initialTab, 0);
    assertTableau(SynthesisData{nqubits, timesteps, reducedCM, qubitChoice, lb, x, z, r, gS, gC}, targetTab, timesteps);

    makeSynthesis(SynthesisData{nqubits, timesteps, reducedCM, qubitChoice, lb, x, z, r, gS, gC});

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
    SynthesisResults                          results{};
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
        results.result               = logicbase::Result::SAT;
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
        results.result = logicbase::Result::UNSAT;
        DEBUG() << "UNSAT" << std::endl;
        return results;
    }
    return results;
}

void CliffordSynthesizer::runMinimizer(
        int timesteps, const CouplingMap& reducedCM,
        const std::vector<std::uint16_t>& qubitChoice) {
    DEBUG() << "Running minimizer" << std::endl;
    SynthesisResults r = mainOptimization(timesteps, reducedCM, qubitChoice,
                                          initialTableau, targetTableau);
    updateResults(r);
}
void CliffordSynthesizer::runStartLow(
        int timesteps, const CouplingMap& reducedCM,
        const std::vector<std::uint16_t>& qubitChoice) {
    DEBUG() << "Running start low" << std::endl;
    SynthesisResults r;
    while (r.result != logicbase::Result::SAT || r.result == logicbase::Result::NDEF) {
        DEBUG() << "Current t=" << timesteps << std::endl;
        r = mainOptimization(timesteps, reducedCM, qubitChoice, initialTableau,
                             targetTableau);
        updateResults(r);
        if (r.result == logicbase::Result::UNSAT) {
            timesteps *= 1.5;
        }
    }
}
void CliffordSynthesizer::runStartHigh(
        int timesteps, const CouplingMap& reducedCM,
        const std::vector<std::uint16_t>& qubitChoice) {
    DEBUG() << "Running start high" << std::endl;
    SynthesisResults r;
    int              oldTimesteps = timesteps;
    while (r.result == logicbase::Result::SAT || r.result == logicbase::Result::NDEF) {
        DEBUG() << "Current t=" << timesteps << std::endl;
        r = mainOptimization(timesteps, reducedCM, qubitChoice, initialTableau,
                             targetTableau);
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
        const std::vector<std::uint16_t>& qubitChoice) {
    DEBUG() << "Running minmax" << std::endl;
    SynthesisResults r;
    int              t     = timesteps;
    int              upper = timesteps;
    int              lower = 0;
    while (std::abs(upper - lower) > 1) {
        DEBUG() << "Current t=" << t << std::endl;
        r = mainOptimization(t, reducedCM, qubitChoice, initialTableau,
                             targetTableau);
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
        qc::QuantumComputation& circuit, SynthesisResults* r,
        CliffordSynthesizer* opt) {
    Tableau targetTableau{};
    Tableau::generateTableau(targetTableau, circuit, 0, (i + 1U) * circuitSplit);
    Tableau initTableau{};
    Tableau::generateTableau(initTableau, circuit, 0, i * circuitSplit);
    (*r) = opt->mainOptimization(split, reducedCM, qubitChoice, targetTableau,
                                 initTableau);
};

void CliffordSynthesizer::runSplitIter(
        const CouplingMap&                reducedCM,
        const std::vector<std::uint16_t>& qubitChoice) {
    if (circuit.size() < 2) {
        return;
    }
    DEBUG() << "Running split iter" << std::endl;
    Tableau                        fullTableau  = targetTableau;
    auto                           circuitSplit = static_cast<unsigned int>(std::log(circuit.getNindividualOps()));
    int                            split        = std::min(5, nqubits / 2);
    std::vector<std::thread*>      threads;
    std::vector<SynthesisResults*> results;
    int                            nThreads = nthreads;
    while (true) {
        results.clear();
        DEBUG() << "Current split size: " << split << std::endl;
        DEBUG() << "Current circuit split size: " << circuitSplit << std::endl;
        auto             start = std::chrono::high_resolution_clock::now();
        SynthesisResults totalResult;
        totalResult.result = logicbase::Result::SAT;
        totalResult.resultCircuit.addQubitRegister(nqubits);
        for (size_t i = 0; i * circuitSplit < circuit.getNindividualOps();
             i += nThreads) {
            threads.clear();
            DEBUG() << "Currently at " << i * circuitSplit << " of "
                    << circuit.getNindividualOps() << std::endl;
            for (int j = 0; j < nThreads; j++) {
                auto* r = new SynthesisResults();
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
    optimalResults.result    = logicbase::Result::SAT;
}

void CliffordSynthesizer::updateResults(SynthesisResults& results) {
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

void CliffordSynthesizer::makeSingleGateConstraints(
        const SynthesisData& data) {
    logicbase::LogicTerm changes = logicbase::LogicTerm(true);
    // CONSISTENCY
    // One gate per qubit, per step
    for (unsigned int gateStep = 1; gateStep < data.timesteps + 1; ++gateStep) {
        std::vector<logicbase::LogicTerm> vars{};
        for (unsigned int a = 0; a < data.nqubits; ++a) {
            for (auto gate: Gates::SINGLE_QUBIT) {
                vars.emplace_back(data.gS[gateStep][Gates::toIndex(gate)][a]);
            }
            for (int b = 0; b < data.nqubits; ++b) {
                if (a == b || data.reducedCM.find({data.qubitChoice.at(a), data.qubitChoice.at(b)}) ==
                                      data.reducedCM.end()) {
                    continue;
                }
                vars.emplace_back(data.gC[gateStep][a][b]);
            }
        }
        data.lb->assertFormula(encodings::exactlyOneCmdr(
                encodings::groupVars(vars, static_cast<std::size_t>(vars.size() / 2U)),
                logicbase::LogicTerm::noneTerm(), data.lb.get()));
    }

    // GATE CONSTRAINTS
    for (unsigned int gateStep = 1; gateStep < data.timesteps + 1; ++gateStep) {
        for (unsigned int a = 0; a < data.nqubits; ++a) {
            // NO GATE
            changes = (data.x[gateStep][a] == data.x[gateStep - 1][a]);
            changes = changes && (data.z[gateStep][a] == data.z[gateStep - 1][a]);

            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (data.x[gateStep][b] == data.x[gateStep - 1][b]);
                changes = changes && (data.z[gateStep][b] == data.z[gateStep - 1][b]);
            }

            changes = changes && (data.r[gateStep] == data.r[gateStep - 1]);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][0][a], changes);
            data.lb->assertFormula(changes);

            // H
            changes = (data.z[gateStep][a] == data.x[gateStep - 1][a]);
            changes = changes && (data.x[gateStep][a] == data.z[gateStep - 1][a]);

            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (data.x[gateStep][b] == data.x[gateStep - 1][b]);
                changes = changes && (data.z[gateStep][b] == data.z[gateStep - 1][b]);
            }

            changes = changes &&
                      (data.r[gateStep] == (data.r[gateStep - 1] ^
                                            (data.x[gateStep - 1][a] & data.z[gateStep - 1][a])));
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][1][a], changes);

            data.lb->assertFormula(changes);

            // S
            changes =
                    (data.z[gateStep][a] == (data.z[gateStep - 1][a] ^ data.x[gateStep - 1][a]));
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);

            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (data.x[gateStep][b] == data.x[gateStep - 1][b]);
                changes = changes && (data.z[gateStep][b] == data.z[gateStep - 1][b]);
            }

            changes = changes &&
                      (data.r[gateStep] == (data.r[gateStep - 1] ^
                                            (data.x[gateStep - 1][a] & data.z[gateStep - 1][a])));
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][2][a], changes);
            data.lb->assertFormula(changes);

            // Sdag
            changes =
                    (data.z[gateStep][a] == (data.z[gateStep - 1][a] ^ data.x[gateStep - 1][a]));
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);

            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (data.x[gateStep][b] == data.x[gateStep - 1][b]);
                changes = changes && (data.z[gateStep][b] == data.z[gateStep - 1][b]);
            }

            changes = changes &&
                      (data.r[gateStep] == (data.r[gateStep - 1] ^
                                            (data.x[gateStep - 1][a] & (data.x[gateStep - 1][a] ^ data.z[gateStep - 1][a]))));
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::Sdag)][a], changes);
            data.lb->assertFormula(changes);

            // Z
            changes =
                    (data.z[gateStep][a] == data.z[gateStep - 1][a]);
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);

            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (data.x[gateStep][b] == data.x[gateStep - 1][b]);
                changes = changes && (data.z[gateStep][b] == data.z[gateStep - 1][b]);
            }

            changes = changes &&
                      (data.r[gateStep] == (data.r[gateStep - 1] ^ data.x[gateStep - 1][a]));
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::Z)][a], changes);
            data.lb->assertFormula(changes);

            // X
            changes =
                    (data.z[gateStep][a] == data.z[gateStep - 1][a]);
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);

            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (data.x[gateStep][b] == data.x[gateStep - 1][b]);
                changes = changes && (data.z[gateStep][b] == data.z[gateStep - 1][b]);
            }

            changes = changes &&
                      (data.r[gateStep] == (data.r[gateStep - 1] ^ data.z[gateStep - 1][a]));
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::X)][a], changes);
            data.lb->assertFormula(changes);

            // Y
            changes =
                    (data.z[gateStep][a] == data.z[gateStep - 1][a]);
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);

            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (data.x[gateStep][b] == data.x[gateStep - 1][b]);
                changes = changes && (data.z[gateStep][b] == data.z[gateStep - 1][b]);
            }

            changes = changes &&
                      (data.r[gateStep] == (data.r[gateStep - 1] ^ (data.z[gateStep - 1][a]) ^ data.x[gateStep - 1][a]));
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::Y)][a], changes);
            data.lb->assertFormula(changes);

            // CNOT
            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (data.reducedCM.find({data.qubitChoice.at(a), data.qubitChoice.at(b)}) ==
                    data.reducedCM.end()) {
                    data.lb->assertFormula(!data.gC[gateStep][a][b]);
                } else {
                    changes =
                            (data.r[gateStep] == (data.r[gateStep - 1] ^
                                                  ((data.x[gateStep - 1][a] & data.z[gateStep - 1][b]) &
                                                   ((data.x[gateStep - 1][b] ^ data.z[gateStep - 1][a]) ^
                                                    logicbase::LogicTerm((1 << data.nqubits) - 1, data.nqubits)))));
                    changes = changes && (data.x[gateStep][b] ==
                                          (data.x[gateStep - 1][b] ^ data.x[gateStep - 1][a]));
                    changes = changes && (data.z[gateStep][a] ==
                                          (data.z[gateStep - 1][a] ^ data.z[gateStep - 1][b]));
                    changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);
                    changes = changes && (data.z[gateStep][b] == data.z[gateStep - 1][b]);

                    for (unsigned int c = 0; c < data.nqubits; ++c) { // All other entries do not change
                        if (a == c || b == c) {
                            continue;
                        }
                        changes = changes && (data.x[gateStep][c] == data.x[gateStep - 1][c]);
                        changes = changes && (data.z[gateStep][c] == data.z[gateStep - 1][c]);
                    }

                    changes = logicbase::LogicTerm::implies(data.gC[gateStep][a][b], changes);
                    data.lb->assertFormula(changes);
                }
            }
        }
    }
}
void CliffordSynthesizer::makeMultipleGateConstraints(
        const SynthesisData& data) {
    logicbase::LogicTerm changes = logicbase::LogicTerm(true);
    // CONSISTENCY
    // One gate per qubit, per step
    for (unsigned int gateStep = 1; gateStep < data.timesteps + 1; ++gateStep) {
        for (unsigned int a = 0; a < data.nqubits; ++a) {
            std::vector<logicbase::LogicTerm> vars{};
            for (auto gate: Gates::SINGLE_QUBIT) {
                vars.emplace_back(data.gS[gateStep][Gates::toIndex(gate)][a]);
            }
            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                vars.emplace_back(data.gC[gateStep][a][b]);
                vars.emplace_back(data.gC[gateStep][b][a]);
            }
            data.lb->assertFormula(encodings::exactlyOneCmdr(
                    encodings::groupVars(vars, static_cast<std::size_t>(vars.size() / 2)),
                    logicbase::LogicTerm::noneTerm(), data.lb.get()));
        }
    }
    // Maximum any combination of 1 and 2 qubit gates adding up to n
    for (unsigned int gateStep = 1; gateStep < data.timesteps + 1; ++gateStep) {
        changes = logicbase::LogicTerm(0);
        for (unsigned int a = 0; a < data.nqubits; ++a) {
            for (auto gate: Gates::SINGLE_QUBIT) {
                changes = changes + data.gS[gateStep][Gates::toIndex(gate)][a];
            }
            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes + data.gC[gateStep][a][b] + data.gC[gateStep][a][b];
            }
        }
        changes = changes < logicbase::LogicTerm(static_cast<int>(data.nqubits + 1));
        data.lb->assertFormula(changes);
    }

    // GATE CONSTRAINTS
    for (unsigned int gateStep = 1; gateStep < data.timesteps + 1; ++gateStep) {
        logicbase::LogicTerm rChanges = data.r[gateStep - 1];
        for (unsigned int a = 0; a < data.nqubits; ++a) {
            // NO GATE
            changes = logicbase::LogicTerm(true);
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);
            changes = changes && (data.z[gateStep][a] == data.z[gateStep - 1][a]);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][0][a], changes);
            data.lb->assertFormula(changes);

            // H
            changes = logicbase::LogicTerm(true);
            changes = changes && (data.z[gateStep][a] == data.x[gateStep - 1][a]);
            changes = changes && (data.x[gateStep][a] == data.z[gateStep - 1][a]);

            rChanges = logicbase::LogicTerm::ite(
                    data.gS[gateStep][1][a],
                    rChanges ^ (data.x[gateStep - 1][a] & data.z[gateStep - 1][a]), rChanges);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][1][a], changes);

            data.lb->assertFormula(changes);

            // S
            changes = logicbase::LogicTerm(true);
            changes = changes && (data.z[gateStep][a] ==
                                  (data.z[gateStep - 1][a] ^ data.x[gateStep - 1][a]));
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);

            rChanges = logicbase::LogicTerm::ite(
                    data.gS[gateStep][2][a],
                    rChanges ^ (data.x[gateStep - 1][a] & data.z[gateStep - 1][a]), rChanges);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][2][a], changes);
            data.lb->assertFormula(changes);

            // Sdag
            changes =
                    (data.z[gateStep][a] == (data.z[gateStep - 1][a] ^ data.x[gateStep - 1][a]));
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);

            rChanges = logicbase::LogicTerm::ite(
                    data.gS[gateStep][Gates::toIndex(Gates::GATES::Sdag)][a],
                    rChanges ^ (data.x[gateStep - 1][a] & (data.x[gateStep - 1][a] ^ data.z[gateStep - 1][a])), rChanges);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::Sdag)][a], changes);
            data.lb->assertFormula(changes);

            // Z
            changes =
                    (data.z[gateStep][a] == data.z[gateStep - 1][a]);
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);

            rChanges = logicbase::LogicTerm::ite(
                    data.gS[gateStep][Gates::toIndex(Gates::GATES::Z)][a],
                    rChanges ^ (data.x[gateStep - 1][a]), rChanges);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::Z)][a], changes);
            data.lb->assertFormula(changes);

            // X
            changes =
                    (data.z[gateStep][a] == data.z[gateStep - 1][a]);
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);

            rChanges = logicbase::LogicTerm::ite(
                    data.gS[gateStep][Gates::toIndex(Gates::GATES::X)][a],
                    rChanges ^ (data.z[gateStep - 1][a]), rChanges);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::X)][a], changes);
            data.lb->assertFormula(changes);

            // Y
            changes =
                    (data.z[gateStep][a] == data.z[gateStep - 1][a]);
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);

            rChanges = logicbase::LogicTerm::ite(
                    data.gS[gateStep][Gates::toIndex(Gates::GATES::Y)][a],
                    rChanges ^ (data.z[gateStep - 1][a] ^ data.x[gateStep - 1][a]), rChanges);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::Y)][a], changes);
            data.lb->assertFormula(changes);

            // CNOT
            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (data.reducedCM.find({data.qubitChoice.at(a), data.qubitChoice.at(b)}) ==
                    data.reducedCM.end()) {
                    data.lb->assertFormula(!data.gC[gateStep][a][b]);
                } else {
                    changes  = logicbase::LogicTerm(true);
                    rChanges = logicbase::LogicTerm::ite(
                            data.gC[gateStep][a][b],
                            (rChanges ^ ((data.x[gateStep - 1][a] & data.z[gateStep - 1][b]) &
                                         ((data.x[gateStep - 1][b] ^ data.z[gateStep - 1][a]) ^
                                          logicbase::LogicTerm((1 << data.nqubits) - 1, data.nqubits)))),
                            rChanges);
                    changes = changes && (data.x[gateStep][b] ==
                                          (data.x[gateStep - 1][b] ^ data.x[gateStep - 1][a]));
                    changes = changes && (data.z[gateStep][a] ==
                                          (data.z[gateStep - 1][a] ^ data.z[gateStep - 1][b]));
                    changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);
                    changes = changes && (data.z[gateStep][b] == data.z[gateStep - 1][b]);
                    changes = logicbase::LogicTerm::implies(data.gC[gateStep][a][b], changes);
                    data.lb->assertFormula(changes);
                }
            }
        }
        data.lb->assertFormula(data.r[gateStep] == rChanges);
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