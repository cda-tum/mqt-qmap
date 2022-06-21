/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
*/

#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#include "LogicTerm/LogicTerm.hpp"
#include "operations/OpType.hpp"
#include "operations/StandardOperation.hpp"
#include "utils.hpp"
#include "utils/logging.hpp"

void CliffordOptimizer::optimize() {
    TRACE() << "Strategy: " << toString(strategy) << std::endl;
    TRACE() << "Target: " << toString(target) << std::endl;
    TRACE() << "Method: " << toString(method) << std::endl;

    auto                     total_start = std::chrono::high_resolution_clock::now();
    std::vector<CouplingMap> reducedMaps;
    architecture.getReducedCouplingMaps(nqubits, reducedMaps);
    auto subsets =
            (choose_best ? highestFidelityMap : reducedMaps);
    for (const auto& subset: subsets) {
        std::vector<unsigned short> qubitMap{Architecture::getQubitList(subset)};

        DEBUG() << "Reduced Coupling Map" << (choose_best ? " (best)" : "") << ": ";
        std::stringstream strings;
        Architecture::printCouplingMap(subset, strings);
        DEBUG() << strings.str();
        DEBUG() << "Qubit Map: " << qubitMap;
        DEBUG() << "Coupling Map Fidelity: "
                << architecture.getAverageArchitectureFidelity(architecture.getCouplingMap(),
                                                               std::set<unsigned short>(qubitMap.begin(), qubitMap.end()),
                                                               architecture.getCalibrationData());
        int timesteps =
                initial_timesteps == 0 ? nqubits * nqubits : initial_timesteps;
        if (strategy == OptimizingStrategy::UseMinimizer) {
            runMinimizer(timesteps, subset, qubitMap);
        }
        if (strategy == OptimizingStrategy::StartLow) {
            runStartLow(timesteps, subset, qubitMap);
        }
        if (strategy == OptimizingStrategy::StartHigh) {
            runStartHigh(timesteps, subset, qubitMap);
        }
        if (strategy == OptimizingStrategy::MinMax) {
            runMinMax(timesteps, subset, qubitMap);
        }
        if (strategy == OptimizingStrategy::SplitIter) {
            runSplitIter(subset, qubitMap);
        }
        if (choose_best && optimal_results.sat)
            break;
    }

    auto                          total_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff      = total_end - total_start;
    INFO() << "Time: " << diff.count() << std::endl;
    optimal_results.total_seconds = diff.count();
}

CliffordOptResults CliffordOptimizer::main_optimization(
        int                                                        timesteps,
        const std::set<std::pair<unsigned short, unsigned short>>& reducedCM,
        const std::vector<unsigned short>& qubitChoice, Tableau& initialTab,
        Tableau& targetTab) {
    LogicBlock* lb;
#ifdef Z3_FOUND
    using namespace z3logic;
    z3::context  c;
    z3::solver   slv(c);
    z3::optimize opt(c);
    z3::params   p(c);
    if (method == OptMethod::Z3) {
        LogicTerm::termType = TermType::BASE;
        if (strategy == OptimizingStrategy::UseMinimizer || strategy == OptimizingStrategy::SplitIter) {
            p.set("pb.compile_equality", true);
            p.set("maxres.hill_climb", true);
            p.set("maxres.pivot_on_correction_set", false);
            // z3::set_param("parallel.enable", true);
            // z3::set_param("parallel.threads.max", nthreads);
            opt.set(p);
            lb = new Z3LogicOptimizer(c, opt, false);
        } else {
            p.set("threads", unsigned(nthreads / 2));
            z3::set_param("parallel.enable", true);
            z3::set_param("parallel.threads.max", nthreads / 2);
            slv.set(p);
            lb = new Z3LogicBlock(c, slv, true);
        }
    } else {
        return CliffordOptResults{};
    }
#endif
    DEBUG() << "lb 1: " << lb << std::endl;
    LogicMatrix   x{};
    LogicMatrix   z{};
    LogicVector   r{};
    LogicMatrix3D g_s{};
    LogicMatrix3D g_c{};

    LogicTerm changes = LogicTerm(true);

    auto start = std::chrono::high_resolution_clock::now();
    /*
   * Tableau Variables x/z
   * k before gate k is applied
   * i column
   */
    std::stringstream x_name{};
    std::stringstream z_name{};
    std::stringstream r_name{};
    for (int k = 0; k < timesteps + 1; ++k) {
        x.emplace_back();
        z.emplace_back();
        for (int i = 0; i < nqubits; ++i) {
            x_name.str("");
            z_name.str("");
            x_name << "x_" << k << "_" << i;
            z_name << "z_" << k << "_" << i;
            x.back().push_back(
                    lb->makeVariable(x_name.str(), CType::BITVECTOR, nqubits));
            z.back().push_back(
                    lb->makeVariable(z_name.str(), CType::BITVECTOR, nqubits));
        }
        r_name.str("");
        r_name << "r_" << k;
        r.push_back(lb->makeVariable(r_name.str(), CType::BITVECTOR, nqubits));
    }

    /*
   * Gate Variables
   * k before gates k are applied
   * i qubit 1
   * j qubit
   */
    std::stringstream g_name{};
    for (int gate_step = 0; gate_step < timesteps + 1; ++gate_step) {
        g_s.emplace_back();
        for (auto gate: Gates::singleQubit) {
            g_s.back().emplace_back();
            for (int j = 0; j < nqubits; ++j) {
                g_name.str("");
                g_name << "g_" << gate_step << "_" << Gates::gateName(gate) << "_" << j;
                g_s.back().back().push_back(lb->makeVariable(g_name.str()));
            }
        }
    }
    for (int gate_step = 0; gate_step < timesteps + 1; ++gate_step) {
        g_c.emplace_back();
        for (int j = 0; j < nqubits; ++j) {
            g_c.back().emplace_back();
            for (int l = 0; l < nqubits; ++l) {
                g_name.str("");
                g_name << "g_" << gate_step << "_CNOT_" << j << "_" << l;
                g_c.back().back().push_back(lb->makeVariable(g_name.str()));
            }
        }
    }

    assertTableau(initialTab, lb, x, z, r, nqubits, 0);
    assertTableau(targetTab, lb, x, z, r, nqubits, timesteps);

    if (target == OptTarget::DEPTH) {
        make_depth_optimizer(timesteps, reducedCM, qubitChoice, lb, x, z, r, g_s,
                             g_c);
    } else if (target == OptTarget::GATES) {
        make_gate_optimizer(timesteps, reducedCM, qubitChoice, lb, x, z, r, g_s,
                            g_c);
    } else if (target == OptTarget::FIDELITY) {
        make_fidelity_optimizer(timesteps, reducedCM, qubitChoice, lb, x, z, r, g_s,
                                g_c);
    } else {
        ERROR() << "Unknown target" << std::endl;
        return CliffordOptResults{};
    }
    auto                          formulation = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff        = formulation - start;
    INFO() << "Time to prduce Formulation: " << diff.count() << std::endl;

    lb->produceInstance();
    // lb->dumpAll(std::cout);

    auto mod_gen = std::chrono::high_resolution_clock::now();
    diff         = mod_gen - formulation;
    INFO() << "Time to generate Model: " << diff.count() << std::endl;

    TRACE() << "Clauses: " << TermImpl::getNextId(lb) << std::endl;
    TRACE() << "None Terms: " << TermImpl::getNextId() << std::endl;

    Result result = lb->solve();

    auto end = std::chrono::high_resolution_clock::now();
    diff     = end - mod_gen;
    INFO() << "Time to solve Model: " << diff.count() << std::endl;
    std::chrono::duration<double, std::milli> elapsed_milliseconds = end - start;
    CliffordOptResults                        results{};
    results.verbose           = verbose;
    results.choose_best       = choose_best;
    results.nqubits           = nqubits;
    results.initial_timesteps = timesteps;
    results.strategy          = strategy;
    results.target            = target;
    results.total_seconds     = elapsed_milliseconds.count();
    results.sat               = result == Result::SAT;
    results.doubleFidelity    = architecture.getCNOTFidelities();
    results.singleFidelity    = architecture.getSingleQubitFidelities();
    results.resultCM          = architecture.getCouplingMap();
    results.resultTableaus.clear();

    if (result == Result::SAT) {
        results.result               = OptResult::SAT;
        Model*                 model = lb->getModel();
        qc::QuantumComputation circuit;
        circuit.addQubitRegister(nqubits);
        results.gate_count = 0;
        results.depth      = 0;
        results.fidelity   = 1;
        for (int gate_step = 0; gate_step < timesteps + 1; ++gate_step) {
            int old_gate_count = results.gate_count;
            TRACE() << "Gate Step: " << gate_step << std::endl
                    << " Actual gate count: " << results.gate_count << std::endl
                    << " Depth: " << results.depth << std::endl
                    << " Fidelity: " << results.fidelity << std::endl;
            if (gate_step > 0) {
                for (int a = 0; a < nqubits; ++a) {
                    for (auto gate: Gates::singleQubitWithoutNOP) {
                        if (model->getBoolValue(g_s[gate_step][Gates::toIndex(gate)][a], lb)) {
                            circuit.emplace_back<qc::StandardOperation>(nqubits, a, Gates::toOpType(gate));
                            if (architecture.isCalibrationDataAvailable())
                                results.fidelity *= (architecture.getSingleQubitFidelities()[a]);
                            TRACE() << Gates::gateName(gate) << "(" << a << ")"
                                    << ") Fidelity: " << architecture.getSingleQubitFidelities()[a]
                                    << std::endl;
                            ++results.gate_count;
                        }
                    }
                    for (int b = 0; b < nqubits; ++b) {
                        if (model->getBoolValue(g_c[gate_step][a][b], lb)) {
                            results.gate_count++;
                            circuit.emplace_back<qc::StandardOperation>(
                                    nqubits, dd::Control{static_cast<dd::Qubit>(a)}, b, qc::X);
                            if (architecture.isCalibrationDataAvailable())
                                results.fidelity *=
                                        (architecture.getCNOTFidelities()[qubitChoice.at(a)]
                                                                         [qubitChoice.at(b)]);
                            TRACE() << "X(" << a << "," << b << ") Fidelity: "
                                    << architecture.getCNOTFidelities()[qubitChoice.at(a)]
                                                                       [qubitChoice.at(b)]
                                    << std::endl;
                        }
                    }
                }
            }
            if (old_gate_count < results.gate_count) {
                results.depth++;
            }
            auto tableau = results.resultTableaus.emplace_back();
            generateTableau(tableau, circuit);
            initTableau(modelTableau);
            for (int i = 0; i < nqubits; ++i) {
                modelTableau.populateTableauFrom(model->getBitvectorValue(x[gate_step][i], lb),
                                                 nqubits, i);
                modelTableau.populateTableauFrom(model->getBitvectorValue(z[gate_step][i], lb),
                                                 nqubits, i + nqubits);
            }
            modelTableau.populateTableauFrom(model->getBitvectorValue(r[gate_step], lb),
                                             nqubits, 2 * nqubits);
            if (verbose >= 5) {
                TRACE() << modelTableau;
            }
        }
        results.resultCircuit = circuit.clone();
    }
    lb->reset();
    delete lb;
    if (result == Result::SAT) {
        DEBUG() << "SAT" << std::endl;
        return results;
    } else {
        results.result = OptResult::UNSAT;
        DEBUG() << "UNSAT" << std::endl;
        return results;
    }
    return results;
}

void CliffordOptimizer::make_depth_optimizer(
        int                                                        timesteps,
        const std::set<std::pair<unsigned short, unsigned short>>& reducedCM,
        const std::vector<unsigned short>& qubitChoice, LogicBlock* lb,
        const LogicMatrix& x, const LogicMatrix& z, const LogicVector& r,
        const LogicMatrix3D& g_s, const LogicMatrix3D& g_c) {
    makeMultipleGateConstraints(lb, x, z, r, nqubits, timesteps, reducedCM,
                                qubitChoice, g_s, g_c);
    // COST
    if (strategy == OptimizingStrategy::UseMinimizer ||
        strategy == OptimizingStrategy::SplitIter) {
        LogicTerm cost = LogicTerm(0);
        for (int gate_step = 1; gate_step < timesteps + 1; ++gate_step) {
            LogicTerm anyGate = LogicTerm(true);
            for (int a = 0; a < nqubits; ++a) {
                for (auto gate: Gates::singleQubitWithoutNOP) {
                    anyGate = anyGate && !g_s[gate_step][Gates::toIndex(gate)][a];
                }
                for (int b = 0; b <= a; ++b) {
                    if (a == b)
                        continue;
                    anyGate = anyGate && !g_c[gate_step][a][b] && !g_c[gate_step][b][a];
                }
            }
            cost = cost + LogicTerm::ite(anyGate, LogicTerm(5), LogicTerm(0));
        }
        dynamic_cast<LogicBlockOptimizer*>(lb)->maximize(cost);
    }
}

void CliffordOptimizer::make_fidelity_optimizer(
        int                                                        timesteps,
        const std::set<std::pair<unsigned short, unsigned short>>& reducedCM,
        const std::vector<unsigned short>& qubitChoice, LogicBlock* lb,
        const LogicMatrix& x, const LogicMatrix& z, const LogicVector& r,
        const LogicMatrix3D& g_s, const LogicMatrix3D& g_c) {
    if (!architecture.isArchitectureAvailable()) {
        util::fatal("No fidelity architecture specified in coupling map.");
    }
    makeMultipleGateConstraints(lb, x, z, r, nqubits, timesteps, reducedCM,
                                qubitChoice, g_s, g_c);
    // COST
    if (strategy == OptimizingStrategy::UseMinimizer ||
        strategy == OptimizingStrategy::SplitIter) {
        LogicTerm cost = LogicTerm(0);
        // For each edge in the coupling map, get the fidelity cost
        for (const auto& edge: reducedCM) {
            LogicTerm fidelity = LogicTerm(
                    architecture.getLogCNOTFidelities()[edge.first][edge.second] * 1000);
            // at each time t if there is a gate on the edge, add the cost
            for (int gate_step = 0; gate_step < timesteps; ++gate_step) {
                cost = cost + (g_c[gate_step][edge.first][edge.second] * fidelity);
            }
        }
        // For each qubit, get the fidelity cost
        for (int a = 0; a < nqubits; ++a) {
            LogicTerm fidelity =
                    LogicTerm(architecture.getLogSingleQubitFidelities()[a] * 1000);
            // at each time t if there is a gate on a, add the cost
            for (int gate_step = 0; gate_step < timesteps; ++gate_step) {
                for (auto gate: Gates::singleQubitWithoutNOP) {
                    cost = cost + (g_s[gate_step][Gates::toIndex(gate)][a] * fidelity);
                }
            }
        }
        dynamic_cast<LogicBlockOptimizer*>(lb)->minimize(cost);
        cost = LogicTerm(0);
        for (int gate_step = 0; gate_step < timesteps; ++gate_step) {
            for (int a = 0; a < nqubits; ++a) {
                cost = cost + g_s[gate_step][1][a] + g_s[gate_step][2][a];
                for (int b = 0; b < nqubits; ++b) {
                    cost = cost + g_c[gate_step][a][b];
                }
            }
        }
        dynamic_cast<LogicBlockOptimizer*>(lb)->maximize(cost);
    }
}

void CliffordOptimizer::make_gate_optimizer(
        int                                                        timesteps,
        const std::set<std::pair<unsigned short, unsigned short>>& reducedCM,
        const std::vector<unsigned short>& qubitChoice, LogicBlock* lb,
        const LogicMatrix& x, const LogicMatrix& z, const LogicVector& r,
        const LogicMatrix3D& g_s, const LogicMatrix3D& g_c) {
    LogicTerm changes = LogicTerm(true);
    makeSingleGateConstraints(lb, x, z, r, nqubits, timesteps, reducedCM,
                              qubitChoice, g_s, g_c);
    // COST
    if (strategy == OptimizingStrategy::UseMinimizer) {
        LogicTerm cost = LogicTerm(0);
        for (int gate_step = 1; gate_step < timesteps + 1; ++gate_step) {
            for (int a = 0; a < nqubits; ++a) {
                for (auto gate: Gates::singleQubitWithoutNOP) {
                    cost = cost + g_s[gate_step][Gates::toIndex(gate)][a];
                }
                for (int b = 0; b <= a; ++b) {
                    if (a == b)
                        continue;
                    cost = cost + g_c[gate_step][a][b] + g_c[gate_step][b][a];
                }
            }
        }
        dynamic_cast<LogicBlockOptimizer*>(lb)->minimize(cost);
    }
}

void CliffordOptimizer::calculateQubitsUsed(qc::QuantumComputation& circ, std::set<signed char>& qubits) {
    for (auto& gate: circ) {
        if (gate->getType() == qc::OpType::Compound) {
            auto compOp = dynamic_cast<qc::CompoundOperation*>(gate.get());
            auto cit    = compOp->begin();
            while (cit != compOp->end()) {
                getGateQubits(*cit, qubits);
                ++cit;
            }
        } else {
            getGateQubits(gate, qubits);
        }
    }
}

int CliffordOptimizer::applyGate(std::unique_ptr<qc::Operation>& gate, Tableau& tableau) const {
    switch (gate->getType()) {
        case qc::OpType::H: // HADAMARD
        {
            if (gate->isControlled()) {
                util::fatal("Expected single-qubit gate");
            }
            const auto a = gate->getTargets().at(0U);
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= (tableau[i][a] & tableau[i][a + nqubits]);
                std::swap(tableau[i][a], tableau[i][a + nqubits]);
            }
            return 1U;
        } break;
        case qc::OpType::S: // PHASE
        {
            if (gate->isControlled()) {
                util::fatal("Expected single-qubit gate");
            }
            const auto a = gate->getTargets().at(0U);
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= (tableau[i][a] & tableau[i][a + nqubits]);
                tableau[i][a + nqubits] ^= tableau[i][a];
            }
            return 1U;
        } break;
        case qc::OpType::X: // CNOT
        {
            if (gate->getNcontrols() != 1U) { // NOT = H x S x S x H
                const auto a = gate->getTargets().at(0U);
                for (auto i = 0U; i < nqubits; i++) {
                    tableau[i][2U * nqubits] ^= (tableau[i][a] & tableau[i][a + nqubits]);
                    std::swap(tableau[i][a], tableau[i][a + nqubits]);
                }
                for (auto i = 0U; i < nqubits; i++) {
                    tableau[i][2U * nqubits] ^= (tableau[i][a] & tableau[i][a + nqubits]);
                    tableau[i][a + nqubits] ^= tableau[i][a];
                }
                for (auto i = 0U; i < nqubits; i++) {
                    tableau[i][2U * nqubits] ^= (tableau[i][a] & tableau[i][a + nqubits]);
                    tableau[i][a + nqubits] ^= tableau[i][a];
                }
                for (auto i = 0U; i < nqubits; i++) {
                    tableau[i][2U * nqubits] ^= (tableau[i][a] & tableau[i][a + nqubits]);
                    std::swap(tableau[i][a], tableau[i][a + nqubits]);
                }
                return 4U;
            } else {
                const auto a = (*gate->getControls().begin()).qubit;
                const auto b = gate->getTargets().at(0);
                if (a == b) {
                    util::fatal("Invalid CNOT with same control and target.");
                }
                for (auto i = 0U; i < nqubits; i++) {
                    const auto xa = tableau[i][a];
                    const auto za = tableau[i][a + nqubits];
                    const auto xb = tableau[i][b];
                    const auto zb = tableau[i][b + nqubits];
                    tableau[i][2 * nqubits] ^= ((xa & zb) & ((xb ^ za) ^ 1));
                    tableau[i][a + nqubits] = za ^ zb;
                    tableau[i][b]           = xb ^ xa;
                }
                return 1U;
            }
        } break;

        case qc::OpType::Sdag: { // Sdag  = S x S x S
            if (gate->isControlled()) {
                util::fatal("Expected single-qubit gate");
            }
            const auto a = gate->getTargets().at(0U);
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= (tableau[i][a] & tableau[i][a + nqubits]);
                tableau[i][a + nqubits] ^= tableau[i][a];
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= (tableau[i][a] & tableau[i][a + nqubits]);
                tableau[i][a + nqubits] ^= tableau[i][a];
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= (tableau[i][a] & tableau[i][a + nqubits]);
                tableau[i][a + nqubits] ^= tableau[i][a];
            }
            return 3U;
        } break;
        case qc::OpType::Z: { // Z = S x S
            if (gate->isControlled()) {
                util::fatal("Expected single-qubit gate");
            }
            const auto a = gate->getTargets().at(0U);
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= (tableau[i][a] & tableau[i][a + nqubits]);
                tableau[i][a + nqubits] ^= tableau[i][a];
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= (tableau[i][a] & tableau[i][a + nqubits]);
                tableau[i][a + nqubits] ^= tableau[i][a];
            }
            return 2U;
        } break;
        case qc::OpType::Y: { // Y = H x S x S x H x S x S
            if (gate->isControlled()) {
                util::fatal("Expected single-qubit gate");
            }
            const auto a = gate->getTargets().at(0U);
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= (tableau[i][a] & tableau[i][a + nqubits]);
                std::swap(tableau[i][a], tableau[i][a + nqubits]);
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= (tableau[i][a] & tableau[i][a + nqubits]);
                tableau[i][a + nqubits] ^= tableau[i][a];
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= (tableau[i][a] & tableau[i][a + nqubits]);
                tableau[i][a + nqubits] ^= tableau[i][a];
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= (tableau[i][a] & tableau[i][a + nqubits]);
                std::swap(tableau[i][a], tableau[i][a + nqubits]);
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= (tableau[i][a] & tableau[i][a + nqubits]);
                tableau[i][a + nqubits] ^= tableau[i][a];
            }
            for (auto i = 0U; i < nqubits; i++) {
                tableau[i][2U * nqubits] ^= (tableau[i][a] & tableau[i][a + nqubits]);
                tableau[i][a + nqubits] ^= tableau[i][a];
            }
            return 6U;
        } break;
        default:
            util::fatal("Unsupported gate encountered: " + std::to_string(gate->getType()));
            break;
    }
    return 0U;
}

void CliffordOptimizer::generateTableau(Tableau& tableau, qc::QuantumComputation& circuit, int begin, int end) {
    initTableau(tableau);
    int current_g = 0;
    for (auto& gate: circuit) {
        if (current_g >= begin && (current_g < end || end < 0)) {
            if (gate->getType() == qc::OpType::Compound) {
                auto compOp = dynamic_cast<qc::CompoundOperation*>(gate.get());
                auto cit    = compOp->begin();
                while (cit != compOp->end() && current_g >= begin &&
                       (current_g < end || end < 0)) {
                    applyGate((*cit), tableau);
                    ++cit;
                    ++current_g;
                }
            } else {
                applyGate(gate, tableau);
                ++current_g;
            }
        }
    }
}

void CliffordOptimizer::initTableau(Tableau& tableau) const {
    // tableau.clear();
    tableau.init(nqubits);
}

void CliffordOptimizer::runMinimizer(
        int timesteps, const CouplingMap& reducedCM,
        const std::vector<unsigned short>& qubitChoice) {
    DEBUG() << "Running minimizer" << std::endl;
    CliffordOptResults r = main_optimization(timesteps, reducedCM, qubitChoice,
                                             initialTableau, targetTableau);
    updateResults(r);
}
void CliffordOptimizer::runStartLow(
        int timesteps, const CouplingMap& reducedCM,
        const std::vector<unsigned short>& qubitChoice) {
    DEBUG() << "Running start low" << std::endl;
    CliffordOptResults r;
    while (r.result != OptResult::SAT) {
        DEBUG() << "Current t=" << timesteps << std::endl;
        r = main_optimization(timesteps, reducedCM, qubitChoice, initialTableau,
                              targetTableau);
        updateResults(r);
        if (r.result == OptResult::UNSAT) {
            timesteps *= 1.5;
        }
    }
}
void CliffordOptimizer::runStartHigh(
        int timesteps, const CouplingMap& reducedCM,
        const std::vector<unsigned short>& qubitChoice) {
    DEBUG() << "Running start high" << std::endl;
    CliffordOptResults r;
    int                old_timesteps = timesteps;
    while (r.result == OptResult::SAT) {
        DEBUG() << "Current t=" << timesteps << std::endl;
        r = main_optimization(timesteps, reducedCM, qubitChoice, initialTableau,
                              targetTableau);
        updateResults(r);
        if (r.result == OptResult::SAT) {
            old_timesteps = timesteps;
            timesteps *= 0.5;
        } else {
            timesteps = old_timesteps;
        }
    }
}
void CliffordOptimizer::runMinMax(
        int timesteps, const CouplingMap& reducedCM,
        const std::vector<unsigned short>& qubitChoice) {
    DEBUG() << "Running minmax" << std::endl;
    CliffordOptResults r;
    int                t     = timesteps;
    int                upper = timesteps, lower = 0;
    while (std::abs(upper - lower) > 1) {
        DEBUG() << "Current t=" << t << std::endl;
        r = main_optimization(t, reducedCM, qubitChoice, initialTableau,
                              targetTableau);
        updateResults(r);
        if (r.result == OptResult::SAT) {
            upper = t;
        } else if (r.result == OptResult::UNSAT) {
            lower = t;
        } else {
            break;
        }
        if (upper - lower < 1 && r.result == OptResult::UNSAT) {
            upper *= 1.5;
        }
        t = lower + std::abs(upper - lower) / 2;
    }
}

void CliffordOptimizer::runSplinter(
        int i, unsigned int circuit_split, unsigned int split,
        const CouplingMap& reducedCM, const std::vector<unsigned short>& qubitChoice,
        qc::QuantumComputation& circuit, CliffordOptResults* r,
        CliffordOptimizer* opt) {
    Tableau targetTableau{};
    opt->generateTableau(targetTableau, circuit, 0, (i + 1U) * circuit_split);
    Tableau initTableau{};
    opt->generateTableau(initTableau, circuit, 0, i * circuit_split);
    (*r) = opt->main_optimization(split, reducedCM, qubitChoice, targetTableau,
                                  initTableau);
};

void CliffordOptimizer::runSplitIter(
        const CouplingMap&                 reducedCM,
        const std::vector<unsigned short>& qubitChoice) {
    DEBUG() << "Running split iter" << std::endl;
    Tableau                          fullTableau   = targetTableau;
    auto                             circuit_split = static_cast<unsigned int>(std::log(circuit.getNindividualOps()));
    int                              split         = std::max(5, nqubits / 2);
    std::vector<std::thread*>        threads;
    std::vector<CliffordOptResults*> results;
    int                              nThreads = nthreads;
    while (true) {
        results.clear();
        DEBUG() << "Current split size: " << split << std::endl;
        DEBUG() << "Current circuit split size: " << circuit_split << std::endl;
        auto               start = std::chrono::high_resolution_clock::now();
        CliffordOptResults total_result;
        total_result.result = OptResult::SAT;
        total_result.resultCircuit.addQubitRegister(nqubits);
        for (int i = 0; i * circuit_split < circuit.getNindividualOps();
             i += nThreads) {
            threads.clear();
            DEBUG() << "Currently at " << i * circuit_split << " of "
                    << circuit.getNindividualOps() << std::endl;
            for (int j = 0; j < nThreads; j++) {
                auto r = new CliffordOptResults();
                auto t = new std::thread(CliffordOptimizer::runSplinter, i, circuit_split,
                                         split, std::ref(reducedCM), std::ref(qubitChoice),
                                         std::ref(circuit), r, this);
                threads.push_back(t);
                results.push_back(r);
            }
            for (auto t: threads) {
                t->join();
            }
            for (auto t: threads) {
                delete t;
            }
            for (auto r: results) {
                if (r->result == OptResult::UNSAT) {
                    total_result.result = OptResult::UNSAT;
                    break;
                }
            }
            if (total_result.result == OptResult::UNSAT) {
                DEBUG() << "UNSAT, increasing split size." << std::endl;
                split += std::max(1.0, split * 0.2);
                break;
            }
        }
        for (auto r: results) {
            for (const auto& gate: r->resultCircuit) {
                total_result.resultCircuit.insert(total_result.resultCircuit.end(),
                                                  gate->clone());
            }
            delete r;
        }
        if (total_result.result == OptResult::SAT) {
            // TRACE() << "Full Tableau" << std::endl;
            // TRACE() << util::pretty_s(fullTableau);
            // TRACE() << "Target Tableau" << std::endl;
            // TRACE() << util::pretty_s(targetTableau);
            // TRACE() << "Result Tableau" << std::endl;
            // TRACE() << util::pretty_s(r.resultTableaus.back());
            Tableau resultingTableau{};
            generateTableau(resultingTableau, total_result.resultCircuit);
            DEBUG() << "Equality (Results): "
                    << (fullTableau == resultingTableau)
                    << std::endl;
            DEBUG() << "Original Circuit size: " << circuit.getNindividualOps()
                    << std::endl;
            DEBUG() << "Optimized Circuit size: "
                    << total_result.resultCircuit.getNindividualOps() << std::endl;
            TRACE() << "Resulting Circuit: " << std::endl;
            std::ostringstream ss;
            total_result.resultCircuit.dump(ss, qc::Format::OpenQASM);
            TRACE() << ss.str() << std::endl;
            // split *= 1.5;
            // circ_split *= 1.5;
            auto                          end  = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> diff = end - start;
            INFO() << "Time for complete run: " << diff.count() << std::endl;
            if (circuit.getNindividualOps() ==
                total_result.resultCircuit.getNindividualOps()) {
                split *= 1.2;
                break;
            }
            circuit = total_result.resultCircuit.clone();
        }
    }
    optimal_results.resultCircuit = circuit.clone();
    optimal_results.resultTableaus.emplace_back(targetTableau);
    optimal_results.gate_count = circuit.getNindividualOps();
}

void CliffordOptimizer::updateResults(CliffordOptResults& results) {
    if (!results.sat)
        return;
    switch (target) {
        case OptTarget::GATES:
            if (results.gate_count < optimal_results.gate_count ||
                optimal_results.gate_count == 0) {
                optimal_results = results;
            }
            break;
        case OptTarget::DEPTH:
            if (results.depth < optimal_results.depth || optimal_results.depth == 0) {
                optimal_results = results;
            }
            break;
        case OptTarget::FIDELITY:
            if (results.fidelity >= optimal_results.fidelity ||
                optimal_results.fidelity == 0) {
                optimal_results = results;
            }
            break;
        default:
            break;
    }
}

void CliffordOptimizer::assertTableau(const Tableau& tableau, LogicBlock* lb,
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

void CliffordOptimizer::makeSingleGateConstraints(
        LogicBlock* lb, const LogicMatrix& x, const LogicMatrix& z,
        const LogicVector& r, int nqubits, int timesteps,
        const std::set<std::pair<unsigned short, unsigned short>>& reducedCM,
        const std::vector<unsigned short>& qubitChoice, const LogicMatrix3D& g_s,
        const LogicMatrix3D& g_c) {
    LogicTerm changes = LogicTerm(true);
    // CONSISTENCY
    // One gate per qubit, per step
    for (int gate_step = 1; gate_step < timesteps + 1; ++gate_step) {
        std::vector<LogicTerm> vars{};
        for (int a = 0; a < nqubits; ++a) {
            for (auto gate: Gates::singleQubit) {
                vars.emplace_back(g_s[gate_step][Gates::toIndex(gate)][a]);
            }
            for (int b = 0; b < nqubits; ++b) {
                if (a == b || reducedCM.find({qubitChoice.at(a), qubitChoice.at(b)}) ==
                                      reducedCM.end())
                    continue;
                vars.emplace_back(g_c[gate_step][a][b]);
            }
        }
        lb->assertFormula(ExactlyOneCMDR(
                groupVars(vars, static_cast<std::size_t>(vars.size() / 2U)),
                LogicTerm::noneTerm(), lb));
    }

    // GATE CONSTRAINTS
    for (int gate_step = 1; gate_step < timesteps + 1; ++gate_step) {
        for (int a = 0; a < nqubits; ++a) {
            // NO GATE
            changes = (x[gate_step][a] == x[gate_step - 1][a]);
            changes = changes && (z[gate_step][a] == z[gate_step - 1][a]);

            for (int b = 0; b < nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (x[gate_step][b] == x[gate_step - 1][b]);
                changes = changes && (z[gate_step][b] == z[gate_step - 1][b]);
            }

            changes = changes && (r[gate_step] == r[gate_step - 1]);
            changes = LogicTerm::implies(g_s[gate_step][0][a], changes);
            lb->assertFormula(changes);

            // H
            changes = (z[gate_step][a] == x[gate_step - 1][a]);
            changes = changes && (x[gate_step][a] == z[gate_step - 1][a]);

            for (int b = 0; b < nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (x[gate_step][b] == x[gate_step - 1][b]);
                changes = changes && (z[gate_step][b] == z[gate_step - 1][b]);
            }

            changes = changes &&
                      (r[gate_step] == (r[gate_step - 1] ^
                                        (x[gate_step - 1][a] & z[gate_step - 1][a])));
            changes = LogicTerm::implies(g_s[gate_step][1][a], changes);

            lb->assertFormula(changes);

            // S
            changes =
                    (z[gate_step][a] == (z[gate_step - 1][a] ^ x[gate_step - 1][a]));
            changes = changes && (x[gate_step][a] == x[gate_step - 1][a]);

            for (int b = 0; b < nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (x[gate_step][b] == x[gate_step - 1][b]);
                changes = changes && (z[gate_step][b] == z[gate_step - 1][b]);
            }

            changes = changes &&
                      (r[gate_step] == (r[gate_step - 1] ^
                                        (x[gate_step - 1][a] & z[gate_step - 1][a])));
            changes = LogicTerm::implies(g_s[gate_step][2][a], changes);
            lb->assertFormula(changes);

            // Sdag
            changes =
                    (z[gate_step][a] == (z[gate_step - 1][a] ^ x[gate_step - 1][a]));
            changes = changes && (x[gate_step][a] == x[gate_step - 1][a]);

            for (int b = 0; b < nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (x[gate_step][b] == x[gate_step - 1][b]);
                changes = changes && (z[gate_step][b] == z[gate_step - 1][b]);
            }

            changes = changes &&
                      (r[gate_step] == (r[gate_step - 1] ^
                                        (x[gate_step - 1][a] & (x[gate_step - 1][a] ^ z[gate_step - 1][a]))));
            changes = LogicTerm::implies(g_s[gate_step][Gates::toIndex(Gates::GATES::Sdag)][a], changes);
            lb->assertFormula(changes);

            // Z
            changes =
                    (z[gate_step][a] == z[gate_step - 1][a]);
            changes = changes && (x[gate_step][a] == x[gate_step - 1][a]);

            for (int b = 0; b < nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (x[gate_step][b] == x[gate_step - 1][b]);
                changes = changes && (z[gate_step][b] == z[gate_step - 1][b]);
            }

            changes = changes &&
                      (r[gate_step] == (r[gate_step - 1] ^ x[gate_step - 1][a]));
            changes = LogicTerm::implies(g_s[gate_step][Gates::toIndex(Gates::GATES::Z)][a], changes);
            lb->assertFormula(changes);

            // X
            changes =
                    (z[gate_step][a] == z[gate_step - 1][a]);
            changes = changes && (x[gate_step][a] == x[gate_step - 1][a]);

            for (int b = 0; b < nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (x[gate_step][b] == x[gate_step - 1][b]);
                changes = changes && (z[gate_step][b] == z[gate_step - 1][b]);
            }

            changes = changes &&
                      (r[gate_step] == (r[gate_step - 1] ^ z[gate_step - 1][a]));
            changes = LogicTerm::implies(g_s[gate_step][Gates::toIndex(Gates::GATES::X)][a], changes);
            lb->assertFormula(changes);

            // Y
            changes =
                    (z[gate_step][a] == z[gate_step - 1][a]);
            changes = changes && (x[gate_step][a] == x[gate_step - 1][a]);

            for (int b = 0; b < nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (x[gate_step][b] == x[gate_step - 1][b]);
                changes = changes && (z[gate_step][b] == z[gate_step - 1][b]);
            }

            changes = changes &&
                      (r[gate_step] == (r[gate_step - 1] ^ (z[gate_step - 1][a]) ^ x[gate_step - 1][a]));
            changes = LogicTerm::implies(g_s[gate_step][Gates::toIndex(Gates::GATES::Y)][a], changes);
            lb->assertFormula(changes);

            // CNOT
            for (int b = 0; b < nqubits; ++b) {
                if (reducedCM.find({qubitChoice.at(a), qubitChoice.at(b)}) ==
                    reducedCM.end()) {
                    lb->assertFormula(!g_c[gate_step][a][b]);
                } else {
                    changes =
                            (r[gate_step] == (r[gate_step - 1] ^
                                              ((x[gate_step - 1][a] & z[gate_step - 1][b]) &
                                               ((x[gate_step - 1][b] ^ z[gate_step - 1][a]) ^
                                                LogicTerm((1 << nqubits) - 1, nqubits)))));
                    changes = changes && (x[gate_step][b] ==
                                          (x[gate_step - 1][b] ^ x[gate_step - 1][a]));
                    changes = changes && (z[gate_step][a] ==
                                          (z[gate_step - 1][a] ^ z[gate_step - 1][b]));
                    changes = changes && (x[gate_step][a] == x[gate_step - 1][a]);
                    changes = changes && (z[gate_step][b] == z[gate_step - 1][b]);

                    for (int c = 0; c < nqubits; ++c) { // All other entries do not change
                        if (a == c || b == c)
                            continue;
                        changes = changes && (x[gate_step][c] == x[gate_step - 1][c]);
                        changes = changes && (z[gate_step][c] == z[gate_step - 1][c]);
                    }

                    changes = LogicTerm::implies(g_c[gate_step][a][b], changes);
                    lb->assertFormula(changes);
                }
            }
        }
    }
}
void CliffordOptimizer::makeMultipleGateConstraints(
        LogicBlock* lb, const LogicMatrix& x, const LogicMatrix& z,
        const LogicVector& r, int nqubits, int timesteps,
        const std::set<std::pair<unsigned short, unsigned short>>& reducedCM,
        const std::vector<unsigned short>& qubitChoice, const LogicMatrix3D& g_s,
        const LogicMatrix3D& g_c) {
    LogicTerm changes = LogicTerm(true);
    // CONSISTENCY
    // One gate per qubit, per step
    for (int gate_step = 1; gate_step < timesteps + 1; ++gate_step) {
        for (int a = 0; a < nqubits; ++a) {
            std::vector<LogicTerm> vars{};
            for (auto gate: Gates::singleQubit) {
                vars.emplace_back(g_s[gate_step][Gates::toIndex(gate)][a]);
            }
            for (int b = 0; b < nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                vars.emplace_back(g_c[gate_step][a][b]);
                vars.emplace_back(g_c[gate_step][b][a]);
            }
            lb->assertFormula(ExactlyOneCMDR(
                    groupVars(vars, static_cast<std::size_t>(vars.size() / 2)),
                    LogicTerm::noneTerm(), lb));
        }
    }
    // Maximum any combination of 1 and 2 qubit gates adding up to n
    for (int gate_step = 1; gate_step < timesteps + 1; ++gate_step) {
        changes = LogicTerm(0);
        for (int a = 0; a < nqubits; ++a) {
            for (auto gate: Gates::singleQubit) {
                changes = changes + g_s[gate_step][Gates::toIndex(gate)][a];
            }
            for (int b = 0; b < nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes + g_c[gate_step][a][b] + g_c[gate_step][a][b];
            }
        }
        changes = changes < LogicTerm(static_cast<int>(nqubits + 1));
        lb->assertFormula(changes);
    }

    // GATE CONSTRAINTS
    for (int gate_step = 1; gate_step < timesteps + 1; ++gate_step) {
        LogicTerm r_changes = r[gate_step - 1];
        for (int a = 0; a < nqubits; ++a) {
            // NO GATE
            changes = LogicTerm(true);
            changes = changes && (x[gate_step][a] == x[gate_step - 1][a]);
            changes = changes && (z[gate_step][a] == z[gate_step - 1][a]);
            changes = LogicTerm::implies(g_s[gate_step][0][a], changes);
            lb->assertFormula(changes);

            // H
            changes = LogicTerm(true);
            changes = changes && (z[gate_step][a] == x[gate_step - 1][a]);
            changes = changes && (x[gate_step][a] == z[gate_step - 1][a]);

            r_changes = LogicTerm::ite(
                    g_s[gate_step][1][a],
                    r_changes ^ (x[gate_step - 1][a] & z[gate_step - 1][a]), r_changes);
            changes = LogicTerm::implies(g_s[gate_step][1][a], changes);

            lb->assertFormula(changes);

            // S
            changes = LogicTerm(true);
            changes = changes && (z[gate_step][a] ==
                                  (z[gate_step - 1][a] ^ x[gate_step - 1][a]));
            changes = changes && (x[gate_step][a] == x[gate_step - 1][a]);

            r_changes = LogicTerm::ite(
                    g_s[gate_step][2][a],
                    r_changes ^ (x[gate_step - 1][a] & z[gate_step - 1][a]), r_changes);
            changes = LogicTerm::implies(g_s[gate_step][2][a], changes);
            lb->assertFormula(changes);

            // Sdag
            changes =
                    (z[gate_step][a] == (z[gate_step - 1][a] ^ x[gate_step - 1][a]));
            changes = changes && (x[gate_step][a] == x[gate_step - 1][a]);

            r_changes = LogicTerm::ite(
                    g_s[gate_step][Gates::toIndex(Gates::GATES::Sdag)][a],
                    r_changes ^ (x[gate_step - 1][a] & (x[gate_step - 1][a] ^ z[gate_step - 1][a])), r_changes);
            changes = LogicTerm::implies(g_s[gate_step][Gates::toIndex(Gates::GATES::Sdag)][a], changes);
            lb->assertFormula(changes);

            // Z
            changes =
                    (z[gate_step][a] == z[gate_step - 1][a]);
            changes = changes && (x[gate_step][a] == x[gate_step - 1][a]);

            r_changes = LogicTerm::ite(
                    g_s[gate_step][Gates::toIndex(Gates::GATES::Z)][a],
                    r_changes ^ (x[gate_step - 1][a]), r_changes);
            changes = LogicTerm::implies(g_s[gate_step][Gates::toIndex(Gates::GATES::Z)][a], changes);
            lb->assertFormula(changes);

            // X
            changes =
                    (z[gate_step][a] == z[gate_step - 1][a]);
            changes = changes && (x[gate_step][a] == x[gate_step - 1][a]);

            r_changes = LogicTerm::ite(
                    g_s[gate_step][Gates::toIndex(Gates::GATES::X)][a],
                    r_changes ^ (z[gate_step - 1][a]), r_changes);
            changes = LogicTerm::implies(g_s[gate_step][Gates::toIndex(Gates::GATES::X)][a], changes);
            lb->assertFormula(changes);

            // Y
            changes =
                    (z[gate_step][a] == z[gate_step - 1][a]);
            changes = changes && (x[gate_step][a] == x[gate_step - 1][a]);

            r_changes = LogicTerm::ite(
                    g_s[gate_step][Gates::toIndex(Gates::GATES::Y)][a],
                    r_changes ^ (z[gate_step - 1][a] ^ x[gate_step - 1][a]), r_changes);
            changes = LogicTerm::implies(g_s[gate_step][Gates::toIndex(Gates::GATES::Y)][a], changes);
            lb->assertFormula(changes);

            // CNOT
            for (int b = 0; b < nqubits; ++b) {
                if (reducedCM.find({qubitChoice.at(a), qubitChoice.at(b)}) ==
                    reducedCM.end()) {
                    lb->assertFormula(!g_c[gate_step][a][b]);
                } else {
                    changes   = LogicTerm(true);
                    r_changes = LogicTerm::ite(
                            g_c[gate_step][a][b],
                            (r_changes ^ ((x[gate_step - 1][a] & z[gate_step - 1][b]) &
                                          ((x[gate_step - 1][b] ^ z[gate_step - 1][a]) ^
                                           LogicTerm((1 << nqubits) - 1, nqubits)))),
                            r_changes);
                    changes = changes && (x[gate_step][b] ==
                                          (x[gate_step - 1][b] ^ x[gate_step - 1][a]));
                    changes = changes && (z[gate_step][a] ==
                                          (z[gate_step - 1][a] ^ z[gate_step - 1][b]));
                    changes = changes && (x[gate_step][a] == x[gate_step - 1][a]);
                    changes = changes && (z[gate_step][b] == z[gate_step - 1][b]);
                    changes = LogicTerm::implies(g_c[gate_step][a][b], changes);
                    lb->assertFormula(changes);
                }
            }
        }
        lb->assertFormula(r[gate_step] == r_changes);
    }
}
