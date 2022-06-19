/*
 * This file is part of the MQT QMAP library which is released under the MIT license.
 * See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
 */

#include "exact/ExactMapper.hpp"

void ExactMapper::map(const Configuration& settings) {
    results.config     = settings;
    const auto& config = results.config;

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    initResults();

    // 0) perform pre-mapping optimizations
    preMappingOptimizations(config);

    // 1) create layers according to different criteria
    createLayers();
    if (config.verbose) {
        printLayering(std::cout);
    }
    unsigned long k = 0;
    for (const auto& layer: layers) {
        bool onlySingleQubit = true;
        for (const auto& gate: layer) {
            if (!gate.singleQubit()) {
                onlySingleQubit = false;
                break;
            }
        }
        if (!onlySingleQubit) {
            reducedLayerIndices.emplace_back(k);
        }
        ++k;
    }

    // quickly terminate if the circuit only contains single-qubit gates
    if (reducedLayerIndices.empty()) {
        qcMapped = qc.clone();
        postMappingOptimizations(config);
        results.output.gates = 0U;
        countGates(qcMapped, results.output);
        finalizeMappedCircuit();

        results.output.layers = results.input.layers;
        results.time          = static_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - start).count();
        results.timeout       = false;
        return;
    }

    unsigned long long maxIndex = factorial(qc.getNqubits()) * reducedLayerIndices.size();
    if (maxIndex > std::numeric_limits<int>::max()) {
        std::cerr << "The exact approach can only be used for up to " << std::numeric_limits<int>::max() << " permutation variables, due to 'layers * nq!' overflowing Z3's expr_vector class (uses 'int' index) when trying to instantiate permutation variables y_k_pi. Try reducing the number of layers or the number of qubits." << std::endl;
        return;
    }

    maxIndex = qc.getNqubits() * qc.getNqubits() * reducedLayerIndices.size();
    if (maxIndex > std::numeric_limits<int>::max()) {
        std::cerr << "The exact approach can only be used for up to " << std::numeric_limits<int>::max() << " X variables, due to nq*nq*nlayers overflowing Z3's expr_vector class (uses 'int' index). Try reducing the number of layers or the number of qubits." << std::endl;
        return;
    }

    // 2) For all possibilities k (=m over n) to pick n qubits from m physical qubits
    std::vector<unsigned short> qubitRange{};
    if (!config.subgraph.empty()) {
        const auto subgraphQubits = config.subgraph.size();
        if (subgraphQubits < qc.getNqubits()) {
            std::cerr << "The subgraph must contain at least as many qubits as the circuit has physical qubits." << std::endl;
            return;
        }

        CouplingMap reducedCouplingMap{};
        architecture.getReducedCouplingMap(config.subgraph, reducedCouplingMap);

        // check if the subgraph is connected
        if (!Architecture::isConnected(config.subgraph, reducedCouplingMap)) {
            std::cerr << "The subgraph is not connected." << std::endl;
            return;
        }
        qubitRange = Architecture::getQubitList(reducedCouplingMap);
    } else {
        qubitRange.clear();
        auto qubitSet = architecture.getQubitSet();
        qubitRange.reserve(qubitSet.size());
        std::copy(qubitSet.begin(), qubitSet.end(), std::back_inserter(qubitRange));
    }
    std::vector<QubitChoice> allPossibleQubitChoices{};
    if (config.useSubsets) {
        allPossibleQubitChoices = architecture.getAllConnectedSubsets(qc.getNqubits());
    } else {
        QubitChoice allQubits(qubitRange.begin(), qubitRange.end());
        allPossibleQubitChoices.push_back(allQubits);
    }
    // 3) determine exact mapping for this qubit choice
    std::vector<Swaps> swaps(reducedLayerIndices.size(), Swaps{});
    mappingSwaps.reserve(reducedLayerIndices.size());
    int runs = 1;
    for (auto& choice: allPossibleQubitChoices) {
        std::size_t limit      = 0U;
        std::size_t maxLimit   = 0U;
        std::size_t upperLimit = config.swapLimit;
        if (config.useSubsets) {
            maxLimit = architecture.getCouplingLimit(choice) - 1U;
        } else {
            maxLimit = architecture.getCouplingLimit() - 1U;
        }
        if (config.swapReduction == SwapReduction::CouplingLimit) {
            limit = maxLimit;
        } else if (config.swapReduction == SwapReduction::Increasing) {
            limit = 0U;
        } else { //CustomLimit
            limit = upperLimit;
        }

        unsigned int timeout = 0U;
        do {
            if (config.swapReduction == SwapReduction::Increasing) {
                timeout += settings.timeout * (static_cast<double>(limit * 0.5) / (maxLimit < upperLimit ? upperLimit : maxLimit));
                if (timeout <= 10000U) {
                    timeout = 10000U;
                }
                if (settings.verbose) {
                    std::cout << "Timeout: " << timeout << "  Max-Timeout: " << settings.timeout << std::endl;
                }
            } else {
                timeout = settings.timeout;
            }

            // reset swaps
            for (auto& layer: swaps) {
                layer.clear();
            }

            MappingResults choiceResults{};
            choiceResults.copyInput(results);
            choiceResults.config.swapLimit        = limit;
            choiceResults.output.swaps            = 0U;
            choiceResults.output.directionReverse = 0U;
            choiceResults.output.gates            = std::numeric_limits<unsigned long>::max();

            // 4) reduce coupling map
            CouplingMap reducedCouplingMap = {};
            architecture.getReducedCouplingMap(choice, reducedCouplingMap);

            if (reducedCouplingMap.empty()) {
                break;
            }

            if (config.verbose) {
                std::cout << "-------- qubit choice: ";
                for (const auto Q: choice) {
                    std::cout << Q << " ";
                }
                std::cout << "---------- ";
                if (config.swapReduction != SwapReduction::None) {
                    std::cout << "SWAP limit: " << limit;
                }
                std::cout << std::endl;
            }

            // 6) call actual mapping routine
            coreMappingRoutine(choice, reducedCouplingMap, choiceResults, swaps, static_cast<long unsigned int>(limit), timeout);

            if (config.verbose) {
                if (!choiceResults.timeout) {
                    std::cout << "Costs: " << choiceResults.output.swaps << " SWAP(s)";
                    if (!architecture.bidirectional()) {
                        std::cout << ", " << choiceResults.output.directionReverse << " direction reverses";
                    }
                    std::cout << std::endl;
                } else {
                    std::cout << "Did not yield a result" << std::endl;
                }
            }

            // 7) Check if new optimum found
            if (!choiceResults.timeout && choiceResults.output.gates < results.output.gates) {
                results      = choiceResults;
                mappingSwaps = swaps;
            }
            if (limit == 0) {
                limit = 1;
            } else {
                limit += runs;
                runs++;
            }
        } while (config.swapReduction == SwapReduction::Increasing && (limit <= upperLimit || config.swapLimit == 0) && limit < architecture.getCouplingLimit());

        // stop if a perfect result has been found
        if (!results.timeout && results.output.swaps == 0U && results.output.directionReverse == 0U) {
            break;
        }
    }

    // return in case no result has been found
    if (results.timeout) {
        return;
    }

    // 8) Write best result and statistics
    auto layerIterator = reducedLayerIndices.begin();
    auto swapsIterator = mappingSwaps.begin();

    if (settings.verbose) {
        auto it = reducedLayerIndices.begin();
        for (const auto& layer: mappingSwaps) {
            std::cout << *it << ": ";
            for (const auto& swap: layer) {
                std::cout << "(" << swap.first << "<->" << swap.second << ") ";
            }
            ++it;
            std::cout << std::endl;
        }
    }

    for (const auto& q: qubitRange) {
        locations.at(q) = static_cast<short>(q);
    }

    for (unsigned long i = 0U; i < layers.size(); ++i) {
        if (i == 0U) {
            qcMapped.initialLayout.clear();
            qcMapped.outputPermutation.clear();

            // no swaps but initial permutation
            for (const auto& [physical, logical]: *swapsIterator) {
                locations.at(logical)                                        = static_cast<short>(physical);
                qubits.at(physical)                                          = static_cast<short>(logical);
                qcMapped.initialLayout[static_cast<dd::Qubit>(physical)]     = static_cast<dd::Qubit>(logical);
                qcMapped.outputPermutation[static_cast<dd::Qubit>(physical)] = static_cast<dd::Qubit>(logical);
            }

            // place remaining architecture qubits
            placeRemainingArchitectureQubits();

            if (settings.verbose) {
                for (auto q = 0U; q < architecture.getNqubits(); ++q) {
                    std::cout << qubits.at(q) << " ";
                }
                std::cout << std::endl;
            }
            ++swapsIterator;
        }

        // apply all gates of layer
        for (const auto& gate: layers.at(i)) {
            auto op = dynamic_cast<qc::StandardOperation*>(gate.op);
            if (!op) {
                throw QMAPException("Cast to StandardOperation not possible during mapping. Check that circuit contains only StandardOperations");
            }

            if (gate.singleQubit()) {
                if (settings.verbose) {
                    std::cout << i << ": Added single qubit gate with target: " << locations.at(gate.target) << std::endl;
                }

                qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(),
                                                             locations.at(gate.target),
                                                             op->getType(),
                                                             op->getParameter().at(0),
                                                             op->getParameter().at(1),
                                                             op->getParameter().at(2));
            } else {
                Edge cnot = {locations.at(gate.control), locations.at(gate.target)};

                if (architecture.getCouplingMap().find(cnot) == architecture.getCouplingMap().end()) {
                    Edge reverse = {cnot.second, cnot.first};
                    if (architecture.getCouplingMap().find(reverse) == architecture.getCouplingMap().end()) {
                        throw QMAPException("Invalid CNOT: " + std::to_string(reverse.first) + "-" + std::to_string(reverse.second));
                    }
                    if (settings.verbose) {
                        std::cout << i << ": Added (direction-reversed) cnot with control and target: " << cnot.first << " " << cnot.second << std::endl;
                    }
                    qcMapped.h(reverse.first);
                    qcMapped.h(reverse.second);
                    qcMapped.x(reverse.second, dd::Control{static_cast<dd::Qubit>(reverse.first)});
                    qcMapped.h(reverse.second);
                    qcMapped.h(reverse.first);
                } else {
                    if (settings.verbose) {
                        std::cout << i << ": Added cnot with control and target: " << cnot.first << " " << cnot.second << std::endl;
                    }
                    qcMapped.x(cnot.second, dd::Control{static_cast<dd::Qubit>(cnot.first)});
                }
            }
        }

        if (!mappingSwaps.empty() && swapsIterator != mappingSwaps.end() && layerIterator != reducedLayerIndices.end() && i == *layerIterator) {
            // apply swaps before layer
            for (auto it = (*swapsIterator).rbegin(); it != (*swapsIterator).rend(); ++it) {
                auto& swap = *it;
                qcMapped.swap(swap.first, swap.second);
                std::swap(qcMapped.outputPermutation.at(swap.first), qcMapped.outputPermutation.at(swap.second));
                std::swap(qubits.at(swap.first), qubits.at(swap.second));
                std::swap(locations.at(qubits.at(swap.first)), locations.at(qubits.at(swap.second)));

                if (settings.verbose) {
                    for (auto q = 0U; q < architecture.getNqubits(); ++q) {
                        std::cout << qubits.at(q) << " ";
                    }
                    std::cout << std::endl;
                }
            }

            ++swapsIterator;
            ++layerIterator;
        }
    }

    // 9) apply post mapping optimizations
    postMappingOptimizations(config);

    // 10) re-count gates
    results.output.singleQubitGates = 0U;
    results.output.cnots            = 0U;
    results.output.gates            = 0U;
    countGates(qcMapped, results.output);

    // 11) final post-processing
    finalizeMappedCircuit();

    auto                          end  = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    results.time                       = diff.count();
}

void ExactMapper::coreMappingRoutine(const std::set<unsigned short>& qubitChoice, const CouplingMap& rcm, MappingResults& choiceResults, std::vector<std::vector<std::pair<unsigned short, unsigned short>>>& swaps, long unsigned int limit, unsigned int timeout) {
    const auto& config = results.config;

    // Z3 context
    context c;

    std::vector<unsigned short>                        pi(qubitChoice.begin(), qubitChoice.end());
    unsigned long long                                 piCount{};
    unsigned long long                                 internalPiCount{};
    std::unordered_set<unsigned long long>             skipped_pi{};
    std::unordered_map<unsigned short, unsigned short> physicalQubitIndex{};
    unsigned short                                     qIdx = 0;
    for (const auto& Q: qubitChoice) {
        physicalQubitIndex[Q] = qIdx;
        ++qIdx;
    }

    //////////////////////////////////////////
    /// 	Check necessary permutations	//
    //////////////////////////////////////////
    if (config.enableSwapLimits && !config.useBDD) {
        do {
            auto picost = architecture.minimumNumberOfSwaps(pi, static_cast<long>(limit));
            if (picost > limit) {
                skipped_pi.insert(piCount);
            }
            ++piCount;
        } while (std::next_permutation(pi.begin(), pi.end()));
    }
    //////////////////////////////////////////
    /// 	Boolean Variable Definitions	//
    //////////////////////////////////////////
    /*
	 auxilary variable declarations
	*/
    expr_vector auxvars(c);

    /*
	 locical/physical qubit variables x_k_i_j
	 k	before layer k
	 i	physical qubit i
	 j	logical qubit j
	 number of variables: (|L|) * m * n
	 */
    std::vector<matrix> x{};
    std::stringstream   x_name{};
    for (unsigned long k = 0; k < reducedLayerIndices.size(); ++k) {
        x.emplace_back();
        for (unsigned short Q: qubitChoice) {
            x.back().emplace_back(c);
            for (unsigned short q = 0; q < qc.getNqubits(); ++q) {
                x_name.str("");
                x_name << "x_" << k << '_' << Q << '_' << q;
                x.back().back().push_back(c.bool_const(x_name.str().c_str()));
            }
        }
    }

    /*
 permutation variables y_k_pi
 k	before layer k
 pi	arbitrary permutation of the m qubits
 number of variables: (|L|-1) * m!
 */
    std::vector<std::vector<expr>> y{};
    std::stringstream              y_name{};
    for (unsigned long k = 1; k < reducedLayerIndices.size(); ++k) {
        y.emplace_back();
        piCount = 0;
        do {
            if (skipped_pi.count(piCount) == 0 || !config.enableSwapLimits) {
                y_name.str("");
                y_name << "y_" << k << '_' << piCount;
                y.back().push_back(c.bool_const(y_name.str().c_str()));
            }
            ++piCount;
        } while (std::next_permutation(pi.begin(), pi.end()));
    }

    // Z3 optimizer
    optimize opt(c);
    params   p(c);
    p.set("timeout", timeout);
    p.set("pb.compile_equality", true);
    p.set("pp.wcnf", true);
    p.set("maxres.hill_climb", true);
    p.set("maxres.pivot_on_correction_set", false);
    opt.set(p);

    //////////////////////////////////////////
    /// 	Consistency Constraints			//
    //////////////////////////////////////////
    if (config.encoding == Encoding::Naive) {
        for (unsigned long k = 0; k < reducedLayerIndices.size(); ++k) {
            for (unsigned long i = 0; i < qubitChoice.size(); ++i) {
                expr rowConsistency = c.int_val(0);
                for (unsigned short j = 0; j < qc.getNqubits(); ++j) {
                    rowConsistency = rowConsistency +
                                     ite(x[k][i][j], c.int_val(1), c.int_val(0));
                }
                opt.add(rowConsistency.simplify() <= 1);
            }

            for (unsigned short j = 0; j < qc.getNqubits(); ++j) {
                expr colConsistency = c.int_val(0);
                for (unsigned long i = 0; i < qubitChoice.size(); ++i) {
                    colConsistency = colConsistency +
                                     ite(x[k][i][j], c.int_val(1), c.int_val(0));
                }
                opt.add(colConsistency.simplify() == 1);
            }
        }
    } else if (config.encoding == Encoding::Commander) {
        for (unsigned long k = 0; k < reducedLayerIndices.size(); ++k) {
            for (unsigned long i = 0; i < qubitChoice.size(); ++i) {
                std::vector<expr> varIDs;
                for (unsigned short j = 0; j < qc.getNqubits(); ++j) {
                    varIDs.push_back(x[k][i][j]);
                }
                if (config.commanderGrouping == CommanderGrouping::Fixed2) {
                    opt.add(AtMostOneCMDR(varIDs, groupVars(varIDs, 2), -1, auxvars, c));
                } else if (config.commanderGrouping == CommanderGrouping::Fixed3) {
                    opt.add(AtMostOneCMDR(varIDs, groupVars(varIDs, 3), -1, auxvars, c));
                } else if (config.commanderGrouping == CommanderGrouping::Logarithm) {
                    opt.add(AtMostOneCMDR(varIDs, groupVars(varIDs, static_cast<std::size_t>(std::log(varIDs.size()))), -1, auxvars, c));
                } else if (config.commanderGrouping == CommanderGrouping::Halves) {
                    opt.add(AtMostOneCMDR(varIDs, groupVars(varIDs, varIDs.size() / 2), -1, auxvars, c));
                }
            }

            for (unsigned short j = 0; j < qc.getNqubits(); ++j) {
                std::vector<expr> varIDs;
                for (unsigned long i = 0; i < qubitChoice.size(); ++i) {
                    varIDs.push_back(x[k][i][j]);
                }

                if (config.commanderGrouping == CommanderGrouping::Fixed2) {
                    opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, 2), -1, auxvars, c));
                } else if (config.commanderGrouping == CommanderGrouping::Fixed3) {
                    opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, 3), -1, auxvars, c));
                } else if (config.commanderGrouping == CommanderGrouping::Logarithm) {
                    opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, static_cast<std::size_t>(std::log(varIDs.size()))), -1, auxvars, c));
                } else if (config.commanderGrouping == CommanderGrouping::Halves) {
                    opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, varIDs.size() / 2), -1, auxvars, c));
                }
            }
        }
    } else if (config.encoding == Encoding::Bimander) {
        for (unsigned long k = 0; k < reducedLayerIndices.size(); ++k) {
            for (unsigned long i = 0; i < qubitChoice.size(); ++i) {
                std::vector<expr>          vars;
                std::vector<unsigned long> varIDs;
                for (unsigned short j = 0; j < qc.getNqubits(); ++j) {
                    vars.emplace_back(x[k][i][j]);
                    varIDs.emplace_back(j);
                }
                opt.add(AtMostOneBiMander(vars, varIDs, auxvars, c));
            }

            for (unsigned short j = 0; j < qc.getNqubits(); ++j) { //There is no exactly one Bimander
                std::vector<expr> varIDs;
                for (unsigned long i = 0; i < qubitChoice.size(); ++i) {
                    varIDs.push_back(x[k][i][j]);
                }

                if (config.commanderGrouping == CommanderGrouping::Fixed2) {
                    opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, 2), -1, auxvars, c));
                } else if (config.commanderGrouping == CommanderGrouping::Fixed3) {
                    opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, 3), -1, auxvars, c));
                } else if (config.commanderGrouping == CommanderGrouping::Logarithm) {
                    opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, static_cast<std::size_t>(std::log(varIDs.size()))), -1, auxvars, c));
                } else if (config.commanderGrouping == CommanderGrouping::Halves) {
                    opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, varIDs.size() / 2), -1, auxvars, c));
                }
            }
        }
    }

    //////////////////////////////////////////
    ///		Coupling Constraints			//
    //////////////////////////////////////////
    for (unsigned long k = 0; k < reducedLayerIndices.size(); ++k) {
        expr allCouplings = c.bool_val(true);
        for (const auto& gate: layers.at(reducedLayerIndices.at(k))) {
            if (gate.singleQubit())
                continue;

            expr coupling = c.bool_val(false);
            if (architecture.bidirectional()) {
                for (const auto& edge: rcm) {
                    auto indexFC = x[k][physicalQubitIndex[edge.first]][gate.control];
                    auto indexST = x[k][physicalQubitIndex[edge.second]][gate.target];
                    coupling     = coupling || (indexFC && indexST);
                }
            } else {
                for (const auto& edge: rcm) {
                    auto indexFC = x[k][physicalQubitIndex[edge.first]][gate.control];
                    auto indexST = x[k][physicalQubitIndex[edge.second]][gate.target];
                    auto indexFT = x[k][physicalQubitIndex[edge.first]][gate.target];
                    auto indexSC = x[k][physicalQubitIndex[edge.second]][gate.control];

                    coupling = coupling || ((indexFC && indexST) || (indexFT && indexSC));
                }
            }
            allCouplings = allCouplings && coupling.simplify();
        }
        opt.add(allCouplings.simplify());
    }

    //////////////////////////////////////////
    /// 	Permutation Constraints			//
    //////////////////////////////////////////
    for (unsigned long k = 1; k < reducedLayerIndices.size(); ++k) {
        piCount         = 0;
        internalPiCount = 0;
        auto& i         = x[k - 1];
        auto& j         = x[k];
        do {
            if (skipped_pi.count(piCount) == 0 || !config.enableSwapLimits) {
                expr equal = c.bool_val(true);
                for (unsigned short Q: qubitChoice) {
                    for (unsigned short q = 0; q < qc.getNqubits(); ++q) {
                        auto before = i[physicalQubitIndex[Q]][q];
                        auto after  = j[physicalQubitIndex[pi[physicalQubitIndex[Q]]]][q];
                        equal       = equal && (before == after);
                    }
                }
                opt.add(implies(y[k - 1][internalPiCount], equal.simplify()).simplify());
                ++internalPiCount;
            }
            ++piCount;
        } while (std::next_permutation(pi.begin(), pi.end()));
    }

    // Allow only 1 y_k_pi to be true
    if (config.encoding == Encoding::Naive) {
        for (unsigned long k = 1; k < reducedLayerIndices.size(); ++k) {
            expr onlyOne    = c.int_val(0);
            piCount         = 0;
            internalPiCount = 0;
            do {
                if (skipped_pi.count(piCount) == 0 || !config.enableSwapLimits) {
                    onlyOne = onlyOne + ite(y[k - 1][internalPiCount], c.int_val(1), c.int_val(0));
                    ++internalPiCount;
                }
                ++piCount;
            } while (std::next_permutation(pi.begin(), pi.end()));
            opt.add(onlyOne.simplify() == 1);
        }
    } else {
        for (unsigned long k = 1; k < reducedLayerIndices.size(); ++k) {
            std::vector<expr> varIDs;
            piCount         = 0;
            internalPiCount = 0;
            do {
                if (skipped_pi.count(piCount) == 0 || !config.enableSwapLimits) {
                    varIDs.push_back(y[k - 1][internalPiCount]);
                    ++internalPiCount;
                }
                ++piCount;
            } while (std::next_permutation(pi.begin(), pi.end()));
            if (config.commanderGrouping == CommanderGrouping::Fixed2) {
                opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, 2), -1, auxvars, c));
            } else if (config.commanderGrouping == CommanderGrouping::Fixed3) {
                opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, 3), -1, auxvars, c));
            } else if (config.commanderGrouping == CommanderGrouping::Logarithm) {
                opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, static_cast<std::size_t>(std::log(varIDs.size()))), -1, auxvars, c));
            } else if (config.commanderGrouping == CommanderGrouping::Halves) {
                opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, varIDs.size() / 2), -1, auxvars, c));
            }
        }
    }
    //////////////////////////////////////////
    /// 	Objective Function				//
    //////////////////////////////////////////
    // cost for permutations
    piCount         = 0;
    internalPiCount = 0;
    std::vector<std::set<WeightedVar>> weightedVars(reducedLayerIndices.size());
    do {
        if (skipped_pi.count(piCount) == 0 || !config.enableSwapLimits) {
            auto picost = architecture.minimumNumberOfSwaps(pi);
            if (architecture.bidirectional()) {
                picost *= GATES_OF_BIDIRECTIONAL_SWAP;
            } else {
                picost *= GATES_OF_UNIDIRECTIONAL_SWAP;
            }
            for (unsigned long k = 1; k < reducedLayerIndices.size(); ++k) {
                opt.add(!y[k - 1][internalPiCount], picost);
                if (config.useBDD)
                    weightedVars[k].insert(WeightedVar(internalPiCount, static_cast<int>(picost)));
            }
            ++internalPiCount;
        }
        ++piCount;
    } while (std::next_permutation(pi.begin(), pi.end()));

    if (config.enableSwapLimits && config.useBDD) {
        for (unsigned long k = 1; k < reducedLayerIndices.size(); ++k) {
            opt.add(BuildBDD(weightedVars[k], y[k - 1], auxvars, static_cast<int>(limit), c));
        }
    }

    // cost for reversed directions
    if (!architecture.bidirectional()) {
        for (unsigned long k = 0; k < reducedLayerIndices.size(); ++k) {
            for (const auto& gate: layers.at(reducedLayerIndices.at(k))) {
                if (gate.singleQubit())
                    continue;

                expr reverse = c.bool_val(true);
                for (const auto& edge: rcm) {
                    auto indexFT = x[k][physicalQubitIndex[edge.first]][gate.target];
                    auto indexSC = x[k][physicalQubitIndex[edge.second]][gate.control];
                    reverse      = reverse && (!indexFT || !indexSC);
                }
                opt.add(reverse.simplify(), GATES_OF_DIRECTION_REVERSE);
            }
        }
    }

    if (config.includeWCNF) {
        std::stringstream ss{};
        ss << opt;
        try {
            c.check_error();
        } catch (const z3::exception& e) {
            std::cerr << "Z3 reported an exception while trying to gather the WCNF formula: " << e.msg() << std::endl;
            std::cerr << "Most likely, this is due to the usage of Z3's atMostOne and exactlyOne constraints." << std::endl;
            std::cerr << "This can be circumvented by using QMAP's `commander` or `bimander` encoding." << std::endl;
        }
        choiceResults.wcnf = ss.str();
    }

    //////////////////////////////////////////
    /// 	Solving							//
    //////////////////////////////////////////
    if (sat == opt.check()) {
        model m               = opt.get_model();
        choiceResults.timeout = results.timeout = false;

        // quickly determine cost
        choiceResults.output.singleQubitGates = choiceResults.input.singleQubitGates;
        choiceResults.output.cnots            = choiceResults.input.cnots;
        choiceResults.output.gates            = choiceResults.output.singleQubitGates + choiceResults.output.cnots;
        assert(choiceResults.output.swaps == 0U);
        assert(choiceResults.output.directionReverse == 0U);

        // swaps
        for (unsigned long k = 1; k < reducedLayerIndices.size(); ++k) {
            auto& i = x[k - 1];
            auto& j = x[k];

            for (unsigned short Q: qubitChoice) {
                for (unsigned short q = 0; q < qc.getNqubits(); ++q) {
                    if (eq(m.eval(i[physicalQubitIndex[Q]][q]), c.bool_val(true))) {
                        // logical qubit q was mapped to physical qubit Q
                        for (unsigned short P: qubitChoice) {
                            // and has been assigned to physical qubit P going forward
                            if (eq(m.eval(j[physicalQubitIndex[P]][q]), c.bool_val(true))) {
                                pi[physicalQubitIndex[Q]] = P;
                            }
                        }
                    }
                }
            }
            architecture.minimumNumberOfSwaps(pi, swaps.at(k));
            choiceResults.output.swaps += swaps.at(k).size();
            if (architecture.bidirectional()) {
                choiceResults.output.gates += GATES_OF_BIDIRECTIONAL_SWAP * swaps.at(k).size();
            } else {
                choiceResults.output.gates += GATES_OF_UNIDIRECTIONAL_SWAP * swaps.at(k).size();
            }
        }

        // direction reverse
        if (!architecture.bidirectional()) {
            for (unsigned long k = 0; k < reducedLayerIndices.size(); ++k) {
                for (const auto& gate: layers.at(reducedLayerIndices.at(k))) {
                    if (gate.singleQubit())
                        continue;
                    for (const auto& edge: rcm) {
                        auto indexFT = x[k][physicalQubitIndex[edge.first]][gate.target];
                        auto indexSC = x[k][physicalQubitIndex[edge.second]][gate.control];
                        if (eq(m.eval(indexFT && indexSC), c.bool_val(true))) {
                            choiceResults.output.directionReverse++;
                            choiceResults.output.gates += GATES_OF_DIRECTION_REVERSE;
                        }
                    }
                }
            }
        }

        // save initial layout for later
        for (const auto& Q: qubitChoice) {
            for (unsigned short q = 0; q < qc.getNqubits(); ++q) {
                bool set = eq(m.eval(x[0][physicalQubitIndex[Q]][q]), c.bool_val(true));
                if (set) {
                    swaps.at(0).emplace_back(std::pair<unsigned short, unsigned short>{Q, q});
                }
            }
        }

    } else {
        results.timeout = true;
    }
}
