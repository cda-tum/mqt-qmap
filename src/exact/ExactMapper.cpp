/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#include "exact/ExactMapper.hpp"

void ExactMapper::map(const MappingSettings& settings) {

	this->settings = settings;
	std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
	qc.stripIdleQubits(true, true);
	initResults();

	// 1) create layers according to different criteria
	createLayers();
	if (settings.verbose) {
		printLayering(std::cout);
	}
	unsigned long k = 0;
	for (const auto& layer: layers) {
		bool onlySingleQubit = true;
		for (const auto& gate: layer) {
			if(!gate.singleQubit()) {
				onlySingleQubit = false;
				break;
			}
		}
		if (!onlySingleQubit) {
			reducedLayerIndices.emplace_back(k);
		}
		++k;
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
	std::vector<unsigned short> qubits{};
	for (unsigned short i = 0; i < architecture.getNqubits(); ++i) {
		qubits.push_back(i);
	}
	std::vector<std::set<unsigned short>> allPossibleQubitChoices{};
	if (this->settings.useQubitSubsets) {
		do {
			std::set<unsigned short> qubitChoice{};
			for (unsigned short i = 0; i < qc.getNqubits(); ++i) {
				qubitChoice.insert(qubits.at(i));
			}
			allPossibleQubitChoices.push_back(qubitChoice);
		} while (next_combination(qubits.begin(), qubits.begin()+qc.getNqubits(), qubits.end()));
	} else {
		std::set<unsigned short> allQubits(qubits.begin(), qubits.end());
		allPossibleQubitChoices.push_back(allQubits);
	}
	// 3) determine exact mapping for this qubit choice
	std::vector<std::vector<std::pair<unsigned short, unsigned short>>> swaps(reducedLayerIndices.size(), std::vector<std::pair<unsigned short, unsigned short>>{});
	mappingSwaps.reserve(reducedLayerIndices.size());
	int runs = 1;
	for (auto& choice: allPossibleQubitChoices) {
		int limit = 0;
		unsigned int upperLimit = this->settings.limit;
		if (this->settings.strategy == Strategy::ArchitectureSwaps){
			limit = this->architecture.getLongestPath()-1;
		} else if (this->settings.strategy == Strategy::SubsetSwaps) {
			limit = this->architecture.getLongestPath(choice)-1;
		} else if (this->settings.strategy == Strategy::Increasing) {
		} else { //CustomLimit
			limit = upperLimit;
		}
		unsigned int maxLimit = this->architecture.getLongestPath();
		unsigned int timeout = 0;
		do {	
			timeout += settings.timeout * (static_cast<double>(limit*0.5)/(maxLimit<upperLimit?upperLimit:maxLimit));
			if (timeout<=10000)
				timeout = 10000;
			if (this->settings.strategy != Strategy::Increasing)
				timeout = settings.timeout;
			if (settings.verbose)
				std::cout << "Timeout: " << timeout << "  Max-Timeout: " << settings.timeout << std::endl;
			// reset swaps
			for(auto& layer: swaps) {
				layer.clear();
			}

			MappingResults choiceResults{};
			choiceResults.copyInput(results);
			choiceResults.limit = limit;

			// 4) reduce coupling map
			CouplingMap reducedCouplingMap = architecture.getCouplingMap();
			for (const auto& edge: architecture.getCouplingMap()) {
				if (!choice.count(edge.first) || !choice.count(edge.second)) {
					reducedCouplingMap.erase(edge);
				}
			}
			if (reducedCouplingMap.empty()) {
				break;
			}

			// 5) Check if E_k is connected. If yes, then possible subset found
			std::set<unsigned short> reachedQubits{};
			reachedQubits.insert(*(choice.begin()));
			dfs(*(choice.begin()), reachedQubits, reducedCouplingMap);
			if (!(reachedQubits == choice)){
				break;
			}

			// 6) call actual mapping routine
			coreMappingRoutine(choice, reducedCouplingMap, choiceResults, swaps, static_cast<long unsigned int>(limit), timeout);

			if (settings.verbose) {
				std::cout << "SWAPs: " << choiceResults.output_swaps << std::endl;
				std::cout << "Direction reverses: " << choiceResults.output_direction_reverse << std::endl;
			}

			// 7) Check if new optimum found
			if (!choiceResults.timeout && choiceResults.output_gates < results.output_gates) {
				results = choiceResults;
				mappingSwaps = swaps;
			}
			if (limit == 0) {
				limit = 1;
			} else {
				limit += runs;
				runs++;
			}
		} while (this->settings.strategy == Strategy::Increasing && (limit <= upperLimit || this->settings.limit==0) && limit < this->architecture.getLongestPath());
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

	for (const auto& q: qubits) {
		locations.at(q) = static_cast<short>(q);
	}

	for(unsigned long i=0; i<layers.size(); ++i) {
		if (i == 0) {
			qc::Permutation inverseInitialLayout{ };
			for (auto& pu: qc.initialLayout)
				inverseInitialLayout.insert({pu.second, pu.first});

			// no swaps but initial permutation
			for(const auto& assignment: *swapsIterator) {
				locations.at(inverseInitialLayout.at(assignment.second)) = static_cast<short>(assignment.first);
				qubits.at(assignment.first) = inverseInitialLayout.at(assignment.second);
				qcMapped.initialLayout.at(assignment.first) = inverseInitialLayout.at(assignment.second);
				qcMapped.outputPermutation.at(assignment.first) = inverseInitialLayout.at(assignment.second);
			}

			if(settings.verbose) {
				for (const auto& q:qubits) {
					std::cout << q << " " ;
				}
				std::cout << std::endl;
			}
			++swapsIterator;
		}



		// apply all gates of layer
		for(const auto& gate: layers.at(i)) {
			auto op = dynamic_cast<qc::StandardOperation*>(gate.op);
			if (!op) {
				throw QMAPException("Cast to StandardOperation not possible during mapping. Check that circuit contains only StandardOperations");
			}

			if (gate.singleQubit()) {
				if(settings.verbose) {
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

		if (!mappingSwaps.empty() && swapsIterator != mappingSwaps.end() && layerIterator !=  reducedLayerIndices.end() && i == *layerIterator) {
			// apply swaps before layer
			for (auto it=(*swapsIterator).rbegin(); it != (*swapsIterator).rend(); ++it) {
				auto& swap = *it;
				qcMapped.x(swap.second, dd::Control{static_cast<dd::Qubit>(swap.first)});
				std::swap(qcMapped.outputPermutation.at(swap.first), qcMapped.outputPermutation.at(swap.second));
				std::swap(qubits.at(swap.first), qubits.at(swap.second));
				std::swap(locations.at(qubits.at(swap.first)), locations.at(qubits.at(swap.second)));

				if (settings.verbose) {
					for (const auto& q:qubits) {
						std::cout << q << " " ;
					}
					std::cout << std::endl;
				}
			}

			++swapsIterator;
			++layerIterator;
		}

	}

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff = end - start;
	results.time = diff.count();
}

void ExactMapper::coreMappingRoutine(const std::set<unsigned short>& qubitChoice, const CouplingMap& rcm, MappingResults& choiceResults, std::vector<std::vector<std::pair<unsigned short, unsigned short>>>& swaps, long unsigned int limit, unsigned int timeout) {
	// Z3 context
	context c;

	std::vector<unsigned short> pi(qubitChoice.begin(), qubitChoice.end());
	unsigned long long piCount;
	std::unordered_set<unsigned long long> skipped_pi{};
	std::vector<std::unordered_set<unsigned long long>> skipped_pi_per_layer{};
	std::unordered_map<unsigned short, unsigned short> physicalQubitIndex{};
	unsigned short qIdx=0;
	for (const auto& Q: qubitChoice){
		physicalQubitIndex[Q] = qIdx;
		++qIdx;
	}

	//////////////////////////////////////////
	/// 	Check necessary permutations	//
	//////////////////////////////////////////
	skipped_pi_per_layer.resize(reducedLayerIndices.size());
	if (this->settings.enableLimits && !this->settings.useBDD) {
		do {
			auto picost = architecture.minimumNumberOfSwaps(pi);
			if (picost > limit) {
				skipped_pi.insert(piCount);
			}
			++piCount;
		} while(std::next_permutation(pi.begin(), pi.end()));
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
	std::stringstream x_name{};
	for (unsigned long k=0; k<reducedLayerIndices.size(); ++k) {
		x.emplace_back();
		for (unsigned short Q: qubitChoice) {
			x.back().emplace_back(c);
			for (unsigned short q=0; q<qc.getNqubits(); ++q) {
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
	std::stringstream y_name{};
	for (unsigned long k=1; k<reducedLayerIndices.size(); ++k) {
		y.emplace_back();
		piCount = 0;
		do {
			if (skipped_pi.count(piCount) == 0 || !this->settings.enableLimits){
				y_name.str("");
				y_name << "y_" << k << '_' << piCount;
				y.back().push_back(c.bool_const(y_name.str().c_str()));
			}
			++piCount;
		} while(std::next_permutation(pi.begin(), pi.end()));
	}

	// Z3 optimizer
	optimize opt(c);
	params p(c);
	p.set("timeout", timeout);
	p.set("pb.compile_equality", true);
	p.set("maxres.hill_climb", true);
	p.set("maxres.pivot_on_correction_set", false);
	opt.set(p);

	//////////////////////////////////////////
	/// 	Consistency Constraints			//
	//////////////////////////////////////////
	if (this->settings.encoding == Encodings::None) {
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
	}
	else if (this->settings.encoding == Encodings::Commander) {
		for (unsigned long k = 0; k < reducedLayerIndices.size(); ++k) {
			for (unsigned long i = 0; i < qubitChoice.size(); ++i) {
				std::vector<expr> varIDs;
				for (unsigned short j = 0; j < qc.getNqubits(); ++j) {
					varIDs.push_back(x[k][i][j]);
				}
				if (this->settings.grouping == Groupings::Fixed2) {
					opt.add(AtMostOneCMDR(varIDs, groupVars(varIDs, 2), static_cast<int>(auxvars.size() - 1), auxvars, c));
				} else if (this->settings.grouping == Groupings::Fixed3) {
					opt.add(AtMostOneCMDR(varIDs, groupVars(varIDs, 3), static_cast<int>(auxvars.size() - 1), auxvars, c));
				} else if (this->settings.grouping == Groupings::Logarithm) {
					opt.add(AtMostOneCMDR(varIDs, groupVars(varIDs, static_cast<std::size_t>(std::log(varIDs.size()))), static_cast<int>(auxvars.size() - 1), auxvars, c));
				} else if (this->settings.grouping == Groupings::Halves) {
					opt.add(AtMostOneCMDR(varIDs, groupVars(varIDs, varIDs.size() / 2), static_cast<int>(auxvars.size() - 1), auxvars, c));
				}
			}

			for (unsigned short j = 0; j < qc.getNqubits(); ++j) {
				std::vector<expr> varIDs;
				for (unsigned long i = 0; i < qubitChoice.size(); ++i) {
					varIDs.push_back(x[k][i][j]);
				}

				if (this->settings.grouping == Groupings::Fixed2) {
					opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, 2), static_cast<int>(auxvars.size() - 1), auxvars, c));
				} else if (this->settings.grouping == Groupings::Fixed3) {
					opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, 3), static_cast<int>(auxvars.size() - 1), auxvars, c));
				} else if (this->settings.grouping == Groupings::Logarithm) {
					opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, static_cast<std::size_t>(std::log(varIDs.size()))), static_cast<int>(auxvars.size() - 1), auxvars, c));
				} else if (this->settings.grouping == Groupings::Halves) {
					opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, varIDs.size() / 2), static_cast<int>(auxvars.size() - 1), auxvars, c));
				}
			}
		}
	}
	else if (this->settings.encoding == Encodings::Bimander) {
		for (unsigned long k = 0; k < reducedLayerIndices.size(); ++k) {
			for (unsigned long i = 0; i < qubitChoice.size(); ++i) {
				std::vector<expr> vars;
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

				if (this->settings.grouping == Groupings::Fixed2) {
					opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, 2), static_cast<int>(auxvars.size() - 1), auxvars, c));
				} else if (this->settings.grouping == Groupings::Fixed3){
					opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, 3), static_cast<int>(auxvars.size() - 1), auxvars, c));
				} else if (this->settings.grouping == Groupings::Logarithm) {
					opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, static_cast<std::size_t>(std::log(varIDs.size()))), static_cast<int>(auxvars.size() - 1), auxvars, c));
				} else if (this->settings.grouping == Groupings::Halves) {
					opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, varIDs.size() / 2), static_cast<int>(auxvars.size() - 1), auxvars, c));
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
					auto indexFC = x[k][physicalQubitIndex[edge.first]][qc.initialLayout.at(gate.control)];
					auto indexST = x[k][physicalQubitIndex[edge.second]][qc.initialLayout.at(gate.target)];
					coupling = coupling || (indexFC && indexST);
				}
			} else {
				for (const auto& edge: rcm) {
					auto indexFC = x[k][physicalQubitIndex[edge.first]][qc.initialLayout.at(gate.control)];
					auto indexST = x[k][physicalQubitIndex[edge.second]][qc.initialLayout.at(gate.target)];
					auto indexFT = x[k][physicalQubitIndex[edge.first]][qc.initialLayout.at(gate.target)];
					auto indexSC = x[k][physicalQubitIndex[edge.second]][qc.initialLayout.at(gate.control)];

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
		piCount=0;
		auto& i = x[k-1];
		auto& j = x[k];
		do {
			if (skipped_pi.count(piCount) == 0 || !this->settings.enableLimits){
				if ( skipped_pi_per_layer[k].count(piCount) == 0) {
					expr equal = c.bool_val(true);
					for (unsigned short Q: qubitChoice) {
						for (unsigned short q=0; q< qc.getNqubits(); ++q) {
							auto before = i[physicalQubitIndex[Q]][q];
							auto after = j[physicalQubitIndex[pi[physicalQubitIndex[Q]]]][q];
							equal = equal && (before == after);
						}
					}
					opt.add(implies(y[k-1][piCount],equal.simplify()).simplify());
				}
				++piCount;
			}
		} while(std::next_permutation(pi.begin(), pi.end()));
	}

	// Allow only 1 y_k_pi to be true
    if (this->settings.encoding == Encodings::None) {
		for (unsigned long k = 1; k < reducedLayerIndices.size(); ++k) {
			expr onlyOne = c.int_val(0);
			piCount = 0;
			do {
				if (skipped_pi.count(piCount) == 0 || !this->settings.enableLimits){
					if (skipped_pi_per_layer[k].count(piCount) == 0)
						onlyOne = onlyOne + ite(y[k - 1][piCount], c.int_val(1), c.int_val(0));
					++piCount;
				}
			} while (std::next_permutation(pi.begin(), pi.end()));
			opt.add(onlyOne.simplify() == 1);
		}
	} else {
		for (unsigned long k = 1; k < reducedLayerIndices.size(); ++k) {
			std::vector<expr> varIDs;
			piCount = 0;
			do {
				if (skipped_pi.count(piCount) == 0 || !this->settings.enableLimits){
					if (skipped_pi_per_layer[k].count(piCount) == 0)
						varIDs.push_back(y[k - 1][piCount]);
					++piCount;
				}
			} while (std::next_permutation(pi.begin(), pi.end()));
			if (this->settings.grouping == Groupings::Fixed2) {
				opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, 2), -1, auxvars, c));
			} else if (this->settings.grouping == Groupings::Fixed3) {
				opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, 3), -1, auxvars, c));
			} else if (this->settings.grouping == Groupings::Logarithm) {
				opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, static_cast<std::size_t>(std::log(varIDs.size()))), -1, auxvars, c));
			} else if (this->settings.grouping == Groupings::Halves) {
				opt.add(ExactlyOneCMDR(varIDs, groupVars(varIDs, varIDs.size() / 2), -1, auxvars, c));
			}

		}
	}
	//////////////////////////////////////////
	/// 	Objective Function				//
	//////////////////////////////////////////
	// cost for permutations
	piCount=0;
	std::vector<std::set<WeightedVar>> weightedVars(reducedLayerIndices.size());
	do {
			if (skipped_pi.count(piCount) == 0 || !this->settings.enableLimits){
				auto picost = architecture.minimumNumberOfSwaps(pi);
				if (architecture.bidirectional()) {
					picost *= GATES_OF_BIDIRECTIONAL_SWAP;
				} else {
					picost *= GATES_OF_UNIDIRECTIONAL_SWAP;
				}
				for (unsigned long k = 1; k < reducedLayerIndices.size(); ++k) {
					if (!this->settings.enableLimits || skipped_pi_per_layer[k].count(piCount) == 0) {
						opt.add(!y[k-1][piCount], picost);
						if (this->settings.useBDD)
							weightedVars[k].insert(WeightedVar(piCount, static_cast<int>(picost)));
					}
				}
				++piCount;
			}
	} while(std::next_permutation(pi.begin(), pi.end()));

    if (this->settings.enableLimits && this->settings.useBDD) {
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
					auto indexFT = x[k][physicalQubitIndex[edge.first]][qc.initialLayout.at(gate.target)];
					auto indexSC = x[k][physicalQubitIndex[edge.second]][qc.initialLayout.at(gate.control)];
					reverse = reverse && (!indexFT || !indexSC);
				}
				opt.add(reverse.simplify(), GATES_OF_DIRECTION_REVERSE);
			}
		}
	}

	//////////////////////////////////////////
	/// 	Solving							//
	//////////////////////////////////////////
	if (sat == opt.check()) {
		model m = opt.get_model();
		choiceResults.timeout = results.timeout = false;

		if (settings.verbose) {
			std::cout << "-------- qubit choice: ";
			for (const auto Q: qubitChoice) {
				std::cout << Q << " ";
			}
			std::cout << "----------" << std::endl;
		}

		// quickly determine cost
		choiceResults.output_singlequbitgates = choiceResults.input_singlequbitgates;
		choiceResults.output_cnots = choiceResults.input_cnots;
		choiceResults.output_gates = choiceResults.output_singlequbitgates + choiceResults.output_cnots;

		// swaps
		for (unsigned long k = 1; k < reducedLayerIndices.size(); ++k) {
			auto& i = x[k-1];
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
			choiceResults.output_swaps += swaps.at(k).size();
			if (architecture.bidirectional()) {
				choiceResults.output_gates += GATES_OF_BIDIRECTIONAL_SWAP * swaps.at(k).size();
			} else {
				choiceResults.output_gates += GATES_OF_UNIDIRECTIONAL_SWAP * swaps.at(k).size();
			}
		}

		// direction reverse
		if (!architecture.bidirectional()) {
			for (unsigned long k = 0; k < reducedLayerIndices.size(); ++k) {
				for (const auto& gate: layers.at(reducedLayerIndices.at(k))) {
					if (gate.singleQubit())
						continue;
					for (const auto& edge: rcm) {
						auto indexFT = x[k][physicalQubitIndex[edge.first]][qc.initialLayout.at(gate.target)];
						auto indexSC = x[k][physicalQubitIndex[edge.second]][qc.initialLayout.at(gate.control)];
						if (eq(m.eval(indexFT && indexSC), c.bool_val(true))) {
							choiceResults.output_direction_reverse++;
							choiceResults.output_gates += GATES_OF_DIRECTION_REVERSE;
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
					swaps.at(0).emplace_back(std::pair<unsigned short, unsigned short>{Q,q});
				}
			}
		}

	}
	else {
		results.timeout = true;
	}
}

void ExactMapper::initResults() {
	Mapper::initResults();
	results.method = Method::Exact;
	results.output_gates = std::numeric_limits<unsigned long>::max();
	results.encoding = this->settings.encoding;
	results.grouping = this->settings.grouping;
	results.strategy = this->settings.strategy;
	results.limit = this->settings.limit;
}


