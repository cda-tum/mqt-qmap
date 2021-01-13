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
	do {
		std::set<unsigned short> qubitChoice{};
		for (unsigned short i = 0; i < qc.getNqubits(); ++i) {
			qubitChoice.insert(qubits.at(i));
		}
		allPossibleQubitChoices.push_back(qubitChoice);
	} while (next_combination(qubits.begin(), qubits.begin()+qc.getNqubits(), qubits.end()));

	// 3) determine exact mapping for this qubit choice
	std::vector<std::vector<std::pair<unsigned short, unsigned short>>> swaps(reducedLayerIndices.size(), std::vector<std::pair<unsigned short, unsigned short>>{});
	mappingSwaps.reserve(reducedLayerIndices.size());
	for (auto& choice: allPossibleQubitChoices) {
		// reset swaps
		for(auto& layer: swaps) {
			layer.clear();
		}

		MappingResults choiceResults{};
		choiceResults.copyInput(results);

		// 4) reduce coupling map
		CouplingMap reducedCouplingMap = architecture.getCouplingMap();
		for (const auto& edge: architecture.getCouplingMap()) {
			if (!choice.count(edge.first) || !choice.count(edge.second)) {
				reducedCouplingMap.erase(edge);
			}
		}
		if (reducedCouplingMap.empty()) continue;

		// 5) Check if E_k is connected. If yes, then possible subset found
		std::set<unsigned short> reachedQubits{};
		reachedQubits.insert(*(choice.begin()));
		dfs(*(choice.begin()), reachedQubits, reducedCouplingMap);
		if (!(reachedQubits == choice)) continue;

		// 6) call actual mapping routine
		coreMappingRoutine(choice, reducedCouplingMap, choiceResults, swaps);

		if (settings.verbose) {
			std::cout << "SWAPs: " << choiceResults.output_swaps << std::endl;
			std::cout << "Direction reverses: " << choiceResults.output_direction_reverse << std::endl;
		}

		// 7) Check if new optimum found
		if (!choiceResults.timeout && choiceResults.output_gates < results.output_gates) {
			results = choiceResults;
			mappingSwaps = swaps;
		}
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
		locations.at(q) = q;
	}

	for(unsigned long i=0; i<layers.size(); ++i) {
		if (i == 0) {
			qc::permutationMap inverseInitialLayout{ };
			for (auto& pu: qc.initialLayout)
				inverseInitialLayout.insert({pu.second, pu.first});

			// no swaps but initial permutation
			for(const auto& assignment: *swapsIterator) {
				locations.at(inverseInitialLayout.at(assignment.second)) = assignment.first;
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
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(), reverse.first, qc::H);
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(), reverse.second, qc::H);
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(), qc::Control(reverse.first), reverse.second, qc::X);
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(), reverse.second, qc::H);
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(), reverse.first, qc::H);

				} else {
					if (settings.verbose) {
						std::cout << i << ": Added cnot with control and target: " << cnot.first << " " << cnot.second << std::endl;
					}
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(),
					                                             qc::Control(cnot.first), cnot.second, qc::X);
				}
			}
		}

		if (swapsIterator != mappingSwaps.end() && layerIterator !=  reducedLayerIndices.end() && i == *layerIterator) {
			// apply swaps before layer
			for (auto it=(*swapsIterator).rbegin(); it != (*swapsIterator).rend(); ++it) {
				auto& swap = *it;
				qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(), std::vector<qc::Control>{ }, swap.first, swap.second, qc::SWAP);
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

void ExactMapper::coreMappingRoutine(const std::set<unsigned short>& qubitChoice, const CouplingMap& rcm, MappingResults& choiceResults, std::vector<std::vector<std::pair<unsigned short, unsigned short>>>& swaps) {
	// Z3 context
	context c;

	std::vector<unsigned short> pi(qubitChoice.begin(), qubitChoice.end());
	unsigned long long piCount;
	std::unordered_map<unsigned short, unsigned short> physicalQubitIndex{};
	unsigned short qIdx=0;
	for (const auto& Q: qubitChoice){
		physicalQubitIndex[Q] = qIdx;
		++qIdx;
	}

	//////////////////////////////////////////
	/// 	Boolean Variable Definitions	//
	//////////////////////////////////////////
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
			y_name.str("");
			y_name << "y_" << k << '_' << piCount;
			y.back().push_back(c.bool_const(y_name.str().c_str()));
			++piCount;
		} while(std::next_permutation(pi.begin(), pi.end()));
	}

	// Z3 optimizer
	optimize opt(c);
	params p(c);
	p.set("timeout", settings.timeout);
	p.set("pb.compile_equality", true);
	p.set("maxres.hill_climb", true);
	p.set("maxres.pivot_on_correction_set", false);
	opt.set(p);

	//////////////////////////////////////////
	/// 	Consistency Constraints			//
	//////////////////////////////////////////
	for (unsigned long k = 0; k < reducedLayerIndices.size(); ++k) {
		for (unsigned long i=0; i<qubitChoice.size(); ++i) {
			expr rowConsistency = c.int_val(0);
			for (unsigned short j=0; j< qc.getNqubits(); ++j) {
				rowConsistency = rowConsistency +
				                 ite(x[k][i][j], c.int_val(1), c.int_val(0));
			}
			opt.add(rowConsistency.simplify() <= 1);
		}

		for (unsigned short j=0; j< qc.getNqubits(); ++j) {
			expr colConsistency = c.int_val(0);
			for (unsigned long i=0; i<qubitChoice.size(); ++i) {
				colConsistency = colConsistency +
				                 ite(x[k][i][j], c.int_val(1), c.int_val(0));
			}
			opt.add(colConsistency.simplify() == 1);
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
			expr equal = c.bool_val(true);
			for (unsigned short Q: qubitChoice) {
				for (unsigned short q=0; q< qc.getNqubits(); ++q) {
					auto before = i[physicalQubitIndex[Q]][q];
					auto after = j[physicalQubitIndex[pi[physicalQubitIndex[Q]]]][q];
					equal = equal && (before == after);
				}
			}
			opt.add(implies(y[k-1][piCount],equal.simplify()).simplify());
			++piCount;
		} while(std::next_permutation(pi.begin(), pi.end()));
	}

	// Allow only 1 y_k_pi to be true
	for (unsigned long k = 1; k < reducedLayerIndices.size(); ++k) {
		expr onlyOne = c.int_val(0);
		piCount=0;
		do {
			onlyOne = onlyOne + ite(y[k-1][piCount], c.int_val(1), c.int_val(0));
			++piCount;
		} while(std::next_permutation(pi.begin(), pi.end()));
		opt.add(onlyOne.simplify() == 1);
	}

	//////////////////////////////////////////
	/// 	Objective Function				//
	//////////////////////////////////////////
	// cost for permutations
	piCount=0;
	do {
		auto picost = architecture.minimumNumberOfSwaps(pi);
		if (architecture.bidirectional()) {
			picost *= GATES_OF_BIDIRECTIONAL_SWAP;
		} else {
			picost *= GATES_OF_UNIDIRECTIONAL_SWAP;
		}
		for (unsigned long k = 1; k < reducedLayerIndices.size(); ++k) {
			opt.add(!y[k-1][piCount], picost);
		}
		++piCount;
	} while(std::next_permutation(pi.begin(), pi.end()));

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
}


