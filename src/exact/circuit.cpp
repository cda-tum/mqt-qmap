/*
Minimal Mapping of Quantum Circuits to IBM QX Architectures by JKU Linz, Austria

Developer: Robert Wille, Lukas Burgholzer, Alwin Zulehner

For more information, please visit http://iic.jku.at/eda/research/ibm_qx_mapping

If you have any questions feel free to contact us using
robert.wille@jku.at, lukas.burgholzer@jku.at or alwin.zulehner@jku.at

If you use the compiler for your research, we would be thankful if you referred to it
by citing the following publication:

@inproceedings{wille2019mapping,
    title={Mapping Quantum Circuits to {IBM QX} Architectures Using the Minimal Number of {SWAP} and {H} Operations},
    author={Wille, Robert and Burgholzer, Lukas and Zulehner, Alwin},
    booktitle={Design Automation Conference},
    year={2019}
}
*/



#if Z3_INCLUDE_PATH

#include "circuit.hpp"

/// Static driver routine
/// \param filename QASM-file to read circuit from
/// \param timeout Timeout threshold
/// \param cm Coupling map
/// \param physicalQubits Available physical qubits
/// \param cost Cost function
/// \param settings Mapping settings
/// \return Mapping results
MappingResults Circuit::run(std::string& filename, unsigned int timeout, const CouplingMap& cm, std::vector<int>& physicalQubits, const std::function< unsigned int(std::vector<int>&)>& cost, MappingSettings settings) {
	MappingResults optimum = MappingResults();

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	unsigned long nrPhysicalQubits = physicalQubits.size();

	// 0) Parse input file
	QASMparser *parser = new QASMparser(filename);
	parser->Parse();

	// 1) Extract gates and layers
	// 1a) Compute list of non-unary gates
	auto layers = parser->getLayers();
	unsigned long nrLayers = layers.size();
	unsigned long ignoredUnaryGates = 0;
	std::vector<QASMparser::gate> gates;

	for (auto& layer: layers) {
		for (auto& gate: layer) {
			// ignore unary gates
			if (gate.control != -1) {
				gates.push_back(gate);
			} else {
				ignoredUnaryGates++;
			}
		}
	}
	unsigned long nrGates = gates.size();

	// 1b) create layers according to different criteria
	std::vector<std::vector<QASMparser::gate>> reducedLayers;

	if (settings.exactStrategy) {
		for (auto& gate: gates) {
			reducedLayers.emplace_back(std::vector<QASMparser::gate>());
			reducedLayers.back().push_back(gate);
		}

	} else if (settings.disjointQubitsStrategy) {
		// recompute layers
		auto qubitsInLayer = std::set<int>();
		reducedLayers.emplace_back(std::vector<QASMparser::gate>());
		for(auto& gate: gates) {
			if (qubitsInLayer.count(gate.control) == 0 && qubitsInLayer.count(gate.target) == 0) {
				qubitsInLayer.insert(gate.control);
				qubitsInLayer.insert(gate.target);
				reducedLayers.back().push_back(gate);
			} else {
				qubitsInLayer.clear();
				qubitsInLayer.insert(gate.control);
				qubitsInLayer.insert(gate.target);
				reducedLayers.emplace_back(std::vector<QASMparser::gate>());
				reducedLayers.back().push_back(gate);
			}
		}
	} else if (settings.oddGatesStrategy) {
		int g = 0;
		for (auto& gate: gates) {
			if (g%2 == 0) {
				reducedLayers.emplace_back(std::vector<QASMparser::gate>());
				reducedLayers[g/2].push_back(gate);
			} else {
				reducedLayers[g/2].push_back(gate);
			}
			g++;
		}
	} else if (settings.qubitTriangleStrategy) {
		auto qubitsInCurrentLayer = std::set<int>();
		reducedLayers.emplace_back(std::vector<QASMparser::gate>());
		for (auto& gate: gates) {
			// consider gate for current layer
			qubitsInCurrentLayer.insert(gate.control);
			qubitsInCurrentLayer.insert(gate.target);

			// add to current layer
			if (qubitsInCurrentLayer.size() <= 3) {
				reducedLayers.back().push_back(gate);
			}
			// need to start new layer
			else {
				reducedLayers.emplace_back(std::vector<QASMparser::gate>());
				reducedLayers.back().push_back(gate);
				qubitsInCurrentLayer.clear();
				qubitsInCurrentLayer.insert(gate.control);
				qubitsInCurrentLayer.insert(gate.target);
			}
		}
	}

    std::cout << "---------------- Layering -------------------" << std::endl;
    for(auto& layer: reducedLayers) {
        for (auto& gate: layer) {
            std::cout << "(" << gate.control << " " << gate.target << ") ";
        }
        std::cout << std::endl;
    }
    std::cout << "---------------------------------------------" << std::endl;

	unsigned long nrLayersReducedCurcuit = reducedLayers.size();
	delete parser;

	// 2) Check used qubits -> logicalQubits and n
	std::set<int> logicalQubits;
	unsigned long nrUsedQubits;
	for (auto gate: gates) {
		logicalQubits.insert(gate.control);
		logicalQubits.insert(gate.target);
	}
	nrUsedQubits = logicalQubits.size();

	optimum.nrLayersReducedCurcuit = nrLayersReducedCurcuit;
	optimum.nrGatesReducedCurcuit = nrGates;
	optimum.nrIgnoredUnaryGates = ignoredUnaryGates;
	optimum.nrGatesOriginalCurcuit = nrGates + ignoredUnaryGates;
	optimum.logicalQubits = logicalQubits;

	if (!settings.useMinimumSetOfQubits) {
		for(auto q: physicalQubits) logicalQubits.insert(q);
		nrUsedQubits = logicalQubits.size();
	}

	// 3) For all possibilities k (=m over n) to pick n qubits from m physical qubits
	// 3a) Obtain usedPhysicalQubits -> V_k
	std::vector<std::set<int>> allPossibleQubitChoices;
	do {
		std::set<int> qubitChoice;
		for (int i = 0; i < nrUsedQubits; ++i) {
			qubitChoice.insert(physicalQubits[i]);
		}
		allPossibleQubitChoices.push_back(qubitChoice);
	} while (next_combination(physicalQubits.begin(), physicalQubits.begin()+nrUsedQubits, physicalQubits.end()));

	for (auto& choice: allPossibleQubitChoices) {

		// 3b) Delete all edges in E with v not in V_k -> E_k
		CouplingMap reducedCouplingMap = cm;
		for (auto edge: reducedCouplingMap) {
			if (!choice.count(edge.first) || !choice.count(edge.second)) {
				reducedCouplingMap.erase(edge);
			}
		}

		// 3c) Check if E_k is connected. If yes, then possible subset found
		if (reducedCouplingMap.empty()) continue;
		std::set<int> reachedQubits;
		reachedQubits.insert(*choice.begin());
		dfs(*choice.begin(), reachedQubits, reducedCouplingMap);

		bool isConnected = (reachedQubits == choice);

		if (!isConnected) continue;

		// 3d) Apply mapping routine

		// only consider logical qubits which are really used
		logicalQubits.clear();
		for (auto gate: gates) {
			logicalQubits.insert(gate.control);
			logicalQubits.insert(gate.target);
		}
		// Construct circuit and run mapping routine
		Circuit reducedCurcuit = Circuit(logicalQubits, choice, physicalQubits, reducedLayers, timeout);
		MappingResults results = reducedCurcuit.mapping(reducedCouplingMap, cost);

		results.nrLayersReducedCurcuit = nrLayersReducedCurcuit;
		results.nrIgnoredUnaryGates = ignoredUnaryGates;
		results.nrGatesOriginalCurcuit = results.nrGatesReducedCurcuit + results.nrIgnoredUnaryGates;
		results.nrAllGatesMappedCurcuit = results.nrGatesMappedCurcuit + results.nrIgnoredUnaryGates;

		results.print();

		// 3e) Check if new optimum found
		if (!results.timeout && results.totalCost < optimum.totalCost)
			optimum = results;

	}

	// 4) Write best result and statistics
	optimum.printOptimum();

	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	auto duration = duration_cast<milliseconds>(t2 - t1).count();

	std::cout << "Computation Time: " << duration << "ms" << std::endl;

	return optimum;
}

/// Core mapping routine
/// \param cm Coupling map
/// \param cost Cost function
/// \return Mapping results
MappingResults Circuit::mapping(const CouplingMap& cm, const std::function< unsigned int(std::vector<int>&)>& cost) {

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	// Z3 context
	context c;

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
	expr_vector x(c);
	for (int k=0; k<nrLayers; k++) {
		for (int usedPhysicalQubit : usedPhysicalQubits) {
			for (int logicalQubit : logicalQubits) {
				std::stringstream x_name;
				x_name << "x_" << k << '_' << usedPhysicalQubit << '_' << logicalQubit;
				x.push_back(c.bool_const(x_name.str().c_str()));
			}
		}
	}

	/*
	 permutation variables y_k_pi
	 k	before layer k
	 pi	arbitrary permutation of the m qubits
	 number of variables: (|L|-1) * m!
	 */
	std::vector<int> pi;
	for (int usedPhysicalQubit: usedPhysicalQubits) {
		pi.push_back(usedPhysicalQubit);
	}

	expr_vector y(c);
	for (int k=1; k<nrLayers; k++) {
		do {
			std::stringstream y_name;
			y_name << "y_" << k << '_' << printPi(pi);
			y.push_back(c.bool_const(y_name.str().c_str()));
		} while(std::next_permutation(pi.begin(), pi.end()));
	}

	/*
	 direction reverse variables z_k_m
	 k	before layer k
	 m	gate m of layer k
	 number of variables: |G|
	 */
	expr_vector z(c);
	for (int k=0; k<nrLayers; k++) {
		for(int m=0; m<layers[k].size(); m++) {
			std::stringstream z_name;
			z_name << "d_" << k << "_" << m;
			z.push_back(c.bool_const(z_name.str().c_str()));
		}
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
	for (int k=0; k<nrLayers; k++) {
		for (int usedPhysicalQubit: usedPhysicalQubits) {
			expr rowConsistency = c.int_val(0);
			for (int logicalQubit: logicalQubits) {
				rowConsistency = rowConsistency +
						ite(x[idx(k, usedPhysicalQubit, logicalQubit)], c.int_val(1), c.int_val(0));
			}
			opt.add(rowConsistency.simplify() <= 1);
		}

		for (int logicalQubit: logicalQubits) {
			expr colConsistency = c.int_val(0);
			for (int usedPhysicalQubit: usedPhysicalQubits) {
				colConsistency = colConsistency +
						ite(x[idx(k, usedPhysicalQubit, logicalQubit)], c.int_val(1), c.int_val(0));
			}
			opt.add(colConsistency.simplify() == 1);
		}
	}


	//////////////////////////////////////////
	///		Coupling Constraints			//
	//////////////////////////////////////////
	for (int k=0; k<nrLayers; k++) {
		expr allCouplings = c.bool_val(true);
		for(int m=0; m<layers[k].size(); m++) {
			expr coupling = c.bool_val(false);
			for (auto edge: cm) {
				coupling = coupling or ((x[idx(k, edge.first, layers[k][m].control)] and
										 x[idx(k, edge.second, layers[k][m].target)]) or
										(x[idx(k, edge.first, layers[k][m].target)] and
										 x[idx(k, edge.second, layers[k][m].control)]));
			}
			allCouplings = allCouplings and coupling.simplify();
		}
		opt.add(allCouplings.simplify());
	}

	//////////////////////////////////////////
	/// 	Permutation Constraints			//
	//////////////////////////////////////////
	unsigned long n = factorial(nrUsedPhysicalQubits);
	for (int k=1; k<nrLayers; k++) {
		int piCnt=0;
		do {
			expr equal = c.bool_val(true);
			for (int usedPhysicalQubit: usedPhysicalQubits) {
				for (int logicalQubit: logicalQubits) {
					// find pi index
					int counti = 0;
					for (int q : usedPhysicalQubits) {
						if (q == usedPhysicalQubit) break;
						counti++;
					}

					equal = equal and (
							x[idx(k-1, usedPhysicalQubit, logicalQubit)]
							==
							x[idx(k, pi[counti], logicalQubit)]);
				}
			}
			opt.add(implies(y[(k-1)*n+piCnt],equal.simplify()).simplify());
			piCnt++;
		} while(std::next_permutation(pi.begin(), pi.end()));
	}

	// Allow only 1 y_k_pi to be true
	for (int k=1; k<nrLayers; k++) {
		expr onlyOne = c.int_val(0);
		int piCnt=0;
		do {
			onlyOne = onlyOne + ite(y[(k-1)*n+piCnt], c.int_val(1), c.int_val(0));
			piCnt++;
		} while(std::next_permutation(pi.begin(), pi.end()));
		opt.add(onlyOne.simplify() == 1);
	}

	//////////////////////////////////////////
	/// 	Direction Reverse Constraints	//
	//////////////////////////////////////////
	int g = 0;
	for (int k=0; k<nrLayers; k++) {
		for (int m = 0; m < layers[k].size(); m++) {
			expr reverse = c.bool_val(false);
			for(auto edge: cm) {
				reverse = reverse or (x[idx(k, edge.first, layers[k][m].target)] and
									  x[idx(k, edge.second, layers[k][m].control)]);
			}
			opt.add(z[g] == reverse.simplify());
			g++;
		}
	}

	//////////////////////////////////////////
	/// 	Objective Function				//
	//////////////////////////////////////////
    // cost for permutations
	for (int k=1; k<nrLayers; k++) {
		int piCnt=0;
		do {
			// augment permutation to all physical qubits to obtain correct costs
			std::vector<int> augmentedPi;
			for(auto q: physicalQubits) {
				// check if q is used and get its index
				int counti = 0;
				for (int usedPhysicalQubit : usedPhysicalQubits) {
					if (q == usedPhysicalQubit) break;
					counti++;
				}
				if (counti == nrUsedPhysicalQubits) {
					augmentedPi.push_back(q);
				} else {
					augmentedPi.push_back(pi[counti]);
				}
			}
			opt.add(not(y[(k-1)*n+piCnt]), cost(augmentedPi));
			piCnt++;
		} while(std::next_permutation(pi.begin(), pi.end()));
	}

    // cost for reversed directions
    int gateIdx = 0;
    for(const auto& layer: layers) {
        for (const auto &gate: layer) {
            opt.add(not(z[gateIdx]), 4);
            gateIdx++;
        }
    }

	//////////////////////////////////////////
	/// 	Solving							//
	//////////////////////////////////////////
	MappingResults results = MappingResults();
	results.logicalQubits = logicalQubits;
	results.usedPhysicalQubits = usedPhysicalQubits;
	results.nrLayersReducedCurcuit = nrLayers;
	results.nrGatesReducedCurcuit = nrGates;

	if (sat == opt.check()) {
		model m = opt.get_model();
		results.timeout = false;

		for (int k=0; k<nrLayers; k++) {
			results.X.emplace_back(std::vector<std::vector<int>>());
			for (int usedPhysicalQubit: usedPhysicalQubits) {
				results.X.back().emplace_back(std::vector<int>());
				for (int logicalQubit: logicalQubits) {
					results.X.back().back().push_back(eq(m.eval(x[idx(k, usedPhysicalQubit, logicalQubit)]), c.bool_val(true))? 1 : 0);
				}
			}
		}

        unsigned long swapCost = 0;
        for (int k=1; k<nrLayers; k++) {
			int piCnt=0;
			do {
                // augment permutation to all physical qubits to obtain correct costs
                std::vector<int> augmentedPi;
                for (auto q: physicalQubits) {
                    // check if q is used and get its index
                    int counti = 0;
                    for (int usedPhysicalQubit : usedPhysicalQubits) {
                        if (q == usedPhysicalQubit) break;
                        counti++;
                    }
                    if (counti == nrUsedPhysicalQubits) {
                        augmentedPi.push_back(q);
                    } else {
                        augmentedPi.push_back(pi[counti]);
                    }
                }
				if (eq(m.eval(y[(k-1)*n+piCnt]), c.bool_val(true))) {
					results.Y.push_back(pi);
                    swapCost += cost(augmentedPi);
				}
				piCnt++;
			} while(std::next_permutation(pi.begin(), pi.end()));
		}

		gateIdx = 0;
        unsigned long reverseCost = 0;
        for (const auto& layer: layers) {
			for (const auto& gate: layer) {
				results.Z.push_back((eq(m.eval(z[gateIdx]), c.bool_val(true))? 1 : 0));
                reverseCost += eq(m.eval(z[gateIdx]), c.bool_val(true)) ? 4 : 0;
				gateIdx++;
			}
		}

		results.reverseCost = reverseCost;
		results.swapCost = swapCost;
		results.totalCost = reverseCost + swapCost;
		results.nrGatesMappedCurcuit = nrGates + results.totalCost;
	}
	else {
		results.timeout = true;
	}

	return results;
}

/// Create a string representation of a given permutation
/// \param pi permutation
/// \return string representation of pi
std::string printPi(std::vector<int>& pi){
	if (std::is_sorted(pi.begin(),pi.end())) {
		return "( )";
	}

	std::stringstream perm;
	perm << '(';
	for (int i=0; i<pi.size()-1; i++) {
		perm << pi[i] << ',';
	}
	perm << pi[pi.size()-1] << ')';

	return perm.str();
}

/// Iterating routine through all combinations
/// \tparam Iterator iterator type
/// \param first iterator to beginning
/// \param k current iterator
/// \param last iterator to end
/// \return true if another combination was found
template<typename Iterator>
inline bool next_combination(const Iterator first, Iterator k, const Iterator last) {
	/* Credits: Thomas Draper */
	if ((first == last) || (first == k) || (last == k))
		return false;
	Iterator itr1 = first;
	Iterator itr2 = last;
	++itr1;
	if (last == itr1)
		return false;
	itr1 = last;
	--itr1;
	itr1 = k;
	--itr2;
	while (first != itr1) {
		if (*--itr1 < *itr2) {
			Iterator j = k;
			while (!(*itr1 < *j)) ++j;
			std::iter_swap(itr1, j);
			++itr1;
			++j;
			itr2 = k;
			std::rotate(itr1, j, last);
			while (last != j) {
				++j;
				++itr2;
			}
			std::rotate(k, itr2, last);
			return true;
		}
	}
	std::rotate(first, k, last);
	return false;
}

/// Simple depth-first-search implementation used to check whether a given subset of qubits is
/// connected on the given architecture
/// \param current index of current qubit
/// \param visited visited qubits
/// \param cm coupling map of architecture
void dfs(int current, std::set<int>& visited, CouplingMap& cm) {
	for (auto edge: cm) {
		if (edge.first == current) {
			if(!visited.count(edge.second)) {
				visited.insert(edge.second);
				dfs(edge.second, visited, cm);
			}
		} else if (edge.second == current) {
			if (!visited.count(edge.first)) {
				visited.insert(edge.first);
				dfs(edge.first, visited, cm);
			}
		}
	}
}


#endif /* Z3_INCLUDE_PATH */