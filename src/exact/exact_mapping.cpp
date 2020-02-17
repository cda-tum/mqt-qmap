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

#include <vector>
#include <iostream>
#include <sstream>
#include "circuit.hpp"

/// Cost of permutations on IBM QX2/QX4 architecture.
/// \param pi permutation
/// \return cost of permutation (=7*numberOfSwaps)
static unsigned int costIBMQX4(std::vector<int>& pi) {
	std::stringstream pi_name;
	for (int i : pi) {
		pi_name << i;
	}
	const int COST_PER_SWAP = 7;
	switch (std::stoi(pi_name.str())) {
		case 1234:
			return 0;
		case 10234: case 21034: case 2134: case 1432: case 1324: case 1243:
			return 1 * COST_PER_SWAP;
		case 31204: case 41230: case 3214: case 4231:
			return 3 * COST_PER_SWAP;
		case 12034: case 1342: case 21304: case 2431: case 21430: case 2314:
		case 20134: case 1423: case 31024: case 4132: case 41032: case 3124:
			return 2 * COST_PER_SWAP;
		case 13204: case 14230: case 31240: case 3241:
		case 30214: case 40231: case 41203: case 4213:
			return 4 * COST_PER_SWAP;
		case 21043: case 2143: case 10243: case 10432: case 10324:
			return 2 * COST_PER_SWAP;
		case 42130: case 41320: case 23014: case 3412: case 32104: case 31402: case 24031: case 4321:
			return 4 * COST_PER_SWAP;
		case 43210: case 34201:
			return 6 * COST_PER_SWAP;
		case 12304: case 12430: case 20314: case 20431: case 31042: case 3142:
		case 41023: case 4123: case 24130: case 3421: case 41302: case 32014:
		case 30124: case 40132: case 13024: case 14032: case 21403: case 2413:
		case 21340: case 2341: case 42031: case 4312: case 31420: case 23104:
			return 3 * COST_PER_SWAP;
		case 13240: case 34210: case 30241:
		case 40213: case 43201: case 14203:
			return 5 * COST_PER_SWAP;
		case 12043: case 10342:
		case 20143: case 10423:
			return 3 * COST_PER_SWAP;
		case 24301: case 32401: case 23410: case 42310: case 13402: case 14320: case 32140: case 23041:
		case 34021: case 34102: case 43012: case 43120: case 30412: case 40321: case 42103: case 24013:
			return 5 * COST_PER_SWAP;
		case 12340: case 12403: case 13042: case 13420: case 14023: case 14302:
		case 20341: case 20413: case 23140: case 23401: case 24103: case 24310:
		case 30142: case 30421: case 32041: case 32410: case 34012: case 34120:
		case 40123: case 40312: case 42013: case 42301: case 43021: case 43102:
			return 4 * COST_PER_SWAP;

		default:
			std::cout << "Permutation with no associated cost: " << std::stoi(pi_name.str()) << std::endl;
			return 0;
	}
}

/// Coupling map of the IBM QX4 architecture
const CouplingMap ibmQX4 = {
		{ 1, 0 },
		{ 2, 0 },
		{ 2, 1 },
		{ 3, 2 },
		{ 3, 4 },
		{ 2, 4 }
};

int exact_mapping(std::string filename) {

	unsigned int timeout = 3600000; // 60min timeout
	std::vector<int> physicalQubits = { 0, 1, 2, 3, 4 };

	MappingSettings settings = MappingSettings();
	//std::string filename = "../testExample.qasm";
	
	std::cout << "################### Basis ###################" << std::endl;
	MappingResults basis = Circuit::run(filename, timeout, ibmQX4, physicalQubits, costIBMQX4, settings);

	std::cout << "################### Basis Reduced ###################" << std::endl;
	settings.considerQubitSubsets();
	MappingResults basisReduced = Circuit::run(filename, timeout, ibmQX4, physicalQubits, costIBMQX4, settings);

	std::cout << "################### Disjoint Qubits ###################" << std::endl;
	settings.useDisjointQubitsStrategy();
	MappingResults disjoint = Circuit::run(filename, timeout, ibmQX4, physicalQubits, costIBMQX4, settings);

	std::cout << "################### Odd Gates ###################" << std::endl;
	settings.useOddGatesStrategy();
	MappingResults otherGate = Circuit::run(filename, timeout, ibmQX4, physicalQubits, costIBMQX4, settings);

	std::cout << "################### Qubit Triangle ###################" << std::endl;
	settings.useQubitTriangleStrategy();
	MappingResults triangle = Circuit::run(filename, timeout, ibmQX4, physicalQubits, costIBMQX4, settings);

    return 0;
}
