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

#include <iostream>
#include <sstream>
#include <locale>

#include "exact/ExactMapper.hpp"

/*
int exact_mapping(const std::string& filename) {

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
 */

int main(int argc, char** argv) {
	std::string circuit = "../../examples/3_17_13.qasm";
	std::string cm = "../../extern/architectures/ibmq_london.arch";
	std::string cal = "../../extern/calibration/ibmq_london.csv";

	ExactMapper mapper{circuit, cm, cal};

	MappingSettings ms{};
	ms.layeringStrategy = DisjointQubits;
	mapper.map(ms);

	mapper.printResult(std::cout);

	mapper.dumpResult("3_17_13_mapped.qasm");

	return 0;

	// TODO: Write command line app using boost program options
}
