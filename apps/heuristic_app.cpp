
#include <regex>
#include <cmath>
#include <iostream>
#include <fstream>

#include "heuristic/HeuristicMapper.hpp"


int main(int argc, char** argv) {
	std::string circuit = "../../examples/3_17_13.qasm";
	std::string cm = "../../extern/architectures/ibmq_london.arch";
	std::string cal = "../../extern/calibration/ibmq_london.csv";

	HeuristicMapper mapper{circuit, cm, cal};

	MappingSettings ms{};
	ms.layeringStrategy = DisjointQubits;
	mapper.map(ms);

	mapper.printResult(std::cout);

	mapper.dumpResult("3_17_13_mapped.qasm");

	return 0;

	// TODO: Write command line app using boost program options
}
