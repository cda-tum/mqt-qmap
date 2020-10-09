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
#include <boost/program_options.hpp>

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
    namespace po = boost::program_options;
    po::options_description description("JKQ QMAP heuristic mapper by https://iic.jku.at/eda/ -- Options");
    description.add_options()
            ("help,h", "produce help message")
            ("in", po::value<std::string>()->required(), "File to read from")
            ("out", po::value<std::string>()->required(), "File to write to")
            ("arch", po::value<std::string>()->required(), "Architecture to use (points to a file)")
            ("calibration", po::value<std::string>(), "Calibration to use (points to a file)")
            ("ps", "print statistics")
            ("verbose", "Increase verbosity and output additional information to stderr")
            ;
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, description), vm);
        if (vm.count("help")) {
            std::cout << description;
            return 0;
        }
        po::notify(vm);
    } catch (const po::error &e) {
        std::cerr << "[ERROR] " << e.what() << "! Try option '--help' for available commandline options.\n";
        std::exit(1);
    }


    const std::string circuit = vm["in"].as<std::string>();
    const std::string cm = vm["arch"].as<std::string>();

    ExactMapper* mapper = nullptr;
    if(vm.count("calibration")) {
        const std::string cal = vm["calibration"].as<std::string>();
        mapper = new ExactMapper{circuit, cm, cal};
    } else {
        mapper = new ExactMapper{circuit, cm};
    }



    MappingSettings ms{};
    ms.layeringStrategy = DisjointQubits;
    ms.verbose = vm.count("verbose") > 0;
    mapper->map(ms);

    mapper->dumpResult(vm["out"].as<std::string>());

    mapper->printResult(std::cout);
}
