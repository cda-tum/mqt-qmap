#include <cmath>
#include <iostream>
#include <boost/program_options.hpp>

#include "heuristic/HeuristicMapper.hpp"


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

    HeuristicMapper* mapper = nullptr;
	if(vm.count("calibration")) {
        const std::string cal = vm["calibration"].as<std::string>();
        mapper = new HeuristicMapper{circuit, cm, cal};
	} else {
        mapper = new HeuristicMapper{circuit, cm};
	}



	MappingSettings ms{};
	ms.layeringStrategy = DisjointQubits;
	ms.verbose = vm.count("verbose") > 0;
	mapper->map(ms);

	mapper->dumpResult(vm["out"].as<std::string>());

	mapper->printResult(std::cout);
}
