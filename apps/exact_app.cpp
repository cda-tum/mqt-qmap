/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#include <iostream>
#include <locale>
#include <boost/program_options.hpp>

#include "exact/ExactMapper.hpp"

int main(int argc, char** argv) {
    namespace po = boost::program_options;
    po::options_description description("JKQ QMAP exact mapper by https://iic.jku.at/eda/quantum -- Options");
    description.add_options()
            ("help,h", "produce help message")
            ("in", po::value<std::string>()->required(), "File to read from")
            ("out", po::value<std::string>()->required(), "File to write to")
            ("arch", po::value<std::string>()->required(), "Architecture to use (points to a file)")
            ("calibration", po::value<std::string>(), "Calibration to use (points to a file)")
            ("layering", po::value<std::string>(), R"(Layering strategy ("individual" | "disjoint" | "odd" | "triangle"))")
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
	qc::QuantumComputation qc{};
	try {
		qc.import(circuit);
	} catch (std::exception const& e) {
		std::stringstream ss{};
		ss << "Could not import circuit: " << e.what();
		std::cerr << ss.str() << std::endl;
		std::exit(1);
	}
    const std::string cm = vm["arch"].as<std::string>();
	Architecture arch{};
	try {
		arch.loadCouplingMap(cm);
	} catch (std::exception const& e) {
		std::stringstream ss{};
		ss << "Could not import coupling map: " << e.what();
		std::cerr << ss.str() << std::endl;
		std::exit(1);
	}

	if (vm.count("calibration")) {
		const std::string cal = vm["calibration"].as<std::string>();
		try {
			arch.loadCalibrationData(cal);
		} catch (std::exception const& e) {
			std::stringstream ss{};
			ss << "Could not import calibration data: " << e.what();
			std::cerr << ss.str() << std::endl;
			std::exit(1);
		}
	}

	ExactMapper mapper(qc, arch);

    MappingSettings ms{};
    ms.initialLayoutStrategy = InitialLayoutStrategy::None;
	if (vm.count("layering")) {
		std::string layering = vm["layering"].as<std::string>();
		if (layering == "individual") {
			ms.layeringStrategy = LayeringStrategy::IndividualGates;
		} else if (layering == "disjoint") {
			ms.layeringStrategy = LayeringStrategy::DisjointQubits;
		} else if (layering == "odd") {
			ms.layeringStrategy = LayeringStrategy::OddGates;
		} else if (layering == "triangle") {
			ms.layeringStrategy = LayeringStrategy::QubitTriangle;
		} else {
			ms.layeringStrategy = LayeringStrategy::None;
		}
	}
    ms.verbose = vm.count("verbose") > 0;
    mapper.map(ms);

    mapper.dumpResult(vm["out"].as<std::string>());

    mapper.printResult(std::cout, vm.count("ps"));
}
