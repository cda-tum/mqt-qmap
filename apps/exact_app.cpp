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
			("encoding", po::value<std::string>(), R"(Choose encoding for AMO and exactly one ("none" | "commander" | "bimander"))")
			("grouping", po::value<std::string>(), R"(Choose method of grouping ("fixed2" | "fixed3" | "logarithm" | "halves"))")
			("bdd", "Enable bdd for limiting swaps per layer")
			("bddStrategy", po::value<std::string>(), R"(Choose method of applying bdd limits ("none" | "custom" | "architectureswaps" | "subsetswaps" | "increasing"))")
			("bdd_limit", po::value<std::string>(), "Set a custom limit for max swaps per layer, for increasing it sets the max swaps")
			("useSubsets", "Use qubit subsets, or consider all available physical qubits at once")
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
	if (vm.count("encoding"))
	{
		const std::string encoding = vm["encoding"].as<std::string>();
		if (encoding == "none")
		{
			ms.encoding = Encodings::None;
		}
		else if (encoding == "commander")
		{
			ms.encoding = Encodings::Commander;
		}
		else if (encoding == "bimander")
		{
			ms.encoding = Encodings::Bimander;
		}
		else
		{
			ms.encoding = Encodings::None;
		}
	}
	if (vm.count("grouping"))
	{
		const std::string grouping = vm["grouping"].as<std::string>();
		if (grouping == "fixed3")
		{
			ms.grouping = Groupings::Fixed3;
		}
		else if (grouping == "fixed2")
		{
			ms.grouping = Groupings::Fixed2;
		}
		else if (grouping == "logarithm")
		{
			ms.grouping = Groupings::Logarithm;
		} else if (grouping == "halves"){
			ms.grouping = Groupings::Halves;
		}
	}
	if (vm.count("bdd"))
	{
		ms.enableBDDLimits = true;
		if (vm.count("bddStrategy")) {
			const std::string bddStrat = vm["bddStrategy"].as<std::string>();
			if (bddStrat == "custom") {
				ms.bddStrategy = BDDStrategy::Custom;
				if (vm.count("bdd_limit")) {
					const std::string bdd_limit = vm["bdd_limit"].as<std::string>();
					ms.bddLimit = std::stoi(bdd_limit.c_str());
				}
			} else if (bddStrat == "architectureswaps") {
				ms.bddStrategy = BDDStrategy::ArchitectureSwaps;
			} else if (bddStrat == "subsetswaps") {
				ms.bddStrategy = BDDStrategy::SubsetSwaps;
			} else if (bddStrat == "increasing") {
				ms.bddStrategy = BDDStrategy::Increasing;
				if (vm.count("bdd_limit")) {
					const std::string bdd_limit = vm["bdd_limit"].as<std::string>();
					ms.bddLimit = std::stoi(bdd_limit.c_str());
				}
			} else {
				ms.bddStrategy = BDDStrategy::None;
				ms.enableBDDLimits = false;
			}
		}
		
	}
	ms.useQubitSubsets = vm.count("useSubsets") > 0;
    ms.verbose = vm.count("verbose") > 0;
    mapper.map(ms);

    mapper.dumpResult(vm["out"].as<std::string>());

    mapper.printResult(std::cout, vm.count("ps"));
}
