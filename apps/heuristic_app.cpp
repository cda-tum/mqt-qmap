/*
 * This file is part of the MQT QMAP library which is released under the MIT license.
 * See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
 */

#include "heuristic/HeuristicMapper.hpp"

#include <boost/program_options.hpp>
#include <cmath>
#include <iostream>

int main(int argc, char** argv) {
    namespace po = boost::program_options;
    // clang-format off
    po::options_description description("MQT QMAP heuristic mapper by https://iic.jku.at/eda/quantum -- Options");
    description.add_options()
            ("help,h", "produce help message")
            ("in", po::value<std::string>()->required(), "File to read from")
            ("out", po::value<std::string>()->required(), "File to write to")
            ("arch", po::value<std::string>()->required(), "Architecture to use (points to a file)")
            ("calibration", po::value<std::string>(), "Calibration to use (points to a file)")
            ("initial_layout", po::value<std::string>(), R"(Initial layout strategy ("identity" | "static" | "dynamic"))")
            ("layering", po::value<std::string>(), R"(Layering strategy ("individual" | "disjoint"))")
            ("teleportation", po::value<unsigned long long int>()->implicit_value(0), "Use teleportation with optionally specifying the seed for the RNG used for initial placement")
            ("teleportation_fake", "Assign qubits as ancillary for teleportation in the initial placement but don't actually use them")
            ("verbose", "Increase verbosity and output additional information to stderr")
            ;
    // clang-format on
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, description), vm);
        if (vm.count("help")) {
            std::cout << description;
            return 0;
        }
        po::notify(vm);
    } catch (const po::error& e) {
        std::cerr << "[ERROR] " << e.what() << "! Try option '--help' for available commandline options.\n";
        std::exit(1);
    }

    const std::string      circuit = vm["in"].as<std::string>();
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
    Architecture      arch{};
    try {
        try {
            auto available = architectureFromString(cm);
            arch.loadCouplingMap(available);
        } catch (std::exception const& e) {
            arch.loadCouplingMap(cm);
        }
    } catch (std::exception const& e) {
        std::stringstream ss{};
        ss << "Could not import coupling map: " << e.what();
        std::cerr << ss.str() << std::endl;
        std::exit(1);
    }

    if (vm.count("calibration")) {
        const std::string cal = vm["calibration"].as<std::string>();
        try {
            arch.loadProperties(cal);
        } catch (std::exception const& e) {
            std::stringstream ss{};
            ss << "Could not import calibration data: " << e.what();
            std::cerr << ss.str() << std::endl;
            std::exit(1);
        }
    }

    HeuristicMapper mapper(qc, arch);

    Configuration ms{};
    ms.layering = Layering::IndividualGates;
    if (vm.count("layering")) {
        std::string layering = vm["layering"].as<std::string>();
        ms.layering          = layeringFromString(layering);
    }

    ms.initialLayout = InitialLayout::Dynamic;
    if (vm.count("initial_layout")) {
        std::string initialLayout = vm["initial_layout"].as<std::string>();
        ms.initialLayout          = initialLayoutFromString(initialLayout);
    }

    ms.verbose = vm.count("verbose") > 0;

    if (vm.count("teleportation")) {
        ms.teleportationQubits = std::min((arch.getNqubits() - qc.getNqubits()) & ~1u, 8u);
        ms.teleportationSeed   = vm["teleportation"].as<unsigned long long int>();
        ms.teleportationFake   = vm.count("teleportationFake") > 0;
    }

    mapper.map(ms);

    mapper.dumpResult(vm["out"].as<std::string>());

    mapper.printResult(std::cout);
}
