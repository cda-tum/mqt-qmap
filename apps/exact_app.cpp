/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#include "exact/ExactMapper.hpp"

#include <boost/program_options.hpp>
#include <iostream>
#include <locale>

int main(int argc, char** argv) {
    namespace po = boost::program_options;
    po::options_description description("JKQ QMAP exact mapper by https://iic.jku.at/eda/quantum -- Options");
    // clang-format off
    description.add_options()
            ("help,h", "produce help message")
            ("in", po::value<std::string>()->required(), "File to read from")
            ("out", po::value<std::string>()->required(), "File to write to")
            ("arch", po::value<std::string>()->required(), "Architecture to use (points to a file)")
            ("calibration", po::value<std::string>(), "Calibration to use (points to a file)")
            ("layering", po::value<std::string>(), R"(Layering strategy ("individual" | "disjoint" | "odd" | "triangle"))")
            ("verbose", "Increase verbosity and output additional information to stderr")
			("encoding", po::value<std::string>(), R"(Choose encoding for AMO and exactly one ("none" | "commander" | "bimander"))")
			("commander_grouping", po::value<std::string>(), R"(Choose method of grouping ("fixed2" | "fixed3" | "logarithm" | "halves"))")
			//("limitswaps", "Enable bdd for limiting swaps per layer")
			("use_bdd", "Choose to use BDDs instead of directly limiting the permutation variables")
			("swap_reduction", po::value<std::string>(), R"(Choose method of limiting the search space ("none" | "custom" | "coupling_limit" | "increasing"))")
			("swap_limit", po::value<std::string>(), "Set a custom limit for max swaps per layer, for increasing it sets the max swaps")
			("use_subsets", "Use qubit subsets, or consider all available physical qubits at once")
			("timeout", po::value<std::string>(), "timeout for the execution")
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
            arch.loadCalibrationData(cal);
        } catch (std::exception const& e) {
            std::stringstream ss{};
            ss << "Could not import calibration data: " << e.what();
            std::cerr << ss.str() << std::endl;
            std::exit(1);
        }
    }

    ExactMapper mapper(qc, arch);

    Configuration ms{};
    ms.initialLayout = InitialLayout::None;
    if (vm.count("layering")) {
        std::string layering = vm["layering"].as<std::string>();
        ms.layering          = layeringFromString(layering);
    }
    if (vm.count("encoding")) {
        const std::string encoding = vm["encoding"].as<std::string>();
        ms.encoding                = encodingFromString(encoding);
    }
    if (vm.count("commander_grouping")) {
        const std::string grouping = vm["commander_grouping"].as<std::string>();
        ms.commanderGrouping       = groupingFromString(grouping);
    }
    if (vm.count("swap_reduction")) {
        ms.enableSwapLimits = true;
        if (vm.count("use_bdd")) {
            ms.useBDD = true;
        }
        const std::string swapReduction = vm["swap_reduction"].as<std::string>();
        if (swapReduction == "custom") {
            ms.swapReduction = SwapReduction::Custom;
            if (vm.count("swap_limit")) {
                const std::string swap_limit = vm["swap_limit"].as<std::string>();
                ms.swapLimit                 = std::stoi(swap_limit);
            }
        } else if (swapReduction == "coupling_limit") {
            ms.swapReduction = SwapReduction::CouplingLimit;
        } else if (swapReduction == "increasing") {
            ms.swapReduction = SwapReduction::Increasing;
            if (vm.count("limit")) {
                const std::string swap_limit = vm["limit"].as<std::string>();
                ms.swapLimit                 = std::stoi(swap_limit);
            }
        } else {
            ms.swapReduction    = SwapReduction::None;
            ms.enableSwapLimits = false;
            ms.useBDD           = false;
        }
    }
    if (vm.count("timeout")) {
        const std::string timeout = vm["timeout"].as<std::string>();
        ms.setTimeout(std::stoi(timeout) * 1000);
    }
    ms.useSubsets = vm.count("use_subsets") > 0;
    ms.verbose    = vm.count("verbose") > 0;
    mapper.map(ms);

    mapper.dumpResult(vm["out"].as<std::string>());

    mapper.printResult(std::cout);
}
