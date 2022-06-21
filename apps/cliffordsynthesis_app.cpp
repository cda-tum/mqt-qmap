/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
*/

#include "LogicBlock/LogicBlock.hpp"
#include "algorithms/RandomCliffordCircuit.hpp"
#include "cliffordsynthesis/CliffordSynthesizer.hpp"
#include "utils/logging.hpp"

#include <boost/program_options.hpp>
#include <iostream>
#include <locale>

int main(int argc, char** argv) {
    namespace po = boost::program_options;
    po::options_description description("Clifford-Optimizer -- Options");
    // clang-format off
    description.add_options()("help,h", "produce help message")(
        "in,i", po::value<std::string>(), "File to read from")(
        "out,o", po::value<std::string>(), "File to write to")(
        "stats", po::value<std::string>(), "File to write statistics to")(
        "arch,a", po::value<std::string>(),
        "Architecture that the circuit should be executed on/mapped to")(
        "fidelity,f", po::value<std::string>(),
        "Fidelities of the architectures")(
        "initial_timesteps,t", po::value<int>(),
        "Initial timesteps for the generated circuit (Depth for "
        "Depth-Optimization, Gates for Gate-Optimization)\n\t\t Sensible Values "
        "are for Depth: nQubit+log(nQubit) For Gates: nQubits*log(nQubits)")(
        "strategy,s", po::value<std::string>(),
        "choose one of use_minimizer, start_high, start_low, minmax, split_iter")(
        "target,r", po::value<std::string>(),
        R"(choose one metric to optimize ("gates" | "depth" | "fidelity"))")(
        "method,m", po::value<std::string>(),
        R"(choose method used to solve ("z3" | "optimath" | "smtlibv2" | "dimacs"))")(
        "verbose,v", po::value<int>(),
        "print more information")("testing", "toggle switch for testing mode")(
        "qubits", po::value<int>(), "qubits for test circuit generation")(
        "seed", po::value<int>(),
        "seed for test circuit generation, default 0 chooses randomly")(
        "circ_depth", po::value<int>(), "Max circ depth for testing mode")(
        "logfile", po::value<std::string>(),
        "path to a file (supports %N for logfile rotation), or 'std'")(
        "nthread", po::value<int>(),
        "maximum number of threads for use in split_iter strategy (default: 1)")(
        "choosebest", "only choose subset of coupling map with best fidelities "
                    "(default: false)")("useembed", "useembed"
                                                    "(default: false)");
    // clang-format on

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, description), vm);
        if (vm.count("help")) {
            std::cout << "";
            return 0;
        }
        po::notify(vm);
    } catch (const po::error& e) {
        ERROR() << e.what()
                << "! Try option '--help' for available commandline options.\n";
        std::exit(1);
    }

    CliffordOptimizer opt{};
    if (vm.count("logfile")) {
        std::string logfile = vm["logfile"].as<std::string>();
        util::init(logfile);
    } else {
        util::init();
    }

    if (vm.count("choosebest")) {
        opt.choose_best = true;
    }

    if (vm.count("useembed")) {
        opt.use_embedding = true;
    }

    if (vm.count("verbose")) {
        opt.verbose = vm["verbose"].as<int>();
        switch (opt.verbose) {
            case 0:
                std::cout << "Verbosity: Error\n";
                break;
            case 1:
                std::cout << "Verbosity: Warning\n";
                break;
            case 2:
                std::cout << "Verbosity: Info\n";
                break;
            case 3:
                std::cout << "Verbosity: Debug\n";
                break;
            case 4:
                std::cout << "Verbosity: Trace\n";
                break;
            default:
                std::cout << "Verbosity: Error\n";
                break;
        }
    }

    if (vm.count("method")) {
        const std::string method = vm["method"].as<std::string>();
        if (method == "z3") {
            opt.method = OptMethod::Z3;
        } else if (method == "optimath") {
            opt.method = OptMethod::MATHSAT;
        } else if (method == "smtlibv2") {
            opt.method = OptMethod::SMTLibV2;
        } else if (method == "dimacs") {
            opt.method = OptMethod::DIMACS;
        } else {
            ERROR() << "[ERROR] Unknown method '" << method << "'!\n";
            std::exit(1);
        }
    }

    Architecture architecture{};
    if (vm.count("arch")) {
        const std::string cm = vm["arch"].as<std::string>();
        try {
            architecture.loadCouplingMap(cm);
        } catch (std::exception const& e) {
            ERROR() << "Could not import coupling map: " << e.what() << std::endl;
            std::exit(1);
        }
    }

    if (vm.count("fidelity")) {
        const std::string fid = vm["fidelity"].as<std::string>();
        architecture.loadCalibrationData(fid);
    }
    opt.setCouplingMap(architecture);

    if (vm.count("target")) {
        const std::string target = vm["target"].as<std::string>();
        if (target == "gates") {
            opt.target = OptTarget::GATES;
        } else if (target == "depth") {
            opt.target = OptTarget::DEPTH;
        } else if (target == "fidelity") {
            opt.target = OptTarget::FIDELITY;
        } else {
            ERROR() << "Unknown target: " << target << std::endl;
            std::exit(1);
        }
    }

    if (!vm.count("testing")) {
        qc::QuantumComputation qc{};
        try {
            const std::string circuit = vm["in"].as<std::string>();
            qc.import(circuit);
            opt.nqubits = qc.getNqubits();
            opt.circuit = qc.clone();
        } catch (std::exception const& e) {
            ERROR() << "Could not import circuit: " << e.what() << std::endl;
            std::exit(1);
        }

        opt.generateTableau(opt.targetTableau, opt.circuit);
        opt.initTableau(opt.initialTableau);
    } else {
        if (vm.count("qubits")) {
            int qubits  = vm["qubits"].as<int>();
            opt.nqubits = qubits;
        } else {
            opt.nqubits = 10;
        }
        int circ_depth = 5;
        if (vm.count("circ_depth")) {
            circ_depth = vm["circ_depth"].as<int>();
        }
        int seed = 0;
        if (vm.count("seed")) {
            seed = vm["seed"].as<int>();
        }
        qc::RandomCliffordCircuit rnd(opt.nqubits, circ_depth, seed);

        if (opt.verbose >= 5) {
            rnd.dumpOpenQASM(std::cout);
        }
        if (opt.verbose >= 2) {
            rnd.printStatistics(std::cout);
        }
        opt.generateTableau(opt.targetTableau, rnd);
        opt.initTableau(opt.initialTableau);
        opt.circuit        = rnd.clone();
    }

    if (vm.count("initial_timesteps")) {
        int initial_timesteps = vm["initial_timesteps"].as<int>();
        opt.initial_timesteps = initial_timesteps;
    } else {
        opt.initial_timesteps = 4 * (opt.nqubits + log(opt.nqubits));
    }

    if (vm.count("strategy")) {
        const std::string strategy = vm["strategy"].as<std::string>();
        if (strategy == "start_high") {
            opt.strategy = OptimizingStrategy::StartHigh;
        } else if (strategy == "start_low") {
            opt.strategy = OptimizingStrategy::StartLow;
        } else if (strategy == "minmax") {
            opt.strategy = OptimizingStrategy::MinMax;
        } else if (strategy == "split_iter") {
            opt.strategy = OptimizingStrategy::SplitIter;

        } else {
            opt.strategy = OptimizingStrategy::UseMinimizer;
        }
    }
    if (vm.count("nthread")) {
        int nthreads = vm["nthread"].as<int>();
        opt.nthreads = nthreads;
    }
    opt.optimize();
    Tableau resultTableau{};
    opt.generateTableau(resultTableau, opt.circuit);
    if (opt.verbose >= 2) {
        DEBUG() << "TargetTableau:" << std::endl
                << opt.targetTableau
                << "ResultTableau:" << std::endl
                << resultTableau << std::endl
                << "Used Gates:" << opt.optimal_results.gate_count << std::endl
                << "Depth: " << opt.optimal_results.depth << std::endl
                << "Fidelity: " << opt.optimal_results.fidelity << std::endl;
    }
    DEBUG() << "ResultTableau-Equality: "
            << (opt.targetTableau ==  resultTableau?
                        "true" :
                        "false")
            << std::endl;
    if ((vm.count("testing") || opt.verbose >= 2) &&
        opt.optimal_results.gate_count > 0) {
        opt.dumpResult(std::cout, qc::Format::OpenQASM);
    }
    if (vm.count("out") && opt.optimal_results.gate_count > 0) {
        const std::string out = vm["out"].as<std::string>();
        opt.dumpResult(out, qc::Format::OpenQASM);
    }
    if (vm.count("stats")) {
        const std::string stats = vm["stats"].as<std::string>();
        std::ofstream     ofs(stats);
        opt.optimal_results.dump(ofs);
    } else {
        opt.optimal_results.dump(std::cout);
    }
}
