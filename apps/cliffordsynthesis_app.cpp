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
        "initialTimesteps,t", po::value<int>(),
        "Initial timesteps for the generated circuit (Depth for "
        "Depth-Synthesis, Gates for Gate-Synthesis)\n\t\t Sensible Values "
        "are for Depth: nQubit+log(nQubit) For Gates: nQubits*log(nQubits)")(
        "strategy,s", po::value<std::string>(),
        "choose one of use_minimizer, start_high, start_low, minmax, split_iter")(
        "target,r", po::value<std::string>(),
        R"(choose one metric to optimize ("gates" | "gates_only_cnot" | "depth" | "fidelity"))")(
        "method,m", po::value<std::string>(),
        R"(choose method used to solve ("z3" | "optimath" | "smtlibv2" | "dimacs"))")(
        "verbosity,v", po::value<int>(),
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
"(default: false)")("string", "Use String representation from Qiskit as input"
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

    CliffordSynthesizer opt = CliffordSynthesizer();
    if (vm.count("logfile")) {
        std::string logfile = vm["logfile"].as<std::string>();
        util::init(logfile);
    } else {
        util::init();
    }

    if (vm.count("choosebest")) {
        opt.chooseBest = true;
    }

    if (vm.count("useembed")) {
        opt.useEmbedding = true;
    }

    if (vm.count("verbosity")) {
        opt.verbosity = vm["verbosity"].as<int>();
        opt.verbosity = 5;
        switch (opt.verbosity) {
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
                opt.verbosity = 5;
                break;
        }
    }

    if (vm.count("method")) {
        const std::string method = vm["method"].as<std::string>();
        if (method == "z3") {
            opt.method = SynthesisMethod::Z3;
        } else if (method == "optimath") {
            opt.method = SynthesisMethod::MATHSAT;
        } else if (method == "smtlibv2") {
            opt.method = SynthesisMethod::SMTLibV2;
        } else if (method == "dimacs") {
            opt.method = SynthesisMethod::DIMACS;
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
        architecture.loadProperties(fid);
    }
    opt.setArchitecture(architecture);

    if (vm.count("target")) {
        const std::string target = vm["target"].as<std::string>();
        if (target == "gates") {
            opt.target = SynthesisTarget::GATES;
        } else if (target == "gates_only_cnot") {
            opt.target = SynthesisTarget::GATES_ONLY_CNOT;
        } else if (target == "depth") {
            opt.target = SynthesisTarget::DEPTH;
        } else if (target == "fidelity") {
            opt.target = SynthesisTarget::FIDELITY;
        } else {
            ERROR() << "Unknown target: " << target << std::endl;
            std::exit(1);
        }
    }

    if (!vm.count("testing")) {
        qc::QuantumComputation qc{};
        if (!vm.count("string")) {
            try {
                const std::string circuit = vm["in"].as<std::string>();
                if (circuit.substr(circuit.find_last_of('.') + 1) == "tabl") {
                    opt.targetTableau.import(circuit);
                    opt.nqubits = opt.targetTableau.getQubitCount();
                } else {
                    qc.import(circuit);
                    opt.nqubits = qc.getNqubits();
                    opt.circuit = qc.clone();
                    Tableau::generateTableau(opt.targetTableau, opt.circuit);
                }
            } catch (std::exception const& e) {
                ERROR() << "Could not import file: " << e.what() << std::endl;
                std::exit(1);
            }
        } else {
            const std::string tableau = vm["in"].as<std::string>();
            opt.targetTableau.fromString(tableau);
            opt.nqubits = opt.targetTableau.getQubitCount();
        }

        Tableau::initTableau(opt.initialTableau, opt.nqubits);
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

        if (opt.verbosity >= 5) {
            rnd.dumpOpenQASM(std::cout);
        }
        if (opt.verbosity >= 2) {
            rnd.printStatistics(std::cout);
        }
        Tableau::generateTableau(opt.targetTableau, rnd);
        Tableau::initTableau(opt.initialTableau, opt.nqubits);
        opt.circuit = rnd.clone();
    }

    if (vm.count("initialTimesteps")) {
        int initialTimesteps = vm["initialTimesteps"].as<int>();
        opt.initialTimesteps = initialTimesteps;
    } else {
        opt.initialTimesteps = 4 * (opt.nqubits + log(opt.nqubits));
    }

    if (vm.count("strategy")) {
        const std::string strategy = vm["strategy"].as<std::string>();
        if (strategy == "start_high") {
            opt.strategy = SynthesisStrategy::StartHigh;
        } else if (strategy == "start_low") {
            opt.strategy = SynthesisStrategy::StartLow;
        } else if (strategy == "minmax") {
            opt.strategy = SynthesisStrategy::MinMax;
        } else if (strategy == "split_iter") {
            opt.strategy = SynthesisStrategy::SplitIter;

        } else {
            opt.strategy = SynthesisStrategy::UseMinimizer;
        }
    }
    if (vm.count("nthread")) {
        int nthreads = vm["nthread"].as<int>();
        opt.nthreads = nthreads;
    }
    opt.optimize();
    Tableau resultTableau{};
    Tableau::generateTableau(resultTableau, opt.optimalResults.resultCircuit);
    if (opt.verbosity >= 2) {
        DEBUG() << "TargetTableau:" << std::endl
                << opt.targetTableau
                << "ResultTableau:" << std::endl
                << resultTableau << std::endl
                << "Used Gates:" << opt.optimalResults.gateCount << std::endl
                << "Depth: " << opt.optimalResults.depth << std::endl
                << "Fidelity: " << opt.optimalResults.fidelity << std::endl;
    }
    DEBUG() << "ResultTableau-Equality: "
            << (opt.targetTableau == resultTableau ?
                        "true" :
                        "false")
            << std::endl;
    if ((vm.count("testing") || opt.verbosity >= 2) &&
        opt.optimalResults.gateCount > 0) {
        opt.dumpResult(std::cout, qc::Format::OpenQASM);
    }
    if (vm.count("out") && opt.optimalResults.gateCount > 0) {
        const std::string out = vm["out"].as<std::string>();
        opt.dumpResult(out, qc::Format::OpenQASM);
    }
    if (vm.count("stats")) {
        const std::string stats = vm["stats"].as<std::string>();
        std::ofstream     ofs(stats);
        opt.optimalResults.dump(ofs);
    } else {
        opt.optimalResults.dump(std::cout);
    }
}
