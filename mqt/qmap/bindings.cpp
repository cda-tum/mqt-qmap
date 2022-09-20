/*
* This file is part of MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/quantum/ for more information.
*/

#include "cliffordsynthesis/CliffordSynthesizer.hpp"
#include "exact/ExactMapper.hpp"
#include "heuristic/HeuristicMapper.hpp"
#include "nlohmann/json.hpp"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11_json/pybind11_json.hpp"
#include "qiskit/QuantumCircuit.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace pybind11::literals;

void loadQC(qc::QuantumComputation& qc, const py::object& circ) {
    try {
        if (py::isinstance<py::str>(circ)) {
            auto&& file = circ.cast<std::string>();
            qc.import(file);
        } else {
            qc::qiskit::QuantumCircuit::import(qc, circ);
        }
    } catch (std::exception const& e) {
        std::stringstream ss{};
        ss << "Could not import circuit: " << e.what();
        throw std::invalid_argument(ss.str());
    }
}

// c++ binding function
MappingResults map(const py::object& circ, Architecture& arch, Configuration& config) {
    qc::QuantumComputation qc{};

    loadQC(qc, circ);

    if (config.useTeleportation) {
        config.teleportationQubits = std::min((arch.getNqubits() - qc.getNqubits()) & ~1, 8);
    }

    std::unique_ptr<Mapper> mapper;
    try {
        if (config.method == Method::Heuristic) {
            mapper = std::make_unique<HeuristicMapper>(qc, arch);
        } else if (config.method == Method::Exact) {
            mapper = std::make_unique<ExactMapper>(qc, arch);
        }
    } catch (std::exception const& e) {
        std::stringstream ss{};
        ss << "Could not construct mapper: " << e.what();
        throw std::invalid_argument(ss.str());
    }

    try {
        mapper->map(config);
    } catch (std::exception const& e) {
        std::stringstream ss{};
        ss << "Error during mapping: " << e.what();
        throw std::invalid_argument(ss.str());
    }

    auto& results = mapper->getResults();

    std::stringstream qasm{};
    mapper->dumpResult(qasm, qc::OpenQASM);
    results.mappedCircuit = qasm.str();

    return results;
}

// c++ binding function
CliffordOptimizationResults optimize(const py::object& circ, Architecture& arch, SynthesisStrategy& strategy) {
    qc::QuantumComputation qc{};

    loadQC(qc, circ);

    std::unique_ptr<CliffordSynthesizer> optimizer;
    try {
        optimizer = std::make_unique<CliffordSynthesizer>(false, false, 0, 0, strategy, SynthesisTarget::GATES);
    } catch (std::exception const& e) {
        std::stringstream ss{};
        ss << "Could not construct optimizer: " << e.what();
        throw std::invalid_argument(ss.str());
    }

    try {
        optimizer->setCircuit(qc);
        optimizer->setArchitecture(arch);
    } catch (std::exception const& e) {
        std::stringstream ss{};
        ss << "Error during initialization: " << e.what();
        throw std::invalid_argument(ss.str());
    }

    try {
        optimizer->optimize();
    } catch (std::exception const& e) {
        std::stringstream ss{};
        ss << "Error during optimization: " << e.what();
        throw std::invalid_argument(ss.str());
    }

    optimizer->optimalResults.generateStringCircuit();

    return optimizer->optimalResults;
}

// c++ binding function
CliffordOptimizationResults synthesize(const std::string& tableau, Architecture& arch, SynthesisStrategy& strategy) {
    auto tab = Tableau();
    try {
        tab.fromString(tableau);
    } catch (std::exception const& e) {
        std::stringstream ss{};
        ss << "Could not parse tableau: " << e.what();
        throw std::invalid_argument(ss.str());
    }

    std::unique_ptr<CliffordSynthesizer> optimizer;
    try {
        optimizer = std::make_unique<CliffordSynthesizer>(false, false, 0, 0, strategy, SynthesisTarget::GATES);
    } catch (std::exception const& e) {
        std::stringstream ss{};
        ss << "Could not construct optimizer: " << e.what();
        throw std::invalid_argument(ss.str());
    }

    try {
        optimizer->setTargetTableau(tab);
        optimizer->setArchitecture(arch);
    } catch (std::exception const& e) {
        std::stringstream ss{};
        ss << "Error during initialization: " << e.what();
        throw std::invalid_argument(ss.str());
    }

    try {
        optimizer->optimize();
    } catch (std::exception const& e) {
        std::stringstream ss{};
        ss << "Error during optimization: " << e.what();
        throw std::invalid_argument(ss.str());
    }

    optimizer->optimalResults.generateStringCircuit();

    return optimizer->optimalResults;
}

PYBIND11_MODULE(pyqmap, m) {
    m.doc() = "pybind11 for the MQT QMAP quantum circuit mapping tool";

    // Pre-defined architecture available within QMAP
    py::enum_<AvailableArchitecture>(m, "Arch")
            .value("IBM_QX4", AvailableArchitecture::IBM_QX4,
                   "5 qubit, directed bow tie layout")
            .value("IBM_QX5", AvailableArchitecture::IBM_QX5,
                   "16 qubit, directed ladder layout")
            .value("IBMQ_Yorktown", AvailableArchitecture::IBMQ_Yorktown,
                   "5 qubit, undirected bow tie layout")
            .value("IBMQ_London", AvailableArchitecture::IBMQ_London,
                   "5 qubit, undirected T-shape layout")
            .value("IBMQ_Bogota", AvailableArchitecture::IBMQ_Bogota,
                   "5 qubit, undirected linear chain layout")
            .value("IBMQ_Casablanca", AvailableArchitecture::IBMQ_Casablanca,
                   "7 qubit, undirected H-shape layout")
            .value("IBMQ_Tokyo", AvailableArchitecture::IBMQ_Tokyo,
                   "20 qubit, undirected brick-like layout")
            .value("Rigetti_Agave", AvailableArchitecture::Rigetti_Agave,
                   "8 qubit, undirected ring layout")
            .value("Rigetti_Aspen", AvailableArchitecture::Rigetti_Aspen,
                   "16 qubit, undirected dumbbell layout")
            .export_values()
            // allow construction from string
            .def(py::init([](const std::string& str) -> AvailableArchitecture { return architectureFromString(str); }));

    // Mapping methodology to use
    py::enum_<Method>(m, "Method")
            .value("heuristic", Method::Heuristic)
            .value("exact", Method::Exact)
            .export_values()
            // allow construction from string
            .def(py::init([](const std::string& str) -> Method { return methodFromString(str); }));

    // Initial layout strategy
    py::enum_<InitialLayout>(m, "InitialLayout")
            .value("identity", InitialLayout::Identity)
            .value("static", InitialLayout::Static)
            .value("dynamic", InitialLayout::Dynamic)
            .export_values()
            // allow construction from string
            .def(py::init([](const std::string& str) -> InitialLayout { return initialLayoutFromString(str); }));

    // Gate clustering / layering strategy
    py::enum_<Layering>(m, "Layering")
            .value("individual_gates", Layering::IndividualGates)
            .value("disjoint_qubits", Layering::DisjointQubits)
            .value("odd_gates", Layering::OddGates)
            .value("qubit_triangle", Layering::QubitTriangle)
            .export_values()
            // allow construction from string
            .def(py::init([](const std::string& str) -> Layering { return layeringFromString(str); }));

    // Encoding settings for at-most-one and exactly-one constraints
    py::enum_<Encoding>(m, "Encoding")
            .value("naive", Encoding::Naive)
            .value("commander", Encoding::Commander)
            .value("bimander", Encoding::Bimander)
            .export_values()
            // allow construction from string
            .def(py::init([](const std::string& str) -> Encoding { return encodingFromString(str); }));

    // Grouping settings if using the commander encoding
    py::enum_<CommanderGrouping>(m, "CommanderGrouping")
            .value("fixed2", CommanderGrouping::Fixed2)
            .value("fixed3", CommanderGrouping::Fixed3)
            .value("halves", CommanderGrouping::Halves)
            .value("logarithm", CommanderGrouping::Logarithm)
            .export_values()
            // allow construction from string
            .def(py::init([](const std::string& str) -> CommanderGrouping { return groupingFromString(str); }));

    // Strategy for reducing the number of permutations/swaps considered in front of every gate
    py::enum_<SwapReduction>(m, "SwapReduction")
            .value("none", SwapReduction::None)
            .value("coupling_limit", SwapReduction::CouplingLimit)
            .value("custom", SwapReduction::Custom)
            .value("increasing", SwapReduction::Increasing)
            .export_values()
            // allow construction from string
            .def(py::init([](const std::string& str) -> SwapReduction { return swapReductionFromString(str); }));

    py::enum_<SynthesisTarget>(m, "SynthesisTarget")
            .value("gates", SynthesisTarget::GATES)
            .value("depth", SynthesisTarget::DEPTH)
            .value("fidelity", SynthesisTarget::FIDELITY)
            .value("gates_only_cnot", SynthesisTarget::GATES_ONLY_CNOT)
            .export_values()
            .def(py::init([](const std::string& str) -> SynthesisTarget { return synthesisTargetFromString(str); }));

    py::enum_<SynthesisStrategy>(m, "SynthesisStrategy")
            .value("use_minimizer", SynthesisStrategy::UseMinimizer)
            .value("minmax", SynthesisStrategy::MinMax)
            .value("start_low", SynthesisStrategy::StartLow)
            .value("start_high", SynthesisStrategy::StartHigh)
            .value("split_iter", SynthesisStrategy::SplitIter)
            .export_values()
            .def(py::init([](const std::string& str) -> SynthesisStrategy { return synthesisStrategyFromString(str); }));

    py::enum_<SynthesisResult>(m, "SynthesisResult")
            .value("sat", SynthesisResult::SAT)
            .value("unsat", SynthesisResult::UNSAT)
            .value("undef", SynthesisResult::UNDEF)
            .export_values()
            .def(py::init([](const std::string& str) -> SynthesisResult { return synthesisResultFromString(str); }));

    py::class_<Configuration>(m, "Configuration", "Configuration options for the MQT QMAP quantum circuit mapping tool")
            .def(py::init<>())
            .def_readwrite("method", &Configuration::method)
            .def_readwrite("verbose", &Configuration::verbose)
            .def_readwrite("layering", &Configuration::layering)
            .def_readwrite("initial_layout", &Configuration::initialLayout)
            .def_readwrite("lookahead", &Configuration::lookahead)
            .def_readwrite("admissible_heuristic", &Configuration::admissibleHeuristic)
            .def_readwrite("lookaheads", &Configuration::nrLookaheads)
            .def_readwrite("first_lookahead_factor", &Configuration::firstLookaheadFactor)
            .def_readwrite("lookahead_factor", &Configuration::lookaheadFactor)
            .def_readwrite("use_teleportation", &Configuration::useTeleportation)
            .def_readwrite("teleportation_qubits", &Configuration::teleportationQubits)
            .def_readwrite("teleportation_seed", &Configuration::teleportationSeed)
            .def_readwrite("teleportation_fake", &Configuration::teleportationFake)
            .def_readwrite("timeout", &Configuration::timeout)
            .def_readwrite("encoding", &Configuration::encoding)
            .def_readwrite("commander_grouping", &Configuration::commanderGrouping)
            .def_readwrite("use_subsets", &Configuration::useSubsets)
            .def_readwrite("include_WCNF", &Configuration::includeWCNF)
            .def_readwrite("enable_limits", &Configuration::enableSwapLimits)
            .def_readwrite("swap_reduction", &Configuration::swapReduction)
            .def_readwrite("swap_limit", &Configuration::swapLimit)
            .def_readwrite("use_bdd", &Configuration::useBDD)
            .def_readwrite("subgraph", &Configuration::subgraph)
            .def_readwrite("pre_mapping_optimizations", &Configuration::preMappingOptimizations)
            .def_readwrite("post_mapping_optimizations", &Configuration::postMappingOptimizations)
            .def_readwrite("add_measurements_to_mapped_circuit", &Configuration::addMeasurementsToMappedCircuit)
            .def("json", &Configuration::json)
            .def("__repr__", &Configuration::toString);

    // Results of the mapping process
    py::class_<MappingResults>(m, "MappingResults", "Results of the MQT QMAP quantum circuit mapping tool")
            .def(py::init<>())
            .def_readwrite("input", &MappingResults::input)
            .def_readwrite("output", &MappingResults::output)
            .def_readwrite("configuration", &MappingResults::config)
            .def_readwrite("time", &MappingResults::time)
            .def_readwrite("timeout", &MappingResults::timeout)
            .def_readwrite("mapped_circuit", &MappingResults::mappedCircuit)
            .def_readwrite("wcnf", &MappingResults::wcnf)
            .def("json", &MappingResults::json)
            .def("csv", &MappingResults::csv)
            .def("__repr__", &MappingResults::toString);

    // Main class for storing circuit information
    py::class_<MappingResults::CircuitInfo>(m, "CircuitInfo", "Circuit information")
            .def(py::init<>())
            .def_readwrite("name", &MappingResults::CircuitInfo::name)
            .def_readwrite("qubits", &MappingResults::CircuitInfo::qubits)
            .def_readwrite("gates", &MappingResults::CircuitInfo::gates)
            .def_readwrite("single_qubit_gates", &MappingResults::CircuitInfo::singleQubitGates)
            .def_readwrite("cnots", &MappingResults::CircuitInfo::cnots)
            .def_readwrite("layers", &MappingResults::CircuitInfo::layers)
            .def_readwrite("swaps", &MappingResults::CircuitInfo::swaps)
            .def_readwrite("direction_reverse", &MappingResults::CircuitInfo::directionReverse)
            .def_readwrite("teleportations", &MappingResults::CircuitInfo::teleportations);

    auto arch       = py::class_<Architecture>(m, "Architecture", "Class representing device/backend information");
    auto properties = py::class_<Architecture::Properties>(arch, "Properties", "Class representing properties of an architecture");

    // Properties of an architecture (e.g. number of qubits, connectivity, error rates, ...)
    properties.def(py::init<>())
            .def_property("name", &Architecture::Properties::getName, &Architecture::Properties::setName)
            .def_property("num_qubits", &Architecture::Properties::getNqubits, &Architecture::Properties::setNqubits)
            .def("get_single_qubit_error", &Architecture::Properties::getSingleQubitErrorRate, "qubit"_a, "operation"_a)
            .def("set_single_qubit_error", &Architecture::Properties::setSingleQubitErrorRate, "qubit"_a, "operation"_a, "error_rate"_a)
            .def("get_two_qubit_error", &Architecture::Properties::getTwoQubitErrorRate, "control"_a, "target"_a, "operation"_a = "cx")
            .def("set_two_qubit_error", &Architecture::Properties::setTwoQubitErrorRate, "control"_a, "target"_a, "error_rate"_a, "operation"_a = "cx")
            .def(
                    "get_readout_error", [](const Architecture::Properties& props, uint16_t qubit) { return props.readoutErrorRate.get(qubit); }, "qubit"_a)
            .def(
                    "set_readout_error", [](Architecture::Properties& props, uint16_t qubit, double rate) { props.readoutErrorRate.set(qubit, rate); }, "qubit"_a, "readout_error_rate"_a)
            .def(
                    "get_t1", [](const Architecture::Properties& props, uint16_t qubit) { return props.t1Time.get(qubit); }, "qubit"_a)
            .def(
                    "set_t1", [](Architecture::Properties& props, uint16_t qubit, double t1) { props.t1Time.set(qubit, t1); }, "qubit"_a, "t1"_a)
            .def(
                    "get_t2", [](const Architecture::Properties& props, uint16_t qubit) { return props.t2Time.get(qubit); }, "qubit"_a)
            .def(
                    "set_t2", [](Architecture::Properties& props, uint16_t qubit, double t2) { props.t2Time.set(qubit, t2); }, "qubit"_a, "t2"_a)
            .def(
                    "get_frequency", [](const Architecture::Properties& props, uint16_t qubit) { return props.qubitFrequency.get(qubit); }, "qubit"_a)
            .def(
                    "set_frequency", [](Architecture::Properties& props, uint16_t qubit, double freq) { props.qubitFrequency.set(qubit, freq); }, "qubit"_a, "qubit_frequency"_a)
            .def(
                    "get_calibration_date", [](const Architecture::Properties& props, uint16_t qubit) { return props.calibrationDate.get(qubit); }, "qubit"_a)
            .def(
                    "set_calibration_date", [](Architecture::Properties& props, uint16_t qubit, const std::string& date) { props.calibrationDate.set(qubit, date); }, "qubit"_a, "calibration_date"_a)
            .def("json", &Architecture::Properties::json,
                 "Returns a JSON-style dictionary of all the information present in the :class:`.Properties`")
            .def("__repr__", &Architecture::Properties::toString,
                 "Prints a JSON-formatted representation of all the information present in the :class:`.Properties`");

    // Interface to the QMAP internal architecture class
    arch.def(py::init<>())
            .def(py::init<uint16_t, const CouplingMap&>(), "num_qubits"_a, "coupling_map"_a)
            .def(py::init<uint16_t, const CouplingMap&, const Architecture::Properties&>(), "num_qubits"_a, "coupling_map"_a, "properties"_a)
            .def_property("name", &Architecture::getName, &Architecture::setName)
            .def_property("num_qubits", &Architecture::getNqubits, &Architecture::setNqubits)
            .def_property("coupling_map", py::overload_cast<>(&Architecture::getCouplingMap), &Architecture::setCouplingMap)
            .def_property("properties", py::overload_cast<>(&Architecture::getProperties), &Architecture::setProperties)
            .def("load_coupling_map", py::overload_cast<AvailableArchitecture>(&Architecture::loadCouplingMap), "available_architecture"_a)
            .def("load_coupling_map", py::overload_cast<const std::string&>(&Architecture::loadCouplingMap), "coupling_map_file"_a)
            .def("load_properties", py::overload_cast<const Architecture::Properties&>(&Architecture::loadProperties), "properties"_a)
            .def("load_properties", py::overload_cast<const std::string&>(&Architecture::loadProperties), "properties"_a);

    py::class_<CliffordOptimizationResults>(m, "CliffordOptimizationResults", "Results of the MQT QMAP Clifford synthesis tool")
            .def(py::init<>())
            .def_readwrite("sat", &CliffordOptimizationResults::result, "Whether the optimization problem was satisfiable")
            .def_readwrite("result_circuit", &CliffordOptimizationResults::resultStringCircuit, "The resulting circuit")
            .def_readwrite("verbose", &CliffordOptimizationResults::verbose, "Verbosity of the debug messages")
            .def_readwrite("choose_best", &CliffordOptimizationResults::chooseBest, "If true, the subgraph of an architecture with the lowest overall fidelity has been chosen, otherwise all possible subgraphs are tried")
            .def_readwrite("strategy", &CliffordOptimizationResults::strategy, "The strategy used to optimize the circuit")
            .def_readwrite("target", &CliffordOptimizationResults::target, "The synthesis target, either 'gates', 'gates_only_cnot', 'depth', or 'fidelity'")
            .def_readwrite("method", &CliffordOptimizationResults::method, "The synthesis method, at the moment only 'z3' is supported")
            .def_readwrite("qubits", &CliffordOptimizationResults::nqubits, "The number of qubits in the resulting circuit")
            .def_readwrite("initial_timesteps", &CliffordOptimizationResults::initialTimesteps, "The number of initial timesteps alloted for synthesis")
            .def_readwrite("gate_count", &CliffordOptimizationResults::gateCount, "The number of gates in the resulting circuit")
            .def_readwrite("depth", &CliffordOptimizationResults::depth, "The depth of the resulting circuit")
            .def_readwrite("fidelity", &CliffordOptimizationResults::fidelity, "The fidelity of the resulting circuit, only available if fidelity data is given")
            .def_readwrite("total_seconds", &CliffordOptimizationResults::totalSeconds, "The total time taken to synthesize the circuit")
            .def("json", &CliffordOptimizationResults::json)
            .def("__repr__", &CliffordOptimizationResults::getStrRepr);

    m.def("map", &map, "map a quantum circuit");
    m.def("synthesize", &synthesize, "synthesize a clifford circuit");
    m.def("optimize", &synthesize, "optimize a clifford circuit");
#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
