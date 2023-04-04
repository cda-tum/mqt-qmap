//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

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
MappingResults map(const py::object& circ, Architecture& arch,
                   Configuration& config) {
  qc::QuantumComputation qc{};

  loadQC(qc, circ);

  if (config.useTeleportation) {
    config.teleportationQubits =
        std::min((arch.getNqubits() - qc.getNqubits()) & ~1U,
                 static_cast<std::size_t>(8));
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
  mapper->dumpResult(qasm, qc::Format::OpenQASM);
  results.mappedCircuit = qasm.str();

  return results;
}

PYBIND11_MODULE(pyqmap, m) {
  m.doc() = "pybind11 for the MQT QMAP quantum circuit mapping tool";

  // Pre-defined architecture available within QMAP
  py::enum_<AvailableArchitecture>(m, "Arch")
      .value("IBM_QX4", AvailableArchitecture::IbmQx4,
             "5 qubit, directed bow tie layout")
      .value("IBM_QX5", AvailableArchitecture::IbmQx5,
             "16 qubit, directed ladder layout")
      .value("IBMQ_Yorktown", AvailableArchitecture::IbmqYorktown,
             "5 qubit, undirected bow tie layout")
      .value("IBMQ_London", AvailableArchitecture::IbmqLondon,
             "5 qubit, undirected T-shape layout")
      .value("IBMQ_Bogota", AvailableArchitecture::IbmqBogota,
             "5 qubit, undirected linear chain layout")
      .value("IBMQ_Casablanca", AvailableArchitecture::IbmqCasablanca,
             "7 qubit, undirected H-shape layout")
      .value("IBMQ_Tokyo", AvailableArchitecture::IbmqTokyo,
             "20 qubit, undirected brick-like layout")
      .value("Rigetti_Agave", AvailableArchitecture::RigettiAgave,
             "8 qubit, undirected ring layout")
      .value("Rigetti_Aspen", AvailableArchitecture::RigettiAspen,
             "16 qubit, undirected dumbbell layout")
      .export_values()
      // allow construction from string
      .def(py::init([](const std::string& str) -> AvailableArchitecture {
        return architectureFromString(str);
      }));

  // Mapping methodology to use
  py::enum_<Method>(m, "Method")
      .value("heuristic", Method::Heuristic)
      .value("exact", Method::Exact)
      .export_values()
      // allow construction from string
      .def(py::init([](const std::string& str) -> Method {
        return methodFromString(str);
      }));

  // Initial layout strategy
  py::enum_<InitialLayout>(m, "InitialLayout")
      .value("identity", InitialLayout::Identity)
      .value("static", InitialLayout::Static)
      .value("dynamic", InitialLayout::Dynamic)
      .export_values()
      // allow construction from string
      .def(py::init([](const std::string& str) -> InitialLayout {
        return initialLayoutFromString(str);
      }));

  // Gate clustering / layering strategy
  py::enum_<Layering>(m, "Layering")
      .value("individual_gates", Layering::IndividualGates)
      .value("disjoint_qubits", Layering::DisjointQubits)
      .value("odd_gates", Layering::OddGates)
      .value("qubit_triangle", Layering::QubitTriangle)
      .value("disjoint_2q_blocks", Layering::Disjoint2qBlocks)
      .export_values()
      // allow construction from string
      .def(py::init([](const std::string& str) -> Layering {
        return layeringFromString(str);
      }));

  // Encoding settings for at-most-one and exactly-one constraints
  py::enum_<Encoding>(m, "Encoding")
      .value("naive", Encoding::Naive)
      .value("commander", Encoding::Commander)
      .value("bimander", Encoding::Bimander)
      .export_values()
      // allow construction from string
      .def(py::init([](const std::string& str) -> Encoding {
        return encodingFromString(str);
      }));

  // Grouping settings if using the commander encoding
  py::enum_<CommanderGrouping>(m, "CommanderGrouping")
      .value("fixed2", CommanderGrouping::Fixed2)
      .value("fixed3", CommanderGrouping::Fixed3)
      .value("halves", CommanderGrouping::Halves)
      .value("logarithm", CommanderGrouping::Logarithm)
      .export_values()
      // allow construction from string
      .def(py::init([](const std::string& str) -> CommanderGrouping {
        return groupingFromString(str);
      }));

  // Strategy for reducing the number of permutations/swaps considered in front
  // of every gate
  py::enum_<SwapReduction>(m, "SwapReduction")
      .value("none", SwapReduction::None)
      .value("coupling_limit", SwapReduction::CouplingLimit)
      .value("custom", SwapReduction::Custom)
      .value("increasing", SwapReduction::Increasing)
      .export_values()
      // allow construction from string
      .def(py::init([](const std::string& str) -> SwapReduction {
        return swapReductionFromString(str);
      }));

  // All configuration options for QMAP
  py::class_<Configuration>(
      m, "Configuration",
      "Configuration options for the MQT QMAP quantum circuit mapping tool")
      .def(py::init<>())
      .def_readwrite("method", &Configuration::method)
      .def_readwrite("verbose", &Configuration::verbose)
      .def_readwrite("debug", &Configuration::debug)
      .def_readwrite("layering", &Configuration::layering)
      .def_readwrite("initial_layout", &Configuration::initialLayout)
      .def_readwrite("lookahead", &Configuration::lookahead)
      .def_readwrite("admissible_heuristic",
                     &Configuration::admissibleHeuristic)
      .def_readwrite("consider_fidelity", &Configuration::considerFidelity)
      .def_readwrite("lookaheads", &Configuration::nrLookaheads)
      .def_readwrite("first_lookahead_factor",
                     &Configuration::firstLookaheadFactor)
      .def_readwrite("lookahead_factor", &Configuration::lookaheadFactor)
      .def_readwrite("use_teleportation", &Configuration::useTeleportation)
      .def_readwrite("teleportation_qubits",
                     &Configuration::teleportationQubits)
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
      .def_readwrite("pre_mapping_optimizations",
                     &Configuration::preMappingOptimizations)
      .def_readwrite("post_mapping_optimizations",
                     &Configuration::postMappingOptimizations)
      .def_readwrite("add_measurements_to_mapped_circuit",
                     &Configuration::addMeasurementsToMappedCircuit)
      .def("json", &Configuration::json)
      .def("__repr__", &Configuration::toString);

  // Results of the mapping process
  py::class_<MappingResults>(
      m, "MappingResults",
      "Results of the MQT QMAP quantum circuit mapping tool")
      .def(py::init<>())
      .def_readwrite("input", &MappingResults::input)
      .def_readwrite("output", &MappingResults::output)
      .def_readwrite("configuration", &MappingResults::config)
      .def_readwrite("time", &MappingResults::time)
      .def_readwrite("timeout", &MappingResults::timeout)
      .def_readwrite("mapped_circuit", &MappingResults::mappedCircuit)
      .def_readwrite("heuristic_benchmark", &MappingResults::heuristicBenchmark)
      .def_readwrite("layer_heuristic_benchmark",
                     &MappingResults::layerHeuristicBenchmark)
      .def_readwrite("wcnf", &MappingResults::wcnf)
      .def("json", &MappingResults::json)
      .def("csv", &MappingResults::csv)
      .def("__repr__", &MappingResults::toString);

  // Main class for storing circuit information
  py::class_<MappingResults::CircuitInfo>(m, "CircuitInfo",
                                          "Circuit information")
      .def(py::init<>())
      .def_readwrite("name", &MappingResults::CircuitInfo::name)
      .def_readwrite("qubits", &MappingResults::CircuitInfo::qubits)
      .def_readwrite("gates", &MappingResults::CircuitInfo::gates)
      .def_readwrite("single_qubit_gates",
                     &MappingResults::CircuitInfo::singleQubitGates)
      .def_readwrite("cnots", &MappingResults::CircuitInfo::cnots)
      .def_readwrite("layers", &MappingResults::CircuitInfo::layers)
      .def_readwrite("swaps", &MappingResults::CircuitInfo::swaps)
      .def_readwrite("direction_reverse",
                     &MappingResults::CircuitInfo::directionReverse)
      .def_readwrite("teleportations",
                     &MappingResults::CircuitInfo::teleportations);

  // Heuristic benchmark information
  py::class_<MappingResults::HeuristicBenchmarkInfo>(
      m, "HeuristicBenchmarkInfo", "Heuristic benchmark information")
      .def(py::init<>())
      .def_readwrite("expanded_nodes",
                     &MappingResults::HeuristicBenchmarkInfo::expandedNodes)
      .def_readwrite("generated_nodes",
                     &MappingResults::HeuristicBenchmarkInfo::generatedNodes)
      .def_readwrite("solution_depth",
                     &MappingResults::HeuristicBenchmarkInfo::solutionDepth)
      .def_readwrite("time_per_node",
                     &MappingResults::HeuristicBenchmarkInfo::timePerNode)
      .def_readwrite(
          "average_branching_factor",
          &MappingResults::HeuristicBenchmarkInfo::averageBranchingFactor)
      .def_readwrite(
          "effective_branching_factor",
          &MappingResults::HeuristicBenchmarkInfo::effectiveBranchingFactor);

  auto arch = py::class_<Architecture>(
      m, "Architecture", "Class representing device/backend information");
  auto properties = py::class_<Architecture::Properties>(
      arch, "Properties", "Class representing properties of an architecture");

  // Properties of an architecture (e.g. number of qubits, connectivity, error
  // rates, ...)
  properties.def(py::init<>())
      .def_property("name", &Architecture::Properties::getName,
                    &Architecture::Properties::setName)
      .def_property("num_qubits", &Architecture::Properties::getNqubits,
                    &Architecture::Properties::setNqubits)
      .def("get_single_qubit_error",
           &Architecture::Properties::getSingleQubitErrorRate, "qubit"_a,
           "operation"_a)
      .def("set_single_qubit_error",
           &Architecture::Properties::setSingleQubitErrorRate, "qubit"_a,
           "operation"_a, "error_rate"_a)
      .def("get_two_qubit_error",
           &Architecture::Properties::getTwoQubitErrorRate, "control"_a,
           "target"_a, "operation"_a = "cx")
      .def("set_two_qubit_error",
           &Architecture::Properties::setTwoQubitErrorRate, "control"_a,
           "target"_a, "error_rate"_a, "operation"_a = "cx")
      .def(
          "get_readout_error",
          [](const Architecture::Properties& props, std::uint16_t qubit) {
            return props.readoutErrorRate.get(qubit);
          },
          "qubit"_a)
      .def(
          "set_readout_error",
          [](Architecture::Properties& props, std::uint16_t qubit,
             double rate) { props.readoutErrorRate.set(qubit, rate); },
          "qubit"_a, "readout_error_rate"_a)
      .def(
          "get_t1",
          [](const Architecture::Properties& props, std::uint16_t qubit) {
            return props.t1Time.get(qubit);
          },
          "qubit"_a)
      .def(
          "set_t1",
          [](Architecture::Properties& props, std::uint16_t qubit, double t1) {
            props.t1Time.set(qubit, t1);
          },
          "qubit"_a, "t1"_a)
      .def(
          "get_t2",
          [](const Architecture::Properties& props, std::uint16_t qubit) {
            return props.t2Time.get(qubit);
          },
          "qubit"_a)
      .def(
          "set_t2",
          [](Architecture::Properties& props, std::uint16_t qubit, double t2) {
            props.t2Time.set(qubit, t2);
          },
          "qubit"_a, "t2"_a)
      .def(
          "get_frequency",
          [](const Architecture::Properties& props, std::uint16_t qubit) {
            return props.qubitFrequency.get(qubit);
          },
          "qubit"_a)
      .def(
          "set_frequency",
          [](Architecture::Properties& props, std::uint16_t qubit,
             double freq) { props.qubitFrequency.set(qubit, freq); },
          "qubit"_a, "qubit_frequency"_a)
      .def(
          "get_calibration_date",
          [](const Architecture::Properties& props, std::uint16_t qubit) {
            return props.calibrationDate.get(qubit);
          },
          "qubit"_a)
      .def(
          "set_calibration_date",
          [](Architecture::Properties& props, std::uint16_t qubit,
             const std::string& date) {
            props.calibrationDate.set(qubit, date);
          },
          "qubit"_a, "calibration_date"_a)
      .def("json", &Architecture::Properties::json,
           "Returns a JSON-style dictionary of all the information present in "
           "the :class:`.Properties`")
      .def("__repr__", &Architecture::Properties::toString,
           "Prints a JSON-formatted representation of all the information "
           "present in the :class:`.Properties`");

  // Interface to the QMAP internal architecture class
  arch.def(py::init<>())
      .def(py::init<std::uint16_t, const CouplingMap&>(), "num_qubits"_a,
           "coupling_map"_a)
      .def(py::init<std::uint16_t, const CouplingMap&,
                    const Architecture::Properties&>(),
           "num_qubits"_a, "coupling_map"_a, "properties"_a)
      .def_property("name", &Architecture::getName, &Architecture::setName)
      .def_property("num_qubits", &Architecture::getNqubits,
                    &Architecture::setNqubits)
      .def_property("coupling_map",
                    py::overload_cast<>(&Architecture::getCouplingMap),
                    &Architecture::setCouplingMap)
      .def_property("properties",
                    py::overload_cast<>(&Architecture::getProperties),
                    &Architecture::setProperties)
      .def("load_coupling_map",
           py::overload_cast<AvailableArchitecture>(
               &Architecture::loadCouplingMap),
           "available_architecture"_a)
      .def(
          "load_coupling_map",
          py::overload_cast<const std::string&>(&Architecture::loadCouplingMap),
          "coupling_map_file"_a)
      .def("load_properties",
           py::overload_cast<const Architecture::Properties&>(
               &Architecture::loadProperties),
           "properties"_a)
      .def("load_properties",
           py::overload_cast<const std::string&>(&Architecture::loadProperties),
           "properties"_a);

  // Main mapping function
  m.def("map", &map, "map a quantum circuit", "circ"_a, "arch"_a, "config"_a);

  // Target metric for the Clifford synthesizer
  py::enum_<cs::TargetMetric>(m, "TargetMetric")
      .value("gates", cs::TargetMetric::Gates, "Optimize gate count.")
      .value("two_qubit_gates", cs::TargetMetric::TwoQubitGates,
             "Optimize two-qubit gate count.")
      .value("depth", cs::TargetMetric::Depth, "Optimize circuit depth.")
      .export_values()
      .def(py::init([](const std::string& name) {
        return cs::targetMetricFromString(name);
      }));
  py::implicitly_convertible<py::str, cs::TargetMetric>();

  py::enum_<plog::Severity>(m, "Verbosity")
      .value("none", plog::Severity::none, "No output.")
      .value("fatal", plog::Severity::fatal, "Only show fatal errors.")
      .value("error", plog::Severity::error, "Show errors.")
      .value("warning", plog::Severity::warning, "Show warnings.")
      .value("info", plog::Severity::info, "Show general information.")
      .value("debug", plog::Severity::debug,
             "Show additional debug information.")
      .value("verbose", plog::Severity::verbose, "Show all information.")
      .export_values()
      .def(py::init([](const std::string& name) {
        return plog::severityFromString(name.c_str());
      }));
  py::implicitly_convertible<py::str, plog::Severity>();

  // Configuration for the synthesis
  py::class_<cs::Configuration>(
      m, "SynthesisConfiguration",
      "Configuration options for the MQT QMAP Clifford synthesis tool.")
      .def(py::init<>())
      .def_readwrite("initial_timestep_limit",
                     &cs::Configuration::initialTimestepLimit,
                     "Initial timestep limit for the Clifford synthesis. "
                     "Defaults to `0`, which implies that the initial timestep "
                     "limit is determined automatically.")
      .def_readwrite(
          "use_maxsat", &cs::Configuration::useMaxSAT,
          "Use MaxSAT to solve the synthesis problem or to really on the "
          "binary search scheme for finding the optimum. Defaults to `false`.")
      .def_readwrite(
          "target_metric", &cs::Configuration::target,
          "Target metric for the Clifford synthesis. Defaults to `gates`.")
      .def_readwrite("use_symmetry_breaking",
                     &cs::Configuration::useSymmetryBreaking,
                     "Use symmetry breaking clauses to speed up the synthesis "
                     "process. Defaults to `true`.")
      .def_readwrite("dump_intermediate_results",
                     &cs::Configuration::dumpIntermediateResults,
                     "Dump intermediate results of the synthesis process. "
                     "Defaults to `false`.")
      .def_readwrite("intermediate_results_path",
                     &cs::Configuration::intermediateResultsPath,
                     "Path to the directory where intermediate results should "
                     "be dumped. Defaults to `./`. The path needs to include a "
                     "path separator at the end.")
      .def_readwrite(
          "verbosity", &cs::Configuration::verbosity,
          "Verbosity level for the synthesis process. Defaults to 'warning'.")
      .def_readwrite(
          "n_threads", &cs::Configuration::nThreads,
          "Number of threads to use for the synthesis. Defaults to `1`.")
      .def_readwrite(
          "minimize_gates_after_depth_optimization",
          &cs::Configuration::minimizeGatesAfterDepthOptimization,
          "Depth optimization might produce a circuit with more gates than "
          "necessary. This option enables an additional run of the synthesizer "
          "to minimize the overall number of gates. Defaults to `false`.")
      .def_readwrite(
          "try_higher_gate_limit_for_two_qubit_gate_optimization",
          &cs::Configuration::tryHigherGateLimitForTwoQubitGateOptimization,
          "When optimizing two-qubit gates, the synthesizer might fail "
          "to find an optimal solution for a certain timestep limit, but there "
          "might be a better solution for some higher timestep limit. This "
          "option enables an additional run of the synthesizer with a higher "
          "gate limit. Defaults to `false`.")
      .def_readwrite("gate_limit_factor", &cs::Configuration::gateLimitFactor,
                     "Factor by which the gate limit is increased when "
                     "trying to find a better solution for the two-qubit "
                     "gate optimization. Defaults to `1.1`.")
      .def_readwrite(
          "minimize_gates_after_two_qubit_gate_optimization",
          &cs::Configuration::minimizeGatesAfterTwoQubitGateOptimization,
          "Two-qubit gate optimization might produce a circuit "
          "with more gates than necessary. This option enables "
          "an additional run of the synthesizer to minimize the "
          "overall number of gates. Defaults to `false`.")
      .def_readwrite("heuristic", &cs::Configuration::heuristic,
                     "Use heuristic to synthesize the circuit. "
                     "This method synthesizes shallow intermediate circuits "
                     "and combines them. Defaults to `false`.")
      .def_readwrite("heuristic", &cs::Configuration::split_size,
                     "Size of subcircuits used in heuristic. "
                     "Defaults to `5`.")
      .def("json", &cs::Configuration::json,
           "Returns a JSON-style dictionary of all the information present in "
           "the :class:`.Configuration`")
      .def(
          "__repr__",
          [](const cs::Configuration& config) { return config.json().dump(2); },
          "Prints a JSON-formatted representation of all the information "
          "present in the :class:`.Configuration`");

  // Results of the synthesis
  py::class_<cs::Results>(m, "SynthesisResults",
                          "Results of the MQT QMAP Clifford synthesis tool.")
      .def(py::init<>())
      .def_property_readonly("gates", &cs::Results::getGates,
                             "Returns the number of gates in the circuit.")
      .def_property_readonly("single_qubit_gates",
                             &cs::Results::getSingleQubitGates,
                             "Returns the number of single-qubit gates in the "
                             "synthesized circuit.")
      .def_property_readonly("two_qubit_gates", &cs::Results::getTwoQubitGates,
                             "Returns the number of two-qubit gates in the "
                             "synthesized circuit.")
      .def_property_readonly("depth", &cs::Results::getDepth,
                             "Returns the depth of the synthesized circuit.")
      .def_property_readonly("runtime", &cs::Results::getRuntime,
                             "Returns the runtime of the synthesis in seconds.")
      .def_property_readonly("solver_calls", &cs::Results::getSolverCalls,
                             "Returns the number of calls to the SAT solver.")
      .def_property_readonly(
          "circuit", &cs::Results::getResultCircuit,
          "Returns the synthesized circuit as a qasm string.")
      .def_property_readonly("tableau", &cs::Results::getResultTableau,
                             "Returns a string representation of the "
                             "synthesized circuit's tableau.")
      .def("sat", &cs::Results::sat,
           "Returns `true` if the synthesis was successful.")
      .def("unsat", &cs::Results::unsat,
           "Returns `true` if the synthesis was unsuccessful.");

  auto tableau = py::class_<cs::Tableau>(
      m, "Tableau", "A class for representing stabilizer tableaus.");
  tableau.def(py::init<std::size_t, bool>(), "n"_a,
              "include_destabilizers"_a = false,
              "Creates a tableau for an n-qubit Clifford.");
  tableau.def(
      py::init<const std::string&>(), "tableau"_a,
      "Constructs a tableau from a string description. This can either be a "
      "semicolon separated binary matrix or a list of Pauli strings.");
  tableau.def(
      py::init<const std::string&, const std::string&>(), "stabilizers"_a,
      "destabilizers"_a,
      "Constructs a tableau from two lists of Pauli strings, the Stabilizers"
      "and Destabilizers.");

  auto quantumComputation = py::class_<qc::QuantumComputation>(
      m, "QuantumComputation",
      "A class for the intermediate representation of quantum circuits in the "
      "Munich Quantum Toolkit.");
  quantumComputation.def_static(
      "from_file",
      [](const std::string& filename) {
        return qc::QuantumComputation(filename);
      },
      "filename"_a, "Reads a quantum circuit from a file.");
  quantumComputation.def_static(
      "from_qasm_str",
      [](const std::string& qasm) {
        std::stringstream      ss(qasm);
        qc::QuantumComputation qc{};
        qc.import(ss, qc::Format::OpenQASM);
        return qc;
      },
      "qasm"_a, "Reads a quantum circuit from a qasm string.");
  quantumComputation.def_static(
      "from_qiskit",
      [](const py::object& circuit) {
        qc::QuantumComputation qc{};
        qc::qiskit::QuantumCircuit::import(qc, circuit);
        return qc;
      },
      "circuit"_a,
      "Reads a quantum circuit from a Qiskit :class:`QuantumCircuit`.");

  auto synthesizer = py::class_<cs::CliffordSynthesizer>(
      m, "CliffordSynthesizer", "A class for synthesizing Clifford circuits.");

  synthesizer.def(py::init<cs::Tableau, cs::Tableau>(), "initial_tableau"_a,
                  "target_tableau"_a,
                  "Constructs a synthesizer for two tableaus representing the "
                  "initial and target state.");
  synthesizer.def(py::init<cs::Tableau>(), "target_tableau"_a,
                  "Constructs a synthesizer for a tableau representing the "
                  "target state.");
  synthesizer.def(py::init<qc::QuantumComputation&, bool>(), "qc"_a,
                  "use_destabilizers"_a
                  "Constructs a synthesizer for a quantum computation "
                  "representing the target state.");
  synthesizer.def(
      py::init<cs::Tableau, qc::QuantumComputation&>(), "initial_tableau"_a,
      "qc"_a,
      "Constructs a synthesizer for a quantum computation representing the "
      "target state that starts in an initial state represented by a tableau.");
  synthesizer.def("synthesize", &cs::CliffordSynthesizer::synthesize,
                  "config"_a = cs::Configuration(),
                  "Runs the synthesis with the given configuration.");
  synthesizer.def_property_readonly("results",
                                    &cs::CliffordSynthesizer::getResults,
                                    "Returns the results of the synthesis.");

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
