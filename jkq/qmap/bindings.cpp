/*
 * This file is part of JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum/ for more information.
 */

#include "exact/ExactMapper.hpp"
#include "heuristic/HeuristicMapper.hpp"
#include "nlohmann/json.hpp"
#include "pybind11/pybind11.h"
#include "pybind11_json/pybind11_json.hpp"
#include "qiskit/QuantumCircuit.hpp"

namespace py = pybind11;
namespace nl = nlohmann;
using namespace pybind11::literals;

// c++ binding function
nl::json map(const py::object& circ, const py::object& arch, const nl::json& jsonConfig) {
    qc::QuantumComputation qc{};
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
        return {{"error", ss.str()}};
    }

    Architecture architecture{};
    try {
        if (py::isinstance<py::str>(arch)) {
            auto&& cm = arch.cast<std::string>();
            architecture.loadCouplingMap(cm);
        } else {
            auto&& cm = arch.cast<AvailableArchitectures>();
            architecture.loadCouplingMap(cm);
        }

        if (jsonConfig.contains("calibration")) {
            auto&& cal = jsonConfig["calibration"].get<std::string>();
            if (!cal.empty())
                architecture.loadCalibrationData(cal);
        }
    } catch (std::exception const& e) {
        std::stringstream ss{};
        ss << "Could not import architecture: " << e.what();
        return {{"error", ss.str()}};
    }

    MappingSettings ms{};
    Method          method = Method::Heuristic;
    if (jsonConfig.contains("method")) {
        method = jsonConfig["method"].get<Method>();
    }

    ms.initialLayoutStrategy = InitialLayoutStrategy::Dynamic;
    if (jsonConfig.contains("initialLayout")) {
        ms.initialLayoutStrategy = jsonConfig["initialLayout"].get<InitialLayoutStrategy>();
    }

    ms.layeringStrategy = LayeringStrategy::IndividualGates;
    if (jsonConfig.contains("layering")) {
        ms.layeringStrategy = jsonConfig["layering"].get<LayeringStrategy>();
    }

    ms.encoding = Encodings::Naive;
    if (jsonConfig.contains("encoding")) {
        ms.encoding = jsonConfig["encoding"].get<Encodings>();
    }

    ms.grouping = CMDRVariableGroupings::Halves;
    if (jsonConfig.contains("grouping")) {
        ms.grouping = jsonConfig["grouping"].get<CMDRVariableGroupings>();
    }

    if (jsonConfig.contains("strategy")) {
        ms.enableLimits = true;
        ms.strategy     = SwapReductionStrategy::CouplingLimit;
        ms.strategy     = jsonConfig["strategy"].get<SwapReductionStrategy>();
        if (jsonConfig.contains("limit")) {
            ms.limit = jsonConfig["limit"].get<int>();
        }
        if (jsonConfig.contains("useBDD")) {
            ms.useBDD = true;
        }
    }

    if (jsonConfig.contains("use_teleportation")) {
        auto useTeleportation = jsonConfig["use_teleportation"].get<bool>();
        if (useTeleportation) {
            ms.teleportationQubits = std::min((architecture.getNqubits() - qc.getNqubits()) & ~1u, 8u);
            ms.teleportationSeed   = jsonConfig["teleportation_seed"].get<unsigned long long>();
            ms.teleportationFake   = jsonConfig["teleportation_fake"].get<bool>();
        }
    }

    if (jsonConfig.contains("verbose")) {
        ms.verbose = jsonConfig["verbose"].get<bool>();
    }

    if (jsonConfig.contains("use_subsets")) {
        ms.useQubitSubsets = jsonConfig["use_subsets"].get<bool>();
    }

    bool printStatistics = false;
    if (jsonConfig.contains("statistics")) {
        printStatistics = jsonConfig["statistics"].get<bool>();
    }

    bool printCSV = false;
    if (jsonConfig.contains("csv")) {
        printCSV = jsonConfig["csv"].get<bool>();
    }

    bool saveMappedCircuit = false;
    if (jsonConfig.contains("saveMappedCircuit")) {
        saveMappedCircuit = jsonConfig["saveMappedCircuit"].get<bool>();
    }

    std::unique_ptr<Mapper> mapper;
    try {
        if (method == Method::Heuristic) {
            mapper = std::make_unique<HeuristicMapper>(qc, architecture);
        } else {
            mapper = std::make_unique<ExactMapper>(qc, architecture);
        }
    } catch (std::exception const& e) {
        std::stringstream ss{};
        ss << "Could not construct mapper: " << e.what();
        return {{"error", ss.str()}};
    }

    try {
        mapper->map(ms);
    } catch (std::exception const& e) {
        std::stringstream ss{};
        ss << "Error during mapping: " << e.what();
        return {{"error", ss.str()}};
    }

    std::stringstream ss{};
    mapper->dumpResult(ss, qc::OpenQASM);
    auto result = mapper->produceJSON(printStatistics);
    if (printCSV)
        result["csv"] = mapper->produceCSVEntry();

    if (saveMappedCircuit) {
        std::stringstream qasm{};
        mapper->dumpResult(qasm, qc::OpenQASM);
        result["mapped_circuit"]["qasm"] = qasm.str();
    }
    return result;
}

PYBIND11_MODULE(pyqmap, m) {
    m.doc() = "pybind11 for the JKQ QMAP quantum circuit mapping tool";
    m.def("map", &map, "map a quantum circuit");

    py::enum_<AvailableArchitectures>(m, "Arch")
            .value("IBM_QX4", AvailableArchitectures::IBM_QX4)
            .value("IBM_QX5", AvailableArchitectures::IBM_QX5)
            .value("IBMQ_Yorktown", AvailableArchitectures::IBMQ_Yorktown)
            .value("IBMQ_London", AvailableArchitectures::IBMQ_London)
            .value("IBMQ_Bogota", AvailableArchitectures::IBMQ_Bogota)
            .value("IBMQ_Tokyo", AvailableArchitectures::IBMQ_Tokyo)
            .export_values();

    py::enum_<Method>(m, "Method")
            .value("heuristic", Method::Heuristic)
            .value("exact", Method::Exact)
            .export_values();

    py::enum_<InitialLayoutStrategy>(m, "InitialLayoutStrategy")
            .value("identity", InitialLayoutStrategy::Identity)
            .value("static", InitialLayoutStrategy::Static)
            .value("dynamic", InitialLayoutStrategy::Dynamic)
            .export_values();

    py::enum_<LayeringStrategy>(m, "LayeringStrategy")
            .value("individual_gates", LayeringStrategy::IndividualGates)
            .value("disjoint_qubits", LayeringStrategy::DisjointQubits)
            .value("odd_gates", LayeringStrategy::OddGates)
            .value("qubit_triangle", LayeringStrategy::QubitTriangle)
            .export_values();

    py::enum_<Encodings>(m, "Encoding")
            .value("none", Encodings::Naive)
            .value("commander", Encodings::Commander)
            .value("bimander", Encodings::Bimander)
            .export_values();

    py::enum_<CMDRVariableGroupings>(m, "Grouping")
            .value("fixed2", CMDRVariableGroupings::Fixed2)
            .value("fixed3", CMDRVariableGroupings::Fixed3)
            .value("halves", CMDRVariableGroupings::Halves)
            .value("logarithm", CMDRVariableGroupings::Logarithm)
            .export_values();

    py::enum_<SwapReductionStrategy>(m, "Strategy")
            .value("none", SwapReductionStrategy::None)
            .value("coupling_limit", SwapReductionStrategy::CouplingLimit)
            .value("custom", SwapReductionStrategy::Custom)
            .value("increasing", SwapReductionStrategy::Increasing)
            .export_values();

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
