/*
 * Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
 * Copyright (c) 2025 Munich Quantum Software Company GmbH
 * All rights reserved.
 *
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License
 */

#include "ir/QuantumComputation.hpp"
#include "na/zoned/Architecture.hpp"
#include "na/zoned/Compiler.hpp"
#include "na/zoned/code_generator/CodeGenerator.hpp"
#include "na/zoned/placer/AStarPlacer.hpp"
#include "na/zoned/placer/VertexMatchingPlacer.hpp"

#include <cstddef>
#include <pybind11/attr.h>
#include <pybind11/cast.h>
#include <pybind11/detail/common.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11_json/pybind11_json.hpp>
#include <spdlog/common.h>
#include <string>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(MQT_QMAP_MODULE_NAME, m, py::mod_gil_not_used()) {
  m.doc() =
      "Python bindings module for MQT QMAP's Zoned Neutral Atom Compiler.";

  py::class_<na::zoned::Architecture> architecture(
      m, "ZonedNeutralAtomArchitecture",
      "Class representing a Zoned Neutral Atom Architecture.");
  architecture.def_static("from_json_file",
                          &na::zoned::Architecture::fromJSONFile, "filename"_a,
                          R"pbdoc(
    Create an architecture from a JSON file.

    Args:
        filename: is the path to the JSON file

    Returns:
        the architecture

    Raises:
        ValueError: if the file does not exist or is not a valid JSON file
)pbdoc");
  architecture.def_static("from_json_string",
                          &na::zoned::Architecture::fromJSONString, "json"_a,
                          R"pbdoc(
    Create an architecture from a JSON string.

    Args:
        json: is the JSON string

    Returns:
        the architecture

    Raises:
        ValueError: if the string is not a valid JSON
)pbdoc");

  py::class_<na::zoned::RoutingAgnosticCompiler> routingAgnosticCompiler(
      m, "RoutingAgnosticCompiler",
      "MQT QMAP's routing-agnostic Zoned Neutral Atom Compiler.");
  routingAgnosticCompiler.def(
      py::init([](const na::zoned::Architecture& arch,
                  const std::string& logLevel, const bool useWindow,
                  const size_t windowSize, const bool dynamicPlacement,
                  const size_t parkingOffset, const bool warnUnsupportedGates)
                   -> na::zoned::RoutingAgnosticCompiler {
        na::zoned::RoutingAgnosticCompiler::Config config;
        config.logLevel = spdlog::level::from_str(logLevel);
        config.placerConfig = {useWindow, windowSize, dynamicPlacement};
        config.codeGeneratorConfig = {parkingOffset, warnUnsupportedGates};
        return {arch, config};
      }),
      py::keep_alive<1, 2>(), "arch"_a, "log_level"_a = "WARN",
      "use_window"_a = true, "window_size"_a = 10, "dynamic_placement"_a = true,
      "parking_offset"_a = 1, "warn_unsupported_gates"_a = true, R"pbdoc(
    Create a routing-agnostic compiler for the given architecture and configurations.

    Args:
        arch: is the zoned neutral atom architecture
        log_level: is the log level for the compiler, possible values are
            "INFO", "WARNING", "ERROR", "CRITICAL"
        use_window: whether to use a window for the placer
        window_size: the size of the window for the placer
        dynamic_placement: whether to use dynamic placement for the placer
        parking_offset: the parking offset of the code generator
        warn_unsupported_gates: whether to warn about unsupported gates in the code generator
)pbdoc");
  routingAgnosticCompiler.def_static(
      "from_json_string",
      [](const na::zoned::Architecture& arch,
         const std::string& json) -> na::zoned::RoutingAgnosticCompiler {
        return {arch, nlohmann::json::parse(json)};
      },
      "arch"_a, "json"_a, R"pbdoc(
    Create a routing-agnostic compiler for the given architecture and with configurations from a JSON string.

    Args:
        arch: is the zoned neutral atom architecture
        json: is the JSON string

    Returns:
        the initialized compiler

    Raises:
        ValueError: if the string is not a valid JSON
)pbdoc");
  routingAgnosticCompiler.def(
      "compile",
      [](na::zoned::RoutingAgnosticCompiler& self,
         const qc::QuantumComputation& qc) -> std::string {
        return self.compile(qc).toString();
      },
      "qc"_a, R"pbdoc(
    Compile a quantum circuit for the zoned neutral atom architecture.

    Args:
        qc: is the quantum circuit

    Returns:
        the compilations result as a string in the .naviz format.
)pbdoc");
  routingAgnosticCompiler.def(
      "stats",
      [](const na::zoned::RoutingAgnosticCompiler& self) -> nlohmann::json {
        return self.getStatistics();
      },
      R"pbdoc(
    Get the statistics of the last compilation.

    Returns:
        the statistics as a dictionary
)pbdoc");

  py::class_<na::zoned::RoutingAwareCompiler> routingAwareCompiler(
      m, "RoutingAwareCompiler",
      "MQT QMAP's routing-aware Zoned Neutral Atom Compiler.");
  routingAwareCompiler.def(
      py::init([](const na::zoned::Architecture& arch,
                  const std::string& logLevel, const bool useWindow,
                  const size_t windowMinWidth, const double windowRatio,
                  const double windowShare, const float deepeningFactor,
                  const float deepeningValue, const float lookaheadFactor,
                  const float reuseLevel, const size_t maxNodes,
                  const size_t parkingOffset, const bool warnUnsupportedGates)
                   -> na::zoned::RoutingAwareCompiler {
        na::zoned::RoutingAwareCompiler::Config config;
        config.logLevel = spdlog::level::from_str(logLevel);
        config.placerConfig = {useWindow,       windowMinWidth,  windowRatio,
                               windowShare,     deepeningFactor, deepeningValue,
                               lookaheadFactor, reuseLevel,      maxNodes};
        config.codeGeneratorConfig = {parkingOffset, warnUnsupportedGates};
        return {arch, config};
      }),
      py::keep_alive<1, 2>(), "arch"_a, "log_level"_a = "WARN",
      "use_window"_a = true, "window_min_width"_a = 8, "window_ratio"_a = 1.0,
      "window_share"_a = 0.6, "deepening_factor"_a = 0.8,
      "deepening_value"_a = 0.2, "lookahead_factor"_a = 0.2,
      "reuse_level"_a = 5.0, "max_nodes"_a = 50000000, "parking_offset"_a = 1,
      "warn_unsupported_gates"_a = true, R"pbdoc(
    Create a routing-aware compiler for the given architecture and configurations.

    Args:
        arch: is the zoned neutral atom architecture
        log_level: is the log level for the compiler, possible values are
            "INFO", "WARNING", "ERROR", "CRITICAL"
        use_window: is a flag whether to use a window for the placer
        window_min_width: is the minimum width of the window for the placer
        window_ratio: is the ratio between the height and the width of the window
        window_share: is the share of free sites in the window in relation to the
            number of atoms to be moved in this step
        deepening_factor: controls the impact of the term in the heuristic of the
            A* search that resembles the standard deviation of the differences
            between the current and target sites of the atoms to be moved in every
            orientation
        deepening_value: is added to the sum of standard deviations before it is
            multiplied with the number of unplaced nodes and :attr:`deepening_factor`
        lookahead_factor: controls the lookahead's influence that considers the
            distance of atoms to their interaction partner in the next layer
        reuse_level: is the reuse level that corresponds to the estimated extra
            fidelity loss due to the extra trap transfers when the atom is not
            reused and instead moved to the storage zone and back to the
            entanglement zone
        max_nodes: is the maximum number of nodes that are considered in the A*
            search. If this number is exceeded, the search is aborted and an error
            is raised. In the current implementation, one node roughly consumes 120
            Byte. Hence, allowing 50,000,000 nodes results in memory consumption of
            about 6 GB plus the size of the rest of the data structures.
        parking_offset: is the parking offset of the code generator
        warn_unsupported_gates: is a flag whether to warn about unsupported gates
            in the code generator
)pbdoc");
  routingAwareCompiler.def_static(
      "from_json_string",
      [](const na::zoned::Architecture& arch,
         const std::string& json) -> na::zoned::RoutingAwareCompiler {
        return {arch, nlohmann::json::parse(json)};
      },
      "arch"_a, "json"_a, R"pbdoc(
    Create a routing-aware compiler for the given architecture and configurations from a JSON string.

    Args:
        arch: is the zoned neutral atom architecture
        json: is the JSON string

    Returns:
        the initialized compiler

    Raises:
        ValueError: if the string is not a valid JSON
)pbdoc");
  routingAwareCompiler.def(
      "compile",
      [](na::zoned::RoutingAwareCompiler& self,
         const qc::QuantumComputation& qc) -> std::string {
        return self.compile(qc).toString();
      },
      "qc"_a, R"pbdoc(
    Compile a quantum circuit for the zoned neutral atom architecture.

    Args:
        qc: is the quantum circuit

    Returns:
        the compilations result as a string in the .naviz format.
)pbdoc");
  routingAwareCompiler.def(
      "stats",
      [](const na::zoned::RoutingAwareCompiler& self) -> nlohmann::json {
        return self.getStatistics();
      },
      R"pbdoc(
    Get the statistics of the last compilation.

    Returns:
        the statistics as a dictionary
)pbdoc");
}
