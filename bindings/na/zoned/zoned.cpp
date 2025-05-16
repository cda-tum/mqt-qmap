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

#include <nlohmann/json.hpp>
#include <pybind11/attr.h>
#include <pybind11/cast.h>
#include <pybind11/pybind11.h>
#include <string>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(MQT_QMAP_MODULE_NAME, m, py::mod_gil_not_used()) {
  m.doc() = "Bindings for mqt.qmap.na.zoned";

  py::class_<na::zoned::Architecture>(m, "ZonedNeutralAtomArchitecture", R"(
The representation of a zoned neutral atom architecture.
)")
      .def(py::init<>(), R"(
Create a plain architecture.
)")
      .def_readwrite("name", &na::zoned::Architecture::name, R"(
Name of the architecture.
)") // todo: add remaining properties
      .def_static("from_json_file", &na::zoned::Architecture::fromJSONFile,
                  "filename"_a, R"(
Create an architecture from a JSON file.

:param filename: is the path to the JSON file
:raises ValueError: if the file does not exist or is not a valid JSON file
:returns: the architecture
)")
      .def_static("from_json_string", &na::zoned::Architecture::fromJSONString,
                  "json"_a, R"(
Create an architecture from a JSON string.

:param json: is the JSON string
:raises ValueError: if the string is not a valid JSON string
:returns: the architecture
)");

  py::class_<na::zoned::RoutingAgnosticCompiler>(m, "RoutingAgnosticCompiler",
                                                 R"(
MQT QMAP's routing-agnostic Zoned Neutral Atom Compiler is a general purpose
compiler for Zoned Neutral Atom Architectures.
)")
      .def(py::init<const na::zoned::Architecture&,
                    const na::zoned::RoutingAgnosticCompiler::Config&>(),
           py::keep_alive<1, 2>(), "architecture"_a, "settings"_a, R"(
Create a routing-agnostic compiler for the given architecture and settings.

:param architecture: is the zoned neutral atom architecture
:param settings: is a dictionary with the settings for the compiler
)")
      .def(py::init<const na::zoned::Architecture&>(), py::keep_alive<1, 2>(),
           "architecture"_a, R"(
Create a routing-agnostic compiler for the given architecture and settings.

:param architecture: is the zoned neutral atom architecture
)")
      .def(
          "compile",
          [](na::zoned::RoutingAgnosticCompiler& self,
             const qc::QuantumComputation& qc) -> std::string {
            return self.compile(qc).toString();
          },
          "qc"_a,
          R"(
Compile a quantum circuit for the zoned neutral atom architecture.

:param qc: is the quantum circuit
:returns: the compilation results as an NAComputation.
)")
      .def(
          "stats",
          [](const na::zoned::RoutingAgnosticCompiler& self) -> nlohmann::json {
            return self.getStatistics();
          },
          R"(
Get the statistics of the last compilation.

:returns: the statistics as a dictionary
)");

  py::class_<na::zoned::RoutingAwareCompiler>(m, "RoutingAwareCompiler", R"(
MQT QMAP's routing-aware Zoned Neutral Atom Compiler is a general purpose
compiler for Zoned Neutral Atom Architectures.
)")
      .def(py::init<const na::zoned::Architecture&,
                    const na::zoned::RoutingAwareCompiler::Config&>(),
           py::keep_alive<1, 2>(), "architecture"_a, "settings"_a, R"(
Create a routing-aware compiler for the given architecture and settings.

:param architecture: is the zoned neutral atom architecture
:param settings: is a dictionary with the settings for the compiler
)")
      .def(py::init<const na::zoned::Architecture&>(), py::keep_alive<1, 2>(),
           "architecture"_a, R"(
Create a routing-aware compiler for the given architecture and settings.

:param architecture: is the zoned neutral atom architecture
)")
      .def(
          "compile",
          [](na::zoned::RoutingAwareCompiler& self,
             const qc::QuantumComputation& qc) -> std::string {
            return self.compile(qc).toString();
          },
          "qc"_a,
          R"(
Compile a quantum circuit for the zoned neutral atom architecture.

:param qc: is the quantum circuit
:returns: the compilation results as an NAComputation.
)")
      .def(
          "stats",
          [](na::zoned::RoutingAwareCompiler& self) -> nlohmann::json {
            return self.getStatistics();
          },
          R"(
Get the statistics of the last compilation.

:returns: the statistics as a dictionary
)");
}
