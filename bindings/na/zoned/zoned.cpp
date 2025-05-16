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
#include "na/zoned/AStarPlacer.hpp"
#include "na/zoned/Architecture.hpp"
#include "na/zoned/CodeGenerator.hpp"
#include "na/zoned/Compiler.hpp"
#include "na/zoned/VMPlacer.hpp"

#include <nlohmann/json_fwd.hpp>
#include <pybind11/attr.h>
#include <pybind11/cast.h>
#include <pybind11/detail/common.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11_json/pybind11_json.hpp>
#include <string>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(MQT_QMAP_MODULE_NAME, m, py::mod_gil_not_used()) {

  py::class_<na::zoned::Architecture> architecture(
      m, "ZonedNeutralAtomArchitecture");
  architecture.def_static("from_json_file",
                          &na::zoned::Architecture::fromJSONFile, "filename"_a);
  architecture.def_static("from_json_string",
                          &na::zoned::Architecture::fromJSONString, "json"_a);

  py::class_<na::zoned::VMPlacer::Config> vmPlacerConfig(m, "VMPlacerConfig");
  vmPlacerConfig.def(py::init<>());
  vmPlacerConfig.def_readwrite("use_window",
                               &na::zoned::VMPlacer::Config::useWindow);
  vmPlacerConfig.def_readwrite("window_size",
                               &na::zoned::VMPlacer::Config::windowSize);
  vmPlacerConfig.def_readwrite("dynamic_placement",
                               &na::zoned::VMPlacer::Config::dynamicPlacement);

  py::class_<na::zoned::CodeGenerator::Config> codeGeneratorConfig(
      m, "CodeGeneratorConfig");
  codeGeneratorConfig.def(py::init<>());
  codeGeneratorConfig.def_readwrite(
      "parking_offset", &na::zoned::CodeGenerator::Config::parkingOffset);
  codeGeneratorConfig.def_readwrite(
      "warn_unsupported_gates",
      &na::zoned::CodeGenerator::Config::warnUnsupportedGates);

  py::class_<na::zoned::AStarPlacer::Config> aStarPlacerConfig(
      m, "AStarPlacerConfig");
  aStarPlacerConfig.def(py::init<>());
  aStarPlacerConfig.def_readwrite("use_window",
                                  &na::zoned::AStarPlacer::Config::useWindow);
  aStarPlacerConfig.def_readwrite(
      "window_min_width", &na::zoned::AStarPlacer::Config::windowMinWidth);
  aStarPlacerConfig.def_readwrite("window_ratio",
                                  &na::zoned::AStarPlacer::Config::windowRatio);
  aStarPlacerConfig.def_readwrite("window_share",
                                  &na::zoned::AStarPlacer::Config::windowShare);
  aStarPlacerConfig.def_readwrite(
      "deepening_factor", &na::zoned::AStarPlacer::Config::deepeningFactor);
  aStarPlacerConfig.def_readwrite(
      "deepening_value", &na::zoned::AStarPlacer::Config::deepeningValue);
  aStarPlacerConfig.def_readwrite(
      "lookahead_factor", &na::zoned::AStarPlacer::Config::lookaheadFactor);
  aStarPlacerConfig.def_readwrite("reuse_level",
                                  &na::zoned::AStarPlacer::Config::reuseLevel);
  aStarPlacerConfig.def_readwrite("max_nodes",
                                  &na::zoned::AStarPlacer::Config::maxNodes);

  py::class_<na::zoned::RoutingAgnosticCompiler> routingAgnosticCompiler(
      m, "RoutingAgnosticCompiler");

  py::class_<na::zoned::RoutingAgnosticCompiler::Config>
      routingAgnosticCompilerConfig(routingAgnosticCompiler, "Config");

  routingAgnosticCompilerConfig.def(py::init<>());
  routingAgnosticCompilerConfig.def_readwrite(
      "placer_config",
      &na::zoned::RoutingAgnosticCompiler::Config::placerConfig);
  routingAgnosticCompilerConfig.def_readwrite(
      "coder_generator_config",
      &na::zoned::RoutingAgnosticCompiler::Config::codeGeneratorConfig);
  routingAgnosticCompilerConfig.def_static(
      "from_json_string",
      [](const std::string& json)
          -> na::zoned::RoutingAgnosticCompiler::Config {
        return nlohmann::json::parse(json);
      },
      "json"_a);

  routingAgnosticCompiler.def(py::init<const na::zoned::Architecture&>(),
                              py::keep_alive<1, 2>(), "arch"_a);
  routingAgnosticCompiler.def(
      py::init<const na::zoned::Architecture&,
               const na::zoned::RoutingAgnosticCompiler::Config&>(),
      py::keep_alive<1, 2>(), "arch"_a, "config"_a);
  routingAgnosticCompiler.def(
      "compile",
      [](na::zoned::RoutingAgnosticCompiler& self,
         const qc::QuantumComputation& qc) -> std::string {
        return self.compile(qc).toString();
      },
      "qc"_a);
  routingAgnosticCompiler.def(
      "stats",
      [](const na::zoned::RoutingAgnosticCompiler& self) -> nlohmann::json {
        return self.getStatistics();
      });

  py::class_<na::zoned::RoutingAwareCompiler> routingAwareCompiler(
      m, "RoutingAwareCompiler");

  py::class_<na::zoned::RoutingAwareCompiler::Config>
      routingAwareCompilerConfig(routingAwareCompiler, "Config");
  routingAwareCompilerConfig.def(py::init<>());
  routingAwareCompilerConfig.def_readwrite(
      "placer_config", &na::zoned::RoutingAwareCompiler::Config::placerConfig);
  routingAwareCompilerConfig.def_readwrite(
      "coder_generator_config",
      &na::zoned::RoutingAwareCompiler::Config::codeGeneratorConfig);
  routingAwareCompilerConfig.def_static(
      "from_json_string",
      [](const std::string& json) -> na::zoned::RoutingAwareCompiler::Config {
        return nlohmann::json::parse(json);
      },
      "json"_a);

  routingAwareCompiler.def(py::init<const na::zoned::Architecture&>(),
                           py::keep_alive<1, 2>(), "arch"_a);
  routingAwareCompiler.def(
      py::init<const na::zoned::Architecture&,
               const na::zoned::RoutingAwareCompiler::Config&>(),
      py::keep_alive<1, 2>(), "arch"_a, "config"_a);
  routingAwareCompiler.def(
      "compile",
      [](na::zoned::RoutingAwareCompiler& self,
         const qc::QuantumComputation& qc) -> std::string {
        return self.compile(qc).toString();
      },
      "qc"_a);
  routingAwareCompiler.def(
      "stats",
      [](const na::zoned::RoutingAwareCompiler& self) -> nlohmann::json {
        return self.getStatistics();
      });
}
