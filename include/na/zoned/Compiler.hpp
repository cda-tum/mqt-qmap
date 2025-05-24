/*
 * Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
 * Copyright (c) 2025 Munich Quantum Software Company GmbH
 * All rights reserved.
 *
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License
 */

#pragma once

#include "Architecture.hpp"
#include "code_generator/CodeGenerator.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/Operation.hpp"
#include "na/NAComputation.hpp"
#include "placer/AStarPlacer.hpp"
#include "placer/VertexMatchingPlacer.hpp"
#include "reuse_analyzer/VertexMatchingReuseAnalyzer.hpp"
#include "router/IndependentSetRouter.hpp"
#include "scheduler/ASAPScheduler.hpp"

#include <cassert>
#include <chrono>
#include <functional>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>

namespace na::zoned {
#define SELF (*static_cast<ConcreteType*>(this))
template <class ConcreteType, class Scheduler, class ReuseAnalyzer,
          class Placer, class Router, class CodeGenerator>
class Compiler : protected Scheduler,
                 protected ReuseAnalyzer,
                 protected Placer,
                 protected Router,
                 protected CodeGenerator {
  friend ConcreteType;

public:
  struct Config {
    typename Scheduler::Config schedulerConfig{};
    typename ReuseAnalyzer::Config reuseAnalyzerConfig{};
    typename Placer::Config placerConfig{};
    typename Router::Config routerConfig{};
    typename CodeGenerator::Config codeGeneratorConfig{};
    spdlog::level::level_enum logLevel = spdlog::level::info;
    NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(Config, schedulerConfig,
                                                reuseAnalyzerConfig,
                                                placerConfig, routerConfig,
                                                codeGeneratorConfig, logLevel);
  };
  struct Statistics {
    std::chrono::microseconds schedulingTime;
    std::chrono::microseconds reuseAnalysisTime;
    std::chrono::microseconds placementTime;
    std::chrono::microseconds routingTime;
    std::chrono::microseconds codeGenerationTime;
    std::chrono::microseconds totalTime;
    NLOHMANN_DEFINE_TYPE_INTRUSIVE_ONLY_SERIALIZE(Statistics, schedulingTime,
                                                  reuseAnalysisTime,
                                                  placementTime, routingTime,
                                                  codeGenerationTime,
                                                  totalTime);
  };

private:
  std::reference_wrapper<const Architecture> architecture_;
  nlohmann::json config_;
  Statistics statistics_;

  Compiler(const Architecture& architecture, const Config& config)
      : Scheduler(architecture, config.schedulerConfig),
        ReuseAnalyzer(architecture, config.reuseAnalyzerConfig),
        Placer(architecture, config.placerConfig),
        Router(architecture, config.routerConfig),
        CodeGenerator(architecture, config.codeGeneratorConfig),
        architecture_(architecture), config_(config) {
    spdlog::set_level(config.logLevel);
  }

  explicit Compiler(const Architecture& architecture)
      : Compiler(architecture, {}) {}

public:
  [[nodiscard]] auto compile(const qc::QuantumComputation& qComp)
      -> NAComputation {
    SPDLOG_INFO("*** MQT QMAP Zoned Neutral Atom Compiler ***");
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_DEBUG
    if (spdlog::should_log(spdlog::level::debug)) {
      SPDLOG_DEBUG("Used compiler settings:");
      std::string jsonStr =
          config_.dump(2); // Pretty-print with 2-space indentation
      std::istringstream iss(jsonStr);
      std::string line;
      while (std::getline(iss, line)) {
        SPDLOG_DEBUG(line);
      }
      SPDLOG_DEBUG("Number of qubits: {}", qComp.getNqubits());
      const auto nTwoQubitGates = static_cast<size_t>(
          std::count_if(qComp.cbegin(), qComp.cend(),
                        [](const std::unique_ptr<qc::Operation>& op) {
                          return op->getNqubits() == 2;
                        }));
      SPDLOG_DEBUG("Number of two-qubit gates: {}", nTwoQubitGates);
      const auto nSingleQubitGates = static_cast<size_t>(
          std::count_if(qComp.cbegin(), qComp.cend(),
                        [](const std::unique_ptr<qc::Operation>& op) {
                          return op->getNqubits() == 1;
                        }));
      SPDLOG_DEBUG("Number of single-qubit gates: {}", nSingleQubitGates);
    }
#endif // SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_DEBUG

    const auto& schedulingStart = std::chrono::system_clock::now();
    // CodeQL was not very happy about the structural binding here, hence I
    // removed it.
    const auto& schedule = SELF.schedule(qComp);
    const auto& singleQubitGateLayers = schedule.first;
    const auto& twoQubitGateLayers = schedule.second;
    statistics_.schedulingTime =
        std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::system_clock::now() - schedulingStart);
    SPDLOG_INFO("Time for scheduling: {}us",
                statistics_.schedulingTime.count());
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_DEBUG
    SPDLOG_DEBUG("Number of single-qubit gate layers: {}",
                 singleQubitGateLayers.size());
    SPDLOG_DEBUG("Number of two-qubit gate layers: {}",
                 twoQubitGateLayers.size());
    if (!twoQubitGateLayers.empty() &&
        spdlog::should_log(spdlog::level::debug)) {
      const auto& [min, sum, max] = std::accumulate(
          twoQubitGateLayers.cbegin(), twoQubitGateLayers.cend(),
          std::array<size_t, 3>{std::numeric_limits<size_t>::max(), 0UL, 0UL},
          [](const auto& acc, const auto& layer) -> std::array<size_t, 3> {
            const auto& [minAcc, sumAcc, maxAcc] = acc;
            const auto n = layer.size();
            return {std::min(minAcc, n), sumAcc + n, std::max(maxAcc, n)};
          });
      const auto avg = static_cast<double>(sum) /
                       static_cast<double>(twoQubitGateLayers.size());
      SPDLOG_DEBUG(
          "Number of two-qubit gates per layer: min: {}, avg: {}, max: {}", min,
          avg, max);
    }
#endif // SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_DEBUG

    const auto& reuseAnalysisStart = std::chrono::system_clock::now();
    const auto& reuseQubits = SELF.analyzeReuse(twoQubitGateLayers);
    statistics_.reuseAnalysisTime =
        std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::system_clock::now() - reuseAnalysisStart);
    SPDLOG_INFO("Time for reuse analysis: {}us",
                statistics_.reuseAnalysisTime.count());

    const auto& placementStart = std::chrono::system_clock::now();
    const auto& placement =
        SELF.place(qComp.getNqubits(), twoQubitGateLayers, reuseQubits);
    statistics_.placementTime =
        std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::system_clock::now() - placementStart);
    SPDLOG_INFO("Time for placement: {}us", statistics_.placementTime.count());

    const auto& routingStart = std::chrono::system_clock::now();
    const auto& routing = SELF.route(placement);
    statistics_.routingTime =
        std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::system_clock::now() - routingStart);
    SPDLOG_INFO("Time for routing: {}us", statistics_.routingTime.count());

    const auto& codeGenerationStart = std::chrono::system_clock::now();
    NAComputation code =
        SELF.generate(singleQubitGateLayers, placement, routing);
    assert(code.validate().first);
    statistics_.codeGenerationTime =
        std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::system_clock::now() - codeGenerationStart);
    SPDLOG_INFO("Time for code generation: {}us",
                statistics_.codeGenerationTime.count());

    statistics_.totalTime =
        std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::system_clock::now() - schedulingStart);
    SPDLOG_INFO("Total time: {}us", statistics_.totalTime.count());
    return code;
  }
  [[nodiscard]] auto getStatistics() const -> const Statistics& {
    return statistics_;
  }
};

class RoutingAgnosticCompiler final
    : public Compiler<RoutingAgnosticCompiler, ASAPScheduler,
                      VertexMatchingReuseAnalyzer, VertexMatchingPlacer,
                      IndependentSetRouter, CodeGenerator> {
public:
  RoutingAgnosticCompiler(const Architecture& architecture,
                          const Config& config)
      : Compiler(architecture, config) {}
  RoutingAgnosticCompiler(const Architecture& architecture)
      : Compiler(architecture) {}
};

class RoutingAwareCompiler final
    : public Compiler<RoutingAwareCompiler, ASAPScheduler,
                      VertexMatchingReuseAnalyzer, AStarPlacer,
                      IndependentSetRouter, CodeGenerator> {
public:
  RoutingAwareCompiler(const Architecture& architecture, const Config& config)
      : Compiler(architecture, config) {}
  RoutingAwareCompiler(const Architecture& architecture)
      : Compiler(architecture) {}
};
} // namespace na::zoned

NLOHMANN_JSON_NAMESPACE_BEGIN
template <> struct adl_serializer<std::chrono::microseconds> {
  static void to_json(json& j, const std::chrono::microseconds& ms) {
    j = ms.count();
  }
};
NLOHMANN_JSON_NAMESPACE_END
