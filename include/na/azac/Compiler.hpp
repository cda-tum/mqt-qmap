#pragma once

#include "Architecture.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/Operation.hpp"
#include "na/NAComputation.hpp"
#include "na/azac/ASAPScheduler.hpp"
#include "na/azac/AStarPlacer.hpp"
#include "na/azac/CodeGenerator.hpp"
#include "na/azac/ISRouter.hpp"
#include "na/azac/VMPlacer.hpp"
#include "na/azac/VMReuseAnalyzer.hpp"

#include <cassert>
#include <chrono>
#include <iostream>
#include <nlohmann/json_fwd.hpp>
#include <utility>

namespace na {
#define self (*static_cast<ConcreteType*>(this))
template <class ConcreteType, class... Mixins>
class Compiler : public Mixins... {
  friend ConcreteType;
  Architecture architecture_;
  nlohmann::json config_;
  std::chrono::microseconds schedulingTime_;
  std::chrono::microseconds reuseAnalysisTime_;
  std::chrono::microseconds placementTime_;
  std::chrono::microseconds routingTime_;
  std::chrono::microseconds codeGenerationTime_;
  std::chrono::microseconds totalTime_;

  Compiler(Architecture architecture, nlohmann::json config)
      : Mixins(architecture_, config_)...,
        architecture_(std::move(architecture)), config_(std::move(config)) {}

public:
  [[nodiscard]] auto compile(const qc::QuantumComputation& qComp)
      -> NAComputation {
    std::cout << "[INFO] AZAC: An advanced compiler for zoned neutral atom "
                 "architecture\n";
    std::cout << "[INFO]           Number of qubits: " << qComp.getNqubits()
              << "\n";
    const auto nTwoQubitGates = std::count_if(
        qComp.cbegin(), qComp.cend(),
        [](const qc::Operation& op) { return op.getNqubits() == 2; });
    const auto nOneQubitGates = std::count_if(
        qComp.cbegin(), qComp.cend(),
        [](const qc::Operation& op) { return op.getNqubits() == 1; });
    std::cout << "[INFO]           Number of two-qubit gates: "
              << nTwoQubitGates << "\n";
    std::cout << "[INFO]           Number of single-qubit gates: "
              << nOneQubitGates << "\n";

    const auto& schedulingStart = std::chrono::system_clock::now();
    const auto& [oneQubitGateLayers, twoQubitGateLayers] = self.schedule(qComp);
    schedulingTime_ = std::chrono::system_clock::now() - schedulingStart;
    std::cout << "[INFO]           Time for scheduling: "
              << schedulingTime_.count() << "µs\n";

    const auto& reuseAnalysisStart = std::chrono::system_clock::now();
    const auto& reuseQubits = self.analyzeReuse(twoQubitGateLayers);
    reuseAnalysisTime_ = std::chrono::system_clock::now() - reuseAnalysisStart;
    std::cout << "[INFO]           Time for reuse analysis: "
              << reuseAnalysisTime_.count() << "µs\n";

    const auto& placementStart = std::chrono::system_clock::now();
    const auto& placement = static_cast<ConcreteType*>(this)->place(
        qComp.getNqubits(), twoQubitGateLayers, reuseQubits);
    placementTime_ = std::chrono::system_clock::now() - placementStart;
    std::cout << "[INFO]           Time for placement: "
              << placementTime_.count() << "µs\n";

    const auto& routingStart = std::chrono::system_clock::now();
    const auto& routing = self.route(placement);
    routingTime_ = std::chrono::system_clock::now() - routingStart;
    std::cout << "[INFO]           Time for routing: " << routingTime_.count()
              << "µs\n";

    const auto& codeGenerationStart = std::chrono::system_clock::now();
    const NAComputation& code =
        self.generateCode(oneQubitGateLayers, placement, routing);
    assert(code.validate().first);
    codeGenerationTime_ =
        std::chrono::system_clock::now() - codeGenerationStart;
    std::cout << "[INFO]           Time for code generation: "
              << codeGenerationTime_.count() << "µs\n";

    totalTime_ = std::chrono::system_clock::now() - routingStart;
    std::cout << "[INFO]           Total time: " << totalTime_.count()
              << "µs\n";
    return code;
  }
};

class ZACompiler final
    : public Compiler<ZACompiler, ASAPScheduler, VMReuseAnalyzer, VMPlacer,
                      ISRouter, CodeGenerator> {
public:
  ZACompiler(const Architecture& architecture, const nlohmann::json& config)
      : Compiler(architecture, config) {}
};

class AZACompiler final
    : public Compiler<AZACompiler, ASAPScheduler, VMReuseAnalyzer, AStarPlacer,
                      ISRouter, CodeGenerator> {
public:
  AZACompiler(const Architecture& architecture, const nlohmann::json& config)
      : Compiler(architecture, config) {}
};
} // namespace na
