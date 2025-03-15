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
#include <functional>
#include <iostream>
#include <memory>
#include <nlohmann/json_fwd.hpp>

namespace na {
#define self (*static_cast<ConcreteType*>(this))
template <class ConcreteType, class... Mixins>
class Compiler : protected Mixins... {
  friend ConcreteType;

public:
  struct Statistics {
    std::chrono::microseconds schedulingTime;
    std::chrono::microseconds reuseAnalysisTime;
    std::chrono::microseconds placementTime;
    std::chrono::microseconds routingTime;
    std::chrono::microseconds codeGenerationTime;
    std::chrono::microseconds totalTime;
    [[nodiscard]] auto asJson() const -> nlohmann::json {
      return {{"scheduling_time", schedulingTime.count()},
              {"reuse_analysis_time", reuseAnalysisTime.count()},
              {"placement_time", placementTime.count()},
              {"routing_time", routingTime.count()},
              {"code_generation_time", codeGenerationTime.count()},
              {"total_time", totalTime.count()}};
    }
  };

private:
  std::reference_wrapper<const Architecture> architecture_;
  nlohmann::json config_;
  Statistics statistics_;

  Compiler(const Architecture& architecture, const nlohmann::json& config)
      : Mixins(architecture, config)..., architecture_(architecture),
        config_(config) {}

public:
  [[nodiscard]] auto compile(const qc::QuantumComputation& qComp)
      -> NAComputation {
    std::cout << "\033[1;32m[INFO]\033[0m AZAC: An advanced compiler for zoned "
                 "neutral atom "
                 "architecture\n";
    std::cout << "\033[1;32m[INFO]\033[0m           Number of qubits: "
              << qComp.getNqubits() << "\n";
    const auto nTwoQubitGates =
        std::count_if(qComp.cbegin(), qComp.cend(),
                      [](const std::unique_ptr<qc::Operation>& op) {
                        return op->getNqubits() == 2;
                      });
    const auto nOneQubitGates =
        std::count_if(qComp.cbegin(), qComp.cend(),
                      [](const std::unique_ptr<qc::Operation>& op) {
                        return op->getNqubits() == 1;
                      });
    std::cout << "\033[1;32m[INFO]\033[0m           Number of two-qubit gates: "
              << nTwoQubitGates << "\n";
    std::cout
        << "\033[1;32m[INFO]\033[0m           Number of single-qubit gates: "
        << nOneQubitGates << "\n";

    const auto& schedulingStart = std::chrono::system_clock::now();
    const auto& [oneQubitGateLayers, twoQubitGateLayers] = self.schedule(qComp);
    statistics_.schedulingTime =
        std::chrono::system_clock::now() - schedulingStart;
    std::cout << "\033[1;32m[INFO]\033[0m           Time for scheduling: "
              << statistics_.schedulingTime.count() << "µs\n";

    const auto& reuseAnalysisStart = std::chrono::system_clock::now();
    const auto& reuseQubits = self.analyzeReuse(twoQubitGateLayers);
    statistics_.reuseAnalysisTime =
        std::chrono::system_clock::now() - reuseAnalysisStart;
    std::cout << "\033[1;32m[INFO]\033[0m           Time for reuse analysis: "
              << statistics_.reuseAnalysisTime.count() << "µs\n";

    const auto& placementStart = std::chrono::system_clock::now();
    const auto& placement = static_cast<ConcreteType*>(this)->place(
        qComp.getNqubits(), twoQubitGateLayers, reuseQubits);
    statistics_.placementTime =
        std::chrono::system_clock::now() - placementStart;
    std::cout << "\033[1;32m[INFO]\033[0m           Time for placement: "
              << statistics_.placementTime.count() << "µs\n";

    const auto& routingStart = std::chrono::system_clock::now();
    const auto& routing = self.route(placement);
    statistics_.routingTime = std::chrono::system_clock::now() - routingStart;
    std::cout << "\033[1;32m[INFO]\033[0m           Time for routing: "
              << statistics_.routingTime.count() << "µs\n";

    const auto& codeGenerationStart = std::chrono::system_clock::now();
    NAComputation code = self.generate(oneQubitGateLayers, placement, routing);
    assert(code.validate().first);
    statistics_.codeGenerationTime =
        std::chrono::system_clock::now() - codeGenerationStart;
    std::cout << "\033[1;32m[INFO]\033[0m           Time for code generation: "
              << statistics_.codeGenerationTime.count() << "µs\n";

    statistics_.totalTime = std::chrono::system_clock::now() - schedulingStart;
    std::cout << "\033[1;32m[INFO]\033[0m           Total time: "
              << statistics_.totalTime.count() << "µs\n";
    return code;
  }
  [[nodiscard]] auto getStatistics() const -> const Statistics& {
    return statistics_;
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
