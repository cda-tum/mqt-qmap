/*
 * This file is part of the MQT QMAP library which is released under the MIT
 * license. See file README.md or go to
 * https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
 */

#ifndef CS_CONFIGURATION_HPP
#define CS_CONFIGURATION_HPP

#include "OptimizationStrategy.hpp"
#include "TargetMetric.hpp"
#include "cliffordsynthesis/Tableau.hpp"
#include "nlohmann/json.hpp"

namespace cs {
struct Configuration {
  Configuration() = default;

  explicit Configuration(bool chooseBest, std::uint8_t nqubits,
                         std::uint16_t initialTimesteps,
                         std::uint32_t fidelityScaling,
                         double        limitFindingFactor,
                         double circuitSplittingIncrease, std::uint8_t nThreads,
                         std::uint8_t verbosity, OptimizationStrategy strategy,
                         TargetMetric target)
      : chooseBest(chooseBest), nqubits(nqubits),
        initialTimestep(initialTimesteps), fidelityScaling(fidelityScaling),
        limitFindingFactor(limitFindingFactor),
        circuitSplittingIncrease(circuitSplittingIncrease), nThreads(nThreads),
        verbosity(verbosity), strategy(strategy), target(target) {}
  explicit Configuration(bool chooseBest, std::uint8_t nqubits,
                         std::uint16_t        initialTimesteps,
                         OptimizationStrategy strategy, TargetMetric target)
      : chooseBest(chooseBest), nqubits(nqubits),
        initialTimestep(initialTimesteps), strategy(strategy), target(target) {}

  bool                 chooseBest               = false;
  std::uint8_t         nqubits                  = 0;
  std::uint16_t        initialTimestep          = 0;
  std::uint32_t        fidelityScaling          = 1000;
  double               limitFindingFactor       = 0.5;
  double               circuitSplittingIncrease = 0.2;
  std::uint8_t         nThreads                 = 1;
  std::uint8_t         verbosity                = 0;
  OptimizationStrategy strategy = OptimizationStrategy::BinarySearch;
  TargetMetric         target   = TargetMetric::TWO_QUBIT_GATES;

  std::shared_ptr<qc::QuantumComputation> targetCircuit{};
  std::shared_ptr<Tableau>                targetTableau{};
  std::shared_ptr<Tableau>                initialTableau{};

  Architecture architecture{};

  [[nodiscard]] nlohmann::json json() const {
    nlohmann::json j;
    j["choose_best"]                = chooseBest;
    j["nqubits"]                    = nqubits;
    j["initial_timestep"]           = initialTimestep;
    j["fidelity_scaling"]           = fidelityScaling;
    j["limit_finding_factor"]       = limitFindingFactor;
    j["circuit_splitting_increase"] = circuitSplittingIncrease;
    j["verbosity"]                  = verbosity;
    j["optimizing_strategy"]        = strategy;
    j["target_metric"]              = target;
    j["target_tableau"]             = targetTableau->toString();
    j["initial_tableau"]            = initialTableau->toString();
    j["architecture"]               = architecture.getName();
    return j;
  }
  [[nodiscard]] std::string toString() const { return json().dump(2); }
};
} // namespace cs

#endif // CS_CONFIGURATION_HPP
