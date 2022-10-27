/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
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

        explicit Configuration(bool chooseBest, std::uint8_t nqubits, std::uint16_t initialTimesteps, std::uint8_t nThreads,
                               std::uint8_t verbosity, OptimizationStrategy strategy, TargetMetric target):
            chooseBest(chooseBest),
            nqubits(nqubits), initialTimestep(initialTimesteps), nThreads(nThreads), verbosity(verbosity), strategy(strategy), target(target){}
        explicit Configuration(bool chooseBest, std::uint8_t nqubits, std::uint16_t initialTimesteps, OptimizationStrategy strategy, TargetMetric target):
            chooseBest(chooseBest), nqubits(nqubits), initialTimestep(initialTimesteps), nThreads(2), verbosity(1), strategy(strategy), target(target){}

        Configuration(const Configuration& other):
            chooseBest(other.chooseBest), nqubits(other.nqubits), initialTimestep(other.initialTimestep), nThreads(other.nThreads), verbosity(other.verbosity), strategy(other.strategy), target(other.target), targetTableau(other.targetTableau), initialTableau(other.initialTableau) {
            this->targetCircuit = other.targetCircuit.clone();
        }

        bool                 chooseBest      = false;
        std::uint8_t         nqubits         = 0;
        std::uint16_t        initialTimestep = 0;
        std::uint8_t         nThreads        = 1;
        std::uint8_t         verbosity       = 0;
        OptimizationStrategy strategy        = OptimizationStrategy::UseMinimizer;
        TargetMetric         target          = TargetMetric::GATES;

        qc::QuantumComputation targetCircuit{};
        Tableau                targetTableau{};
        Tableau                initialTableau{};

        Architecture architecture{};

        [[nodiscard]] nlohmann::json json() const {
            nlohmann::json j;
            j["choose_best"]         = chooseBest;
            j["nqubits"]             = nqubits;
            j["initial_timestep"]    = initialTimestep;
            j["verbosity"]           = verbosity;
            j["optimizing_strategy"] = strategy;
            j["target_metric"]       = target;
            j["target_tableau"]      = targetTableau.toString();
            j["initial_tableau"]     = initialTableau.toString();
            j["architecture"]        = architecture.getName();
            return j;
        }
        [[nodiscard]] std::string toString() const {
            return json().dump(2);
        }
    };
} // namespace cs

#endif //CS_CONFIGURATION_HPP
