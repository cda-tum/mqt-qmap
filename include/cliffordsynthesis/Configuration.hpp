/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
*/

#ifndef CS_CONFIGURATION_HPP
#define CS_CONFIGURATION_HPP

#include "OptimizationStrategy.hpp"
#include "ReasoningEngine.hpp"
#include "TargetMetric.hpp"
#include "cliffordsynthesis/Tableau.hpp"
#include "nlohmann/json.hpp"

namespace cs {
    struct Configuration {
        Configuration() = default;

        explicit Configuration(bool chooseBest, bool useEmbedding, std::uint8_t nqubits, std::uint16_t initialTimesteps, std::uint8_t nThreads,
                               std::uint8_t verbosity, OptimizationStrategy strategy, TargetMetric target, ReasoningEngine method):
            chooseBest(chooseBest),
            useEmbedding(useEmbedding), nqubits(nqubits), initialTimesteps(initialTimesteps), nThreads(nThreads), verbosity(verbosity), strategy(strategy), target(target), method(method) {}
        explicit Configuration(bool chooseBest, bool useEmbedding, std::uint8_t nqubits, std::uint16_t initialTimesteps, OptimizationStrategy strategy, TargetMetric target):
            chooseBest(chooseBest), useEmbedding(useEmbedding), nqubits(nqubits), initialTimesteps(initialTimesteps), nThreads(2), verbosity(1), strategy(strategy), target(target), method(ReasoningEngine::Z3) {}

        Configuration(const Configuration& other):
            chooseBest(other.chooseBest), useEmbedding(other.useEmbedding), nqubits(other.nqubits), initialTimesteps(other.initialTimesteps), nThreads(other.nThreads), verbosity(other.verbosity), strategy(other.strategy), target(other.target), method(other.method), targetTableau(other.targetTableau), initialTableau(other.initialTableau) {
            this->targetCircuit = other.targetCircuit.clone();
        }

        bool                 chooseBest       = false;
        bool                 useEmbedding     = false;
        std::uint8_t         nqubits          = 0;
        std::uint16_t        initialTimesteps = 0;
        std::uint8_t         nThreads         = 1;
        std::uint8_t         verbosity        = 0;
        OptimizationStrategy strategy         = OptimizationStrategy::UseMinimizer;
        TargetMetric         target           = TargetMetric::GATES;
        ReasoningEngine      method           = ReasoningEngine::Z3;

        qc::QuantumComputation targetCircuit{};
        Tableau                targetTableau{};
        Tableau                initialTableau{};

        Architecture architecture{};

        [[nodiscard]] nlohmann::json json() const {
            nlohmann::json j;
            j["choose_best"]       = chooseBest;
            j["use_embedding"]     = useEmbedding;
            j["nqubits"]          = nqubits;
            j["initial_timesteps"] = initialTimesteps;
            j["verbosity"]        = verbosity;
            j["strategy"]         = strategy;
            j["target"]           = target;
            j["method"]           = method;
            j["target_tableau"]    = targetTableau.toString();
            j["initial_tableau"]   = initialTableau.toString();
            j["architecture"]     = architecture.getName();
            return j;
        }
        [[nodiscard]] std::string toString() const {
            return json().dump(2);
        }
    };
} // namespace cs

#endif //CS_CONFIGURATION_HPP
