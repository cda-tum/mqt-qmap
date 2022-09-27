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

        Configuration(const Configuration& other):
            chooseBest(other.chooseBest), nqubits(other.nqubits), useEmbedding(other.useEmbedding), verbosity(other.verbosity), strategy(other.strategy), nThreads(other.nThreads), target(other.target), method(other.method), initialTableau(other.initialTableau), targetTableau(other.targetTableau), initialTimesteps(other.initialTimesteps) {
            this->targetCircuit = other.targetCircuit.clone();
        }

        bool                 chooseBest       = false;
        bool                 useEmbedding     = false;
        std::uint8_t         nqubits          = 0;
        std::uint16_t        initialTimesteps = 0;
        std::uint8_t         nThreads         = 1;
        std::uint8_t         verbosity        = 0;
        OptimizationStrategy strategy         = OptimizationStrategy::UseMinimizer;
        TargetMetric target           = TargetMetric::GATES;
        ReasoningEngine      method           = ReasoningEngine::Z3;

        qc::QuantumComputation targetCircuit{};
        Tableau                targetTableau{};
        Tableau                initialTableau{};

        [[nodiscard]] nlohmann::json json() const {
            nlohmann::json j;
            j["chooseBest"]       = chooseBest;
            j["useEmbedding"]     = useEmbedding;
            j["nqubits"]          = nqubits;
            j["initialTimesteps"] = initialTimesteps;
            j["verbosity"]        = verbosity;
            j["strategy"]         = strategy;
            j["target"]           = TargetMetric::toString(target);
            j["method"]           = method;
            j["targetTableau"]    = targetTableau.toString();
            j["initialTableau"]   = initialTableau.toString();
            return j;
        }
        [[nodiscard]] std::string toString() const {
            return json().dump(2);
        }
    };
} // namespace cs

#endif //CS_CONFIGURATION_HPP
