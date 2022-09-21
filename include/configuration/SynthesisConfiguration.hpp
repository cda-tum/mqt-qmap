//
// Created by Sarah on 21.09.2022.
//

#ifndef QMAP_SYNTHESISCONFIGURATION_HPP
#define QMAP_SYNTHESISCONFIGURATION_HPP

#include "SynthesisMethod.hpp"
#include "SynthesisStrategy.hpp"
#include "SynthesisTarget.hpp"
#include "Tableau.hpp"
#include "nlohmann/json.hpp"

struct SynthesisConfiguration {
    SynthesisConfiguration() = default;

    bool            chooseBest       = false;
    bool useEmbedding = false;
    std::uint8_t   nqubits          = 0;
    std::uint16_t initialTimesteps = 0;
    std::uint8_t verbosity = 0;
    SynthesisStrategy strategy = SynthesisStrategy::UseMinimizer;
    SynthesisTarget   target   = SynthesisTarget::GATES;
    SynthesisMethod   method   = SynthesisMethod::Z3;

    qc::QuantumComputation targetCircuit{};
    Tableau              targetTableau{};
    Tableau initialTableau{};


    [[nodiscard]] nlohmann::json json() const {
        nlohmann::json j;
        j["chooseBest"] = chooseBest;
        j["useEmbedding"] = useEmbedding;
        j["nqubits"] = nqubits;
        j["initialTimesteps"] = initialTimesteps;
        j["verbosity"] = verbosity;
        j["strategy"] = strategy;
        j["target"] = target;
        j["method"] = method;
        j["targetTableau"] = targetTableau.toString();
        j["initialTableau"] = initialTableau.toString();
        return j;
    }
    [[nodiscard]] std::string toString() const {
        return json().dump(2);
    }
};

#endif //QMAP_SYNTHESISCONFIGURATION_HPP
