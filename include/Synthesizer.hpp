//
// Created by Sarah on 21.09.2022.
//

#ifndef QMAP_SYNTHESIZER_HPP
#define QMAP_SYNTHESIZER_HPP

#include "Architecture.hpp"
#include "QuantumComputation.hpp"
#include "SynthesisResults.hpp"
#include "configuration/SynthesisConfiguration.hpp"

class Synthesizer {
protected:
    Architecture architecture;

    SynthesisResults results{};

    qc::QuantumComputation resultCircuit{};

    std::vector<CouplingMap> highestFidelityCouplingMap;

    virtual void initResults();

    virtual void initCouplingMap(std::uint32_t nqubits);

    virtual void updateResults(SynthesisResults& results) = 0;

public:
    Synthesizer(Architecture architecture);
    Synthesizer()                                                        = default;
    virtual void synthesize(const SynthesisConfiguration& configuration) = 0;
};

#endif //QMAP_SYNTHESIZER_HPP
