/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
*/

#ifndef CS_HEURISTICSTRATEGY_HPP
#define CS_HEURISTICSTRATEGY_HPP
#include "Architecture.hpp"
#include "CliffordSynthesizer.hpp"
#include "Configuration.hpp"
namespace cs {
    class HeuristicStrategy {
        static void runSplitIter(const CouplingMap&                reducedCM,
                                 const std::vector<std::uint16_t>& qubitChoice, const Configuration& configuration, CliffordSynthesizer& synthesizer);
        static void runSplinter(int i, unsigned int circSplit, unsigned int split,
                                const CouplingMap&                reducedCM,
                                const std::vector<std::uint16_t>& qubitChoice,
                                qc::QuantumComputation&           circuit,
                                Results* r, CliffordSynthesizer* opt, const Configuration& configuration);

    public:
        static void runHeuristicStrategy(const CouplingMap&                reducedCM,
                                         const std::vector<std::uint16_t>& qubitChoice, const Configuration& configuration, CliffordSynthesizer& synthesizer);
    };

} // namespace cs

#endif //CS_HEURISTICSTRATEGY_HPP
