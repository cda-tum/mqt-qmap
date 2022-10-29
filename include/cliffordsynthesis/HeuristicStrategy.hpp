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
        static void runSplitIter(const CouplingMap& reducedCM, const QubitSubset& qubitChoice,
                                 const Configuration& configuration, CliffordSynthesizer& synthesizer);
        static void runSplinter(int i, std::size_t circSplit, std::size_t split,
                                const CouplingMap&       reducedCM,
                                const QubitSubset&       qubitChoice,
                                qc::QuantumComputation&  circuit,
                                std::shared_ptr<Results> r, CliffordSynthesizer* opt, const Configuration& configuration);

    public:
        static void runHeuristicStrategy(const CouplingMap& reducedCM, const QubitSubset& qubitChoice,
                                         const Configuration& configuration, CliffordSynthesizer& synthesizer);
    };

} // namespace cs

#endif //CS_HEURISTICSTRATEGY_HPP
