/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
*/

#ifndef CS_EXACTSTRATEGY_HPP
#define CS_EXACTSTRATEGY_HPP
#include "Architecture.hpp"
#include "CliffordSynthesizer.hpp"
#include "Configuration.hpp"
namespace cs {
    class ExactStrategy {
        static void runMaxSat(int timesteps, const CouplingMap& reducedCM,
                              const std::vector<std::uint16_t>& qubitChoice, const Configuration& configuration, CliffordSynthesizer& synthesizer);
        static void runStartLow(int timesteps, const CouplingMap& reducedCM,
                                const std::vector<std::uint16_t>& qubitChoice, const Configuration& configuration, CliffordSynthesizer& synthesizer);
        static void runStartHigh(int timesteps, const CouplingMap& reducedCM,
                                 const std::vector<std::uint16_t>& qubitChoice, const Configuration& configuration, CliffordSynthesizer& synthesizer);
        static void runBinarySearch(int timesteps, const CouplingMap& reducedCM,
                                    const std::vector<std::uint16_t>& qubitChoice, const Configuration& configuration, CliffordSynthesizer& synthesizer);

    public:
        static void runExactStrategy(int timesteps, const CouplingMap& reducedCM,
                                     const std::vector<std::uint16_t>& qubitChoice, const Configuration& configuration, CliffordSynthesizer& synthesizer);
    };

} // namespace cs

#endif //CS_EXACTSTRATEGY_HPP
