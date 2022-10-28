/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
*/

#ifndef CS_TARGETMETRICHANDLER_HPP
#define CS_TARGETMETRICHANDLER_HPP
#include "Architecture.hpp"
#include "Configuration.hpp"
#include "Results.hpp"
#include "SynthesisData.hpp"
namespace cs {

    class TargetMetricHandler {
        static void makeGateMetric(const SynthesisData& data, bool onlyCNOT);
        static void makeDepthMetric(const SynthesisData& data);
        static void makeFidelityMetric(const SynthesisData& data, const Architecture& architecture, const std::uint32_t fidelityScaling);

    public:
        static void makeTargetMetric(const SynthesisData& data, const Configuration& configuration);

        static void updateResults(const Configuration& configuration, Results& results, Results& currentResults);
    };
} // namespace cs
#endif //CS_TARGETMETRICHANDLER_HPP
