/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
*/

#ifndef CS_GATEENCODING_HPP
#define CS_GATEENCODING_HPP

#include "CliffordSynthesizer.hpp"

namespace cs {
    class GateEncoding {
        static void makeSingleGateEncoding(
                const SynthesisData& data);
        static void makeMultiGateEncoding(
                const SynthesisData& data);

        static void makeNotChangingSet(const SynthesisData& data, unsigned int a, unsigned int gateStep, encodings::LogicTerm& changes);
        static void makeIConstraints(const SynthesisData& data, unsigned int a, unsigned int gateStep, encodings::LogicTerm& changes);
        static void makeHConstraints(const SynthesisData& data, unsigned int a, unsigned int gateStep, encodings::LogicTerm& changes);
        static void makeSConstraints(const SynthesisData& data, unsigned int a, unsigned int gateStep, encodings::LogicTerm& changes);
        static void makeSdagConstraints(const SynthesisData& data, unsigned int a, unsigned int gateStep, encodings::LogicTerm& changes);
        static void makeXConstraints(const SynthesisData& data, unsigned int a, unsigned int gateStep, encodings::LogicTerm& changes);
        static void makeYConstraints(const SynthesisData& data, unsigned int a, unsigned int gateStep, encodings::LogicTerm& changes);
        static void makeZConstraints(const SynthesisData& data, unsigned int a, unsigned int gateStep, encodings::LogicTerm& changes);
        static void makeCNOTConstraints(const SynthesisData& data, unsigned int a, unsigned int b, unsigned int gateStep, encodings::LogicTerm& changes);

        static void makeSingleGateConsistency(const SynthesisData&);
        static void makeMultiGateConsistency(const SynthesisData&);

    public:
        static void makeGateEncoding(const SynthesisData& data, const Configuration& configuration);
    };
} // namespace cs
#endif //CS_GATEENCODING_HPP
