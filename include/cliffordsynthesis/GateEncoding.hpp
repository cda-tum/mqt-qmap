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

    public:
        static void makeGateEncoding(const SynthesisData& data, const Configuration& configuration);
    };
} // namespace cs
#endif //CS_GATEENCODING_HPP
