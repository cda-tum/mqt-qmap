/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
*/

#ifndef CS_RESONINGENGINE_HPP
#define CS_RESONINGENGINE_HPP
#include <iostream>

namespace cs {
    enum class ReasoningEngine { Z3,
                                 MATHSAT,
                                 SMTLibV2,
                                 DIMACS };

    [[maybe_unused]] static ReasoningEngine reasoningEngineFromString(const std::string& method) {
        if (method == "Z3") {
            return ReasoningEngine::Z3;
        }
        if (method == "MATHSAT") {
            return ReasoningEngine::MATHSAT;
        }
        if (method == "SMTLibV2") {
            return ReasoningEngine::SMTLibV2;
        }
        if (method == "DIMACS") {
            return ReasoningEngine::DIMACS;
        }
        return ReasoningEngine::Z3;
    }
    static std::string toString(const ReasoningEngine method) {
        switch (method) {
            case ReasoningEngine::Z3:
                return "Z3";
            case ReasoningEngine::MATHSAT:
                return "MATHSAT";
            case ReasoningEngine::SMTLibV2:
                return "SMTLibV2";
            case ReasoningEngine::DIMACS:
                return "DIMACS";
        }
        return "Error";
    }
} // namespace cs
#endif //CS_RESONINGENGINE_HPP
