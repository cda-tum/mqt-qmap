/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
*/

#ifndef CS_TARGETMETRIC_HPP
#define CS_TARGETMETRIC_HPP
#include <iostream>

namespace cs {
    enum class TargetMetric { GATES,
                              GATES_ONLY_CNOT,
                              DEPTH,
                              FIDELITY };

    static std::string toString(const TargetMetric target) {
        switch (target) {
            case TargetMetric::GATES:
                return "gates";
            case TargetMetric::GATES_ONLY_CNOT:
                return "gates_only_cnot";
            case TargetMetric::DEPTH:
                return "depth";
            case TargetMetric::FIDELITY:
                return "fidelity";
        }
        return "Error";
    }

    [[maybe_unused]] static TargetMetric synthesisTargetFromString(const std::string& target) {
        if (target == "gates") {
            return TargetMetric::GATES;
        }
        if (target == "gates_only_cnot") {
            return TargetMetric::GATES_ONLY_CNOT;
        }
        if (target == "depth") {
            return TargetMetric::DEPTH;
        }
        if (target == "fidelity") {
            return TargetMetric::FIDELITY;
        }
        return TargetMetric::GATES;
    }
} // namespace cs
#endif //CS_TARGETMETRIC_HPP
