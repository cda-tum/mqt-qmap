/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
*/

#ifndef CS_TARGETMETRIC_HPP
#define CS_TARGETMETRIC_HPP

#include "Configuration.hpp"

#include <iostream>

namespace cs {
    class TargetMetric {
    public:
        enum class Metric { GATES,
                            GATES_ONLY_CNOT,
                            DEPTH,
                            FIDELITY };

        static std::string toString(const TargetMetric target) {
            switch (target.metric) {
                case TargetMetric::Metric::GATES:
                    return "gates";
                case TargetMetric::Metric::GATES_ONLY_CNOT:
                    return "gates_only_cnot";
                case TargetMetric::Metric::DEPTH:
                    return "depth";
                case TargetMetric::Metric::FIDELITY:
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

        explicit TargetMetric(Metric lMetric):
            metric(lMetric) {}

        const static TargetMetric GATES;
        const static TargetMetric GATES_ONLY_CNOT;
        const static TargetMetric DEPTH;
        const static TargetMetric FIDELITY;

    private:
        Metric metric;
    };

    const TargetMetric TargetMetric::GATES         = TargetMetric(Metric::GATES);
    const TargetMetric TargetMetric::GATES_ONLY_CNOT = TargetMetric(Metric::GATES_ONLY_CNOT);
    const TargetMetric TargetMetric::DEPTH         = TargetMetric(Metric::DEPTH);
    const TargetMetric TargetMetric::FIDELITY      = TargetMetric(Metric::FIDELITY);
} // namespace cs
#endif //CS_TARGETMETRIC_HPP
