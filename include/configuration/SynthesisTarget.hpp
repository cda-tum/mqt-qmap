//
// Created by Sarah on 19.09.2022.
//

#ifndef QMAP_SYNTHESISTARGET_HPP
#define QMAP_SYNTHESISTARGET_HPP
#include <iostream>

enum class SynthesisTarget { GATES,
                             GATES_ONLY_CNOT,
                             DEPTH,
                             FIDELITY };

static std::string toString(const SynthesisTarget target) {
    switch (target) {
        case SynthesisTarget::GATES:
            return "gates";
        case SynthesisTarget::GATES_ONLY_CNOT:
            return "gates_only_cnot";
        case SynthesisTarget::DEPTH:
            return "depth";
        case SynthesisTarget::FIDELITY:
            return "fidelity";
    }
    return "Error";
}

[[maybe_unused]] static SynthesisTarget synthesisTargetFromString(const std::string& target) {
    if (target == "gates")
        return SynthesisTarget::GATES;
    if (target == "gates_only_cnot")
        return SynthesisTarget::GATES_ONLY_CNOT;
    if (target == "depth")
        return SynthesisTarget::DEPTH;
    if (target == "fidelity")
        return SynthesisTarget::FIDELITY;
    return SynthesisTarget::GATES;
}
#endif //QMAP_SYNTHESISTARGET_HPP
