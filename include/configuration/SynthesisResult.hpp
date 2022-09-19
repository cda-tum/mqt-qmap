//
// Created by Sarah on 19.09.2022.
//

#ifndef QMAP_SYNTHESISRESULT_HPP
#define QMAP_SYNTHESISRESULT_HPP

#include <iostream>

enum class SynthesisResult { SAT,
                             UNSAT,
                             UNDEF };

static std::string toString(const SynthesisResult result) {
    switch (result) {
        case SynthesisResult::SAT:
            return "sat";
        case SynthesisResult::UNSAT:
            return "unsat";
        case SynthesisResult::UNDEF:
            return "undef";
    }
    return "Error";
}

[[maybe_unused]] static SynthesisResult synthesisResultFromString(const std::string& result) {
    if (result == "sat")
        return SynthesisResult::SAT;
    if (result == "unsat")
        return SynthesisResult::UNSAT;
    if (result == "undef")
        return SynthesisResult::UNDEF;
    return SynthesisResult::UNDEF;
}
#endif //QMAP_SYNTHESISRESULT_HPP
