//
// Created by Sarah on 19.09.2022.
//

#ifndef QMAP_SYNTHESISMETHOD_HPP
#define QMAP_SYNTHESISMETHOD_HPP
#include <iostream>

enum class SynthesisMethod { Z3,
                             MATHSAT,
                             SMTLibV2,
                             DIMACS };

[[maybe_unused]] static SynthesisMethod synthesisMethodfromString(const std::string& method) {
    if (method == "Z3")
        return SynthesisMethod::Z3;
    if (method == "MATHSAT")
        return SynthesisMethod::MATHSAT;
    if (method == "SMTLibV2")
        return SynthesisMethod::SMTLibV2;
    if (method == "DIMACS")
        return SynthesisMethod::DIMACS;
    return SynthesisMethod::Z3;
}
static std::string toString(const SynthesisMethod method) {
    switch (method) {
        case SynthesisMethod::Z3:
            return "Z3";
        case SynthesisMethod::MATHSAT:
            return "MATHSAT";
        case SynthesisMethod::SMTLibV2:
            return "SMTLibV2";
        case SynthesisMethod::DIMACS:
            return "DIMACS";
    }
    return "Error";
}
#endif //QMAP_SYNTHESISMETHOD_HPP
