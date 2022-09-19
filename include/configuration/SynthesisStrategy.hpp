//
// Created by Sarah on 19.09.2022.
//

#ifndef QMAP_SYNTHESISSTRATEGY_HPP
#define QMAP_SYNTHESISSTRATEGY_HPP

#include <string>

enum class SynthesisStrategy {
    StartLow,
    StartHigh,
    UseMinimizer,
    MinMax,
    SplitIter
};

static std::string toString(SynthesisStrategy strategy) {
    switch (strategy) {
        case SynthesisStrategy::MinMax:
            return "minmax";
        case SynthesisStrategy::StartHigh:
            return "start_high";
        case SynthesisStrategy::StartLow:
            return "start_low";
        case SynthesisStrategy::UseMinimizer:
            return "useminimizer";
        case SynthesisStrategy::SplitIter:
            return "split_iterative";
    }
    return "Error";
}

[[maybe_unused]] static SynthesisStrategy synthesisStrategyFromString(const std::string& strategy) {
    if (strategy == "minmax")
        return SynthesisStrategy::MinMax;
    if (strategy == "start_high")
        return SynthesisStrategy::StartHigh;
    if (strategy == "start_low")
        return SynthesisStrategy::StartLow;
    if (strategy == "useminimizer")
        return SynthesisStrategy::UseMinimizer;
    if (strategy == "split_iterative")
        return SynthesisStrategy::SplitIter;
    return SynthesisStrategy::MinMax;
}
#endif //QMAP_SYNTHESISSTRATEGY_HPP
