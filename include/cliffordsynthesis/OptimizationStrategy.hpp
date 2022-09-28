/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
*/

#ifndef CS_OPTIMIZATIONSTRATEGY_HPP
#define CS_OPTIMIZATIONSTRATEGY_HPP

#include <string>

namespace cs {
    enum class OptimizationStrategy {
        StartLow,
        StartHigh,
        UseMinimizer,
        MinMax,
        SplitIter
    };

    static std::string toString(OptimizationStrategy strategy) {
        switch (strategy) {
            case OptimizationStrategy::MinMax:
                return "minmax";
            case OptimizationStrategy::StartHigh:
                return "start_high";
            case OptimizationStrategy::StartLow:
                return "start_low";
            case OptimizationStrategy::UseMinimizer:
                return "use_minimizer";
            case OptimizationStrategy::SplitIter:
                return "split_iterative";
        }
        return "Error";
    }

    [[maybe_unused]] static OptimizationStrategy optimizationStrategyFromString(const std::string& strategy) {
        if (strategy == "minmax") {
            return OptimizationStrategy::MinMax;
        }
        if (strategy == "start_high") {
            return OptimizationStrategy::StartHigh;
        }
        if (strategy == "start_low") {
            return OptimizationStrategy::StartLow;
        }
        if (strategy == "use_minimizer") {
            return OptimizationStrategy::UseMinimizer;
        }
        if (strategy == "split_iterative") {
            return OptimizationStrategy::SplitIter;
        }
        return OptimizationStrategy::MinMax;
    }

    [[maybe_unused]] static bool isExact(OptimizationStrategy strategy) {
        switch (strategy) {
            case OptimizationStrategy::MinMax:
            case OptimizationStrategy::StartHigh:
            case OptimizationStrategy::StartLow:
            case OptimizationStrategy::UseMinimizer:
                return true;
            case OptimizationStrategy::SplitIter:
                return false;
        }
        return false;
    }
} // namespace cs
#endif //CS_OPTIMIZATIONSTRATEGY_HPP
