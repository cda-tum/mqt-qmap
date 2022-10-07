/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
*/

#include "cliffordsynthesis/ExactStrategy.hpp"

#include "cliffordsynthesis/TargetMetricHandler.hpp"
namespace cs {
    void cs::ExactStrategy::runExactStrategy(std::size_t timesteps, const CouplingMap& reducedCM,
                                             const QubitSubset& qubitChoice, const Configuration& configuration, CliffordSynthesizer& synthesizer) {
        switch (configuration.strategy) {
            case OptimizationStrategy::UseMinimizer:
                runMaxSat(timesteps, reducedCM, qubitChoice, configuration, synthesizer);
                break;
            case OptimizationStrategy::StartLow:
                runStartLow(timesteps, reducedCM, qubitChoice, configuration, synthesizer);
                break;
            case OptimizationStrategy::StartHigh:
                runStartHigh(timesteps, reducedCM, qubitChoice, configuration, synthesizer);
                break;
            case OptimizationStrategy::MinMax:
                runBinarySearch(timesteps, reducedCM, qubitChoice, configuration, synthesizer);
                break;
            default:
                throw std::runtime_error("Unknown optimization strategy");
        }
    }
    void ExactStrategy::runMaxSat(std::size_t timesteps, const CouplingMap& reducedCM, const QubitSubset& qubitChoice, const Configuration& configuration, CliffordSynthesizer& synthesizer) {
        DEBUG() << "Running minimizer" << std::endl;
        Results r = synthesizer.mainOptimization(timesteps, reducedCM, qubitChoice, configuration.targetTableau, configuration.initialTableau,
                                                 configuration);
        TargetMetricHandler::updateResults(configuration, r, synthesizer.optimalResults);
    }
    void ExactStrategy::runStartLow(std::size_t timesteps, const CouplingMap& reducedCM, const QubitSubset& qubitChoice, const Configuration& configuration, CliffordSynthesizer& synthesizer) {
        DEBUG() << "Running start low" << std::endl;
        Results r;
        while (r.result != logicbase::Result::SAT || r.result == logicbase::Result::NDEF) {
            DEBUG() << "Current t=" << timesteps << std::endl;
            r = synthesizer.mainOptimization(timesteps, reducedCM, qubitChoice, configuration.targetTableau, configuration.initialTableau, configuration);
            TargetMetricHandler::updateResults(configuration, r, synthesizer.optimalResults);
            if (r.result == logicbase::Result::UNSAT) {
                timesteps *= 1.5;
            }
        }
    }
    void ExactStrategy::runStartHigh(std::size_t timesteps, const CouplingMap& reducedCM, const QubitSubset& qubitChoice, const Configuration& configuration, CliffordSynthesizer& synthesizer) {
        DEBUG() << "Running start high" << std::endl;
        Results r;
        auto    oldTimesteps = timesteps;
        while (r.result == logicbase::Result::SAT || r.result == logicbase::Result::NDEF) {
            DEBUG() << "Current t=" << timesteps << std::endl;
            r = synthesizer.mainOptimization(timesteps, reducedCM, qubitChoice, configuration.targetTableau, configuration.initialTableau, configuration);
            TargetMetricHandler::updateResults(configuration, r, synthesizer.optimalResults);
            if (r.result == logicbase::Result::SAT) {
                oldTimesteps = timesteps;
                timesteps *= 0.5;
            } else {
                timesteps = oldTimesteps;
            }
        }
    }
    void ExactStrategy::runBinarySearch(std::size_t timesteps, const CouplingMap& reducedCM, const QubitSubset& qubitChoice, const Configuration& configuration, CliffordSynthesizer& synthesizer) {
        DEBUG() << "Running minmax" << std::endl;
        Results r;
        auto    t     = static_cast<int>(timesteps);
        auto    upper = static_cast<int>(timesteps);
        int     lower = 0;
        while (std::abs(upper - lower) > 1) {
            DEBUG() << "Current t=" << t << std::endl;
            r = synthesizer.mainOptimization(t, reducedCM, qubitChoice, configuration.targetTableau, configuration.initialTableau, configuration);
            TargetMetricHandler::updateResults(configuration, r, synthesizer.optimalResults);
            if (r.result == logicbase::Result::SAT) {
                upper = t;
            } else if (r.result == logicbase::Result::UNSAT) {
                lower = t;
            } else {
                break;
            }
            if (upper - lower < 1 && r.result == logicbase::Result::UNSAT) {
                upper *= 1.5;
            }
            t = lower + std::abs(upper - lower) / 2;
        }
    }
} // namespace cs
