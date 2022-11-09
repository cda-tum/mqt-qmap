/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
*/

#ifndef CS_CLIFFORDSYNTHESIS_H
#define CS_CLIFFORDSYNTHESIS_H

#include "Architecture.hpp"
#include "Configuration.hpp"
#include "Encodings/Encodings.hpp"
#include "Gates.hpp"
#include "LogicBlock/LogicBlock.hpp"
#include "QuantumComputation.hpp"
#include "Results.hpp"
#include "SynthesisData.hpp"
#include "operations/OpType.hpp"
#include "operations/StandardOperation.hpp"

#include <algorithm>
#include <bitset>
#include <cassert>
#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <istream>
#include <iterator>
#include <list>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <thread>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>
namespace cs {
    class CliffordSynthesizer {
    public:
        Tableau modelTableau{};

        virtual ~CliffordSynthesizer() = default;
        CliffordSynthesizer()          = default;

        void synthesize(Configuration& configuration);

        void dumpResult(const std::string& outputFilename, qc::Format format) {
            qc::QuantumComputation splitResult;
            std::istringstream     iss(optimalResults.resultStringCircuit);
            splitResult.import(iss, qc::OpenQASM);
            splitResult.dump(outputFilename, format);
        }

        Results optimalResults{};

        Results mainOptimization(
                const std::size_t    timesteps,
                const CouplingMap&   reducedCM,
                const QubitSubset&   qubitChoice,
                const Tableau&       targetTableau,
                const Tableau&       initialTableau,
                const Configuration& configuration);

    protected:
        qc::QuantumComputation resultCircuit{};

        std::vector<CouplingMap> couplingMaps{};

        void initConfiguration(Configuration& configuration);

        virtual void initResults();

        virtual void initCouplingMaps(const Configuration& configuration);

        static void assertTableau(const SynthesisData& data, const Tableau& tableau, std::size_t position);
    };

} // namespace cs
#endif //CS_CLIFFORDSYNTHESIS_H
