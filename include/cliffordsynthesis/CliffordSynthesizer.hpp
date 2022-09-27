/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
*/

#ifndef CS_CLIFFORDSYNTHESIS_H
#define CS_CLIFFORDSYNTHESIS_H

#include "Architecture.hpp"
#include "Configuration.hpp"
#include "Encodings/Encodings.hpp"
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

        void synthesize(const Configuration& configuration);
        void optimize(Configuration& configuration);

        void dumpResult(const std::string& outputFilename) {
            if (optimalResults.resultCircuit.empty()) {
                std::cerr << "Circuit is empty." << std::endl;
                return;
            }

            size_t      dot       = outputFilename.find_last_of('.');
            std::string extension = outputFilename.substr(dot + 1U);
            std::transform(extension.begin(), extension.end(), extension.begin(), [](unsigned char c) { return ::tolower(c); });
            dumpResult(outputFilename, qc::OpenQASM);
        }

        void dumpResult(const std::string& outputFilename, qc::Format format) {
            optimalResults.resultCircuit.dump(outputFilename, format);
        }

        void dumpResult(std::ostream& os, qc::Format format) {
            optimalResults.resultCircuit.dump(os, format);
        }

        Results optimalResults{};

        Results mainOptimization(
                std::uint32_t                                            timesteps,
                const std::set<std::pair<std::uint16_t, std::uint16_t>>& reducedCM,
                const std::vector<std::uint16_t>&                        qubitChoice,
                const Tableau& targetTableau, const Tableau& initialTableau,
                const Configuration& configuration);

    protected:
        qc::QuantumComputation resultCircuit{};

        std::vector<CouplingMap> highestFidelityCouplingMap;

        virtual void initResults();

        virtual void initCouplingMap(const Configuration& configuration);

        void        runSplitIter(const CouplingMap&                reducedCM,
                                 const std::vector<std::uint16_t>& qubitChoice, const Configuration& configuration);
        static void runSplinter(int i, unsigned int circSplit, unsigned int split,
                                const CouplingMap&                reducedCM,
                                const std::vector<std::uint16_t>& qubitChoice,
                                qc::QuantumComputation&           circuit,
                                Results* r, CliffordSynthesizer* opt, const Configuration& configuration);

        static void assertTableau(const SynthesisData& data, const Tableau& tableau, std::uint32_t position);
    };

    class Gates {
    public:
        enum GATES {
            NOP,
            H,
            S,
            X,
            Y,
            Z,
            Sdag,
            CX,
        };

        static std::string gateName(GATES gate) {
            switch (gate) {
                case NOP:
                    return "NOP";
                case H:
                    return "H";
                case S:
                    return "S";
                case X:
                    return "X";
                case Y:
                    return "Y";
                case Z:
                    return "Z";
                case Sdag:
                    return "Sdag";
                case CX:
                    return "CX";
                default:
                    return "";
            }
        }

        static int toIndex(GATES gate) {
            switch (gate) {
                case NOP:
                    return 0;
                case H:
                    return 1;
                case S:
                    return 2;
                case X:
                    return 3;
                case Y:
                    return 4;
                case Z:
                    return 5;
                case Sdag:
                    return 6;
                case CX:
                    return 7;
                default:
                    return -1;
            }
        }

        static qc::OpType toOpType(GATES gate) {
            switch (gate) {
                case NOP:
                    return qc::OpType::None;
                case H:
                    return qc::OpType::H;
                case S:
                    return qc::OpType::S;
                case X:
                    return qc::OpType::X;
                case Y:
                    return qc::OpType::Y;
                case Z:
                    return qc::OpType::Z;
                case Sdag:
                    return qc::OpType::Sdag;
                case CX:
                    return qc::OpType::X;
                default:
                    return qc::OpType::None;
            }
        }

        static constexpr auto SINGLE_QUBIT = std::array<GATES, 7>{
                GATES::NOP,
                GATES::H,
                GATES::S,
                GATES::X,
                GATES::Y,
                GATES::Z,
                GATES::Sdag};

        static constexpr auto TWO_QUBIT = std::array<GATES, 3>{
                GATES::CX};

        static constexpr auto SINGLE_QUBIT_WITHOUT_NOP = std::array<GATES, 6>{
                GATES::H,
                GATES::S,
                GATES::X,
                GATES::Y,
                GATES::Z,
                GATES::Sdag};

        static constexpr auto ALL_GATES = std::array<GATES, 10>{
                GATES::H,
                GATES::S,
                GATES::X,
                GATES::Y,
                GATES::Z,
                GATES::Sdag,
                GATES::CX};
    };
} // namespace cs
#endif //CS_CLIFFORDSYNTHESIS_H
