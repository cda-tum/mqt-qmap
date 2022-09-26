/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
*/

#ifndef QMAP_CLIFFORDSYNTHESIS_H
#define QMAP_CLIFFORDSYNTHESIS_H

#include "Architecture.hpp"
#include "Encodings/Encodings.hpp"
#include "LogicBlock/LogicBlock.hpp"
#include "QuantumComputation.hpp"
#include "SynthesisResults.hpp"
#include "Synthesizer.hpp"
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

class CliffordSynthesizer: public Synthesizer {
public:
    struct SynthesisData {
        std::uint32_t                                            nqubits;
        std::uint32_t                                            timesteps;
        const std::set<std::pair<std::uint16_t, std::uint16_t>>& reducedCM;
        const std::vector<std::uint16_t>&                        qubitChoice;
        std::unique_ptr<logicbase::LogicBlock>&                  lb;
        const logicbase::LogicMatrix&                            x;
        const logicbase::LogicMatrix&                            z;
        const logicbase::LogicVector&                            r;
        const logicbase::LogicMatrix3D&                          gS;
        const logicbase::LogicMatrix3D&                          gC;
    };

    Tableau modelTableau{};

    virtual ~CliffordSynthesizer() = default;
    CliffordSynthesizer()          = default;

    void synthesize(const SynthesisConfiguration& configuration) override;

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

    SynthesisResults optimalResults{};

protected:
    Architecture     architecture{};
    SynthesisResults mainOptimization(
            std::uint32_t timesteps,
            const std::set<std::pair<std::uint16_t, std::uint16_t>>& reducedCM,
            const std::vector<std::uint16_t>& qubitChoice, const SynthesisConfiguration& configuration);

    virtual void makeSpecificEncoding(const SynthesisData& data) = 0;

    void        runMinimizer(int timesteps, const CouplingMap& reducedCM,
                             const std::vector<std::uint16_t>& qubitChoice, const SynthesisConfiguration& configuration);
    void        runStartLow(int timesteps, const CouplingMap& reducedCM,
                            const std::vector<std::uint16_t>& qubitChoice, const SynthesisConfiguration& configuration);
    void        runStartHigh(int timesteps, const CouplingMap& reducedCM,
                             const std::vector<std::uint16_t>& qubitChoice, const SynthesisConfiguration& configuration);
    void        runMinMax(int timesteps, const CouplingMap& reducedCM,
                          const std::vector<std::uint16_t>& qubitChoice, const SynthesisConfiguration& configuration);
    void        runSplitIter(const CouplingMap&                reducedCM,
                             const std::vector<std::uint16_t>& qubitChoice, const SynthesisConfiguration& configuration);
    static void runSplinter(int i, unsigned int circSplit, unsigned int split,
                            const CouplingMap&                reducedCM,
                            const std::vector<std::uint16_t>& qubitChoice,
                            qc::QuantumComputation&           circuit,
                            SynthesisResults* r, CliffordSynthesizer* opt, const SynthesisConfiguration& configuration);

    static void assertTableau(const SynthesisData& data, const Tableau& tableau, std::uint32_t position);

    static void makeSingleGateConstraints(
            const SynthesisData& data);
    static void makeMultipleGateConstraints(
            const SynthesisData& data);
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
        CY,
        CZ,
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
            case CY:
                return "CY";
            case CZ:
                return "CZ";
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
            case CY:
                return 8;
            case CZ:
                return 9;
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
            case CY:
                return qc::OpType::Y;
            case CZ:
                return qc::OpType::Z;
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
            GATES::CX,
            GATES::CY,
            GATES::CZ};

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
            GATES::CX,
            GATES::CY,
            GATES::CZ};
};
#endif //QMAP_CLIFFORDSYNTHESIS_H
