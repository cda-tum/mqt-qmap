/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
*/

#ifndef QMAP_CLIFFORDSYNTHESIS_H
#define QMAP_CLIFFORDSYNTHESIS_H

#include "Architecture.hpp"
#include "CliffordOptimizationResult.hpp"
#include "Encodings/Encodings.hpp"
#include "LogicBlock/LogicBlock.hpp"
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

#ifdef Z3_FOUND
    #include "LogicBlock/Z3Logic.hpp"
#endif

#include "QuantumComputation.hpp"

class CliffordOptimizer {
public:
    CliffordOptimizer() {
    }
    void optimize();

    void setArchitecture(const Architecture& arch) {
        architecture = arch;
        if (nqubits == 0) {
            nqubits = architecture.getNqubits();
        }
        auto map = highestFidelityMap.emplace_back();
        architecture.getHighestFidelityCouplingMap(nqubits, map);
    }
    bool                   choose_best   = false;
    bool                   use_embedding = false;
    unsigned char          nqubits       = 0U;
    std::set<signed char>  used_qubits{};
    unsigned short         initial_timesteps = 0U;
    int                    verbose           = 0;
    int                    nthreads          = 1;
    OptimizingStrategy     strategy          = OptimizingStrategy::UseMinimizer;
    OptTarget              target            = OptTarget::GATES;
    OptMethod              method            = OptMethod::Z3;
    qc::QuantumComputation circuit;

    std::vector<CouplingMap> highestFidelityMap;

    Tableau initialTableau{};
    Tableau targetTableau{};
    Tableau modelTableau{};

    void dumpResult(const std::string& outputFilename) {
        if (optimal_results.resultCircuit.empty()) {
            std::cerr << "Circuit is empty." << std::endl;
            return;
        }

        size_t      dot       = outputFilename.find_last_of('.');
        std::string extension = outputFilename.substr(dot + 1U);
        std::transform(extension.begin(), extension.end(), extension.begin(), [](unsigned char c) { return ::tolower(c); });
        if (extension == "real") {
            dumpResult(outputFilename, qc::Real);
        } else if (extension == "qasm") {
            dumpResult(outputFilename, qc::OpenQASM);
        }
    }

    void dumpResult(const std::string& outputFilename, qc::Format format) {
        optimal_results.resultCircuit.dump(outputFilename, format);
    }

    void dumpResult(std::ostream& os, qc::Format format) {
        optimal_results.resultCircuit.dump(os, format);
    }

    CliffordOptResults optimal_results{};

protected:
    Architecture       architecture{};
    CliffordOptResults main_optimization(
            int                                                        timesteps,
            const std::set<std::pair<unsigned short, unsigned short>>& reducedCM,
            const std::vector<unsigned short>& qubitChoice, Tableau& initialTableau,
            Tableau& targetTableau);

    void make_depth_optimizer(
            int                                                        timesteps,
            const std::set<std::pair<unsigned short, unsigned short>>& reducedCM,
            const std::vector<unsigned short>& qubitChoice, std::unique_ptr<LogicBlock>& lb,
            const LogicMatrix& x, const LogicMatrix& z, const LogicVector& r,
            const LogicMatrix3D& g_s, const LogicMatrix3D& g_c) const;

    void make_gate_optimizer(
            int                                                        timesteps,
            const std::set<std::pair<unsigned short, unsigned short>>& reducedCM,
            const std::vector<unsigned short>& qubitChoice, std::unique_ptr<LogicBlock>& lb,
            const LogicMatrix& x, const LogicMatrix& z, const LogicVector& r,
            const LogicMatrix3D& g_s, const LogicMatrix3D& g_c) const;
    void make_fidelity_optimizer(
            int                                                        timesteps,
            const std::set<std::pair<unsigned short, unsigned short>>& reducedCM,
            const std::vector<unsigned short>& qubitChoice, std::unique_ptr<LogicBlock>& lb,
            const LogicMatrix& x, const LogicMatrix& z, const LogicVector& r,
            const LogicMatrix3D& g_s, const LogicMatrix3D& g_c) const;

    void        runMinimizer(int timesteps, const CouplingMap& reducedCM,
                             const std::vector<unsigned short>& qubitChoice);
    void        runStartLow(int timesteps, const CouplingMap& reducedCM,
                            const std::vector<unsigned short>& qubitChoice);
    void        runStartHigh(int timesteps, const CouplingMap& reducedCM,
                             const std::vector<unsigned short>& qubitChoice);
    void        runMinMax(int timesteps, const CouplingMap& reducedCM,
                          const std::vector<unsigned short>& qubitChoice);
    void        runSplitIter(const CouplingMap&                 reducedCM,
                             const std::vector<unsigned short>& qubitChoice);
    static void runSplinter(int i, unsigned int circ_split, unsigned int split,
                            const CouplingMap&                 reducedCM,
                            const std::vector<unsigned short>& qubitChoice,
                            qc::QuantumComputation&            circuit,
                            CliffordOptResults* r, CliffordOptimizer* opt);
    void        updateResults(CliffordOptResults& r);

    static void assertTableau(const Tableau& tableau, std::unique_ptr<LogicBlock>& lb,
                              const LogicMatrix& x, const LogicMatrix& z,
                              const LogicVector& r, int nqubits, int position);

    static void makeSingleGateConstraints(
            std::unique_ptr<LogicBlock>& lb, const LogicMatrix& x, const LogicMatrix& z,
            const LogicVector& r, int nqubits, int timesteps,
            const std::set<std::pair<unsigned short, unsigned short>>& reducedCM,
            const std::vector<unsigned short>& qubitChoice, const LogicMatrix3D& g_s,
            const LogicMatrix3D& g_c);
    static void makeMultipleGateConstraints(
            std::unique_ptr<LogicBlock>& lb, const LogicMatrix& x, const LogicMatrix& z,
            const LogicVector& r, int nqubits, int timesteps,
            const std::set<std::pair<unsigned short, unsigned short>>& reducedCM,
            const std::vector<unsigned short>& qubitChoice, const LogicMatrix3D& g_s,
            const LogicMatrix3D& g_c);
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
        CX
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

    static constexpr GATES singleQubit[] = {
            GATES::NOP,
            GATES::H,
            GATES::S,
            GATES::X,
            GATES::Y,
            GATES::Z,
            GATES::Sdag};

    static constexpr GATES twoQubit[] = {
            GATES::CX};

    static constexpr GATES singleQubitWithoutNOP[] = {
            GATES::H,
            GATES::S,
            GATES::X,
            GATES::Y,
            GATES::Z,
            GATES::Sdag};

    static constexpr GATES allGates[] = {
            GATES::H,
            GATES::S,
            GATES::X,
            GATES::Y,
            GATES::Z,
            GATES::Sdag,
            GATES::CX};
};
#endif //QMAP_CLIFFORDSYNTHESIS_H
