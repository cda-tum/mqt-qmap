#pragma once

#include "Architecture.hpp"
#include "Definitions.hpp"
#include "QuantumComputation.hpp"
#include "Tableau.hpp"
#include "configuration/SynthesisMethod.hpp"
#include "configuration/SynthesisResult.hpp"
#include "configuration/SynthesisStrategy.hpp"
#include "configuration/SynthesisTarget.hpp"
#include "operations/Operation.hpp"
#include "utils.hpp"

#include <ostream>
#include <sstream>
#include <string>

class CliffordOptimizationResults {
public:
    int               verbose          = 0;
    bool              chooseBest       = false;
    SynthesisStrategy strategy         = SynthesisStrategy::UseMinimizer;
    SynthesisTarget   target           = SynthesisTarget::GATES;
    SynthesisMethod   method           = SynthesisMethod::Z3;
    SynthesisResult   result           = SynthesisResult::UNDEF;
    unsigned char     nqubits          = 0;
    int               initialTimesteps = 0;
    int               gateCount        = 0;
    int               depth            = 0;
    bool              sat              = false;
    double            fidelity         = 0.0;

    double totalSeconds = 0;
    double finalRunTime = 0;

    qc::QuantumComputation resultCircuit{};
    std::string            resultStringCircuit{};
    std::vector<Tableau>   resultTableaus{};

    CouplingMap                      resultCM{};
    std::vector<double>              singleFidelity{};
    std::vector<std::vector<double>> doubleFidelity{};

    CliffordOptimizationResults()          = default;
    virtual ~CliffordOptimizationResults() = default;

    CliffordOptimizationResults(CliffordOptimizationResults& other):
        verbose(other.verbose), chooseBest(other.chooseBest), strategy(other.strategy), target(other.target), method(other.method), nqubits(other.nqubits), initialTimesteps(other.initialTimesteps), gateCount(other.gateCount), depth(other.depth), sat(other.sat), totalSeconds(other.totalSeconds), finalRunTime(other.finalRunTime) {
        resultCircuit = other.resultCircuit.clone();
        if (other.resultStringCircuit.empty()) {
            std::stringstream ss;
            other.resultCircuit.dumpOpenQASM(ss);
            resultStringCircuit = ss.str();
        } else {
            resultStringCircuit = other.resultStringCircuit;
        }
        resultTableaus = other.resultTableaus;
        resultCM       = other.resultCM;
        singleFidelity = other.singleFidelity;
        doubleFidelity = other.doubleFidelity;
        fidelity       = other.fidelity;
        result         = other.result;
    };

    CliffordOptimizationResults& operator=(CliffordOptimizationResults other) {
        verbose          = other.verbose;
        chooseBest       = other.chooseBest;
        strategy         = other.strategy;
        target           = other.target;
        method           = other.method;
        nqubits          = other.nqubits;
        initialTimesteps = other.initialTimesteps;
        gateCount        = other.gateCount;
        depth            = other.depth;
        sat              = other.sat;
        totalSeconds     = other.totalSeconds;
        finalRunTime     = other.finalRunTime;
        resultCircuit    = other.resultCircuit.clone();
        if (other.resultStringCircuit.empty()) {
            std::stringstream ss;
            other.resultCircuit.dumpOpenQASM(ss);
            resultStringCircuit = ss.str();
        } else {
            resultStringCircuit = other.resultStringCircuit;
        }
        resultTableaus = other.resultTableaus;
        resultCM       = other.resultCM;
        singleFidelity = other.singleFidelity;
        doubleFidelity = other.doubleFidelity;
        fidelity       = other.fidelity;
        result         = other.result;
        return *this;
    };

    void dump(std::ostream& os) {
        os << "{\"CliffordOptimizationResult\":{" << std::endl;
        os << R"("verbose":")" << verbose << "\"," << std::endl;
        os << R"("choose_best":")" << chooseBest << "\"," << std::endl;
        os << R"("result":")" << toString(result) << "\"," << std::endl;
        os << R"("strategy":")" << toString(strategy) << "\"," << std::endl;
        os << R"("target":")" << toString(target) << "\"," << std::endl;
        os << R"("method":")" << toString(method) << "\"," << std::endl;
        os << R"("qubits":")" << std::to_string(nqubits) << "\"," << std::endl;
        os << R"("initial_timesteps":")" << std::to_string(initialTimesteps)
           << "\"," << std::endl;
        os << R"("gate_count":")" << std::to_string(gateCount) << "\","
           << std::endl;
        os << R"("depth":")" << std::to_string(depth) << "\"," << std::endl;
        os << R"("fidelity":")" << std::to_string(fidelity) << "\"," << std::endl;
        os << R"("sat":")" << (sat ? "SAT" : "UNSAT") << "\"," << std::endl;
        os << R"("total_seconds":")" << std::to_string(totalSeconds) << "\","
           << std::endl;
        os << R"("resultCircuit":")";
        std::stringstream ss;
        resultCircuit.dump(ss, qc::Format::OpenQASM);
        os << escapeChars(ss.str(), "\"") << "\"," << std::endl;
        os << "\"resultTableaus\":[" << std::endl;
        bool skipfirst = true;
        for (const auto& tableau: resultTableaus) {
            if (!skipfirst) {
                os << "," << std::endl;
            }
            os << "\"";
            os << escapeChars(tableau.toString(), "\"");
            os << "\"";
            skipfirst = false;
        }
        os << "]," << std::endl;
        std::stringstream strings;
        Architecture::printCouplingMap(resultCM, strings);
        os << R"("CouplingMap":")" << strings.str() << "\","
           << std::endl;
        os << "\"singleFidelity\":[";
        skipfirst = true;
        for (const auto& f: singleFidelity) {
            if (!skipfirst) {
                os << ",";
            }
            os << "\"" << std::to_string(f) << "\"";
            skipfirst = false;
        }
        os << "]," << std::endl;
        os << "\"doubleFidelity\":[";
        skipfirst = true;
        for (const auto& f: doubleFidelity) {
            if (!skipfirst) {
                os << ",";
            }
            os << "[";
            bool skipfirst2 = true;
            for (const auto& f2: f) {
                if (!skipfirst2) {
                    os << ",";
                }
                os << "\"" << std::to_string(f2) << "\"";
                skipfirst2 = false;
            }
            os << "]";
            skipfirst = false;
        }
        os << "]" << std::endl;
        os << "}}" << std::endl;
    }

    [[nodiscard]] virtual nlohmann::json json() const {
        nlohmann::json resultJSON{};
        resultJSON["verbose"]           = verbose;
        resultJSON["choose_best"]       = chooseBest;
        resultJSON["result"]            = toString(result);
        resultJSON["strategy"]          = toString(strategy);
        resultJSON["target"]            = toString(target);
        resultJSON["method"]            = toString(method);
        resultJSON["qubits"]            = nqubits;
        resultJSON["initial_timesteps"] = initialTimesteps;
        resultJSON["gate_count"]        = gateCount;
        resultJSON["depth"]             = depth;
        resultJSON["fidelity"]          = fidelity;
        resultJSON["sat"]               = sat;
        resultJSON["total_seconds"]     = totalSeconds;
        resultJSON["resultCircuit"]     = resultStringCircuit;

        return resultJSON;
    }

    std::string getStrRepr() {
        std::stringstream ss;
        dump(ss);
        return ss.str();
    }

    void generateStringCircuit() {
        std::stringstream ss;
        resultCircuit.dumpOpenQASM(ss);
        resultStringCircuit = ss.str();
    }
};
