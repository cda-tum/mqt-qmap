#pragma once

#include "Architecture.hpp"
#include "Definitions.hpp"
#include "QuantumComputation.hpp"
#include "Tableau.hpp"
#include "operations/Operation.hpp"
#include "utils.hpp"

#include <ostream>
#include <sstream>
#include <string>

enum class OptimizingStrategy {
    StartLow,
    StartHigh,
    UseMinimizer,
    MinMax,
    SplitIter
};
enum class OptimizationResult { SAT,
                                UNSAT,
                                UNDEF };
enum class OptimizationTarget { GATES,
                                GATES_ONLY_CNOT,
                                DEPTH,
                                FIDELITY };
enum class OptimizationMethod { Z3,
                                MATHSAT,
                                SMTLibV2,
                                DIMACS };

inline std::string toString(OptimizationMethod method) {
    switch (method) {
        case OptimizationMethod::Z3:
            return "Z3";
        case OptimizationMethod::MATHSAT:
            return "MATHSAT";
        case OptimizationMethod::SMTLibV2:
            return "SMTLibV2";
        case OptimizationMethod::DIMACS:
            return "DIMACS";
    }
    return "Error";
}
inline OptimizationMethod optMethodFromString(const std::string& method) {
    if (method == "Z3")
        return OptimizationMethod::Z3;
    if (method == "MATHSAT")
        return OptimizationMethod::MATHSAT;
    if (method == "SMTLibV2")
        return OptimizationMethod::SMTLibV2;
    if (method == "DIMACS")
        return OptimizationMethod::DIMACS;
    return OptimizationMethod::Z3;
}
inline std::string toString(OptimizationTarget target) {
    switch (target) {
        case OptimizationTarget::GATES:
            return "gates";
        case OptimizationTarget::GATES_ONLY_CNOT:
            return "gates_only_cnot";
        case OptimizationTarget::DEPTH:
            return "depth";
        case OptimizationTarget::FIDELITY:
            return "fidelity";
    }
    return "Error";
}

inline OptimizationTarget optTargetFromString(const std::string& target) {
    if (target == "gates")
        return OptimizationTarget::GATES;
    if (target == "gates_only_cnot")
        return OptimizationTarget::GATES_ONLY_CNOT;
    if (target == "depth")
        return OptimizationTarget::DEPTH;
    if (target == "fidelity")
        return OptimizationTarget::FIDELITY;
    return OptimizationTarget::GATES;
}

inline std::string toString(OptimizingStrategy strategy) {
    switch (strategy) {
        case OptimizingStrategy::MinMax:
            return "minmax";
        case OptimizingStrategy::StartHigh:
            return "start_high";
        case OptimizingStrategy::StartLow:
            return "start_low";
        case OptimizingStrategy::UseMinimizer:
            return "useminimizer";
        case OptimizingStrategy::SplitIter:
            return "split_iterative";
    }
    return "Error";
}

inline OptimizingStrategy optStrategyFromString(const std::string& strategy) {
    if (strategy == "minmax")
        return OptimizingStrategy::MinMax;
    if (strategy == "start_high")
        return OptimizingStrategy::StartHigh;
    if (strategy == "start_low")
        return OptimizingStrategy::StartLow;
    if (strategy == "useminimizer")
        return OptimizingStrategy::UseMinimizer;
    if (strategy == "split_iterative")
        return OptimizingStrategy::SplitIter;
    return OptimizingStrategy::MinMax;
}

class CliffordOptResults {
public:
    int                verbose           = 0;
    bool               choose_best       = false;
    OptimizingStrategy strategy          = OptimizingStrategy::UseMinimizer;
    OptimizationTarget target            = OptimizationTarget::GATES;
    OptimizationMethod method            = OptimizationMethod::Z3;
    OptimizationResult result            = OptimizationResult::UNDEF;
    unsigned char      nqubits           = 0;
    int                initial_timesteps = 0;
    int                gate_count        = 0;
    int                depth             = 0;
    bool               sat               = false;
    double             fidelity          = 0.0;

    double total_seconds  = 0;
    double final_run_time = 0;

    qc::QuantumComputation resultCircuit{};
    std::string            resultStringCircuit{};
    std::vector<Tableau>   resultTableaus{};

    CouplingMap                      resultCM{};
    std::vector<double>              singleFidelity{};
    std::vector<std::vector<double>> doubleFidelity{};

    CliffordOptResults()          = default;
    virtual ~CliffordOptResults() = default;

    CliffordOptResults(CliffordOptResults& other) {
        verbose           = other.verbose;
        choose_best       = other.choose_best;
        strategy          = other.strategy;
        target            = other.target;
        method            = other.method;
        nqubits           = other.nqubits;
        initial_timesteps = other.initial_timesteps;
        gate_count        = other.gate_count;
        depth             = other.depth;
        sat               = other.sat;
        total_seconds     = other.total_seconds;
        final_run_time    = other.final_run_time;
        resultCircuit     = other.resultCircuit.clone();
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

    CliffordOptResults& operator=(CliffordOptResults other) {
        verbose           = other.verbose;
        choose_best       = other.choose_best;
        strategy          = other.strategy;
        target            = other.target;
        method            = other.method;
        nqubits           = other.nqubits;
        initial_timesteps = other.initial_timesteps;
        gate_count        = other.gate_count;
        depth             = other.depth;
        sat               = other.sat;
        total_seconds     = other.total_seconds;
        final_run_time    = other.final_run_time;
        resultCircuit     = other.resultCircuit.clone();
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
        os << R"("choose_best":")" << choose_best << "\"," << std::endl;
        os << R"("strategy":")" << toString(strategy) << "\"," << std::endl;
        os << R"("target":")" << toString(target) << "\"," << std::endl;
        os << R"("method":")" << toString(method) << "\"," << std::endl;
        os << R"("qubits":")" << std::to_string(nqubits) << "\"," << std::endl;
        os << R"("initial_timesteps":")" << std::to_string(initial_timesteps)
           << "\"," << std::endl;
        os << R"("gate_count":")" << std::to_string(gate_count) << "\","
           << std::endl;
        os << R"("depth":")" << std::to_string(depth) << "\"," << std::endl;
        os << R"("fidelity":")" << std::to_string(fidelity) << "\"," << std::endl;
        os << R"("sat":")" << (sat ? "SAT" : "UNSAT") << "\"," << std::endl;
        os << R"("total_seconds":")" << std::to_string(total_seconds) << "\","
           << std::endl;
        os << R"("resultCircuit":")";
        std::stringstream ss;
        resultCircuit.dump(ss, qc::Format::OpenQASM);
        os << escapeChars(ss.str(), "\"") << "\"," << std::endl;
        os << "\"resultTableaus\":[" << std::endl;
        bool skipfirst = true;
        for (const auto& tableau: resultTableaus) {
            if (!skipfirst)
                os << "," << std::endl;
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
            if (!skipfirst)
                os << ",";
            os << "\"" << std::to_string(f) << "\"";
            skipfirst = false;
        }
        os << "]," << std::endl;
        os << "\"doubleFidelity\":[";
        skipfirst = true;
        for (const auto& f: doubleFidelity) {
            if (!skipfirst)
                os << ",";
            os << "[";
            bool skipfirst2 = true;
            for (const auto& f2: f) {
                if (!skipfirst2)
                    os << ",";
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
        resultJSON["choose_best"]       = choose_best;
        resultJSON["strategy"]          = toString(strategy);
        resultJSON["target"]            = toString(target);
        resultJSON["method"]            = toString(method);
        resultJSON["qubits"]            = nqubits;
        resultJSON["initial_timesteps"] = initial_timesteps;
        resultJSON["gate_count"]        = gate_count;
        resultJSON["depth"]             = depth;
        resultJSON["fidelity"]          = fidelity;
        resultJSON["sat"]               = sat;
        resultJSON["total_seconds"]     = total_seconds;
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
