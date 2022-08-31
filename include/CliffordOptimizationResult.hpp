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
enum class OptResult { SAT,
                       UNSAT,
                       UNDEF };
enum class OptTarget { GATES,
                       GATES_ONLY_CNOT,
                       DEPTH,
                       FIDELITY };
enum class OptMethod { Z3,
                       MATHSAT,
                       SMTLibV2,
                       DIMACS };

inline std::string toString(OptMethod method) {
    switch (method) {
        case OptMethod::Z3:
            return "Z3";
        case OptMethod::MATHSAT:
            return "MATHSAT";
        case OptMethod::SMTLibV2:
            return "SMTLibV2";
        case OptMethod::DIMACS:
            return "DIMACS";
    }
    return "Error";
}
inline OptMethod optMethodFromString(const std::string& method) {
    if (method == "Z3")
        return OptMethod::Z3;
    if (method == "MATHSAT")
        return OptMethod::MATHSAT;
    if (method == "SMTLibV2")
        return OptMethod::SMTLibV2;
    if (method == "DIMACS")
        return OptMethod::DIMACS;
    return OptMethod::Z3;
}
inline std::string toString(OptTarget target) {
    switch (target) {
        case OptTarget::GATES:
            return "Gates";
        case OptTarget::GATES_ONLY_CNOT:
            return "Gates (only CNOT)";
        case OptTarget::DEPTH:
            return "Depth";
        case OptTarget::FIDELITY:
            return "Fidelity";
    }
    return "Error";
}

inline OptTarget optTargetFromString(const std::string& target) {
    if (target == "Gates")
        return OptTarget::GATES;
    if (target == "Gates (only CNOT)")
        return OptTarget::GATES_ONLY_CNOT;
    if (target == "Depth")
        return OptTarget::DEPTH;
    if (target == "Fidelity")
        return OptTarget::FIDELITY;
    return OptTarget::GATES;
}

inline std::string toString(OptimizingStrategy strategy) {
    switch (strategy) {
        case OptimizingStrategy::MinMax:
            return "MinMax";
        case OptimizingStrategy::StartHigh:
            return "Start High";
        case OptimizingStrategy::StartLow:
            return "Start Low";
        case OptimizingStrategy::UseMinimizer:
            return "Minimizer";
        case OptimizingStrategy::SplitIter:
            return "Split Iterative";
    }
    return "Error";
}

inline OptimizingStrategy optStrategyFromString(const std::string& strategy) {
    if (strategy == "MinMax")
        return OptimizingStrategy::MinMax;
    if (strategy == "Start High")
        return OptimizingStrategy::StartHigh;
    if (strategy == "Start Low")
        return OptimizingStrategy::StartLow;
    if (strategy == "Minimizer")
        return OptimizingStrategy::UseMinimizer;
    if (strategy == "Split Iterative")
        return OptimizingStrategy::SplitIter;
    return OptimizingStrategy::MinMax;
}

class CliffordOptResults {
public:
    int                verbose           = 0;
    bool               choose_best       = false;
    OptimizingStrategy strategy          = OptimizingStrategy::UseMinimizer;
    OptTarget          target            = OptTarget::GATES;
    OptMethod          method            = OptMethod::Z3;
    OptResult          result            = OptResult::UNDEF;
    unsigned char      nqubits           = 0;
    int                initial_timesteps = 0;
    int                gate_count        = 0;
    int                depth             = 0;
    bool               sat               = false;
    double             fidelity          = 0.0;

    double total_seconds  = 0;
    double final_run_time = 0;

    qc::QuantumComputation resultCircuit{};
    std::vector<Tableau>   resultTableaus{};

    CouplingMap                      resultCM{};
    std::vector<double>              singleFidelity{};
    std::vector<std::vector<double>> doubleFidelity{};

    CliffordOptResults(){};

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
        resultTableaus    = other.resultTableaus;
        resultCM          = other.resultCM;
        singleFidelity    = other.singleFidelity;
        doubleFidelity    = other.doubleFidelity;
        fidelity          = other.fidelity;
        result            = other.result;
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
        resultTableaus    = other.resultTableaus;
        resultCM          = other.resultCM;
        singleFidelity    = other.singleFidelity;
        doubleFidelity    = other.doubleFidelity;
        fidelity          = other.fidelity;
        result            = other.result;
        return *this;
    };

    CliffordOptResults operator+(CliffordOptResults& other) {
        CliffordOptResults result = *this;
        result.verbose += other.verbose;
        result.choose_best = other.choose_best;
        result.strategy    = other.strategy;
        result.target      = other.target;
        result.method      = other.method;
        result.nqubits     = other.nqubits;
        result.initial_timesteps += other.initial_timesteps;
        result.gate_count += other.gate_count;
        result.depth += other.depth;
        result.sat = other.sat;
        result.total_seconds += other.total_seconds;
        result.final_run_time += other.final_run_time;
        result.resultCircuit  = resultCircuit.clone();
        result.resultTableaus = other.resultTableaus;
        result.resultCM       = other.resultCM;
        result.singleFidelity = other.singleFidelity;
        result.doubleFidelity = other.doubleFidelity;
        result.fidelity += other.fidelity;
        result.result = other.result;
        return result;
    };

    void dump(std::ostream& os) {
        os << "{\"CliffordOptimizationResult\":{" << std::endl;
        os << "\"verbose\":\"" << verbose << "\"," << std::endl;
        os << "\"choose_best\":\"" << choose_best << "\"," << std::endl;
        os << "\"strategy\":\"" << toString(strategy) << "\"," << std::endl;
        os << "\"target\":\"" << toString(target) << "\"," << std::endl;
        os << "\"method\":\"" << toString(method) << "\"," << std::endl;
        os << "\"nqubits\":\"" << std::to_string(nqubits) << "\"," << std::endl;
        os << "\"initial_timesteps\":\"" << std::to_string(initial_timesteps)
           << "\"," << std::endl;
        os << "\"gate_count\":\"" << std::to_string(gate_count) << "\","
           << std::endl;
        os << "\"depth\":\"" << std::to_string(depth) << "\"," << std::endl;
        os << "\"fidelity\":\"" << std::to_string(fidelity) << "\"," << std::endl;
        os << "\"sat\":\"" << (sat ? "SAT" : "UNSAT") << "\"," << std::endl;
        os << "\"total_seconds\":\"" << std::to_string(total_seconds) << "\","
           << std::endl;
        os << "\"resultCircuit\":\"";
        std::stringstream ss;
        resultCircuit.dump(ss, qc::Format::OpenQASM);
        os << escapeChars(ss.str(), "\"") << "\"," << std::endl;
        os << "\"resultTableaus\":[" << std::endl;
        bool skipfirst = true;
        for (const auto& tableau: resultTableaus) {
            if (!skipfirst)
                os << "," << std::endl;
            os << "\"";
            os << escapeChars(tableau.getRepresentation(), "\"");
            os << "\"";
            skipfirst = false;
        }
        os << "]," << std::endl;
        std::stringstream strings;
        Architecture::printCouplingMap(resultCM, strings);
        os << "\"CouplingMap\":\"" << strings.str() << "\","
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
};
