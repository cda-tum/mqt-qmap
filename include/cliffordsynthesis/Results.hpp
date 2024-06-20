//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "CircuitOptimizer.hpp"
#include "cliffordsynthesis/Tableau.hpp"
#include "logicblocks/Logic.hpp"

#include <nlohmann/json.hpp>
#include <sstream>
#include <string>

namespace cs {
class Results {
public:
  Results() = default;
  Results(qc::QuantumComputation& qc, const Tableau& tableau) {
    // SWAP gates are not natively supported in the encoding, so we need to
    // decompose them into sequences of three CNOTs.
    qc::CircuitOptimizer::decomposeSWAP(qc, false);

    setResultCircuit(qc);
    setResultTableau(tableau);
    setDepth(qc.getDepth());
    setSingleQubitGates(qc.getNsingleQubitOps());
    setTwoQubitGates(qc.getNindividualOps() - singleQubitGates);
    setSolverResult(logicbase::Result::SAT);
  }

  virtual ~Results() = default;

  [[nodiscard]] std::size_t getGates() const {
    return getSingleQubitGates() + getTwoQubitGates();
  }
  [[nodiscard]] std::size_t getTwoQubitGates() const { return twoQubitGates; }
  [[nodiscard]] std::size_t getSingleQubitGates() const {
    return singleQubitGates;
  }
  [[nodiscard]] std::size_t getDepth() const { return depth; }
  [[nodiscard]] double getRuntime() const { return runtime; }
  [[nodiscard]] logicbase::Result getSolverResult() const {
    return solverResult;
  }
  [[nodiscard]] std::size_t getSolverCalls() const { return solverCalls; }

  [[nodiscard]] std::string getResultCircuit() const { return resultCircuit; }
  [[nodiscard]] std::string getResultTableau() const { return resultTableau; }

  void setSingleQubitGates(const std::size_t g) { singleQubitGates = g; }
  void setTwoQubitGates(const std::size_t g) { twoQubitGates = g; }
  void setDepth(const std::size_t d) { depth = d; }
  void setRuntime(const double t) { runtime = t; }
  void setSolverResult(const logicbase::Result r) { solverResult = r; }
  void setSolverCalls(const std::size_t c) { solverCalls = c; }

  void setResultCircuit(qc::QuantumComputation& qc) {
    std::stringstream ss;
    qc.dumpOpenQASM3(ss);
    resultCircuit = ss.str();
  }
  void setResultTableau(const Tableau& tableau) {
    std::stringstream ss;
    ss << tableau;
    resultTableau = ss.str();
  }

  [[nodiscard]] bool sat() const {
    return getSolverResult() == logicbase::Result::SAT;
  }
  [[nodiscard]] bool unsat() const {
    return getSolverResult() == logicbase::Result::UNSAT;
  }

  [[nodiscard]] virtual nlohmann::json json() const {
    nlohmann::json resultJSON{};
    resultJSON["solver_result"] = toString(solverResult);
    resultJSON["single_qubit_gates"] = singleQubitGates;
    resultJSON["two_qubit_gates"] = twoQubitGates;
    resultJSON["depth"] = depth;
    resultJSON["runtime"] = runtime;
    resultJSON["solver_calls"] = solverCalls;

    return resultJSON;
  }

  friend std::ostream& operator<<(std::ostream& os, const Results& config) {
    os << config.json().dump(2);
    return os;
  }

protected:
  logicbase::Result solverResult = logicbase::Result::NDEF;
  std::size_t singleQubitGates = std::numeric_limits<std::size_t>::max();
  std::size_t twoQubitGates = std::numeric_limits<std::size_t>::max();
  std::size_t depth = std::numeric_limits<std::size_t>::max();
  double runtime = 0.0;
  std::size_t solverCalls = 0U;

  std::string resultTableau{};
  std::string resultCircuit{};
};

} // namespace cs
