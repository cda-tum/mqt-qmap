/*
 * This file is part of the MQT QMAP library which is released under the MIT
 * license. See file README.md or go to
 * https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
 */
#ifndef CS_RESULTS_HPP
#define CS_RESULTS_HPP

#include "Architecture.hpp"
#include "Definitions.hpp"
#include "LogicTerm/Logic.hpp"
#include "QuantumComputation.hpp"
#include "cliffordsynthesis/OptimizationStrategy.hpp"
#include "cliffordsynthesis/Tableau.hpp"
#include "cliffordsynthesis/TargetMetric.hpp"
#include "operations/Operation.hpp"
#include "utils.hpp"

#include <ostream>
#include <sstream>
#include <string>
namespace cs {
struct Results {
  int                  verbose            = 0;
  bool                 chooseBest         = false;
  OptimizationStrategy strategy           = OptimizationStrategy::UseMinimizer;
  TargetMetric         target             = TargetMetric::GATES;
  logicbase::Result    result             = logicbase::Result::NDEF;
  std::uint16_t        nqubits            = 0U;
  std::uint16_t        architectureQubits = 0U;
  std::size_t          initialTimesteps   = 0U;
  std::size_t          singleQubitGates   = 0U;
  std::size_t          twoQubitGates      = 0U;
  std::size_t          depth              = 0U;
  bool                 sat                = false;
  double               fidelity           = 0.0;
  std::string          architectureName;

  double totalSeconds = 0;
  double finalRunTime = 0;

  std::string          resultStringCircuit{};
  std::vector<Tableau> resultTableaus{};

  CouplingMap                      resultCM{};
  std::vector<double>              singleFidelity{};
  std::vector<std::vector<double>> doubleFidelity{};

  Results()          = default;
  virtual ~Results() = default;

  Results(const Results& other)            = default;
  Results(Results&& other)                 = default;
  Results& operator=(const Results& other) = default;
  Results& operator=(Results&& other)      = default;

  void dump(std::ostream& os) const { os << json(); }

  [[nodiscard]] virtual nlohmann::json json() const {
    nlohmann::json resultJSON{};
    resultJSON["verbosity"]          = verbose;
    resultJSON["choose_best"]        = chooseBest;
    resultJSON["result"]             = toString(result);
    resultJSON["strategy"]           = toString(strategy);
    resultJSON["target"]             = toString(target);
    resultJSON["qubits"]             = nqubits;
    resultJSON["initial_timesteps"]  = initialTimesteps;
    resultJSON["single_qubit_gates"] = singleQubitGates;
    resultJSON["two_qubit_gates"]    = twoQubitGates;
    resultJSON["depth"]              = depth;
    resultJSON["fidelity"]           = fidelity;
    resultJSON["sat"]                = sat;
    resultJSON["total_seconds"]      = totalSeconds;
    resultJSON["resultCircuit"]      = resultStringCircuit;
    std::stringstream strings;
    Architecture::printCouplingMap(resultCM, strings);
    resultJSON["CouplingMap"]    = strings.str();
    resultJSON["singleFidelity"] = nlohmann::json::array();
    for (const auto& f : singleFidelity) {
      resultJSON["singleFidelity"].push_back(f);
    }
    resultJSON["doubleFidelity"] = nlohmann::json::array();
    for (const auto& f : doubleFidelity) {
      resultJSON["doubleFidelity"].push_back(f);
    }
    resultJSON["resultTableaus"] = nlohmann::json::array();
    for (const auto& tableau : resultTableaus) {
      resultJSON["resultTableaus"].push_back(tableau.toString());
    }

    return resultJSON;
  }

  [[nodiscard]] std::string getStrRepr() const {
    std::stringstream ss;
    dump(ss);
    return ss.str();
  }

  void generateStringCircuit(qc::QuantumComputation& resultCircuit) {
    std::stringstream ss;
    resultCircuit.dumpOpenQASM(ss);
    resultStringCircuit = ss.str();
  }
};

} // namespace cs

#endif // CSIMULATOR_RESULTS_HPP
