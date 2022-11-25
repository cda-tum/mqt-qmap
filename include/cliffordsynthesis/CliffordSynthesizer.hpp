//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "QuantumComputation.hpp"
#include "cliffordsynthesis/Configuration.hpp"
#include "cliffordsynthesis/Results.hpp"
#include "cliffordsynthesis/Tableau.hpp"

#include <cstddef>
#include <optional>
#include <sstream>

namespace cs {
class CliffordSynthesizer {
public:
  CliffordSynthesizer() = default;
  CliffordSynthesizer(Tableau initialTableau, Tableau targetTableau)
      : initialTableau(std::move(initialTableau)),
        targetTableau(std::move(targetTableau)) {}
  explicit CliffordSynthesizer(Tableau targetTableau)
      : initialTableau(targetTableau.getQubitCount()),
        targetTableau(std::move(targetTableau)) {}
  CliffordSynthesizer(Tableau initialTableau, const qc::QuantumComputation& qc)
      : initialTableau(std::move(initialTableau)), targetTableau(qc),
        initialGates(qc.getNindividualOps()) {}
  explicit CliffordSynthesizer(const qc::QuantumComputation& qc)
      : initialTableau(qc.getNqubits()), targetTableau(qc) {}

  virtual ~CliffordSynthesizer() = default;

  void synthesize(const Configuration& config = {});

  [[nodiscard]] Results&                getResults() { return results; };
  [[nodiscard]] qc::QuantumComputation& getResultCircuit() {
    std::stringstream ss;
    ss << results.getResultCircuit();
    resultCircuit = std::make_unique<qc::QuantumComputation>();
    resultCircuit->import(ss, qc::OpenQASM);
    return *resultCircuit;
  };
  [[nodiscard]] Tableau& getResultTableau() {
    std::stringstream ss;
    ss << results.getResultTableau();
    resultTableau.fromString(ss.str());
    return resultTableau;
  }

protected:
  Tableau initialTableau{};
  Tableau targetTableau{};

  std::optional<std::size_t> initialGates{};

  std::size_t timestepLimit{};

  Configuration configuration{};

  Results                                 results{};
  std::shared_ptr<qc::QuantumComputation> resultCircuit{};
  Tableau                                 resultTableau{};
  std::size_t                             solverCalls{};

  void    runMaxSAT();
  void    runBinarySearch();
  Results mainOptimization();

  static void updateResults(const Configuration& config,
                            const Results& newResults, Results& currentResults);
};

} // namespace cs
