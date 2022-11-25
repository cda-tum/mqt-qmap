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
  std::size_t lowerTimestepLimit{};

  Configuration configuration{};

  Results                                 results{};
  std::shared_ptr<qc::QuantumComputation> resultCircuit{};
  Tableau                                 resultTableau{};
  std::size_t                             solverCalls{};

  void    determineUpperBound();
  void    minimizeGatesFixedDepth();
  void    runMaxSAT();
  Results mainOptimization(bool useMaxSAT);

  template <typename T>
  void runBinarySearch(T& value, T lowerBound, T upperBound) {
    INFO() << "Running binary search in range [" << lowerBound << ", "
           << upperBound << ")";

    while (lowerBound != upperBound) {
      value = (lowerBound + upperBound) / 2;
      INFO() << "Trying value " << value << " in range [" << lowerBound << ", "
             << upperBound << ")";
      const auto r = mainOptimization(false);
      updateResults(configuration, r, results);
      if (r.sat()) {
        upperBound = value;
        INFO() << "Found solution. New upper bound is " << upperBound;
      } else {
        lowerBound = value + 1;
        INFO() << "No solution found. New lower bound is " << lowerBound;
      }
    }
    INFO() << "Found optimum: " << lowerBound;
  }

  static void updateResults(const Configuration& config,
                            const Results& newResults, Results& currentResults);
};

} // namespace cs
