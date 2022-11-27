//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "QuantumComputation.hpp"
#include "cliffordsynthesis/Configuration.hpp"
#include "cliffordsynthesis/Results.hpp"
#include "cliffordsynthesis/Tableau.hpp"
#include "cliffordsynthesis/encoding/SATEncoder.hpp"

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
      : initialTableau(std::move(initialTableau)), targetTableau(qc) {}
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

  Configuration configuration{};

  Results                                 results{};
  std::shared_ptr<qc::QuantumComputation> resultCircuit{};
  Tableau                                 resultTableau{};
  std::size_t                             solverCalls{};

  static bool requiresMultiGateEncoding(const TargetMetric metric) {
    return metric == TargetMetric::DEPTH;
  }

  static void
  determineInitialTimestepLimit(encoding::SATEncoder::Configuration& config);
  std::pair<std::size_t, std::size_t>
          determineUpperBound(encoding::SATEncoder::Configuration config);
  void    minimizeGatesFixedDepth(encoding::SATEncoder::Configuration config);
  void    runMaxSAT(const encoding::SATEncoder::Configuration& config);
  Results callSolver(const encoding::SATEncoder::Configuration& config);

  template <typename T>
  void runBinarySearch(T& value, T lowerBound, T upperBound,
                       const encoding::SATEncoder::Configuration& config) {
    INFO() << "Running binary search in range [" << lowerBound << ", "
           << upperBound << ")";

    while (lowerBound != upperBound) {
      value = (lowerBound + upperBound) / 2;
      INFO() << "Trying value " << value << " in range [" << lowerBound << ", "
             << upperBound << ")";
      const auto r = callSolver(config);
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
