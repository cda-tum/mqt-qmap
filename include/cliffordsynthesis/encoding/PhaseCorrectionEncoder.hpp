//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include <utility>

#include "LogicTerm/Logic.hpp"
#include "cliffordsynthesis/Tableau.hpp"
namespace cs {
  class PhaseCorrectionEncoder {
  public:
    PhaseCorrectionEncoder(std::size_t N, Tableau  uncorrected, Tableau  target)
      : uncorrected(std::move(uncorrected)), target(std::move(target)) {}

    [[nodiscard]] Tableau getCorrected() const;

    void createVariables();
  protected:
    Tableau uncorrected{};
    Tableau target{};
    // number of qubits N
    std::size_t N{}; // NOLINT (readability-identifier-naming)
    logicbase::LogicMatrix paulis{};
  }; 
} // namespace cs;
