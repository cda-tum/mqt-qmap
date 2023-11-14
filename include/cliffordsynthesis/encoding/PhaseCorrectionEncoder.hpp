//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include <utility>

#include "LogicTerm/Logic.hpp"
#include "LogicBlock/LogicBlock.hpp"
#include "LogicUtil/util_logicblock.hpp"
#include "cliffordsynthesis/Tableau.hpp"
#include "operations/OpType.hpp"
#include <memory>
#include <vector>

namespace cs::encoding {
  using namespace logicbase;
  class PhaseCorrectionEncoder {
  public:
    PhaseCorrectionEncoder(const std::size_t nQubits, const std::size_t tableauSize,
                           Tableau  uncorrectedTableau, Tableau targetTableau)
      : N{nQubits}, S{tableauSize}, uncorrected(std::move(uncorrectedTableau)), target(std::move(targetTableau)), paulis{N} {
      bool success = true;
      lb = logicutil::getZ3LogicBlock(success, true);
        if (!success) {
          FATAL() << "Could not initialize solver engine.";
  }
    }

    std::vector<qc::OpType> phaseCorrection();

  protected:
        // number of qubits N
    std::size_t N{}; // NOLINT (readability-identifier-naming)
    // number of rows in the tableau S
    std::size_t S{}; // NOLINT (readability-identifier-naming)
    // the logic block to use
    std::shared_ptr<LogicBlock> lb{}; 
    Tableau uncorrected{};
    Tableau target{};
    
    logicbase::LogicMatrix paulis{}; //order: I, X, Y, Z
    LogicVector        initialPhase{};
    LogicVector         targetPhase{};

    LogicMatrix xorHelpers{};

    void splitXorR(const LogicVector& changes);

    [[nodiscard]] LogicVector vectorFromBitset(const std::bitset<64>& bs) const;

    void encodePauliConstraints();

    void createVariables();

    std::vector<qc::OpType> extractResult();
  }; 
} // namespace cs;
