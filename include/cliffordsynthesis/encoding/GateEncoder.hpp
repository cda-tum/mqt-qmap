//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "Definitions.hpp"
#include "LogicBlock/LogicBlock.hpp"
#include "cliffordsynthesis/Results.hpp"
#include "cliffordsynthesis/encoding/TableauEncoder.hpp"
#include "operations/OpType.hpp"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

namespace cs::encoding {
class GateSet : public std::vector<qc::OpType> {
  static constexpr std::array<qc::OpType, 9> SINGLE_QUBIT_CLIFFORDS = {
      qc::OpType::I,   qc::OpType::H,  qc::OpType::X,
      qc::OpType::Y,   qc::OpType::Z,  qc::OpType::S,
      qc::OpType::Sdg, qc::OpType::SX, qc::OpType::SXdg};

public:
  using std::vector<qc::OpType>::vector;

  template <qc::OpType Gate> [[nodiscard]] bool containsGate() const {
    for (const auto& g : // NOLINT(readability-use-anyofallof)
         *this) {
      if (g == Gate) {
        return true;
      }
    }
    return false;
  }
  [[nodiscard]] bool containsX() const { return containsGate<qc::OpType::X>(); }
  [[nodiscard]] bool containsY() const { return containsGate<qc::OpType::Y>(); }
  [[nodiscard]] bool containsZ() const { return containsGate<qc::OpType::Z>(); }
  [[nodiscard]] bool containsH() const { return containsGate<qc::OpType::H>(); }
  [[nodiscard]] bool containsS() const { return containsGate<qc::OpType::S>(); }
  [[nodiscard]] bool containsSdg() const {
    return containsGate<qc::OpType::Sdg>();
  }
  [[nodiscard]] bool containsSX() const {
    return containsGate<qc::OpType::SX>();
  }
  [[nodiscard]] bool containsSXdg() const {
    return containsGate<qc::OpType::SXdg>();
  }
  [[nodiscard]] std::size_t gateToIndex(const qc::OpType type) const {
    for (std::size_t i = 0; i < this->size(); ++i) {
      if (this->at(i) == type) {
        return i;
      }
    }
    return 0;
  }

  [[nodiscard]] bool isValidGateSet() const {
    return std::all_of(this->begin(), this->end(), [](const auto& g) {
      return std::find(SINGLE_QUBIT_CLIFFORDS.begin(),
                       SINGLE_QUBIT_CLIFFORDS.end(),
                       g) != SINGLE_QUBIT_CLIFFORDS.end();
    });
  }
};
class GateEncoder {
public:
  GateEncoder(const std::size_t nQubits, const std::size_t tableauSize,
              const std::size_t                      timestepLimit,
              TableauEncoder::Variables*             tableauVars,
              std::shared_ptr<logicbase::LogicBlock> logicBlock,
              GateSet singleQGates = {qc::OpType::None, qc::OpType::S,
                                      qc::OpType::Sdg, qc::OpType::H,
                                      qc::OpType::X, qc::OpType::Y,
                                      qc::OpType::Z})
      : N(nQubits), S(tableauSize), T(timestepLimit), tvars(tableauVars),
        lb(std::move(logicBlock)), singleQubitGates(std::move(singleQGates)) {
    if (!singleQGates.isValidGateSet()) {
      throw qc::QFRException("Invalid gate set");
    }
  }
  virtual ~GateEncoder() = default;

  struct Variables {
    // variables for the single-qubit gates
    logicbase::LogicMatrix3D gS{};
    // variables for the two-qubit gates
    logicbase::LogicMatrix3D gC{};

    void
         collectSingleQubitGateVariables(std::size_t pos, std::size_t qubit,
                                         logicbase::LogicVector& variables) const;
    void collectTwoQubitGateVariables(std::size_t pos, std::size_t qubit,
                                      bool                    target,
                                      logicbase::LogicVector& variables) const;
  };

  // variable creation
  void createSingleQubitGateVariables();
  void createTwoQubitGateVariables();

  // encode the relation between the tableaus and the gates
  virtual void encodeGates() {
    assertConsistency();
    assertGateConstraints();
  }

  virtual void encodeSymmetryBreakingConstraints();

  // extracting the circuit
  void extractCircuitFromModel(Results& res, logicbase::Model& model);

  [[nodiscard]] auto* getVariables() { return &vars; }

protected:
  // number of qubits N
  std::size_t N{}; // NOLINT (readability-identifier-naming)
  // number of rows in the tableau S
  std::size_t S{}; // NOLINT (readability-identifier-naming)
  // timestep limit T
  std::size_t T{}; // NOLINT (readability-identifier-naming)

  // the gate variables
  Variables vars{};

  // the tableau variables
  TableauEncoder::Variables* tvars{};

  // the logic block to use
  std::shared_ptr<logicbase::LogicBlock> lb{};

  // the gates that are used
  GateSet singleQubitGates;

  using TransformationFamily =
      std::pair<logicbase::LogicTerm, std::vector<qc::OpType>>;
  using GateToTransformation =
      std::function<logicbase::LogicTerm(std::size_t, std::size_t, qc::OpType)>;

  void assertExactlyOne(const logicbase::LogicVector& variables) const;

  virtual void assertConsistency() const = 0;

  virtual void assertGateConstraints()                           = 0;
  virtual void assertSingleQubitGateConstraints(std::size_t pos) = 0;
  virtual void assertTwoQubitGateConstraints(std::size_t pos)    = 0;
  [[nodiscard]] std::vector<TransformationFamily>
       collectGateTransformations(std::size_t pos, std::size_t qubit,
                                  const GateToTransformation& gateToTransformation);
  void assertGatesImplyTransform(
      std::size_t pos, std::size_t qubit,
      const std::vector<TransformationFamily>& transformations);
  virtual void assertZConstraints(std::size_t pos, std::size_t qubit);
  virtual void assertXConstraints(std::size_t pos, std::size_t qubit);
  virtual void assertRConstraints(std::size_t pos, std::size_t qubit);
  [[nodiscard]] virtual logicbase::LogicTerm
  createTwoQubitGateConstraint(std::size_t pos, std::size_t ctrl,
                               std::size_t trgt) = 0;

  void extractSingleQubitGatesFromModel(std::size_t             pos,
                                        logicbase::Model&       model,
                                        qc::QuantumComputation& qc,
                                        std::size_t& nSingleQubitGates);
  void extractTwoQubitGatesFromModel(std::size_t pos, logicbase::Model& model,
                                     qc::QuantumComputation& qc,
                                     std::size_t&            nTwoQubitGates);

  virtual void
  assertSingleQubitGateSymmetryBreakingConstraints(std::size_t pos);
  virtual void assertTwoQubitGateSymmetryBreakingConstraints(std::size_t pos);

  virtual void assertSingleQubitGateOrderConstraints(std::size_t pos,
                                                     std::size_t qubit) = 0;
  virtual void assertSingleQubitGateCancellationConstraints(std::size_t pos,
                                                            std::size_t qubit);

  virtual void assertTwoQubitGateOrderConstraints(std::size_t pos,
                                                  std::size_t ctrl,
                                                  std::size_t trgt) = 0;
};
} // namespace cs::encoding
