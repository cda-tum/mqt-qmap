//
// Created by Velsh Aleksei on 16.06.23.
//

#pragma once

#include "cliffordsynthesis/encoding/GateEncoder.hpp"

#include <cstddef>
#include <optional>

namespace cs::encoding {

class TwoQubitEncoder : public GateEncoder {
public:
  using GateEncoder::GateEncoder;

protected:
  // define variables
  logicbase::LogicTerm   rChanges{};
  logicbase::LogicMatrix xorHelpers{};
  logicbase::LogicMatrix   gP{};

  static constexpr std::array<qc::OpType, 4> PAULI_GATES = {
      qc::OpType::None, qc::OpType::X, qc::OpType::Y, qc::OpType::Z};

  // variable creation
  virtual void createSingleQubitGateVariables() override;
  virtual void createTwoQubitGateVariables() override;
  void createPauliGateVariables();

  // create constrains
  [[nodiscard]] logicbase::LogicTerm
  createTwoQubitGateConstraint(std::size_t pos, std::size_t ctrl,
                               std::size_t trgt) override;

  // assert constrains
  void assertConsistency() const override;
  void assertGateConstraints() override;
  void assertRConstraints(std::size_t pos, std::size_t qubit) override;
  void assertSingleQubitGateConstraints(std::size_t pos) override;
  void assertPauliGateConstraints(std::size_t pos);
  void assertTwoQubitGateConstraints(std::size_t pos) override;
  void assertSingleQubitGateOrderConstraints(std::size_t pos,
                                             std::size_t qubit) override;
  void assertTwoQubitGateOrderConstraints(std::size_t pos, std::size_t ctrl,
                                          std::size_t trgt) override;

  // collect TQG variables
  void
  collectPauliGateVariables(std::size_t qubit,
                                 logicbase::LogicVector& variables) const;

  // extracting the circuit
  void extractCircuitFromModel(Results& res, logicbase::Model& model) override;

  // helpers
  void splitXorR(const logicbase::LogicTerm& changes, std::size_t pos);
  virtual void assertGatesImplyTransform(
      std::size_t pos, std::size_t qubit,
      const std::vector<TransformationFamily>& transformations) override;
  virtual void encodeSymmetryBreakingConstraints() override;
};

} // namespace cs::encoding
