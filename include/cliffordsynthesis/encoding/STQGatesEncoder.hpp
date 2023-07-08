//
// Created by Velsh Aleksei on 16.06.23.
//

#pragma once

#include "cliffordsynthesis/encoding/GateEncoder.hpp"
#include <cstddef>
#include <optional>

namespace cs::encoding {

class STQGatesEncoder : public GateEncoder {
public:
  using GateEncoder::GateEncoder;

protected:
  logicbase::LogicTerm rChanges{};

  static constexpr std::array<qc::OpType, 4>
      SINGLE_QUBIT_GATES_FOR_STQ_ENCODING = {
      qc::OpType::None, qc::OpType::H, qc::OpType::S, qc::OpType::SX};

  void assertConsistency() const override;
  void assertGateConstraints() override;
  void assertRConstraints(std::size_t pos, std::size_t qubit) override;
  void assertSingleQubitGateConstraints(std::size_t pos) override;
  void assertTwoQubitGateConstraints(std::size_t pos) override;
  [[nodiscard]] logicbase::LogicTerm
  createTwoQubitGateConstraint(std::size_t pos, std::size_t ctrl,
                               std::size_t trgt) override;

  // variable creation
  void createSingleQubitGateVariables() override;
  void createTwoQubitGateVariables() override;

  void assertSingleQubitGateOrderConstraints(std::size_t pos,
                                             std::size_t qubit) override;
  void assertTwoQubitGateOrderConstraints(std::size_t pos, std::size_t ctrl,
                                          std::size_t trgt) override;

  void collectTwoQubitGateVariables(std::size_t pos, std::size_t qubit,
                                    bool                    target,
                                    logicbase::LogicVector& variables) const;

  //TODO: remove later
  [[nodiscard]] std::vector<TransformationFamily>
  collectGateTransformations(std::size_t pos, std::size_t qubit,
                             const GateToTransformation& gateToTransformation) override;

  void extractSingleQubitGatesFromModel(std::size_t             pos,
                                        logicbase::Model&       model,
                                        qc::QuantumComputation& qc,
                                        std::size_t& nSingleQubitGates) override;

};

} // namespace cs::encoding
