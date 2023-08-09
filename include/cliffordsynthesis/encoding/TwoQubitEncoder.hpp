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
  logicbase::LogicTerm   rChanges{};
  logicbase::LogicMatrix xorHelpers{};

  void assertConsistency() const override;
  void assertGateConstraints() override;
  void assertRConstraints(std::size_t pos, std::size_t qubit) override;
  void assertSingleQubitGateConstraints(std::size_t pos) override;
  void assertTwoQubitGateConstraints(std::size_t pos) override;
  [[nodiscard]] logicbase::LogicTerm
  createTwoQubitGateConstraint(std::size_t pos, std::size_t ctrl,
                               std::size_t trgt) override;
  [[nodiscard]] logicbase::LogicTerm
  createIdentityConstraintOnTQG(std::size_t pos, std::size_t ctrl);

  // assert constrains
  void assertSingleQubitGateOrderConstraints(std::size_t pos,
                                             std::size_t qubit) override;
  void assertTwoQubitGateOrderConstraints(std::size_t pos, std::size_t ctrl,
                                          std::size_t trgt) override;

  // collect TQG variables
  void collectTwoQubitGateVariables(std::size_t pos, std::size_t qubit,
                                    bool                    target,
                                    logicbase::LogicVector& variables) const;
  // extracting the circuit
  virtual void extractCircuitFromModel(Results&          res,
                                       logicbase::Model& model) override;

  void splitXorR(const logicbase::LogicTerm& changes, std::size_t pos);
};

} // namespace cs::encoding
