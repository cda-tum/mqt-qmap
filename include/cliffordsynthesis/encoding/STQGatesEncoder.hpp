//
// Created by Velsh Aleksei on 17.06.23.
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

  void assertConsistency() const override;
  void assertGateConstraints() override;
  void assertRConstraints(std::size_t pos, std::size_t qubit) override;
  void assertSingleQubitGateConstraints(std::size_t pos) override;
  void assertTwoQubitGateConstraints(std::size_t pos) override;
  [[nodiscard]] logicbase::LogicTerm
  createTwoQubitGateConstraint(std::size_t pos, std::size_t ctrl,
                               std::size_t trgt) override;

  void assertSingleQubitGateOrderConstraints(std::size_t pos,
                                             std::size_t qubit) override;
  void assertTwoQubitGateOrderConstraints(std::size_t pos, std::size_t ctrl,
                                          std::size_t trgt) override;
};

} // namespace cs::encoding
