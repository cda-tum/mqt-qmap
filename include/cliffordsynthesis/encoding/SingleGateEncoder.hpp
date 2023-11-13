//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "cliffordsynthesis/encoding/GateEncoder.hpp"

#include <cstddef>
#include <functional>
#include <optional>

namespace cs::encoding {

class SingleGateEncoder : public GateEncoder {
public:
  using GateEncoder::GateEncoder;

protected:
  void assertConsistency() const override;
  void assertGateConstraints() override;
  void assertNoGateNoChangeConstraint(std::size_t pos);
  [[nodiscard]] logicbase::LogicTerm
  createTwoQubitGateConstraint(std::size_t pos, std::size_t ctrl,
                               std::size_t trgt) override;
  [[nodiscard]] logicbase::LogicTerm
  createTwoQubitRConstraint(std::size_t pos, std::size_t ctrl,
                            std::size_t trgt) override;

  [[nodiscard]] logicbase::LogicTerm createNoChangeOnQubit(std::size_t pos,
                                                           std::size_t q);

  [[nodiscard]] logicbase::LogicTerm createNoGateOnQubit(std::size_t pos,
                                                         std::size_t q);

  void assertSingleQubitGateOrderConstraints(std::size_t pos,
                                             std::size_t qubit) override;
  void assertTwoQubitGateOrderConstraints(std::size_t pos, std::size_t ctrl,
                                          std::size_t trgt) override;
};

} // namespace cs::encoding
