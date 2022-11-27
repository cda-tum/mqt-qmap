//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "cliffordsynthesis/encoding/GateEncoder.hpp"

#include <cstddef>
#include <optional>

namespace cs::encoding {

class SingleGateEncoder : public GateEncoder {
public:
  using GateEncoder::GateEncoder;

protected:
  void assertConsistency() const override;
  void assertGateConstraints() override;
  void assertSingleQubitGateConstraints(std::size_t pos) override;
  [[nodiscard]] logicbase::LogicTerm
       createSingleQubitGateConstraint(std::size_t pos, std::size_t qubit,
                                       qc::OpType gate) override;
  void assertTwoQubitGateConstraints(std::size_t pos) override;
  [[nodiscard]] logicbase::LogicTerm
  createTwoQubitGateConstraint(std::size_t pos, std::size_t ctrl,
                               std::size_t trgt) override;

  [[nodiscard]] logicbase::LogicTerm
  createNoChange(std::size_t pos, std::size_t except,
                 std::optional<std::size_t> except2) const;
};

} // namespace cs::encoding
