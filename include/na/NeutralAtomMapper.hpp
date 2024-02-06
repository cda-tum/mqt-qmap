//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "QuantumComputation.hpp"
#include "na/Architecture.hpp"

namespace na {

class NeutralAtomMapper {
protected:
  qc::QuantumComputation initialQc{};
  qc::QuantumComputation mappedQc{};
  na::Architecture arch;
public:
  explicit NeutralAtomMapper(const Architecture& arch) : arch(arch) {}
  virtual ~NeutralAtomMapper() = default;

  virtual void map(const qc::QuantumComputation& qc) = 0;
  [[nodiscard]] const qc::QuantumComputation& getResult() const { return mappedQc; }
};
} // namespace na