//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "Architecture.hpp"
#include "Configuration.hpp"
#include "QuantumComputation.hpp"

#include <cstdint>

namespace na {

class NeutralAtomMapper {
public:
  struct Statistics {
    std::size_t numInitialGates = 0;
    std::size_t numMappedGates  = 0;
    qc::fp      mappingTime     = 0.0; // [ms]
  };

protected:
  qc::QuantumComputation            initialQc{};
  qc::QuantumComputation            mappedQc{};
  na::Architecture                  arch;
  na::Configuration                 config;
  na::NeutralAtomMapper::Statistics stats{};

public:
  explicit NeutralAtomMapper(const Architecture&  arch,
                             const Configuration& config)
      : arch(arch), config(config) {}
  virtual ~NeutralAtomMapper() = default;

  auto               map(const qc::QuantumComputation& qc) -> void;
  [[nodiscard]] auto getResult() const -> const qc::QuantumComputation& {
    return mappedQc;
  }
  [[nodiscard]] auto getStats() const
      -> const na::NeutralAtomMapper::Statistics& {
    return stats;
  }
};
} // namespace na