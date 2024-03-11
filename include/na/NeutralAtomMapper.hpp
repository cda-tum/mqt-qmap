//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/mqt-qmap for more
// information.
//

#pragma once

#include "Architecture.hpp"
#include "Configuration.hpp"
#include "QuantumComputation.hpp"
#include "operations/NAQuantumComputation.hpp"
#include "operations/Operation.hpp"

#include <cstdint>
#include <memory>
#include <sstream>
#include <utility>
#include <vector>

namespace na {

class NeutralAtomMapper {
public:
  struct Statistics {
    std::size_t        numInitialGates = 0;
    std::size_t        initialDepth    = 0;
    std::size_t        numMappedGates  = 0;
    std::size_t        numQubits       = 0;
    std::size_t        maxSeqWidth     = 0;
    qc::fp             preprocessTime  = 0.0; // [ms]
    qc::fp             mappingTime     = 0.0; // [ms]
    qc::fp             postprocessTime = 0.0; // [ms]
    [[nodiscard]] static auto header() -> std::string {
      return "numInitialGates,initialDepth,numMappedGates,numQubits,"
             "maxSeqWidth,preprocessTime,mappingTime,postprocessTime\n";
    }
    [[nodiscard]] auto toString() const -> std::string {
      std::stringstream ss;
      ss << numInitialGates << "," << initialDepth << "," << numMappedGates
         << "," << numQubits << "," << maxSeqWidth << "," << preprocessTime
         << "," << mappingTime << "," << postprocessTime << std::endl;
      return ss.str();
    }
    friend auto operator<< (std::ostream& os, const Statistics& s) -> std::ostream& {
      os << s.toString();
      return os;
    }
  };

protected:
  qc::QuantumComputation            initialQc{};
  NAQuantumComputation              mappedQc{};
  na::Architecture                  initialArch;
  na::Architecture                  arch;
  na::Configuration                 config;
  na::NeutralAtomMapper::Statistics stats{};
  bool                              done = false;

  class Atom {
  public:
    enum class PositionStatus { UNDEFINED, DEFINED };
    PositionStatus         positionStatus  = PositionStatus::UNDEFINED;
    std::shared_ptr<Point> initialPosition = std::make_shared<Point>(0, 0);
    std::shared_ptr<Point> currentPosition = initialPosition;
    std::vector<Zone>      zones;
    explicit Atom(const std::vector<Zone>& zones) : zones(zones){};
    explicit Atom() : Atom({}){};
  };
  auto preprocess() -> void;
  auto postprocess() -> void;
  [[nodiscard]] auto
       checkApplicability(const std::unique_ptr<qc::Operation>& op,
                          const std::vector<Atom>& placement) const -> bool;
  auto updatePlacement(const std::unique_ptr<qc::Operation>& op,
                       std::vector<Atom>& placement) const -> void;
  [[nodiscard]] auto
  getMisplacement(const std::vector<Atom>&                           initial,
                  const std::unordered_map<qc::Qubit, std::int64_t>& target,
                  const qc::Qubit& q) const -> std::int64_t;

public:
  explicit NeutralAtomMapper(Architecture arch, const Configuration& config)
      : initialArch(std::move(arch)), arch(initialArch.withConfig(config)),
        config(config) {}
  virtual ~NeutralAtomMapper() = default;
  auto               map(const qc::QuantumComputation& qc) -> void;
  [[nodiscard]] auto getResult() const -> const NAQuantumComputation& {
    if (!done) {
      throw std::logic_error("No result available.");
    }
    return mappedQc;
  }
  [[nodiscard]] auto getStats() const
      -> const na::NeutralAtomMapper::Statistics& {
    if (!done) {
      throw std::logic_error("No statistics available.");
    }
    return stats;
  }
};
} // namespace na