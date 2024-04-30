//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/mqt-qmap for more
// information.
//

#pragma once

#include "Architecture.hpp"
#include "Configuration.hpp"
#include "QuantumComputation.hpp"
#include "na/NAComputation.hpp"
#include "operations/Operation.hpp"

#include <cstdint>
#include <sstream>
#include <vector>

namespace na {

class NeutralAtomMapper {
public:
  struct Statistics {
    std::size_t               numInitialGates    = 0;
    std::size_t               numEntanglingGates = 0;
    std::size_t               initialDepth       = 0;
    std::size_t               numMappedGates     = 0;
    std::size_t               numQubits          = 0;
    std::size_t               maxSeqWidth        = 0;
    qc::fp                    preprocessTime     = 0.0; // [ms]
    qc::fp                    mappingTime        = 0.0; // [ms]
    qc::fp                    postprocessTime    = 0.0; // [ms]
    [[nodiscard]] static auto header() -> std::string {
      return "numInitialGates,numEntanglingGates,initialDepth,numMappedGates,"
             "numQubits,maxSeqWidth,preprocessTime,mappingTime,"
             "postprocessTime\n";
    }
    [[nodiscard]] auto toString() const -> std::string {
      std::stringstream ss;
      ss << numInitialGates << "," << numEntanglingGates << "," << initialDepth
         << "," << numMappedGates << "," << numQubits << "," << maxSeqWidth
         << "," << preprocessTime << "," << mappingTime << ","
         << postprocessTime << '\n';
      return ss.str();
    }
    friend auto operator<<(std::ostream&     os,
                           const Statistics& s) -> std::ostream& {
      os << s.toString();
      return os;
    }
  };

protected:
  qc::QuantumComputation initialQc;
  NAComputation          mappedQc;
  Architecture           initialArch;
  Architecture           arch;
  Configuration          config;
  Statistics             stats{};
  bool                   done = false;

  class Atom {
  public:
    enum class PositionStatus : std::uint8_t { UNDEFINED, DEFINED };
    PositionStatus         positionStatus  = PositionStatus::UNDEFINED;
    std::shared_ptr<Point> initialPosition = std::make_shared<Point>(0, 0);
    std::shared_ptr<Point> currentPosition = initialPosition;
    std::vector<Zone>      zones;
    explicit Atom(const std::vector<Zone>& zones) : zones(zones) {};
    explicit Atom() : Atom({}) {};
  };
  auto preprocess() -> void { validateCircuit(); }
  auto validateCircuit() -> void;
  auto postprocess() -> void {
    makeLogicalArrays();
    calculateMovements();
  }
  auto makeLogicalArrays() -> void;
  auto calculateMovements() -> void;
  [[nodiscard]] auto
       checkApplicability(const qc::Operation* op,
                          const std::vector<Atom>& placement) const -> bool;
  auto updatePlacement(const qc::Operation* op,
                       std::vector<Atom>& placement) const -> void;
  [[nodiscard]] static auto
  getMisplacement(const std::vector<Atom>&      initial,
                  const std::vector<qc::Qubit>& target,
                  const qc::Qubit&              q) -> std::int64_t;
  /**
   *
   * @param initialFreeSites The sites that are not yet occupied from the start
   * of the circuit until the current execution step.
   * @param currentFreeSites The currently not occupied sites.
   * @param placement The current placement of the atoms. Will be modified
   * throughout the course of executing the function.
   * @param currentlyShuttling The qubits that are currently being shuttled.
   * @param qubits The qubits to be shuttled.
   * @param destination The destination zone.
   */
  auto store(std::vector<bool>&             initialFreeSites,
               std::vector<bool>&             currentFreeSites,
               std::vector<Atom>&             placement,
               std::unordered_set<qc::Qubit>& currentlyShuttling,
               const std::vector<qc::Qubit>& qubits, Zone destination) -> void;
  /**
   *
   * @param initialFreeSites The sites that are not yet occupied from the start
   * of the circuit until the current execution step.
   * @param currentFreeSites The currently not occupied sites.
   * @param placement The current placement of the atoms. Will be modified
   * throughout the course of executing the function.
   * @param currentlyShuttling The qubits that are currently being shuttled.
   * @param qubits The qubits to be shuttled mapped to their final position.
   */
  auto pickUp(std::vector<bool>&                          initialFreeSites,
               std::vector<bool>&                          currentFreeSites,
               std::vector<Atom>&                          placement,
               std::unordered_set<qc::Qubit>&              currentlyShuttling,
               const std::vector<qc::Qubit>& qubits) -> void;

public:
  explicit NeutralAtomMapper(Architecture arch, const Configuration& config)
      : initialArch(std::move(arch)), arch(initialArch.withConfig(config)),
        config(config) {}
  virtual ~NeutralAtomMapper() = default;
  auto               map(const qc::QuantumComputation& qc) -> void;
  [[nodiscard]] auto getResult() const -> const NAComputation& {
    if (!done) {
      throw std::logic_error("No result available.");
    }
    return mappedQc;
  }
  [[nodiscard]] auto
  getStats() const -> const na::NeutralAtomMapper::Statistics& {
    if (!done) {
      throw std::logic_error("No statistics available.");
    }
    return stats;
  }
};
} // namespace na
