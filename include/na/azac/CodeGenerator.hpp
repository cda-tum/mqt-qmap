#pragma once

#include "Types.hpp"
#include "na/NAComputation.hpp"
#include "na/azac/Architecture.hpp"
#include "na/entities/Atom.hpp"
#include "na/entities/Zone.hpp"

#include <cstddef>
#include <nlohmann/json_fwd.hpp>
#include <vector>

namespace na::azac {
/**
 * The class CodeGenerator implements the code generation for the zoned neutral
 * atom compiler.
 */
class CodeGenerator {
  /// The architecture of the neutral atom system
  std::reference_wrapper<const Architecture> architecture_;
  /// The offset for parking spots
  size_t parkingOffset_ = 1;

public:
  /**
   * Create a new CodeGenerator.
   * @details The code generation is based on the given architecture and the
   * placement and routing of the qubits. It generates a neutral atom
   * computation. The second parameter of the constructor is unused.
   * @param architecture is the architecture of the neutral atom system
   * @param config is the configuration of the code generator
   */
  CodeGenerator(const Architecture& architecture, const nlohmann::json& config);
  /**
   * Generate the neutral atom computation based on the results of the previous
   * steps in the compiler.
   * @param oneQubitGateLayers is a list of layers of single-qubit gates
   * @param placement is the placement of the qubits. The very first entry is
   * the initial placement of the atoms. Every consecutive pair of entries
   * encloses one layer of two-qubit gates.
   * @param routing is the routing of the qubits. It consists of groups of atoms
   * that can be moved together to establish the next placement.
   * @return the neutral atom computation
   */
  [[nodiscard]] auto
  generate(const std::vector<OneQubitGateLayer>& oneQubitGateLayers,
           const std::vector<Placement>& placement,
           const std::vector<Routing>& routing) const -> NAComputation;

private:
  /// Append all one-qubit gates of a layer to the code
  auto appendOneQubitGates(
      size_t nQubits, const OneQubitGateLayer& oneQubitGates,
      const std::vector<std::reference_wrapper<const Atom>>& atoms,
      const Zone& globalZone, NAComputation& code) const -> void;

  /// Append all necessary operations to perform the next set of two-qubit gates
  auto appendTwoQubitGates(
      const Placement& currentPlacement, const Routing& executionRouting,
      const Placement& executionPlacement, const Routing& targetRouting,
      const Placement& targetPlacement,
      const std::vector<std::reference_wrapper<const Atom>>& atoms,
      const std::vector<std::reference_wrapper<const Zone>>& zones,
      NAComputation& code) const -> void;

  /// Append all necessary operations to rearrange the atoms
  auto appendRearrangement(
      const Placement& startPlacement, const Routing& routing,
      const Placement& targetPlacement,
      const std::vector<std::reference_wrapper<const Atom>>& atoms,
      NAComputation& code) const -> void;
};
} // namespace na::azac
