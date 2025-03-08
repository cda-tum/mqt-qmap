#pragma once

#include "ir/operations/Operation.hpp"
#include "na/NAComputation.hpp"
#include "na/azac/Architecture.hpp"
#include "na/entities/Atom.hpp"
#include "na/operations/LoadOp.hpp"

#include <cassert>
#include <functional>
#include <string>
#include <unordered_set>

namespace na {
/**
 * The class CodeGenerator implements the code generation for the zoned neutral
 * atom compiler.
 */
class CodeGenerator {
  std::reference_wrapper<const Architecture> architecture_;

protected:
  /**
   * Create a new CodeGenerator.
   * @details The code generation is based on the given architecture and the
   * placement and routing of the qubits. It generates a neutral atom
   * computation. The second parameter of the constructor is unused.
   * @param architecture is the architecture of the neutral atom system
   */
  CodeGenerator(const Architecture& architecture,
                const nlohmann::json& /* unused */)
      : architecture_(architecture) {}
  [[nodiscard]] auto generateCode(
      const std::vector<std::vector<
          std::reference_wrapper<const qc::Operation>>>& oneQubitGateLayers,
      const std::vector<std::vector<std::tuple<const SLM&, size_t, size_t>>>&
          placement,
      const std::vector<std::vector<std::vector<qc::Qubit>>>& routing) const
      -> NAComputation {
    NAComputation code;
    const auto& initialPlacement = placement.front();
    std::vector<std::reference_wrapper<const Atom>> atoms;
    atoms.reserve(initialPlacement.size());
    for (const auto& [slm, r, c] : initialPlacement) {
      atoms.emplace_back(
          code.emplaceBackAtom("atom" + std::to_string(atoms.size())));
      const auto& [x, y] = architecture_.get().exactSlmLocation(slm, r, c);
      code.emplaceInitialLocation(atoms.back(), x, y);
    }
    assert(oneQubitGateLayers.size() == placement.size() + 1);
    assert(placement.size() == routing.size());
    std::unordered_set<std::reference_wrapper<const Atom>> loadedAtoms;
    for (size_t layer = 0; true; ++layer) {
      const auto& oneQubitGates = oneQubitGateLayers[layer];
      for (const auto& op : oneQubitGates) {
        // todo
      }
      if (layer == placement.size()) {
        break;
      }
      const auto& atomLocations = placement[layer];
      const auto& routingSteps = routing[layer];
      std::vector<const Atom*> atomsToLoad;
      for (const auto& routingStep : routingSteps) {
      }
      code.emplaceBack<LoadOp>(atomsToLoad);
    }
    return code;
  }
};
} // namespace na
