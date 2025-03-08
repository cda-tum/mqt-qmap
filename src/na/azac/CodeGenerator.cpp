#include "na/azac/CodeGenerator.hpp"

#include "Definitions.hpp"
#include "ir/operations/OpType.hpp"
#include "ir/operations/Operation.hpp"
#include "na/NAComputation.hpp"
#include "na/azac/Architecture.hpp"
#include "na/entities/Atom.hpp"
#include "na/entities/Zone.hpp"
#include "na/operations/GlobalCZOp.hpp"
#include "na/operations/LoadOp.hpp"
#include "na/operations/LocalRZOp.hpp"
#include "na/operations/MoveOp.hpp"
#include "na/operations/StoreOp.hpp"

#include <cassert>
#include <cstddef>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <nlohmann/json_fwd.hpp>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

auto na::CodeGenerator::appendOneQubitGates(
    const std::vector<std::reference_wrapper<const qc::Operation>>&
        oneQubitGates,
    const std::vector<std::tuple<const SLM&, size_t, size_t>>& atomLocations,
    const std::vector<std::reference_wrapper<const Atom>>& atoms,
    NAComputation& code) const -> void {
  for (const auto& op : oneQubitGates) {
    assert(op.get().getNqubits() == 1);
    const qc::Qubit qubit = op.get().getTargets().front();
    const auto& [slm, r, c] = atomLocations[qubit];
    const auto& [x, y] = architecture_.get().exactSlmLocation(slm, r, c);
    assert(op.get().getType() == qc::Z);
    code.emplaceBack<LocalRZOp>(atoms[qubit], x, y);
  }
}
auto na::CodeGenerator::appendTwoQubitGates(
    const std::vector<std::tuple<const SLM&, size_t, size_t>>& currentLocations,
    const std::vector<std::vector<qc::Qubit>>& executionRouting,
    const std::vector<std::tuple<const SLM&, size_t, size_t>>&
        executionLocations,
    const std::vector<std::vector<qc::Qubit>>& targetRouting,
    const std::vector<std::tuple<const SLM&, size_t, size_t>>& targetLocations,
    const std::vector<std::reference_wrapper<const Atom>>& atoms,
    const Zone& zone, NAComputation& code) const -> void {
  appendRearrangement(currentLocations, executionRouting, executionLocations,
                      atoms, code);
  code.emplaceBack<GlobalCZOp>(zone);
  appendRearrangement(executionLocations, targetRouting, targetLocations, atoms,
                      code);
}
auto na::CodeGenerator::appendRearrangement(
    const std::vector<std::tuple<const SLM&, size_t, size_t>>& startSites,
    const std::vector<std::vector<qc::Qubit>>& routing,
    const std::vector<std::tuple<const SLM&, size_t, size_t>>& targetSites,
    const std::vector<std::reference_wrapper<const Atom>>& atoms,
    NAComputation& code) const -> void {
  for (const auto& qubits : routing) {
    std::map<size_t, std::map<size_t, qc::Qubit>> rowsWithQubits;
    std::vector<const Atom*> atomsToMove;
    std::vector<Location> targetLocations;
    for (const auto& qubit : qubits) {
      // get the current location of the qubit
      const auto& [slm, r, c] = startSites[qubit];
      const auto& [x, y] = architecture_.get().exactSlmLocation(slm, r, c);
      rowsWithQubits.try_emplace(y).first->second.emplace(x, qubit);
      atomsToMove.emplace_back(&atoms[qubit].get());
      // get the target location of the qubit
      const auto& [targetSlm, targetR, targetC] = targetSites[qubit];
      const auto& [targetX, targetY] =
          architecture_.get().exactSlmLocation(targetSlm, targetR, targetC);
      targetLocations.emplace_back(
          Location{static_cast<double>(targetX), static_cast<double>(targetY)});
    }
    std::vector<std::pair<qc::Qubit, std::pair<size_t, size_t>>>
        alreadyLoadedQubits;
    const auto& [minY, firstRow] = *rowsWithQubits.cbegin();
    std::vector<const Atom*> firstAtomsToLoad;
    firstAtomsToLoad.reserve(firstRow.size());
    for (const auto& [x, qubit] : firstRow) {
      alreadyLoadedQubits.emplace_back(qubit, std::make_pair(x, minY));
      firstAtomsToLoad.emplace_back(&atoms[qubit].get());
    }
    code.emplaceBack<LoadOp>(firstAtomsToLoad);
    // if there are more than one row with atoms to move, we pick them up
    // row-by-row as a simple strategy to avoid ghost-spots
    for (auto it = std::next(rowsWithQubits.cbegin());
         it != rowsWithQubits.cend(); ++it) {
      const auto& [y, row] = *it;
      // perform an offset move to avoid ghost-spots
      std::vector<const Atom*> atomsToOffset;
      std::vector<Location> offsetTargetLocations;
      atomsToOffset.reserve(alreadyLoadedQubits.size());
      offsetTargetLocations.reserve(alreadyLoadedQubits.size());
      for (const auto& [qubit, location] : alreadyLoadedQubits) {
        atomsToOffset.emplace_back(&atoms[qubit].get());
        const auto& [x, y] = location;
        if (row.find(x) != row.end()) {
          // new atoms get picked up in the column at x, i.e., only do a
          // vertical offset
          offsetTargetLocations.emplace_back(Location{
              static_cast<double>(x), static_cast<double>(y + parkingOffset_)});
        } else {
          // no new atoms get picked up in the column at x, i.e., do a
          // diagonal offset to avoid any ghost-spots
          offsetTargetLocations.emplace_back(
              Location{static_cast<double>(x + parkingOffset_),
                       static_cast<double>(y + parkingOffset_)});
        }
      }
      code.emplaceBack<MoveOp>(atomsToOffset, offsetTargetLocations);
      // load the new atoms
      std::vector<const Atom*> atomsToLoad;
      atomsToLoad.reserve(row.size());
      for (const auto& [x, qubit] : row) {
        alreadyLoadedQubits.emplace_back(qubit, std::make_pair(x, y));
        atomsToLoad.emplace_back(&atoms[qubit].get());
      }
      code.emplaceBack<LoadOp>(atomsToLoad);
    }
    // all atoms are loaded, now move them to their target locations
    code.emplaceBack<MoveOp>(atomsToMove, targetLocations);
    code.emplaceBack<StoreOp>(atomsToMove);
  }
}
na::CodeGenerator::CodeGenerator(const Architecture& architecture,
                                 const nlohmann::json& config)
    : architecture_(architecture) {
  if (const auto& configIt = config.find("parking_offset");
      configIt != config.end()) {
    const auto& generatorConfig = configIt.value();
    if (const auto& it = generatorConfig.find("parking_offset");
        it != generatorConfig.end() && it->is_number_unsigned()) {
      parkingOffset_ = it.value();
    } else {
      std::cout << "[WARN] Configuration for CodeGenerator does not contain "
                   "a valid parking offset. Using default offset.\n";
    }
  } else {
    std::cout << "[WARN] Configuration does not contain settings for "
                 "CodeGenerator. Using default settings.\n";
  }
}
auto na::CodeGenerator::generateCode(
    const std::vector<std::vector<std::reference_wrapper<const qc::Operation>>>&
        oneQubitGateLayers,
    const std::vector<std::vector<std::tuple<const SLM&, size_t, size_t>>>&
        placement,
    const std::vector<std::vector<std::vector<qc::Qubit>>>& routing) const
    -> NAComputation {
  NAComputation code;
  const auto& rydbergZone = code.emplaceBackZone("zone_cz0");
  const auto& initialPlacement = placement.front();
  std::vector<std::reference_wrapper<const Atom>> atoms;
  atoms.reserve(initialPlacement.size());
  for (const auto& [slm, r, c] : initialPlacement) {
    atoms.emplace_back(
        code.emplaceBackAtom("atom" + std::to_string(atoms.size())));
    const auto& [x, y] = architecture_.get().exactSlmLocation(slm, r, c);
    code.emplaceInitialLocation(atoms.back(), x, y);
  }
  assert(2 * oneQubitGateLayers.size() == placement.size() + 1);
  assert(placement.size() == routing.size() + 1);
  for (size_t layer = 0; true; ++layer) {
    const auto& oneQubitGates = oneQubitGateLayers[layer];
    const auto& atomLocations = placement[2 * layer];
    appendOneQubitGates(oneQubitGates, atomLocations, atoms, code);
    if (layer == oneQubitGateLayers.size() - 1) {
      break;
    }
    appendTwoQubitGates(atomLocations, routing[2 * layer],
                        placement[(2 * layer) + 1], routing[(2 * layer) + 1],
                        placement[2 * (layer + 1)], atoms, rydbergZone, code);
  }
  return code;
}
