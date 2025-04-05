#include "na/azac/CodeGenerator.hpp"

#include "ir/Definitions.hpp"
#include "ir/operations/CompoundOperation.hpp"
#include "ir/operations/OpType.hpp"
#include "ir/operations/Operation.hpp"
#include "na/NAComputation.hpp"
#include "na/azac/Architecture.hpp"
#include "na/azac/Types.hpp"
#include "na/entities/Atom.hpp"
#include "na/entities/Location.hpp"
#include "na/entities/Zone.hpp"
#include "na/operations/GlobalCZOp.hpp"
#include "na/operations/GlobalRYOp.hpp"
#include "na/operations/LoadOp.hpp"
#include "na/operations/LocalRZOp.hpp"
#include "na/operations/LocalUOp.hpp"
#include "na/operations/MoveOp.hpp"
#include "na/operations/StoreOp.hpp"

#include <cassert>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <map>
#include <nlohmann/json_fwd.hpp>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace na::azac {
auto CodeGenerator::appendOneQubitGates(
    const size_t nQubits, const OneQubitGateLayer& oneQubitGates,
    const std::vector<std::reference_wrapper<const Atom>>& atoms,
    const Zone& globalZone, NAComputation& code) const -> void {
  for (const auto& op : oneQubitGates) {
    bool oneQubitGate = false;
    if (op.get().isGlobal(nQubits)) {
      if (op.get().isCompoundOperation()) {
        const auto& compOp =
            dynamic_cast<const qc::CompoundOperation&>(op.get());
        const auto opType = compOp.front()->getType();
        if (opType == qc::RY) {
          code.emplaceBack<GlobalRYOp>(globalZone,
                                       compOp.front()->getParameter().front());
        } else if (opType == qc::Y) {
          code.emplaceBack<GlobalRYOp>(globalZone, qc::PI);
        } else if (nQubits == 1) {
          oneQubitGate =
              true; // special case for one qubit, fall back to local gate
        } else {
          assert(false);
        }
      } else {
        const auto opType = op.get().getType();
        if (opType == qc::RY) {
          code.emplaceBack<GlobalRYOp>(globalZone,
                                       op.get().getParameter().front());
        } else if (opType == qc::Y) {
          code.emplaceBack<GlobalRYOp>(globalZone, qc::PI);
        } else if (nQubits == 1) {
          oneQubitGate =
              true; // special case for one qubit, fall through to local gate
        } else {
          assert(false);
        }
      }
    } else {
      oneQubitGate = true;
    }
    if (oneQubitGate) {
      assert(op.get().getNqubits() == 1);
      const qc::Qubit qubit = op.get().getTargets().front();
      if (op.get().getType() == qc::RZ) {
        code.emplaceBack<LocalRZOp>(atoms[qubit],
                                    op.get().getParameter().front());
      } else if (op.get().getType() == qc::Z) {
        code.emplaceBack<LocalRZOp>(atoms[qubit], qc::PI);
      } else if (op.get().getType() == qc::S) {
        code.emplaceBack<LocalRZOp>(atoms[qubit], qc::PI_2);
      } else if (op.get().getType() == qc::Sdg) {
        code.emplaceBack<LocalRZOp>(atoms[qubit], -qc::PI_2);
      } else if (op.get().getType() == qc::T) {
        code.emplaceBack<LocalRZOp>(atoms[qubit], qc::PI_4);
      } else if (op.get().getType() == qc::Tdg) {
        code.emplaceBack<LocalRZOp>(atoms[qubit], -qc::PI_4);
      } else if (op.get().getType() == qc::P) {
        code.emplaceBack<LocalRZOp>(atoms[qubit],
                                    op.get().getParameter().front());
      } else {
        std::ostringstream oss;
        oss << "\033[1;35m[WARN]\033[0m Gate not part of basis gates will be "
               "inserted as U3 gate: "
            << op.get().getType() << "\n";
        std::cout << oss.str();
        if (op.get().getType() == qc::U) {
          code.emplaceBack<LocalUOp>(
              atoms[qubit], op.get().getParameter().front(),
              op.get().getParameter().at(1), op.get().getParameter().at(2));
        } else if (op.get().getType() == qc::U2) {
          code.emplaceBack<LocalUOp>(atoms[qubit], qc::PI_2,
                                     op.get().getParameter().front(),
                                     op.get().getParameter().at(1));
        } else if (op.get().getType() == qc::RX) {
          code.emplaceBack<LocalUOp>(atoms[qubit],
                                     op.get().getParameter().front(), -qc::PI_2,
                                     qc::PI_2);
        } else if (op.get().getType() == qc::RY) {
          code.emplaceBack<LocalUOp>(atoms[qubit],
                                     op.get().getParameter().front(), 0, 0);
        } else if (op.get().getType() == qc::H) {
          code.emplaceBack<LocalUOp>(atoms[qubit], qc::PI_2, 0, qc::PI);
        } else if (op.get().getType() == qc::X) {
          code.emplaceBack<LocalUOp>(atoms[qubit], qc::PI, 0, qc::PI);
        } else if (op.get().getType() == qc::Y) {
          code.emplaceBack<LocalUOp>(atoms[qubit], qc::PI, qc::PI_2, qc::PI_2);
        } else if (op.get().getType() == qc::V) {
          code.emplaceBack<LocalUOp>(atoms[qubit], -qc::PI_2, -qc::PI_2,
                                     qc::PI_2);
        } else if (op.get().getType() == qc::Vdg) {
          code.emplaceBack<LocalUOp>(atoms[qubit], -qc::PI_2, qc::PI_2,
                                     -qc::PI_2);
        } else if (op.get().getType() == qc::SX) {
          code.emplaceBack<LocalUOp>(atoms[qubit], qc::PI_2, -qc::PI_2,
                                     qc::PI_2);
        } else if (op.get().getType() == qc::SXdg) {
          code.emplaceBack<LocalUOp>(atoms[qubit], -qc::PI_2, -qc::PI_2,
                                     qc::PI_2);
        } else {
          oss.clear();
          oss << "\033[1;31m[ERROR]\033[0m Unsupported one qubit gate will be "
                 "dropped: "
              << op.get().getType() << "\n";
          std::cout << oss.str();
          assert(false);
        }
      }
    }
  }
}
auto CodeGenerator::appendTwoQubitGates(
    const Placement& currentPlacement, const Routing& executionRouting,
    const Placement& executionPlacement, const Routing& targetRouting,
    const Placement& targetPlacement,
    const std::vector<std::reference_wrapper<const Atom>>& atoms,
    const std::vector<std::reference_wrapper<const Zone>>& zones,
    NAComputation& code) const -> void {
  appendRearrangement(currentPlacement, executionRouting, executionPlacement,
                      atoms, code);
  std::vector<const Zone*> zonePtrs;
  zonePtrs.reserve(zones.size());
  std::transform(zones.begin(), zones.end(), std::back_inserter(zonePtrs),
                 [](const auto& zone) { return &zone.get(); });
  code.emplaceBack<GlobalCZOp>(zonePtrs);
  appendRearrangement(executionPlacement, targetRouting, targetPlacement, atoms,
                      code);
}
auto CodeGenerator::appendRearrangement(
    const Placement& startPlacement, const Routing& routing,
    const Placement& targetPlacement,
    const std::vector<std::reference_wrapper<const Atom>>& atoms,
    NAComputation& code) const -> void {
  for (const auto& qubits : routing) {
    std::map<size_t, std::map<size_t, qc::Qubit>> rowsWithQubits;
    std::vector<const Atom*> atomsToMove;
    std::vector<Location> targetLocations;
    for (const auto& qubit : qubits) {
      // get the current location of the qubit
      const auto& [slm, r, c] = startPlacement[qubit];
      const auto& [x, y] = architecture_.get().exactSlmLocation(slm, r, c);
      rowsWithQubits.try_emplace(y).first->second.emplace(x, qubit);
      atomsToMove.emplace_back(&atoms[qubit].get());
      // get the target location of the qubit
      const auto& [targetSlm, targetR, targetC] = targetPlacement[qubit];
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
      alreadyLoadedQubits.emplace_back(qubit, std::pair{x, minY});
      firstAtomsToLoad.emplace_back(&atoms[qubit].get());
    }
    code.emplaceBack<LoadOp>(firstAtomsToLoad);
    // if there are more than one row with atoms to move, we pick them up
    // row-by-row as a simple strategy to avoid ghost-spots
    for (auto it = std::next(rowsWithQubits.cbegin());
         it != rowsWithQubits.cend(); ++it) {
      const auto& [yCoordinateOfRow, row] = *it;
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
        alreadyLoadedQubits.emplace_back(qubit, std::pair{x, yCoordinateOfRow});
        atomsToLoad.emplace_back(&atoms[qubit].get());
      }
      code.emplaceBack<LoadOp>(atomsToLoad);
    }
    // all atoms are loaded, now move them to their target locations
    code.emplaceBack<MoveOp>(atomsToMove, targetLocations);
    code.emplaceBack<StoreOp>(atomsToMove);
  }
}
CodeGenerator::CodeGenerator(const Architecture& architecture,
                             const nlohmann::json& config)
    : architecture_(architecture) {
  if (const auto& configIt = config.find("code_generator");
      configIt != config.end() && configIt->is_object()) {
    bool parkingOffsetSet = false;
    for (const auto& [key, value] : configIt.value().items()) {
      if (key == "parking_offset") {
        if (value.is_number_unsigned()) {
          parkingOffset_ = value;
          parkingOffsetSet = true;
        } else {
          std::ostringstream oss;
          oss << "\033[1;35m[WARN]\033[0m Configuration for CodeGenerator "
                 "contains an invalid value for parking_offset. Using "
                 "default.\n";
          std::cout << oss.str();
        }
      } else {
        std::ostringstream oss;
        oss << "\033[1;35m[WARN]\033[0m Configuration for CodeGenerator "
               "contains an unknown key: "
            << key << ". Ignoring.\n";
        std::cout << oss.str();
      }
    }
    if (!parkingOffsetSet) {
      std::cout
          << "\033[1;35m[WARN]\033[0m Configuration for CodeGenerator "
             "does not contain a value for parking_offset. Using default.\n";
    }
  } else {
    std::cout << "\033[1;35m[WARN]\033[0m Configuration does not contain "
                 "settings for CodeGenerator or is malformed. Using default "
                 "settings.\n";
  }
}
auto CodeGenerator::generate(
    const std::vector<OneQubitGateLayer>& oneQubitGateLayers,
    const std::vector<Placement>& placement,
    const std::vector<Routing>& routing) const -> NAComputation {
  NAComputation code;
  std::vector<std::reference_wrapper<const Zone>> rydbergZones;
  for (size_t i = 0; i < architecture_.get().rydbergRangeMinX.size(); ++i) {
    rydbergZones.emplace_back(code.emplaceBackZone(
        "zone_cz" + std::to_string(i),
        Zone::Extent{
            static_cast<double>(architecture_.get().rydbergRangeMinX.at(i)),
            static_cast<double>(architecture_.get().rydbergRangeMinY.at(i)),
            static_cast<double>(architecture_.get().rydbergRangeMaxX.at(i)),
            static_cast<double>(architecture_.get().rydbergRangeMaxY.at(i))}));
  }
  const auto& globalZone = code.emplaceBackZone(
      "global",
      Zone::Extent{static_cast<double>(architecture_.get().archRangeMinX),
                   static_cast<double>(architecture_.get().archRangeMinY),
                   static_cast<double>(architecture_.get().archRangeMaxX),
                   static_cast<double>(architecture_.get().archRangeMaxY)});
  const auto& initialPlacement = placement.front();
  std::vector<std::reference_wrapper<const Atom>> atoms;
  atoms.reserve(initialPlacement.size());
  for (const auto& [slm, r, c] : initialPlacement) {
    atoms.emplace_back(
        code.emplaceBackAtom("atom" + std::to_string(atoms.size())));
    const auto& [x, y] = architecture_.get().exactSlmLocation(slm, r, c);
    code.emplaceInitialLocation(atoms.back(), x, y);
  }
  // early return if no one-qubit gates are given
  if (oneQubitGateLayers.empty()) {
    return code;
  }
  assert(2 * oneQubitGateLayers.size() == placement.size() + 1);
  assert(placement.size() == routing.size() + 1);
  appendOneQubitGates(atoms.size(), oneQubitGateLayers.front(), atoms,
                      globalZone, code);
  for (size_t layer = 0; layer + 1 < oneQubitGateLayers.size(); ++layer) {
    appendTwoQubitGates(placement[2 * layer], routing[2 * layer],
                        placement[(2 * layer) + 1], routing[(2 * layer) + 1],
                        placement[2 * (layer + 1)], atoms, rydbergZones, code);
    appendOneQubitGates(atoms.size(), oneQubitGateLayers[layer + 1], atoms,
                        globalZone, code);
  }
  return code;
}
} // namespace na::azac
