//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "NeutralAtomMapper.hpp"

#include "Layer.hpp"
#include "QuantumComputation.hpp"
#include "operations/Operation.hpp"
#include "operations/StandardOperation.hpp"

#include <algorithm>
#include <array>
#include <memory>
#include <stdexcept>
#include <vector>

namespace na {

/// Construct an adjacency matrix for the interaction graph for the given set of
/// operations
auto constructInteractionGraph(const std::set<Layer::DAGVertex*>& ops)
    -> std::vector<std::vector<bool>> {
  // get the number of qubits
  const auto nInvolvedQubits = ops.size();
  // create the adjacency matrix
  std::vector<std::vector<bool>> adjacencyMatrix(nInvolvedQubits - 1);
  // fill the adjacency matrix
  for (const auto& op : ops) {
    const std::unique_ptr<qc::Operation>& gate = *op->getOperation();
    if (gate->getNcontrols() + gate->getNtargets() != 2) {
      throw std::invalid_argument(
          "The interaction graph can only be constructed for two-qubit gates.");
    }
    const auto&      qubits = gate->getUsedQubits();
    auto             it = qubits.begin(); // Change the type of 'it' to 'auto'
    const qc::Qubit& q1 = *it;
    const qc::Qubit& q2 = *(++it); // Increment the iterator
  }
  return adjacencyMatrix;
}

auto NeutralAtomMapper::map(const qc::QuantumComputation& qc) -> void {
  initialQc    = qc;
  auto nqubits = initialQc.getNqubits();
  mappedQc     = qc::QuantumComputation(nqubits);
  // store the placement of atoms, both the initial one (needed later) and the
  // current one leave atoms unmapped as long as possible. This mighty induce
  // some restrictions on the mapping in which zone the atom may be placed.
  std::vector<Atom>                          initialPlacement(nqubits);
  std::vector<Atom>                          currentPlacement(nqubits);
  std::vector<std::unordered_set<qc::Qubit>> atomsInZone(arch.getNZones());
  // get start time
  auto start = std::chrono::high_resolution_clock::now();
  //============================ START MAPPING =================================
  Layer const layer         = Layer(initialQc);
  const auto& executableSet = layer.getExecutableSet();
continue_with_updated_set:
  while (!(*executableSet)->empty()) {
    // 1. execute all gates that are directly applicable and do not need
    //    shuttling
    for (const auto& vertex : **executableSet) {
      // the following check is necessary because the executable set may change
      // due to the execution of another gate
      if (!vertex->isExecutable()) {
        // break the inner loop and continue the outer one
        goto continue_with_updated_set; // NOLINT(cppcoreguidelines-avoid-goto)
      } else {
        const std::unique_ptr<qc::Operation>& gate = *vertex->getOperation();
        auto const nQubits = gate->getNcontrols() + gate->getNtargets();
        if (nQubits > 2) {
          throw std::invalid_argument("Global multi-qubit operations acting "
                                      "on more than two qubits are not "
                                      "supported yet.");
        }
        if (gate->isGlobalOperation()) {
          const std::unique_ptr<GlobalOperation>& ggate =
              dynamic_cast<std::unique_ptr<GlobalOperation>&>(*gate);
          // check availability of the global operation on the architecture
          if (!arch.isAllowedGlobally(ggate->getOpsType(),
                                      ggate->getNcontrols())) {
            std::stringstream ss;
            ss << "The chosen architecture does not support the global "
                  "operation ";
            for (std::size_t i = 0; i < ggate->getNcontrols(); ++i) {
              ss << "c";
            }
            ss << qc::toString(ggate->getOpsType());
            throw std::runtime_error(ss.str());
          }
          // check whether the gate is allowed in the current zone
          for (const auto& qubit : ggate->getUsedQubits()) {
            switch (currentPlacement[qubit].pos) {
            case Atom::Position::UNDEFINED:
              break;
            case Atom::Position::ZONE:
              for (const auto& z : currentPlacement[qubit].zone) {
                if (!arch.isAllowedGlobally(ggate->getOpsType(), z,
                                            ggate->getNcontrols())) {
                  goto continue_applicable_gates; // NOLINT(cppcoreguidelines-avoid-goto)
                }
              }
              break;
            case Atom::Position::DEFINED:
              const auto& z = arch.getZone(currentPlacement[qubit].qubit);
              if (!arch.isAllowedGlobally(ggate->getOpsType(), z,
                                          ggate->getNcontrols())) {
                goto continue_applicable_gates; // NOLINT(cppcoreguidelines-avoid-goto)
                break;
              }
              break;
            }
          }
          // check whether the global operation affects all atoms in the zone;
          // first figure out the zone of each affected atom
          for (const auto& qubit : ggate->getUsedQubits()) {
            switch (currentPlacement[qubit].pos) {
            case Atom::Position::UNDEFINED:
              currentPlacement[qubit].pos = Atom::Position::ZONE;
              for (Zone z = 0; z < arch.getNZones(); ++z) {
                if (arch.isAllowedGlobally(ggate->getOpsType(), z,
                                           ggate->getNcontrols())) {
                  initialPlacement[qubit].zone.push_back(z);
                  currentPlacement[qubit].zone.push_back(z);
                  atomsInZone[z].emplace(qubit);
                }
              }
              break;
            case Atom::Position::ZONE:
            case Atom::Position::DEFINED:
              break;
            }
          }
          // for every zone the global operations acts on, check whether
          // the set of atoms in this zone are a subset usedQubits
          for (Zone z = 0; z < arch.getNZones(); ++z) {
            std::vector<qc::Qubit> diff;
            std::set_difference(ggate->getUsedQubits().begin(),
                                ggate->getUsedQubits().end(),
                                atomsInZone[z].begin(), atomsInZone[z].end(),
                                std::back_inserter(diff));
            if (!diff.empty()) {
              if (nQubits == 2) {
                for (const auto& q1 : diff) {
                  for (const auto q2 : atomsInZone[z]) {
                    if (arch.getDistance(q1, q2) <=
                        arch.getInteractionRadius()) {
                      // TODO move uninvolved atoms out of the zone if there is
                      // one
                      throw std::invalid_argument(
                          "The global operation does not act on all atoms in "
                          "the zone. The current implementation does not "
                          "support to move away those atoms.");
                    }
                    // if q1 has no interaction partner, it will not experience
                    // any effect and it does not matter that it is in the zone;
                    // however, the fidelity might be worse
                  }
                }
              }
              // TODO move uninvolved atoms out of the zone if there is one
              throw std::invalid_argument(
                  "The global operation does not act on all atoms in the "
                  "zone. The current implementation does not support to move "
                  "away those atoms.");
            }
          }
          vertex->execute();
          mappedQc.emplace_back<std::unique_ptr<qc::Operation>>(gate->clone());
        } else {
          if (!arch.isAllowedLocally(gate->getType(), gate->getNcontrols())) {
            std::stringstream ss;
            ss << "The chosen architecture does not support the local "
                  "operation ";
            for (std::size_t i = 0; i < gate->getNcontrols(); ++i) {
              ss << "c";
            }
            ss << qc::toString(gate->getType());
            throw std::runtime_error(ss.str());
          }
          auto const nUsedQubits = gate->getNcontrols() + gate->getNtargets();
          if (nUsedQubits == 1) {
            // TODO check whether the gate is allowed in the current zone
            vertex->execute();
            mappedQc.emplace_back<std::unique_ptr<qc::Operation>>(
                gate->clone());
          } else {
            throw std::invalid_argument("Local multi-qubit operations acting "
                                        "on more than one qubits are not "
                                        "supported yet.");
          }
        }
      }
    continue_applicable_gates:
      continue;
    }
    // 2. when no such gates are left, extract an interaction graph of gates
    //    of the same type and two targets, i.e. cz gates
    // 3. move the atoms accordingly and execute the gates
    // 4. move the atoms back to their storage zone
  }
  //============================= END MAPPING ==================================
  // get end time
  auto end = std::chrono::high_resolution_clock::now();
  // build remaining statistics
  stats.numInitialGates = qc.getNops();
  stats.numMappedGates  = qc.getNops();
  // get the mapping time in milliseconds
  stats.mappingTime =
      std::chrono::duration<qc::fp, std::milli>(end - start).count();
  }
} // namespace na