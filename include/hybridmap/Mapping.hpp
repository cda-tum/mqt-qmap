//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "hybridmap/NeutralAtomDefinitions.hpp"
#include "hybridmap/NeutralAtomUtils.hpp"
#include "ir/Definitions.hpp"
#include "ir/Permutation.hpp"
#include "ir/operations/Operation.hpp"

#include <algorithm>
#include <cstddef>
#include <set>
#include <stdexcept>
#include <string>

namespace na {

/**
 * @brief Class to manage the mapping between circuit qubits and hardware qubits
 * in a bijective manner.
 */
class Mapping {
protected:
  // std::map<qc::Qubit, HwQubit>
  qc::Permutation circToHw;

public:
  Mapping() = default;
  Mapping(size_t nQubits, InitialMapping initialMapping) {
    switch (initialMapping) {
    case Identity:
      for (size_t i = 0; i < nQubits; ++i) {
        circToHw.emplace(i, i);
      }
      break;
    default:
      qc::unreachable();
    }
  }
  /**
   * @brief Assigns a circuit qubit to a hardware qubit.
   * @param qubit The circuit qubit to be assigned
   * @param hwQubit The hardware qubit to be assigned
   */
  void setCircuitQubit(qc::Qubit qubit, HwQubit hwQubit) {
    circToHw[qubit] = hwQubit;
  }

  /**
   * @brief Returns the hardware qubit assigned to the given circuit qubit.
   * @param qubit The circuit qubit to be queried
   * @return The hardware qubit assigned to the given circuit qubit
   */
  [[nodiscard]] HwQubit getHwQubit(qc::Qubit qubit) const {
    return circToHw.at(qubit);
  }

  /**
   * @brief Returns the hardware qubits assigned to the given circuit qubits.
   * @param qubits The circuit qubits to be queried
   * @return The hardware qubits assigned to the given circuit qubits
   */
  [[nodiscard]] std::set<HwQubit>
  getHwQubits(std::set<qc::Qubit>& qubits) const {
    std::set<HwQubit> hwQubits;
    for (const auto& qubit : qubits) {
      hwQubits.emplace(this->getHwQubit(qubit));
    }
    return hwQubits;
  }

  /**
   * @brief Returns the circuit qubit assigned to the given hardware qubit.
   * @details Throws an exception if the hardware qubit is not assigned to any
   * circuit qubit.
   * @param qubit The hardware qubit to be queried
   * @return The circuit qubit assigned to the given hardware qubit
   */
  [[nodiscard]] qc::Qubit getCircQubit(HwQubit qubit) const {
    for (const auto& [circQubit, hwQubit] : circToHw) {
      if (hwQubit == qubit) {
        return circQubit;
      }
    }
    throw std::runtime_error("Hardware qubit: " + std::to_string(qubit) +
                             " not found in mapping");
  }

  /**
   * @brief Indicates if any circuit qubit is assigned to the given hardware
   * qubit.
   * @param qubit The hardware qubit to be queried
   * @return True if any circuit qubit is assigned to the given hardware qubit,
   * false otherwise
   */
  [[nodiscard]] bool isMapped(HwQubit qubit) const {
    return std::any_of(
        circToHw.begin(), circToHw.end(),
        [qubit](const auto& pair) { return pair.second == qubit; });
  }

  /**
   * @brief Converts the qubits of an operation from circuit qubits to hardware
   * qubits.
   * @param op The operation to be converted
   */
  void mapToHwQubits(qc::Operation* op) const {
    op->setTargets(circToHw.apply(op->getTargets()));
    if (op->isControlled()) {
      op->setControls(circToHw.apply(op->getControls()));
    }
  }

  /**
   * @brief Interchanges the mapping of two hardware qubits. At least one of it
   * must be mapped to a circuit qubit.
   * @param swap The two circuit qubits to be swapped
   * @throws std::runtime_error if hardware qubits are not mapped
   */
  void applySwap(Swap swap);
};

} // namespace na
