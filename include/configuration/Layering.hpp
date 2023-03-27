//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include <iostream>

enum class Layering {
  None,
  IndividualGates,
  DisjointQubits,
  OddGates,
  QubitTriangle,
  Disjoint2qBlocks
};

[[maybe_unused]] static inline std::string toString(const Layering strategy) {
  switch (strategy) {
  case Layering::IndividualGates:
    return "individual_gates";
  case Layering::DisjointQubits:
    return "disjoint_qubits";
  case Layering::OddGates:
    return "odd_gates";
  case Layering::QubitTriangle:
    return "qubit_triangle";
  case Layering::Disjoint2qBlocks:
    return "disjoint_2q_blocks";
  case Layering::None:
    return "none";
  }
  return " ";
}

[[maybe_unused]] static Layering
layeringFromString(const std::string& layering) {
  if (layering == "none" || layering == "0") {
    return Layering::None;
  }
  if (layering == "individual_gates" || layering == "1") {
    return Layering::IndividualGates;
  }
  if (layering == "disjoint_qubits" || layering == "2") {
    return Layering::DisjointQubits;
  }
  if (layering == "odd_gates" || layering == "3") {
    return Layering::OddGates;
  }
  if (layering == "qubit_triangle" || layering == "4") {
    return Layering::QubitTriangle;
  }
  if (layering == "disjoint_2q_blocks" || layering == "5") {
    return Layering::Disjoint2qBlocks;
  }
  throw std::invalid_argument("Invalid layering value: " + layering);
}
