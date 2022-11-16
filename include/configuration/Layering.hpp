/*
 * This file is part of MQT QMAP library which is released under the MIT
 * license. See file README.md or go to http://iic.jku.at/eda/research/quantum/
 * for more information.
 */

#ifndef QMAP_LAYERING_HPP
#define QMAP_LAYERING_HPP

#include <iostream>

enum class Layering {
  None,
  IndividualGates,
  DisjointQubits,
  OddGates,
  QubitTriangle
};

static std::string toString(const Layering strategy) {
  switch (strategy) {
  case Layering::IndividualGates:
    return "individual_gates";
  case Layering::DisjointQubits:
    return "disjoint_qubits";
  case Layering::OddGates:
    return "odd_gates";
  case Layering::QubitTriangle:
    return "qubit_triangle";
  case Layering::None:
    return "none";
  }
  return " ";
}

[[maybe_unused]] static Layering
layeringFromString(const std::string& layering) {
  if (layering == "none" || layering == "0") {
    return Layering::None;
  } else if (layering == "individual_gates" || layering == "1") {
    return Layering::IndividualGates;
  } else if (layering == "disjoint_qubits" || layering == "2") {
    return Layering::DisjointQubits;
  } else if (layering == "odd_gates" || layering == "3") {
    return Layering::OddGates;
  } else if (layering == "qubit_triangle" || layering == "4") {
    return Layering::QubitTriangle;
  } else {
    throw std::invalid_argument("Invalid layering value: " + layering);
  }
}

#endif // QMAP_LAYERING_HPP
