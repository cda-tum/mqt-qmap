/*
 * Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
 * Copyright (c) 2025 Munich Quantum Software Company GmbH
 * All rights reserved.
 *
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License
 */

#include "na/nasp/CodeGenerator.hpp"
#include "na/nasp/Solver.hpp"

#include <cstdint>
#include <nlohmann/json.hpp>
#include <pybind11/cast.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(MQT_QMAP_MODULE_NAME, m, py::mod_gil_not_used()) {
  m.doc() = "Bindings for mqt.qmap.na.state_preparation";

  // Neutral Atom State Preparation
  py::class_<na::NASolver>(m, "NAStatePreparationSolver", R"(
The neutral atom state preparation solver generates an optimal sequence of
neutral atom operations for a given state preparation circuit.
)")
      .def(py::init<uint16_t, uint16_t, uint16_t, uint16_t, uint16_t, uint16_t,
                    uint16_t, uint16_t, uint16_t, uint16_t>(),
           "max_x"_a, "max_y"_a, "max_c"_a, "max_r"_a, "max_h_offset"_a,
           "max_v_offset"_a, "max_h_dist"_a, "max_v_dist"_a,
           "min_entangling_y"_a, "max_entangling_y"_a, R"(
Create a solver instance for the neutral atom state preparation problem.

The solver is based on a 2D grid abstraction of the neutral atom quantum
computer. The 2D plane is divided into interaction sites. Each interaction site
is denoted by abstract x- and y-coordinates. The parameter `max_x` specifies the
maximum x-coordinate, and `max_y` specifies the maximum y-coordinate. In the
center of an interaction site, sits an SLM trap. Around that trap there are
several possible discrete AOD positions arranged as a grid. The specific
position of an atom within an interaction site is determined by the x- and
y-offset from the SLM trap, i.e., those can be positive and negative. The
maximum absolute value of the x- and y-offset is specified by `max_h_offset` and
`max_v_offset`, respectively. Then, the parameter `max_c` specifies the maximum
number of AOD columns, and `max_r` specifies the maximum number of AOD rows.
Finally, in order to interact during a Rydberg stage, atoms must be located
within a certain distance. The maximum horizontal and vertical distance between
two atoms is specified by `max_h_dist` and `max_v_dist`, respectively. The
parameter `min_entangling_y` specifies the minimum y-coordinate for entangling
operations, and `max_entangling_y` specifies the maximum y-coordinate for
entangling operations. Hence, y-coordinates outside of this range are located in
the storage zone.

.. note::
    The solver can only handle a single storage zone below the entangling zone,
    i.e., in this case `min_entangling_y` must be zero and `max_entangling_y`
    must be less than `max_y`.

:param max_x: is the maximum discrete x-coordinate of the interaction sites
:param max_y: is the maximum discrete y-coordinate of the interaction sites
:param max_c: is the maximum number of AOD columns
:param max_r: is the maximum number of AOD rows
:param max_h_offset: is the maximum horizontal offset of the atoms
:param max_v_offset: is the maximum vertical offset of the atoms
:param max_h_dist: is the maximum horizontal distance between two atoms
:param max_v_dist: is the maximum vertical distance between two atoms
:param min_entangling_y: is the minimum y-coordinate for entangling operations
:param max_entangling_y: is the maximum y-coordinate for entangling operations
:raises ValueError: if one of the parameters is invalid, e.g., is a negative
value
)")
      .def("solve", &na::NASolver::solve, "ops"_a, "num_qubits"_a,
           "num_stages"_a, "num_transfers"_a, "mind_ops_order"_a,
           "shield_idle_qubits"_a, R"(
Solve the neutral atom state preparation problem.

The solver generates an optimal sequence of neutral atom operations for a given
state preparation circuit. The circuit is given as a list of operations, where
each operation is a pair of qubits. The sequence is divided into stages. Each
stage is either a Rydberg stage or a transfer stage. In a Rydberg stage,
adjacent qubits in the entangling zone undergo an entangling gate. In a transfer
stage, atoms can be stored from AOD into SLM traps and loaded from SLM traps
into AOD. At the end of each stage, the atoms are shuttled to their next
position. The number of stages is specified by `num_stages`. The number of
transfers is fixed by `num_transfers` if give. If this parameter is not
specified, then the solver will determine the optimal number of transfers. The
parameter `mind_ops_order` specifies whether the order of the operations in the
circuit should be preserved. The parameter `shield_idle_qubits` specifies
whether idle qubits should be shielded from the entangling operations.

.. note::
    To retrieve the list of qubit pairs from a quantum circuit, use the function
    :func:`get_ops_for_solver`.

.. note::
    The returned solver's result can either directly exported to the JSON format
    by calling the method :func:`json` on the result object or the result object
    can be passed to the function :func:`generate_code` to generate code
    consisting of neutral atom operations.

:param ops: is the list of operations in the circuit
:param num_qubits: is the number of qubits in the circuit
:param num_stages: is the number of stages in the sequence
:param num_transfers: (optional) is the number of transfers in the sequence
:param mind_ops_order: is True if the order of the operations should be
  preserved
:param shield_idle_qubits: is True if idle qubits should be shielded
:returns: the result of the solver
:raises ValueError: if one of the numeral parameters is invalid, e.g., is a
  negative value
)");

  py::class_<na::NASolver::Result>(
      m, "NAStatePreparationSolver.Result",
      "Neutral Atom State Preparation Solver Result")
      .def(py::init<>(), "Create a result object")
      .def(
          "json",
          [](const na::NASolver::Result& result) { return result.json(); }, R"(
Returns the result as a JSON string.

:returns: the result as a JSON string
)");

  m.def(
      "generate_code",
      [](const qc::QuantumComputation& qc, const na::NASolver::Result& result,
         const uint16_t minAtomDist, const uint16_t noInteractionRadius,
         const uint16_t zoneDist) {
        return na::CodeGenerator::generate(qc, result, minAtomDist,
                                           noInteractionRadius, zoneDist)
            .toString();
      },
      "qc"_a, "result"_a, "min_atom_dist"_a = 1, "no_interaction_radius"_a = 10,
      "zone_dist"_a = 24, R"(
Generate code for the given circuit using the solver's result. Some parameters
of the abstraction from the 2D grid used for the solver must be provided again.

:param qc: is the quantum circuit
:param result: is the result of the solver
:param min_atom_dist: is the minimum distance between atoms
:param no_interaction_radius: is the radius around an atom where no other atom
  can be placed during an entangling operation that should not interact with the
  atom
:param zone_dist: is the distance between zones, i.e., the minimal distance
  between two atoms in different zones
:raises ValueError: if one of the numeral parameters is invalid, e.g., is a
  negative value)");

  m.def(
      "get_ops_for_solver",
      [](const qc::QuantumComputation& qc, const std::string& operationType,
         const uint64_t numControls, const bool quiet) {
        auto opTypeLowerStr = operationType;
        std::transform(opTypeLowerStr.begin(), opTypeLowerStr.end(),
                       opTypeLowerStr.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        return na::NASolver::getOpsForSolver(
            qc, qc::opTypeFromString(operationType), numControls, quiet);
      },
      "qc"_a, "operation_type"_a = "Z", "num_operands"_a = 1, "quiet"_a = true,
      R"(
Extract entangling operations as list of qubit pairs from the circuit.

.. warning::
    This function can only extract qubit pairs of two-qubit operations.
    I.e., the operands of the operation plus the controls must be equal to two.

:param qc: is the quantum circuit
:param operation_type: is the type of operation to extract, e.g., "z" for CZ
  gates
:param num_controls: is the number of controls the operation acts on, e.g., 1
  for CZ gates
:param quiet: if True, suppresses warning when the circuit contains operations
  other than the specified operation type
:returns: list of qubit pairs
:raises ValueError: if the circuit contains operations other than the specified
  operation type and quiet is False
:raises ValueError: if the operation has more than two operands including
  controls
)");
}
