# Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
# Copyright (c) 2025 Munich Quantum Software Company GmbH
# All rights reserved.
#
# SPDX-License-Identifier: MIT
#
# Licensed under the MIT License

"""Python bindings module for MQT QMAP's Neutral Atom State Preparation."""

from typing import Any

from mqt.core.ir import QuantumComputation

class NAStatePreparationSolver:
    """The main class for the Neutral Atom State Preparation."""
    def __init__(
        self,
        new_max_x: int,
        new_max_y: int,
        new_max_c: int,
        new_max_r: int,
        new_max_h_offset: int,
        new_max_v_offset: int,
        new_max_h_dist: int,
        new_max_vdist: int,
        new_min_entangling_y: int,
        new_max_entangling_y: int,
    ) -> None: ...

    class Result:
        def __init__(self) -> None: ...
        def json(self) -> dict[str, Any]: ...

    def solve(
        self,
        ops: list[tuple[int, int]],
        new_num_qubits: int,
        new_num_stages: int,
        new_num_transfers: int | None = ...,
        min_ops_order: bool = ...,
        shield_idle_qubits: bool = ...,
    ) -> Result: ...

def get_ops_for_solver(
    circ: QuantumComputation,
    operation_type: str,
    num_controls: int,
    quiet: bool = ...,
) -> list[tuple[int, int]]: ...
def generate_code(
    circ: QuantumComputation,
    result: NAStatePreparationSolver.Result,
    min_atom_dist: int = ...,
    no_interaction_radius: int = ...,
    zone_dist: int = ...,
) -> str: ...
