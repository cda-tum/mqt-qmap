"""Test the state preparation for zoned neutral atom architectures."""

from __future__ import annotations

from pathlib import Path

import pytest

from mqt.qmap.na import NAStatePreparationSolver, get_ops_for_solver

# get root directory of the project
circ_dir = Path(__file__).resolve() / "../../na/nasp/circuits"


@pytest.fixture
def solver() -> NAStatePreparationSolver:
    """Return a NAStatePreparationSolver instance with both sided storage zone."""
    return NAStatePreparationSolver()  # 3, 7, 2, 3, 2, 2, 2, 2, 2, 4)


@pytest.mark.parametrize(
    "circuit_filename",
    [
        "shor.qasm",
        "steane.qasm",
        "surface_3.qasm",
    ],
)
def test_na_state_prep(solver: NAStatePreparationSolver, circuit_filename: str) -> None:
    """Test the state preparation for the zoned neutral atom architecture."""
    ops = get_ops_for_solver(circ_dir / circuit_filename)
    solver.solve(ops, 7, 4, None, False, True)
    # todo check result
