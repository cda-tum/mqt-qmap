"""Test the state preparation for zoned neutral atom architectures."""

from mqt.qmap.na import Solver

import pytest
from pathlib import Path

# get root directory of the project
circ_dir = Path(__file__).resolve() / "circuits"

@pytest.fixture
def solver() -> Solver:
@pytest.mark.parametrize(
    "circuit_filename",
    [
        "shor.qasm", "steane.qasm", "surface_3.qasm",
    ],
)
def test_na_state_prep(solver: Solver, circuit_filename: str) -> None:
    """Test the state preparation for the zoned neutral atom architecture."""
    solver.solve(circ_dir / circuit_filename)
