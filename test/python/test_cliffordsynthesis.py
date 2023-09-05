"""Test cases for Clifford synthesis."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import pytest
from qiskit import QuantumCircuit
from qiskit.quantum_info import Clifford, PauliList

from mqt import qcec, qmap


@dataclass
class Configuration:
    """Configuration for a test case."""

    expected_minimal_gates: int
    expected_minimal_depth: int
    expected_minimal_two_qubit_depth: int
    expected_minimal_gates_at_minimal_depth: int
    expected_minimal_two_qubit_gates: int
    expected_minimal_gates_at_minimal_two_qubit_gates: int

    description: str
    initial_tableau: str | None = None
    target_tableau: str | None = None
    initial_circuit: str | None = None


def create_circuit_tests() -> list[Configuration]:
    """Create test cases for Clifford synthesis."""
    path = Path(__file__).resolve().parent.parent / "cliffordsynthesis" / "circuits.json"
    with path.open() as f:
        circuits = json.load(f)
    return [Configuration(**c) for c in circuits]


def create_tableau_tests() -> list[Configuration]:
    """Create test cases for tableau synthesis."""
    path = Path(__file__).resolve().parent.parent / "cliffordsynthesis" / "tableaus.json"
    with path.open() as f:
        tableaus = json.load(f)
    return [Configuration(**t) for t in tableaus]


@pytest.mark.parametrize("test_config", create_tableau_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "linear_search"])
def test_synthesize_clifford_tqdepth(test_config: Configuration, use_maxsat: bool) -> None:
    """Test TQDepth-optimal tableau synthesis."""
    circ, results = qmap.synthesize_clifford(
        target_tableau=test_config.target_tableau,
        initial_tableau=test_config.initial_tableau,
        use_maxsat=False,
        target_metric="TQDepth",
    )

    assert results.depth == test_config.expected_minimal_two_qubit_depth
    print("\n", circ)
