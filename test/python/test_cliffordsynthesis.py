from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import pytest
from mqt import qmap


@dataclass
class Configuration:
    """Configuration for a test case."""

    expected_minimal_gates: int
    expected_minimal_depth: int
    expected_minimal_gates_at_minimal_depth: int
    expected_minimal_two_qubit_gates: int
    expected_minimal_gates_at_minimal_two_qubit_gates: int

    description: str
    initial_tableau: str | None = None
    target_tableau: str | None = None
    initial_circuit: str | None = None


def create_circuit_tests() -> list[Configuration]:
    path = Path(__file__).resolve().parent.parent / "cliffordsynthesis" / "circuits.json"
    with path.open() as f:
        circuits = json.load(f)
    return [Configuration(**c) for c in circuits]


def create_tableau_tests() -> list[Configuration]:
    path = Path(__file__).resolve().parent.parent / "cliffordsynthesis" / "tableaus.json"
    with path.open() as f:
        tableaus = json.load(f)
    return [Configuration(**t) for t in tableaus]


@pytest.mark.parametrize("test_config", create_circuit_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "binary_search"])
def test_optimize_clifford_gates(test_config: Configuration, use_maxsat: bool) -> None:
    circ, results = qmap.optimize_clifford(
        circuit=test_config.initial_circuit, use_maxsat=use_maxsat, target_metric="gates"
    )

    assert results.gates == test_config.expected_minimal_gates
    print("\n", circ)


@pytest.mark.parametrize("test_config", create_circuit_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "binary_search"])
def test_optimize_clifford_depth(test_config: Configuration, use_maxsat: bool) -> None:
    circ, results = qmap.optimize_clifford(
        circuit=test_config.initial_circuit, use_maxsat=use_maxsat, target_metric="depth"
    )

    assert results.depth == test_config.expected_minimal_depth
    print("\n", circ)


@pytest.mark.parametrize("test_config", create_circuit_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "binary_search"])
def test_optimize_clifford_gates_at_minimal_depth(test_config: Configuration, use_maxsat: bool) -> None:
    circ, results = qmap.optimize_clifford(
        circuit=test_config.initial_circuit,
        use_maxsat=use_maxsat,
        target_metric="depth",
        minimize_gates_after_depth_optimization=True,
    )

    assert results.gates == test_config.expected_minimal_gates_at_minimal_depth
    print("\n", circ)


@pytest.mark.parametrize("test_config", create_circuit_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "binary_search"])
def test_optimize_clifford_two_qubit_gates(test_config: Configuration, use_maxsat: bool) -> None:
    circ, results = qmap.optimize_clifford(
        circuit=test_config.initial_circuit,
        use_maxsat=use_maxsat,
        target_metric="two_qubit_gates",
        try_higher_gate_limit_for_two_qubit_gate_optimization=True,
    )

    assert results.two_qubit_gates == test_config.expected_minimal_two_qubit_gates
    print("\n", circ)


@pytest.mark.parametrize("test_config", create_circuit_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "binary_search"])
def test_optimize_clifford_gates_at_minimal_two_qubit_gates(test_config: Configuration, use_maxsat: bool) -> None:
    circ, results = qmap.optimize_clifford(
        circuit=test_config.initial_circuit,
        use_maxsat=use_maxsat,
        target_metric="two_qubit_gates",
        try_higher_gate_limit_for_two_qubit_gate_optimization=True,
        minimize_gates_after_two_qubit_gate_optimization=True,
    )

    assert results.gates == test_config.expected_minimal_gates_at_minimal_two_qubit_gates
    print("\n", circ)


@pytest.mark.parametrize("test_config", create_tableau_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "binary_search"])
def test_synthesize_clifford_gates(test_config: Configuration, use_maxsat: bool) -> None:
    circ, results = qmap.synthesize_clifford(
        target_tableau=test_config.target_tableau,
        initial_tableau=test_config.initial_tableau,
        use_maxsat=use_maxsat,
        target_metric="gates",
    )

    assert results.gates == test_config.expected_minimal_gates
    print("\n", circ)


@pytest.mark.parametrize("test_config", create_tableau_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "binary_search"])
def test_synthesize_clifford_depth(test_config: Configuration, use_maxsat: bool) -> None:
    circ, results = qmap.synthesize_clifford(
        target_tableau=test_config.target_tableau,
        initial_tableau=test_config.initial_tableau,
        use_maxsat=use_maxsat,
        target_metric="depth",
    )

    assert results.depth == test_config.expected_minimal_depth
    print("\n", circ)


@pytest.mark.parametrize("test_config", create_tableau_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "binary_search"])
def test_synthesize_clifford_gates_at_minimal_depth(test_config: Configuration, use_maxsat: bool) -> None:
    circ, results = qmap.synthesize_clifford(
        target_tableau=test_config.target_tableau,
        initial_tableau=test_config.initial_tableau,
        use_maxsat=use_maxsat,
        target_metric="depth",
        minimize_gates_after_depth_optimization=True,
    )

    assert results.gates == test_config.expected_minimal_gates_at_minimal_depth
    print("\n", circ)


@pytest.mark.parametrize("test_config", create_tableau_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "binary_search"])
def test_synthesize_clifford_two_qubit_gates(test_config: Configuration, use_maxsat: bool) -> None:
    circ, results = qmap.synthesize_clifford(
        target_tableau=test_config.target_tableau,
        initial_tableau=test_config.initial_tableau,
        use_maxsat=use_maxsat,
        target_metric="two_qubit_gates",
        try_higher_gate_limit_for_two_qubit_gate_optimization=True,
    )

    assert results.two_qubit_gates == test_config.expected_minimal_two_qubit_gates
    print("\n", circ)


@pytest.mark.parametrize("test_config", create_tableau_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "binary_search"])
def test_synthesize_clifford_gates_at_minimal_two_qubit_gates(test_config: Configuration, use_maxsat: bool) -> None:
    circ, results = qmap.synthesize_clifford(
        target_tableau=test_config.target_tableau,
        initial_tableau=test_config.initial_tableau,
        use_maxsat=use_maxsat,
        target_metric="two_qubit_gates",
        try_higher_gate_limit_for_two_qubit_gate_optimization=True,
        minimize_gates_after_two_qubit_gate_optimization=True,
    )

    assert results.gates == test_config.expected_minimal_gates_at_minimal_two_qubit_gates
    print("\n", circ)
