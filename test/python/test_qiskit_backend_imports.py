"""Test the Qiskit backend imports."""

from __future__ import annotations

import pytest
from qiskit import QuantumCircuit
from qiskit.providers.fake_provider import Fake5QV1, GenericBackendV2

from mqt import qmap


@pytest.fixture()
def example_circuit() -> QuantumCircuit:
    """Return a simple circuit."""
    qc = QuantumCircuit(3)
    qc.h(0)
    qc.cx(0, 1)
    qc.cx(1, 2)
    qc.measure_all()
    return qc


@pytest.fixture()
def backend() -> GenericBackendV2:
    """Return a test backend."""
    return GenericBackendV2(num_qubits=5, coupling_map=[[0, 1], [1, 0], [1, 2], [2, 1], [1, 3], [3, 1], [3, 4], [4, 3]])


def test_backend_v1(example_circuit: QuantumCircuit) -> None:
    """Test that circuits can be mapped to Qiskit BackendV1 instances providing the new basis_gates."""
    _, results = qmap.compile(example_circuit, arch=Fake5QV1())
    assert results.timeout is False
    assert results.mapped_circuit


def test_architecture_from_v1_backend_properties(example_circuit: QuantumCircuit) -> None:
    """Test that circuits can be mapped by simply providing the backend properties (the BackendV1 way)."""
    _, results = qmap.compile(example_circuit, arch=None, calibration=Fake5QV1().properties())
    assert results.timeout is False
    assert results.mapped_circuit


def test_backend_v2(example_circuit: QuantumCircuit, backend: GenericBackendV2) -> None:
    """Test that circuits can be mapped to Qiskit BackendV1 instances providing the old basis_gates."""
    _, results = qmap.compile(example_circuit, arch=backend)
    assert results.timeout is False
    assert results.mapped_circuit


def test_architecture_from_v2_target(example_circuit: QuantumCircuit, backend: GenericBackendV2) -> None:
    """Test that circuits can be mapped by simply providing the target (the BackendV2 way)."""
    _, results = qmap.compile(example_circuit, arch=None, calibration=backend.target)
    assert results.timeout is False
    assert results.mapped_circuit
