"""Test the exact mapper."""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
import qiskit.circuit.library as qcl
from qiskit import QuantumCircuit
from qiskit.providers.fake_provider import FakeLondon

if TYPE_CHECKING:
    from qiskit.circuit import ControlledGate

from mqt import qmap
from mqt.qcec import verify


def test_exact_no_swaps_trivial_layout() -> None:
    """Verify that the exact mapper works on a simple circuit that requires no swaps on a trivial initial layout."""
    qc = QuantumCircuit(3)
    qc.h(0)
    qc.cx(0, 1)
    qc.cx(1, 2)
    qc.measure_all()

    qc_mapped, results = qmap.compile(qc, arch=FakeLondon(), method="exact")
    assert results.timeout is False
    assert results.mapped_circuit != ""
    assert results.output.swaps == 0

    result = verify(qc, qc_mapped)
    assert result.considered_equivalent() is True


def test_exact_no_swaps_non_trivial_layout() -> None:
    """Verify that the exact mapper works on a simple circuit that requires a non-trivial layout to achieve no swaps."""
    qc = QuantumCircuit(4)
    qc.h(0)
    qc.cx(0, 1)
    qc.cx(0, 2)
    qc.cx(0, 3)
    qc.measure_all()

    qc_mapped, results = qmap.compile(qc, arch=FakeLondon(), method="exact")

    assert results.timeout is False
    assert results.mapped_circuit != ""
    assert results.output.swaps == 0

    result = verify(qc, qc_mapped)
    assert result.considered_equivalent() is True


def test_exact_non_trivial_swaps() -> None:
    """Verify that the exact mapper works on a simple circuit that requires at least a single SWAP."""
    qc = QuantumCircuit(3)
    qc.h(0)
    qc.cx(0, 1)
    qc.cx(1, 2)
    qc.cx(2, 0)
    qc.measure_all()

    qc_mapped, results = qmap.compile(qc, arch=FakeLondon(), method="exact")

    assert results.timeout is False
    assert results.mapped_circuit != ""
    assert results.output.swaps == 1

    result = verify(qc, qc_mapped)
    assert result.considered_equivalent() is True


@pytest.fixture()
def one_way_arch() -> qmap.Architecture:
    """Return a simple one way architecture to test direction reversal."""
    return qmap.Architecture(
        2,
        {
            (0, 1),
        },
    )


@pytest.mark.parametrize(
    "gate",
    [qcl.CXGate(), qcl.CSXGate()],
)
def test_direction_reverse_hadamard(one_way_arch: qmap.Architecture, gate: ControlledGate) -> None:
    """Verify that control and target are flipped using four hadamard gates for some gates where this is possible."""
    qc = QuantumCircuit(2)
    qc.append(gate, [0, 1])
    qc.append(gate, [1, 0])
    qc.measure_all()

    qc_mapped, results = qmap.compile(qc, arch=one_way_arch, method="exact")
    assert qc_mapped.count_ops()["h"] == 4
    assert "swap" not in qc_mapped.count_ops()

    result = verify(qc, qc_mapped)
    assert result.considered_equivalent() is True


@pytest.mark.parametrize(
    "gate",
    [
        qcl.CYGate(),
        qcl.CRXGate(0.5),
        qcl.CRYGate(0.5),
        qcl.CHGate(),
        qcl.CU3Gate(0.25, 0.5, 0.75),
    ],
)
def test_direction_reverse_swap(one_way_arch: qmap.Architecture, gate: ControlledGate) -> None:
    """Verify that control and target are flipped using a swap permutation."""
    qc = QuantumCircuit(2)
    qc.append(gate, [0, 1])
    qc.append(gate, [1, 0])
    qc.measure_all()

    qc_mapped, results = qmap.compile(qc, arch=one_way_arch, method="exact", swap_limit=10, swap_reduction="custom")
    assert "h" not in qc_mapped.count_ops()
    assert qc_mapped.count_ops()["swap"] == 1

    result = verify(qc, qc_mapped)
    assert result.considered_equivalent() is True


@pytest.mark.parametrize(
    "gate",
    [
        qcl.CZGate(),
        qcl.CPhaseGate(0.5),
        qcl.CRZGate(0.5),
        qcl.CPhaseGate(np.pi / 2),
        qcl.CPhaseGate(-np.pi / 2),
        qcl.CPhaseGate(np.pi / 4),
        qcl.CPhaseGate(-np.pi / 4),
        qcl.CU1Gate(0.25),
    ],
)
def test_direction_reverse_identity(one_way_arch: qmap.Architecture, gate: ControlledGate) -> None:
    """Verify that control and target are flipped without adding gates for some gates where this is possible."""
    qc = QuantumCircuit(2)
    qc.append(gate, [0, 1])
    qc.append(gate, [1, 0])
    qc.measure_all()

    qc_mapped, results = qmap.compile(qc, arch=one_way_arch, method="exact")
    assert "h" not in qc_mapped.count_ops()
    assert "swap" not in qc_mapped.count_ops()

    result = verify(qc, qc_mapped)
    assert result.considered_equivalent() is True
