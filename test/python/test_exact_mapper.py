from typing import Callable

import pytest
from mqt import qmap
from mqt.qcec import verify

from qiskit import QuantumCircuit
from qiskit.circuit import InstructionSet
from qiskit.providers.fake_provider import FakeLondon


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


@pytest.fixture
def one_way_arch() -> qmap.Architecture:
    arch = qmap.Architecture(
        2,
        {
            (0, 1),
        },
    )
    return arch


@pytest.mark.parametrize(
    "gate",
    [
        QuantumCircuit.cx,
        QuantumCircuit.csx,
    ],
)
def test_direction_reverse_hadamard(one_way_arch: qmap.Architecture, gate: Callable[..., InstructionSet]) -> None:
    """Verify that control and target are flipped using four hadamard gates for some gates where this is possible."""
    qc = QuantumCircuit(2)
    gate(qc, 0, 1)
    gate(qc, 1, 0)
    qc.measure_all()

    qc_mapped, results = qmap.compile(qc, arch=one_way_arch, method="exact")
    assert qc_mapped.count_ops()["h"] == 4
    assert "swap" not in qc_mapped.count_ops()


@pytest.mark.parametrize(
    "gate",
    [
        QuantumCircuit.cy,
        lambda qc, ctrl, target: QuantumCircuit.crx(qc, 0.5, ctrl, target),
        lambda qc, ctrl, target: QuantumCircuit.cry(qc, 0.5, ctrl, target),
    ],
)
def test_direction_reverse_swap(one_way_arch: qmap.Architecture, gate: Callable[..., InstructionSet]) -> None:
    """Verify that control and target are flipped using two swap gates for some gates where this is possible."""
    qc = QuantumCircuit(2)
    gate(qc, 0, 1)
    gate(qc, 1, 0)
    qc.measure_all()

    qc_mapped, results = qmap.compile(qc, arch=one_way_arch, method="exact")
    assert "h" not in qc_mapped.count_ops()
    assert qc_mapped.count_ops()["swap"] == 2


@pytest.mark.parametrize(
    "gate",
    [
        QuantumCircuit.cz,
        lambda qc, ctrl, target: QuantumCircuit.cp(qc, 0.5, ctrl, target),
        lambda qc, ctrl, target: QuantumCircuit.crz(qc, 0.5, ctrl, target),
    ],
)
def test_direction_reverse_identity(one_way_arch: qmap.Architecture, gate: Callable[..., InstructionSet]) -> None:
    """Verify that control and target are flipped without adding additional for some gates where this is possible."""
    qc = QuantumCircuit(2)
    gate(qc, 0, 1)
    gate(qc, 1, 0)
    qc.measure_all()

    qc_mapped, results = qmap.compile(qc, arch=one_way_arch, method="exact")
    assert "h" not in qc_mapped.count_ops()
    assert "swap" not in qc_mapped.count_ops()


@pytest.mark.parametrize(
    "gate",
    [
        QuantumCircuit.cx,
        QuantumCircuit.cy,
        QuantumCircuit.cz,
        QuantumCircuit.csx,
        QuantumCircuit.ch,
        QuantumCircuit.cnot,
        QuantumCircuit.cs,
        QuantumCircuit.csdg,
    ],
)
def test_direction_reverse_params_0(one_way_arch: qmap.Architecture, gate: Callable[..., InstructionSet]) -> None:
    """Verify that control and target is flipped correctly for all two qubit controlled gates with zero parameters."""
    qc = QuantumCircuit(2)
    gate(qc, 0, 1)
    gate(qc, 1, 0)
    qc.measure_all()
    qc_mapped, results = qmap.compile(qc, arch=one_way_arch, method="exact")

    result = verify(qc, qc_mapped)
    assert result.considered_equivalent() is True


@pytest.mark.parametrize(
    "gate",
    [
        QuantumCircuit.cp,
        QuantumCircuit.crx,
        QuantumCircuit.cry,
        QuantumCircuit.crz,
    ],
)
def test_direction_reverse_params_1(one_way_arch: qmap.Architecture, gate: Callable[..., InstructionSet]) -> None:
    """Verify that control and target is flipped correctly for all two qubit controlled gates with one parameter."""
    qc = QuantumCircuit(2)
    gate(qc, 0.25, 0, 1)
    gate(qc, 0.75, 1, 0)
    qc.measure_all()
    qc_mapped, results = qmap.compile(qc, arch=one_way_arch, method="exact")

    result = verify(qc, qc_mapped)
    assert result.considered_equivalent() is True


@pytest.mark.parametrize(
    "gate",
    [QuantumCircuit.cu],
)
def test_direction_reverse_params_4(one_way_arch: qmap.Architecture, gate: Callable[..., InstructionSet]) -> None:
    """Verify that control and target is flipped correctly for all two qubit controlled gates with four parameters."""
    qc = QuantumCircuit(2)
    gate(qc, 0.1, 0.2, 0.3, 0.4, 0, 1)
    gate(qc, 0.9, 0.8, 0.7, 0.6, 1, 0)
    qc.measure_all()
    qc_mapped, results = qmap.compile(qc, arch=one_way_arch, method="exact")

    result = verify(qc, qc_mapped)
    assert result.considered_equivalent() is True
