"""Test the exact mapper."""
import pytest
import qiskit
from qiskit import QuantumCircuit
from qiskit.providers.fake_provider import FakeLondon

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
    [qiskit.circuit.library.CXGate(), qiskit.circuit.library.CSXGate()],
)
def test_direction_reverse_hadamard(one_way_arch: qmap.Architecture, gate: qiskit.circuit.ControlledGate) -> None:
    """Verify that control and target are flipped using four hadamard gates for some gates where this is possible."""
    qc = QuantumCircuit(2)
    qc.append(gate, [0, 1])
    qc.append(gate, [1, 0])
    qc.measure_all()

    qc_mapped, results = qmap.compile(qc, arch=one_way_arch, method="exact")
    assert qc_mapped.count_ops()["h"] == 4
    assert "swap" not in qc_mapped.count_ops()


@pytest.mark.parametrize(
    "gate",
    [qiskit.circuit.library.CYGate(), qiskit.circuit.library.CRXGate(0.5), qiskit.circuit.library.CRYGate(0.5)],
)
def test_direction_reverse_swap(one_way_arch: qmap.Architecture, gate: qiskit.circuit.ControlledGate) -> None:
    """Verify that control and target are flipped using two swap gates for some gates where this is possible."""
    qc = QuantumCircuit(2)
    qc.append(gate, [0, 1])
    qc.append(gate, [1, 0])
    qc.measure_all()

    qc_mapped, results = qmap.compile(qc, arch=one_way_arch, method="exact")
    assert "h" not in qc_mapped.count_ops()
    assert qc_mapped.count_ops()["swap"] == 2


@pytest.mark.parametrize(
    "gate",
    [qiskit.circuit.library.CZGate(), qiskit.circuit.library.CPhaseGate(0.5), qiskit.circuit.library.CRZGate(0.5)],
)
def test_direction_reverse_identity(one_way_arch: qmap.Architecture, gate: qiskit.circuit.ControlledGate) -> None:
    """Verify that control and target are flipped without adding additional for some gates where this is possible."""
    qc = QuantumCircuit(2)
    qc.append(gate, [0, 1])
    qc.append(gate, [1, 0])
    qc.measure_all()

    qc_mapped, results = qmap.compile(qc, arch=one_way_arch, method="exact")
    assert "h" not in qc_mapped.count_ops()
    assert "swap" not in qc_mapped.count_ops()


@pytest.mark.parametrize(
    "gate",
    [
        qiskit.circuit.library.CHGate(),
        qiskit.circuit.library.CPhaseGate(0.25),
        qiskit.circuit.library.CSdgGate(),
        qiskit.circuit.library.CRXGate(0.25),
        qiskit.circuit.library.CRYGate(0.25),
        qiskit.circuit.library.CRZGate(0.25),
        qiskit.circuit.library.CSGate(),
        qiskit.circuit.library.CSXGate(),
        qiskit.circuit.library.CU1Gate(0.25),
        qiskit.circuit.library.CU3Gate(0.25, 0.5, 0.75),
        qiskit.circuit.library.CUGate(0.2, 0.4, 0.6, 0.8),
        qiskit.circuit.library.CXGate(),
        qiskit.circuit.library.CYGate(),
        qiskit.circuit.library.CZGate(),
    ],
)
def test_direction_reverse_verify(one_way_arch: qmap.Architecture, gate: qiskit.circuit.ControlledGate) -> None:
    """Verify that control and target is flipped correctly for two qubit controlled gates."""
    qc = QuantumCircuit(2)
    qc.append(gate, [0, 1])
    qc.append(gate, [1, 0])
    qc.measure_all()
    qc_mapped, results = qmap.compile(qc, arch=one_way_arch, method="exact")

    result = verify(qc, qc_mapped)
    assert result.considered_equivalent() is True
