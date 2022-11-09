from mqt import qmap

from qiskit import QuantumCircuit
from qiskit.providers.fake_provider import FakeLondon


def test_cliffordsynthesis_sanity_check() -> None:
    """Verify that the clifford synthesis works on a simple circuit that requires no actual optimization"""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.h(0)
    qc.h(0)
    qc.h(0)
    qc.h(0)

    qc_mapped, results = qmap.optimize_clifford(qc, strategy="minmax")

    assert results.single_qubit_gates == 1


def test_cliffordsynthesis_sanity_check_minimizer() -> None:
    """Verify that the clifford synthesis works on a simple circuit that requires no actual optimization"""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.h(0)
    qc.h(0)
    qc.h(0)
    qc.h(0)

    qc_mapped, results = qmap.optimize_clifford(qc, strategy="use_minimizer")

    assert results.single_qubit_gates == 1


def test_cliffordsynthesis_sanity_check_depth() -> None:
    """Verify that the clifford synthesis works on a simple circuit that requires no actual optimization"""
    qc = QuantumCircuit(2)
    qc.cx(0, 1)
    qc.cx(0, 1)
    qc.cx(0, 1)

    qc_mapped, results = qmap.optimize_clifford(qc, target="depth")

    assert results.depth == 1


def test_cliffordsynthesis_sanity_check_depth_minmax() -> None:
    """Verify that the clifford synthesis works on a simple circuit that requires no actual optimization"""
    qc = QuantumCircuit(2)
    qc.cx(0, 1)
    qc.cx(0, 1)
    qc.cx(0, 1)

    qc_mapped, results = qmap.optimize_clifford(qc, target="depth", strategy="minmax")

    assert results.depth == 1


def test_cliffordsynthesis_sanity_check_fidelity() -> None:
    """Verify that fidelity is 0 if none is given"""

    qc = QuantumCircuit(2)
    qc.cx(0, 1)
    qc.cx(0, 1)
    qc.cx(0, 1)

    qc_mapped, results = qmap.optimize_clifford(qc, arch=FakeLondon(), target="fidelity")

    assert results.fidelity != 0.0
