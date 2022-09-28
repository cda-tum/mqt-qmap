from mqt import qmap

from qiskit import QuantumCircuit


def test_cliffordsynthesis_simple_binary_search() -> None:
    """Verify that the clifford synthesis works on a simple circuit that requires no actual optimization"""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cx(0, 1)

    qc_mapped, results = qmap.optimize_clifford(qc, strategy="minmax")

    assert results.mapped_circuit != ""


def test_cliffordsynthesis_simple_minimizer() -> None:
    """Verify that the clifford synthesis works on a simple circuit that requires no actual optimization"""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cx(0, 1)

    qc_mapped, results = qmap.optimize_clifford(qc, strategy="use_minimizer")

    assert results.mapped_circuit != ""
