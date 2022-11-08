from mqt import qmap

from qiskit import QuantumCircuit


def test_cliffordsynthesis_sanity_check() -> None:
    """Verify that the clifford synthesis works on a simple circuit that requires no actual optimization"""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.h(0)
    qc.h(0)
    qc.h(0)
    qc.h(0)

    qc_mapped, results = qmap.optimize_clifford(qc, strategy="minmax")

    assert results.gate_count == 1


def test_cliffordsynthesis_sanity_check_minimizer() -> None:
    """Verify that the clifford synthesis works on a simple circuit that requires no actual optimization"""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.h(0)
    qc.h(0)
    qc.h(0)
    qc.h(0)

    qc_mapped, results = qmap.optimize_clifford(qc, strategy="use_minimizer")

    assert results.result_circuit != ""
