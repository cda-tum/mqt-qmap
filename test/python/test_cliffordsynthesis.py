from mqt import qmap

from qiskit import QuantumCircuit


def test_cliffordsynthesis_simple_binary_search() -> None:
    """Verify that the exact mapper works on a simple circuit that requires no swaps on a trivial initial layout."""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cx(0, 1)

    qc_mapped, results = qmap.synthesize_clifford(qc, use_binarysearch=True)

    assert results.mapped_circuit != ""


def test_cliffordsynthesis_simple_minimizer() -> None:
    """Verify that the exact mapper works on a simple circuit that requires no swaps on a trivial initial layout."""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cx(0, 1)

    qc_mapped, results = qmap.synthesize_clifford(qc, use_binarysearch=False)

    assert results.mapped_circuit != ""
