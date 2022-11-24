import pytest
from mqt import qmap

from qiskit import QuantumCircuit
from qiskit.providers.fake_provider import FakeLondon


@pytest.fixture
def example_2_qubit_circuit() -> QuantumCircuit:
    qc = QuantumCircuit(2)
    qc.cx(0, 1)
    qc.cx(0, 1)
    qc.cx(0, 1)
    return qc


@pytest.fixture
def example_h_qubit_circuit() -> QuantumCircuit:
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.h(0)
    qc.h(0)
    qc.h(0)
    qc.h(0)
    return qc


def test_cliffordsynthesis_sanity_check_binary_search(example_h_qubit_circuit) -> None:
    """Sanity check Synthesis for 1 gate"""

    qc_mapped, results = qmap.optimize_clifford(example_h_qubit_circuit, strategy="binary_search")

    assert results.single_qubit_gates == 1


def test_cliffordsynthesis_sanity_check_minimizer(example_h_qubit_circuit) -> None:
    """Sanity check maxsat minimizer"""

    qc_mapped, results = qmap.optimize_clifford(example_h_qubit_circuit, strategy="use_minimizer")

    assert results.single_qubit_gates == 1


def test_cliffordsynthesis_sanity_check_depth_binary_search(example_2_qubit_circuit) -> None:
    """Verify that the clifford synthesis works on a simple circuit that requires no actual optimization"""

    qc_mapped, results = qmap.optimize_clifford(example_2_qubit_circuit, target="depth")

    assert results.depth == 1


def test_cliffordsynthesis_sanity_check_depth_minimizer(example_2_qubit_circuit) -> None:
    """Verify that the clifford synthesis works on a simple circuit that requires no actual optimization"""

    qc_mapped, results = qmap.optimize_clifford(example_2_qubit_circuit, target="depth", strategy="use_minimizer")

    assert results.depth == 1


def test_cliffordsynthesis_sanity_check_fidelity(example_2_qubit_circuit) -> None:
    """Verify that fidelity is not 0 if none is given"""

    qc_mapped, results = qmap.optimize_clifford(example_2_qubit_circuit, arch=FakeLondon(), target="fidelity")

    assert results.fidelity != 0.0


def test_cliffordsynthesis_sanity_check_empty_tableau() -> None:
    """Verify that fidelity is 1 if none is given"""

    stabilizers = '["+ZI", "+IZ"]'

    qc_mapped, results = qmap.synthesize_clifford(stabilizers, target="gates")

    assert results.fidelity == 1.0


@pytest.mark.parametrize("strategy", ["use_minimizer", "binary_search", "start_low", "start_high", "split_iter"])
@pytest.mark.parametrize("target", ["gates", "depth"])
def test_cliffordsynthesis_sanity_check_all(example_2_qubit_circuit, strategy, target) -> None:
    """Sanity check Synthesis for 1 gate"""

    qc_mapped, results = qmap.optimize_clifford(example_2_qubit_circuit, strategy=strategy, target=target)

    assert results.sat == qmap.SatSolverResult.sat
