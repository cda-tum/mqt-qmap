import pytest
from qiskit import QuantumCircuit
from qiskit.providers.fake_provider import FakeLondon, FakeLondonV2, FakeAthens, FakeAthensV2

from mqt import qmap


@pytest.fixture
def example_circuit() -> QuantumCircuit:
    qc = QuantumCircuit(3)
    qc.h(0)
    qc.cx(0, 1)
    qc.cx(1, 2)
    qc.measure_all()
    return qc


def test_old_backend_v1(example_circuit: QuantumCircuit):
    """Test that circuits can be mapped to Qiskit BackendV1 instances providing the old basis_gates."""
    _, results = qmap.compile(example_circuit, arch=FakeLondon())
    assert results.timeout is False
    assert results.mapped_circuit != ""


def test_old_backend_v2(example_circuit: QuantumCircuit):
    """Test that circuits can be mapped to Qiskit BackendV2 instances providing the old basis_gates."""
    _, results = qmap.compile(example_circuit, arch=FakeLondonV2())
    assert results.timeout is False
    assert results.mapped_circuit != ""


def test_new_backend_v1(example_circuit: QuantumCircuit):
    """Test that circuits can be mapped to Qiskit BackendV1 instances providing the new basis_gates."""
    _, results = qmap.compile(example_circuit, arch=FakeAthens())
    assert results.timeout is False
    assert results.mapped_circuit != ""


def test_new_backend_v2(example_circuit: QuantumCircuit):
    """Test that circuits can be mapped to Qiskit BackendV2 instances providing the new basis_gates."""
    _, results = qmap.compile(example_circuit, arch=FakeAthensV2())
    assert results.timeout is False
    assert results.mapped_circuit != ""


def test_architecture_from_v1_backend_properties(example_circuit: QuantumCircuit):
    """Test that circuits can be mapped by simply providing the backend properties (the BackendV1 way)."""
    _, results = qmap.compile(example_circuit, arch=None, calibration=FakeLondon().properties())
    assert results.timeout is False
    assert results.mapped_circuit != ""


def test_architecture_from_v2_target(example_circuit: QuantumCircuit):
    """Test that circuits can be mapped by simply providing the target (the BackendV2 way)."""
    _, results = qmap.compile(example_circuit, arch=None, calibration=FakeLondonV2().target)
    assert results.timeout is False
    assert results.mapped_circuit != ""
