"""Test the compilation of circuits."""
from __future__ import annotations

from pathlib import Path

import pytest
from qiskit import QuantumCircuit

from mqt import qmap
from mqt.qcec import verify


@pytest.fixture()
def example_circuit() -> QuantumCircuit:
    """Return a simple example circuit."""
    qc = QuantumCircuit(3)
    qc.h(0)
    qc.cx(0, 1)
    qc.cx(1, 2)
    qc.measure_all()
    return qc


def test_either_arch_or_calibration(example_circuit: QuantumCircuit) -> None:
    """Test that either arch or calibration must be provided."""
    with pytest.raises(ValueError, match="Either arch or calibration must be specified"):
        qmap.compile(example_circuit, arch=None, calibration=None)


@pytest.mark.parametrize(
    "arch",
    [
        "IBM_QX4",
        "IBM_QX5",
        "IBMQ_Yorktown",
        "IBMQ_London",
        "IBMQ_Bogota",
        "IBMQ_Tokyo",
        "Rigetti_Agave",
        "Rigetti_Aspen",
    ],
)
def test_available_architectures_str(example_circuit: QuantumCircuit, arch: str) -> None:
    """Test that the available architectures can be properly used."""
    example_circuit_mapped, results = qmap.compile(example_circuit, arch=arch)
    assert results.timeout is False
    assert results.mapped_circuit

    result = verify(example_circuit, example_circuit_mapped)
    assert result.considered_equivalent() is True


# test that all available architecture enumerations can be properly used
@pytest.mark.parametrize(
    "arch",
    [
        qmap.Arch.IBM_QX4,
        qmap.Arch.IBM_QX5,
        qmap.Arch.IBMQ_Yorktown,
        qmap.Arch.IBMQ_London,
        qmap.Arch.IBMQ_Bogota,
        qmap.Arch.IBMQ_Tokyo,
        qmap.Arch.Rigetti_Agave,
        qmap.Arch.Rigetti_Aspen,
    ],
)
def test_available_architectures_enum(example_circuit: QuantumCircuit, arch: qmap.Arch) -> None:
    """Test that the available architecture enums can be properly used."""
    example_circuit_mapped, results = qmap.compile(example_circuit, arch=arch)
    assert results.timeout is False
    assert results.mapped_circuit

    result = verify(example_circuit, example_circuit_mapped)
    assert result.considered_equivalent() is True


def test_architecture_from_file(example_circuit: QuantumCircuit) -> None:
    """Test that architectures from files can be properly used."""
    with Path("test_architecture.arch").open("w+") as f:
        f.write("3\n0 1\n0 2\n1 2\n")

    example_circuit_mapped, results = qmap.compile(example_circuit, arch="test_architecture.arch")
    assert results.timeout is False
    assert results.mapped_circuit

    result = verify(example_circuit, example_circuit_mapped)
    assert result.considered_equivalent() is True


def test_architecture_from_python(example_circuit: QuantumCircuit) -> None:
    """Test that architectures from python can be properly used."""
    arch = qmap.Architecture(3, {(0, 1), (0, 2), (1, 2)})
    example_circuit_mapped, results = qmap.compile(example_circuit, arch=arch)
    assert results.timeout is False
    assert results.mapped_circuit

    result = verify(example_circuit, example_circuit_mapped)
    assert result.considered_equivalent() is True


def test_calibration_from_file(example_circuit: QuantumCircuit) -> None:
    """Test that calibrations from files can be properly used."""
    with Path("test_calibration.cal").open("w+") as f:
        f.write("Header\n")
        f.write('Q0,0,0,0,1e-2,1e-4,"0_1: 1e-2, 0_2: 1e-2"\n')
        f.write('Q1,0,0,0,1e-2,1e-4,"1_2: 1e-2"\n')
        f.write("Q2,0,0,0,1e-2,1e-4, \n")

    example_circuit_mapped, results = qmap.compile(example_circuit, arch=None, calibration="test_calibration.cal")
    assert results.timeout is False
    assert results.mapped_circuit

    result = verify(example_circuit, example_circuit_mapped)
    assert result.considered_equivalent() is True
