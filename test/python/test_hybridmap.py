"""Test the hybrid Neutral Atom mapper."""

from __future__ import annotations

from pathlib import Path

import pytest
from qiskit import QuantumCircuit

from mqt.qmap import HybridMapperParameters, HybridNAMapper, NeutralAtomHybridArchitecture

arch_dir = Path(__file__).parent.parent / "hybridmap" / "architectures"
circuit_dir = Path(__file__).parent.parent / "hybridmap" / "circuits"


@pytest.mark.parametrize(
    "circuit_filename",
    [
        "dj_nativegates_rigetti_qiskit_opt3_10.qasm",
        "modulo_2.qasm",
        "multiply_2.qasm",
        "qft_nativegates_rigetti_qiskit_opt3_10.qasm",
        "random_nativegates_rigetti_qiskit_opt3_10.qasm",
    ],
)
@pytest.mark.parametrize(
    "arch_filename",
    [
        "rubidium.json",
        "rubidium_hybrid.json",
        "rubidium_shuttling.json",
    ],
)
@pytest.mark.parametrize(
    ("lookahead_weight", "decay", "gate_shuttling_weight"), [(0.0, 0.1, 0.1), (0.0, 0.0, 0.1), (0.0, 1.0, 10)]
)
def test_hybrid_na_mapper(
    circuit_filename: str, arch_filename: str, lookahead_weight: float, decay: float, gate_shuttling_weight: float
) -> None:
    """Test the hybrid Neutral Atom mapper."""
    arch = NeutralAtomHybridArchitecture(str(arch_dir / arch_filename))

    params = HybridMapperParameters(
        lookahead_weight_moves=lookahead_weight,
        lookahead_weight_swaps=lookahead_weight,
        decay=decay,
        gate_weight=gate_shuttling_weight,
    )
    mapper = HybridNAMapper(arch, params=params)

    # Create a simple circuit
    qc = QuantumCircuit.from_qasm_file(str(circuit_dir / circuit_filename))

    mapper.map(qc)
    results = mapper.schedule(create_animation_csv=False)

    assert results["totalExecutionTime"] > 0
    assert results["totalIdleTime"] > 0
    assert results["totalFidelities"] > 0
