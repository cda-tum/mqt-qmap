"""Test the hybrid Neutral Atom mapper."""

from __future__ import annotations

import pytest

from mqt.qmap import HybridMapperParameters, HybridNAMapper, NeutralAtomHybridArchitecture

arch_file = "../hybridmap/architectures/rubidium_hybrid.json"


@pytest.mark.parametrize(
    ("lookahead_weight", "decay", "gate_shuttling_weight"), [(0.0, 0.1, 0.1), (0.0, 0.0, 0.1), (0.0, 1.0, 10)]
)
@pytest.mark.parametrize(
    "circuit_file",
    [
        "../hybridmap/circuits/dj_nativegates_rigetti_qiskit_opt3_10.qasm",
        "../hybridmap/circuits/modulo_2.qasm",
    ],
)
def test_hybrid_na_mapper(circuit_file, lookahead_weight, decay, gate_shuttling_weight):
    """Test the hybrid Neutral Atom mapper."""
    arch = NeutralAtomHybridArchitecture(arch_file)
    params = HybridMapperParameters(
        lookahead_weight_moves=lookahead_weight,
        lookahead_weight_swaps=lookahead_weight,
        decay=decay,
        gate_weight=gate_shuttling_weight,
    )
    mapper = HybridNAMapper(arch, params=params)
    mapper.map_qasm_file(circuit_file)
    mapper.schedule(create_animation_csv=False)
