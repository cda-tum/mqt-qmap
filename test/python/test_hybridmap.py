"""Test the hybrid Neutral Atom mapper."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from qiskit import QuantumCircuit

from mqt.qmap import HybridMapperParameters, HybridNAMapper, NeutralAtomHybridArchitecture

json_data = {
    "name": "rubidium",
    "properties": {
        "nRows": 5,
        "nColumns": 6,
        "nAods": 1,
        "nAodCoordinates": 5,
        "interQubitDistance": 3,
        "minimalAodDistance": 0.1,
        "interactionRadius": 3,
        "blockingFactor": 1,
    },
    "parameters": {
        "nQubits": 20,
        "gateTimes": {
            "none": 0.5,
            "rx": 0.5,
            "rz": 0.5,
            "h": 0.5,
            "x": 0.5,
            "p": 0.001,
            "cx": 0.7,
            "cz": 0.2,
            "ccz": 0.4,
            "cccz": 0.6,
            "ccccz": 0.8,
            "mcx": 2.5,
        },
        "gateAverageFidelities": {
            "none": 0.999,
            "h": 0.999,
            "x": 0.999,
            "rx": 0.999,
            "rz": 0.999,
            "p": 0.999,
            "cz": 0.995,
            "ccz": 0.95,
            "cccz": 0.95,
            "ccccz": 0.95,
            "cccccz": 0.95,
            "mcx": 0.95,
        },
        "decoherenceTimes": {"t1": 100000000, "t2": 1500000},
        "shuttlingTimes": {"move": 0.55, "aod_move": 0.55, "aod_activate": 20, "aod_deactivate": 20},
        "shuttlingAverageFidelities": {"move": 1, "aod_move": 1, "aod_activate": 1, "aod_deactivate": 1},
    },
}


@pytest.mark.parametrize(
    ("lookahead_weight", "decay", "gate_shuttling_weight"), [(0.0, 0.1, 0.1), (0.0, 0.0, 0.1), (0.0, 1.0, 10)]
)
def test_hybrid_na_mapper(lookahead_weight, decay, gate_shuttling_weight):
    """Test the hybrid Neutral Atom mapper."""
    # Create a temporary file
    temp_file = Path("arch.json")
    temp_file.write_text(json.dumps(json_data))

    arch = NeutralAtomHybridArchitecture(str(temp_file))

    # Delete the temporary file
    temp_file.unlink()

    params = HybridMapperParameters(
        lookahead_weight_moves=lookahead_weight,
        lookahead_weight_swaps=lookahead_weight,
        decay=decay,
        gate_weight=gate_shuttling_weight,
    )
    mapper = HybridNAMapper(arch, params=params)

    # Create a simple circuit
    qc = QuantumCircuit(6)
    qc.h(0)
    qc.cx(0, 1)
    qc.cx(0, 2)
    qc.cx(0, 3)
    qc.cx(0, 4)
    qc.cx(0, 5)

    mapper.map(qc)
    results = mapper.schedule(create_animation_csv=False)

    assert results["totalExecutionTime"] > 0
    assert results["totalIdleTime"] > 0
    assert results["totalFidelities"] > 0
