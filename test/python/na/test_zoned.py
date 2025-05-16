# Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
# Copyright (c) 2025 Munich Quantum Software Company GmbH
# All rights reserved.
#
# SPDX-License-Identifier: MIT
#
# Licensed under the MIT License

"""Test MQT QMAP's Zoned Neutral Atom Compiler."""

from __future__ import annotations

from pathlib import Path

import pytest
from mqt.core import load
from mqt.qmap.na.zoned import RoutingAwareCompiler, ZonedNeutralAtomArchitecture

# get the circuit directory of the project
circ_dir = Path(__file__).resolve().parent.parent.parent / "na" / "zoned" / "circuits"
# make list of contained .qasm files
circuits = list(circ_dir.glob("*.qasm"))

architecture_specification = """{
    "name": "compiler_architecture",
    "storage_zones": [{
        "zone_id": 0,
        "slms": [{"id": 0, "site_separation": [3, 3], "r": 20, "c": 20, "location": [0, 0]}],
        "offset": [0, 0],
        "dimension": [60, 60]
    }],
    "entanglement_zones": [{
        "zone_id": 0,
        "slms": [
            {"id": 1, "site_separation": [12, 10], "r": 4, "c": 4, "location": [5, 70]},
            {"id": 2, "site_separation": [12, 10], "r": 4, "c": 4, "location": [7, 70]}
        ],
        "offset": [5, 70],
        "dimension": [50, 40]
    }],
    "aods":[{"id": 0, "site_separation": 2, "r": 20, "c": 20}],
    "rydberg_range": [[[5, 70], [55, 110]]]
}"""
routing_aware_configuration = """{
    "placerConfig" : {
        "useWindow" : true,
        "windowMinWidth" : 4,
        "windowRatio" : 1.5,
        "windowShare" : 0.6,
        "deepeningFactor" : 0.6,
        "deepeningValue" : 0.2,
        "lookaheadFactor": 0.2,
        "reuseLevel": 5.0
    },
    "codeGeneratorConfig" : {
        "parkingOffset" : 1
    }
}"""


@pytest.fixture
def compiler() -> RoutingAwareCompiler:
    """Return an MQT QMAP's Zoned Neutral Atom Compiler initialized with the above architecture and settings."""
    architecture = ZonedNeutralAtomArchitecture.from_json_string(architecture_specification)
    configuration = RoutingAwareCompiler.Config.from_json_string(routing_aware_configuration)
    return RoutingAwareCompiler(architecture, configuration)


@pytest.mark.parametrize("circuit_filename", circuits)
def test_na_routing_aware_compiler(compiler: RoutingAwareCompiler, circuit_filename: str) -> None:
    """Test the MQT QMAP's Zoned Neutral Atom Compiler."""
    qc = load(circuit_filename)
    result = compiler.compile(qc)
    assert result is not None
    stats = compiler.stats()
    assert "totalTime" in stats
    assert stats["totalTime"] > 0
