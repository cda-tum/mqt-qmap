"""Test the advanced zoned atom compiler."""

from __future__ import annotations

from json import loads
from pathlib import Path

import pytest

from mqt.core import load
from mqt.qmap.na.azac import AZACArchitecture, AZACompiler

# get circuit directory of the project
circ_dir = Path(__file__).resolve().parent.parent.parent / "na/azac/circuits"
# make list of contained .qasm files
circuits = list(circ_dir.glob("*.qasm"))

settings = """{
    "architecture": {
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
        "arch_range": [[0, 0], [60, 110]],
        "rydberg_range": [[[5, 70], [55, 110]]]
    },
    "vm_placer" : {
        "use_window" : true,
        "window_size" : 10,
        "dynamic_placement" : true
    },
    "code_generator" : {
        "parking_offset" : 1
    },
    "a_star_placer" : {
        "use_window" : true,
        "window_min_width" : 4,
        "window_ratio" : 1.5,
        "window_share" : 0.6,
        "deepening_factor" : 0.6,
        "deepening_value" : 0.2,
        "lookahead_factor": 0.2,
        "reuse_level": 5.0
    }
}"""


@pytest.fixture
def compiler() -> AZACompiler:
    """Return an advanced zoned atom compiler initialized with the above architecture and settings."""
    # get dict from json string settings
    settings_dict = loads(settings)
    architecture = AZACArchitecture(settings_dict["architecture"])
    return AZACompiler(architecture, settings_dict)


@pytest.mark.parametrize("circuit_filename", circuits)
def test_na_azacompiler(compiler: AZACompiler, circuit_filename: str) -> None:
    """Test the advanced zoned atom compiler."""
    qc = load(circuit_filename)
    result = compiler.compile(qc)
    assert result is not None
