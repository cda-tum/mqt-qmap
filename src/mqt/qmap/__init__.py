"""MQT QMAP library.

This file is part of the MQT QMAP library released under the MIT license.
See README.md or go to https://github.com/cda-tum/qmap for more information.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

if sys.platform == "win32" and sys.version_info > (3, 8, 0) and "Z3_ROOT" in os.environ:
    lib_path = Path(os.environ["Z3_ROOT"]) / "lib"
    if lib_path.exists():
        os.add_dll_directory(str(lib_path))
    bin_path = Path(os.environ["Z3_ROOT"]) / "bin"
    if bin_path.exists():
        os.add_dll_directory(str(bin_path))

from ._version import version as __version__
from .clifford_synthesis import optimize_clifford, synthesize_clifford
from .compile import compile
from .pyqmap import (
    Arch,
    Architecture,
    CliffordSynthesizer,
    CommanderGrouping,
    Configuration,
    Encoding,
    Heuristic,
    InitialLayout,
    Layering,
    LookaheadHeuristic,
    MappingResults,
    Method,
    QuantumComputation,
    SwapReduction,
    SynthesisConfiguration,
    SynthesisResults,
    Tableau,
    TargetMetric,
    Verbosity,
)
from .subarchitectures import SubarchitectureOrder

__all__ = [
    "Arch",
    "Architecture",
    "CliffordSynthesizer",
    "CommanderGrouping",
    "Configuration",
    "Encoding",
    "Heuristic",
    "InitialLayout",
    "Layering",
    "LookaheadHeuristic",
    "MappingResults",
    "Method",
    "QuantumComputation",
    "SubarchitectureOrder",
    "SwapReduction",
    "SynthesisConfiguration",
    "SynthesisResults",
    "Tableau",
    "TargetMetric",
    "Verbosity",
    "__version__",
    "compile",
    "optimize_clifford",
    "synthesize_clifford",
]
