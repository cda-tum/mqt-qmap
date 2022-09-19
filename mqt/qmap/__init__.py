#
# This file is part of MQT QMAP library which is released under the MIT license.
# See file README.md or go to http://iic.jku.at/eda/research/quantum_verification/ for more information.
#

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

from mqt.qmap.compile import compile, synthesize_clifford, optimize_clifford
from mqt.qmap.pyqmap import (
    Arch,
    Architecture,
    CliffordOptResults,
    CommanderGrouping,
    Configuration,
    Encoding,
    InitialLayout,
    Layering,
    MappingResults,
    Method,
    OptimizationTarget,
    OptimizingStrategy,
    SwapReduction,
)

__all__ = [
    "compile",
    "Method",
    "InitialLayout",
    "Layering",
    "Arch",
    "CommanderGrouping",
    "SwapReduction",
    "Encoding",
    "Configuration",
    "MappingResults",
    "Architecture",
    "OptimizationTarget",
    "OptimizingStrategy",
    "CliffordOptResults",
    "synthesize_clifford",
    "optimize_clifford",
]
