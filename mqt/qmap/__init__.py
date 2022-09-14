#
# This file is part of MQT QMAP library which is released under the MIT license.
# See file README.md or go to http://iic.jku.at/eda/research/quantum_verification/ for more information.
#

import os
import sys

if sys.platform == "win32" and sys.version_info > (3, 8, 0) and "Z3_ROOT" in os.environ:
    lib_path = os.path.join(os.environ["Z3_ROOT"], "lib")
    if os.path.exists(lib_path):
        os.add_dll_directory(lib_path)
    bin_path = os.path.join(os.environ["Z3_ROOT"], "bin")
    if os.path.exists(bin_path):
        os.add_dll_directory(bin_path)

from mqt.qmap.compile import compile
from mqt.qmap.pyqmap import (
    Arch,
    Architecture,
    CommanderGrouping,
    Configuration,
    Encoding,
    InitialLayout,
    Layering,
    MappingResults,
    Method,
    SwapReduction,
)
from mqt.qmap.subarchitectures import SubarchitectureOrder

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
    "SubarchitectureOrder",
]
