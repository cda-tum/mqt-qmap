# Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
# Copyright (c) 2025 Munich Quantum Software Company GmbH
# All rights reserved.
#
# SPDX-License-Identifier: MIT
#
# Licensed under the MIT License

"""Python bindings module for MQT QMAP's Zoned Neutral Atom Compiler."""

from typing import Any, overload

from mqt.core.ir import QuantumComputation

class ZonedNeutralAtomArchitecture:
    """Class representing a Zoned Neutral Atom Architecture."""

    name: str

    def __init__(self) -> None: ...
    @classmethod
    def from_json_file(cls, filename: str) -> ZonedNeutralAtomArchitecture: ...
    @classmethod
    def from_json_string(cls, json: str) -> ZonedNeutralAtomArchitecture: ...

class RoutingAgnosticCompiler:
    """Class representing the MQT QMAP's Zoned Neutral Atom Compiler."""

    @overload
    def __init__(self, arch: ZonedNeutralAtomArchitecture) -> None: ...
    @overload
    def __init__(self, arch: ZonedNeutralAtomArchitecture, settings: Any) -> None: ...
    def compile(self, circ: QuantumComputation) -> str: ...
    def stats(self) -> dict[str, float]: ...

class RoutingAwareCompiler:
    """Class representing the MQT QMAP's Zoned Neutral Atom Compiler."""

    @overload
    def __init__(self, arch: ZonedNeutralAtomArchitecture) -> None: ...
    @overload
    def __init__(self, arch: ZonedNeutralAtomArchitecture, settings: Any) -> None: ...
    def compile(self, circ: QuantumComputation) -> str: ...
    def stats(self) -> dict[str, float]: ...
