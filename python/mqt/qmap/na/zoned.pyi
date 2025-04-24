# Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
# Copyright (c) 2025 Munich Quantum Software Company GmbH
# All rights reserved.
#
# SPDX-License-Identifier: MIT
#
# Licensed under the MIT License

from mqt.core.ir import QuantumComputation

class ZonedNeutralAtomArchitecture:
    @classmethod
    def from_json_file(cls, filename: str) -> ZonedNeutralAtomArchitecture: ...
    @classmethod
    def from_json_string(cls, json: str) -> ZonedNeutralAtomArchitecture: ...

class RoutingAgnosticCompiler:
    def __init__(
        self,
        arch: ZonedNeutralAtomArchitecture,
        log_level: str = ...,
        use_window: bool = ...,
        window_size: int = ...,
        dynamic_placement: bool = ...,
        parking_offset: int = ...,
        warn_unsupported_gates: bool = ...,
    ) -> None: ...
    @classmethod
    def from_json_string(cls, arch: ZonedNeutralAtomArchitecture, json: str) -> RoutingAgnosticCompiler: ...
    def compile(self, qc: QuantumComputation) -> str: ...
    def stats(self) -> dict[str, float]: ...

class RoutingAwareCompiler:
    def __init__(
        self,
        arch: ZonedNeutralAtomArchitecture,
        log_level: str = ...,
        use_window: bool = ...,
        window_min_width: int = ...,
        window_ratio: float = ...,
        window_share: float = ...,
        deepening_factor: float = ...,
        deepening_value: float = ...,
        lookahead_factor: float = ...,
        reuse_level: float = ...,
        max_nodes: int = ...,
        parking_offset: int = ...,
        warn_unsupported_gates: bool = ...,
    ) -> None: ...
    @classmethod
    def from_json_string(cls, arch: ZonedNeutralAtomArchitecture, json: str) -> RoutingAwareCompiler: ...
    def compile(self, qc: QuantumComputation) -> str: ...
    def stats(self) -> dict[str, float]: ...
