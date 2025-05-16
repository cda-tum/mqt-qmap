# Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
# Copyright (c) 2025 Munich Quantum Software Company GmbH
# All rights reserved.
#
# SPDX-License-Identifier: MIT
#
# Licensed under the MIT License

"""Python bindings module for MQT QMAP's Zoned Neutral Atom Compiler."""

from typing import overload

from mqt.core.ir import QuantumComputation

class ZonedNeutralAtomArchitecture:
    """Class representing a Zoned Neutral Atom Architecture."""

    @classmethod
    def from_json_file(cls, filename: str) -> ZonedNeutralAtomArchitecture:
        """Create an architecture from a JSON file.

        Args:
            filename: is the path to the JSON file

        Returns:
            the architecture

        Raises:
            ValueError: if the file does not exist or is not a valid JSON file
        """

    @classmethod
    def from_json_string(cls, json: str) -> ZonedNeutralAtomArchitecture:
        """Create an architecture from a JSON string.

        Args:
            json: is the JSON string

        Returns:
            the architecture

        Raises:
            ValueError: if the string is not a valid JSON
        """

class VMPlacerConfig:
    """Class representing the configuration of the placer."""

    use_window: bool
    """Whether to use a window for the placer."""
    window_size: int
    """The size of the window for the placer."""
    dynamic_placement: bool
    """Whether to use dynamic placement for the placer."""
    def __init__(self) -> None:
        """Initialize the placer config with default settings."""

class CodeGeneratorConfig:
    """Class representing the configuration of the code generator."""

    parking_offset: int
    """The parking offset of the code generator."""
    warn_unsupported_gates: bool
    """Whether to warn about unsupported gates in the code generator."""
    def __init__(self) -> None:
        """Initialize the code generator config with default settings."""

class AStarPlacerConfig:
    """Class representing the configuration of the placer."""

    use_window: bool
    """Whether to use a window for the placer."""
    window_min_width: int
    """The minimum width of the window for the placer."""
    window_ratio: float
    """The ratio between the height and the width of the window."""
    window_share: float
    """The the share of free sites in the window in relation to the number of atoms to be moved in this step."""
    deepening_factor: float
    """The heuristic used in the A* search contains a term that resembles
    the standard deviation of the differences between the current and target
    sites of the atoms to be moved in every orientation.
    """
    deepening_value: float
    """Before the sum of standard deviations is multiplied with the
    number of unplaced nodes and :attr:`deepening_factor`, this value is added
    to the sum to amplify the influence of the unplaced nodes count.
    """
    lookahead_factor: float
    """The cost function can consider the distance of atoms to their
    interaction partner in the next layer.
    """
    reuse_level: float
    """The reuse level corresponds to the estimated extra fidelity loss
    due to the extra trap transfers when the atom is not reused and instead
    moved to the storage zone and back to the entanglement zone.
    """
    max_nodes: int
    def __init__(self) -> None:
        """Initialize the placer config with default settings."""

class RoutingAgnosticCompiler:
    """MQT QMAP's routing-agnostic Zoned Neutral Atom Compiler."""

    class Config:
        """Class representing the configuration of the compiler."""

        placer_config: VMPlacerConfig
        """The configuration of the placer."""
        code_generator_config: CodeGeneratorConfig
        """The configuration of the code generator."""

        def __init__(self) -> None:
            """Initialize the compiler config with default settings."""

        @classmethod
        def from_json_string(cls, json: str) -> RoutingAgnosticCompiler.Config:
            """Create a configuration from a JSON string.

            Args:
                json: is the JSON string

            Returns:
                the configuration

            Raises:
                ValueError: if the string is not a valid JSON
            """

    @overload
    def __init__(self, arch: ZonedNeutralAtomArchitecture) -> None:
        """Create a routing-agnostic compiler for the given architecture and settings.

        Args:
            arch: is the zoned neutral atom architecture
        """

    @overload
    def __init__(self, arch: ZonedNeutralAtomArchitecture, config: Config) -> None:
        """Create a routing-agnostic compiler for the given architecture and settings.

        Args:
            arch: is the zoned neutral atom architecture
            config: is a dictionary with the settings for the compiler
        """

    def compile(self, qc: QuantumComputation) -> str:
        """Compile a quantum circuit for the zoned neutral atom architecture.

        Args:
            qc: is the quantum circuit

        Returns:
            the compilations result as a string in the .naviz format.
        """

    def stats(self) -> dict[str, float]:
        """Get the statistics of the last compilation.

        Returns:
            the statistics as a dictionary
        """

class RoutingAwareCompiler:
    """MQT QMAP's routing-aware Zoned Neutral Atom Compiler."""

    class Config:
        """Class representing the configuration of the compiler."""

        placer_config: AStarPlacerConfig
        """The configuration of the placer."""
        code_generator_config: CodeGeneratorConfig
        """The configuration of the code generator."""

        def __init__(self) -> None:
            """Initialize the compiler config with default settings."""

        @classmethod
        def from_json_string(cls, json: str) -> RoutingAwareCompiler.Config:
            """Create a configuration from a JSON string.

            Args:
                json: is the JSON string

            Returns:
                the configuration

            Raises:
                ValueError: if the string is not a valid JSON
            """

    @overload
    def __init__(self, arch: ZonedNeutralAtomArchitecture) -> None:
        """Create a routing-aware compiler for the given architecture and settings.

        Args:
            arch: is the zoned neutral atom architecture
        """

    @overload
    def __init__(self, arch: ZonedNeutralAtomArchitecture, config: Config) -> None:
        """Create a routing-aware compiler for the given architecture and settings.

        Args:
            arch: is the zoned neutral atom architecture
            config: is a dictionary with the settings for the compiler
        """

    def compile(self, qc: QuantumComputation) -> str:
        """Compile a quantum circuit for the zoned neutral atom architecture.

        Args:
            qc: is the quantum circuit

        Returns:
            the compilations result as a string in the .naviz format.
        """

    def stats(self) -> dict[str, float]:
        """Get the statistics of the last compilation.

        Returns:
            the statistics as a dictionary
        """
