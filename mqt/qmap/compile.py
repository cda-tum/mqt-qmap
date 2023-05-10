"""Main entry point for the compilation module."""
from __future__ import annotations

from typing import TYPE_CHECKING

from qiskit import QuantumCircuit, QuantumRegister
from qiskit.transpiler import Layout

if TYPE_CHECKING:  # pragma: no cover
    from qiskit.providers import Backend
    from qiskit.providers.models import BackendProperties
    from qiskit.transpiler.target import Target

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
    map,
)

from .load_architecture import load_architecture
from .load_calibration import load_calibration


def extract_initial_layout_from_qasm(qasm: str, qregs: list[QuantumRegister]) -> Layout:
    """Extract the initial layout resulting from compiling a circuit from a QASM file.

    Args:
        qasm: The QASM file to extract the initial layout from.
        qregs: The quantum registers of the circuit.

    Returns:
        The initial layout.
    """
    for line in qasm.split("\n"):
        if line.startswith("// i "):
            # strip away initial part of line
            stripped_line = line[5:]
            # split line into tokens
            tokens = stripped_line.split(" ")
            # convert tokens to integers
            int_tokens = [int(token) for token in tokens]
            # create an empty layout
            layout = Layout().from_intlist(int_tokens, *qregs)
            return layout
    msg = "No initial layout found in QASM file."
    raise ValueError(msg)


def compile(  # noqa: A001
    circ: QuantumCircuit | str,
    arch: str | Arch | Architecture | Backend | None,
    calibration: str | BackendProperties | Target | None = None,
    method: str | Method = "heuristic",
    initial_layout: str | InitialLayout = "dynamic",
    layering: str | Layering = "individual_gates",
    use_teleportation: bool = False,
    teleportation_fake: bool = False,
    teleportation_seed: int = 0,
    encoding: str | Encoding = "commander",
    commander_grouping: str | CommanderGrouping = "fixed3",
    use_bdd: bool = False,
    swap_reduction: str | SwapReduction = "coupling_limit",
    swap_limit: int = 0,
    include_WCNF: bool = False,  # noqa: N803
    use_subsets: bool = True,
    subgraph: set[int] | None = None,
    pre_mapping_optimizations: bool = True,
    post_mapping_optimizations: bool = True,
    add_measurements_to_mapped_circuit: bool = True,
    verbose: bool = False,
    debug: bool = False,
) -> tuple[QuantumCircuit, MappingResults]:
    """Interface to the MQT QMAP tool for mapping quantum circuits.

    Args:
        circ: The circuit to map.
        arch: The architecture to map to.
        calibration: The calibration to use.
        method: The mapping method to use. Either "heuristic" or "exact". Defaults to "heuristic".
        initial_layout: The initial layout to use. Defaults to "dynamic".
        layering: The layering strategy to use. Defaults to "individual_gates".
        encoding: The encoding to use for the AMO and exactly one constraints. Defaults to "naive".
        commander_grouping: The grouping strategy to use for the commander and bimander encoding. Defaults to "halves".
        use_bdd: Whether to use BDDs to limit the search space. Defaults to False. Use with caution.
        swap_reduction: The swap reduction strategy to use. Defaults to "coupling_limit".
        swap_limit: Set a custom limit for max swaps per layer, for the increasing reduction strategy it sets the max swaps per layer. Defaults to 0.
        include_WCNF: Include WCNF file in the results. Defaults to False.
        use_subsets: Use qubit subsets, or consider all available physical qubits at once. Defaults to True.
        subgraph: List of qubits to consider for mapping (in exact mapper), if None all qubits are considered. Defaults to None.
        use_teleportation: Use teleportation in addition to swaps. Defaults to False.
        teleportation_fake: Assign qubits as ancillary for teleportation in the initial placement but don't actually use them (used for comparisons). Defaults to False.
        teleportation_seed: Fix a seed for the RNG in the initial ancilla placement (0 means the RNG will be seeded from /dev/urandom/ or similar). Defaults to 0.
        pre_mapping_optimizations: Run pre-mapping optimizations. Defaults to True.
        post_mapping_optimizations: Run post-mapping optimizations. Defaults to True.
        add_measurements_to_mapped_circuit: Whether to add measurements at the end of the mapped circuit. Defaults to True.
        verbose: Print more detailed information during the mapping process. Defaults to False.
        debug: Gather additional information during the mapping process (e.g. number of generated nodes, branching factors, ...). Defaults to False.

    Returns:
        The mapped circuit and the mapping results.
    """
    if subgraph is None:
        subgraph = set()

    if arch is None and calibration is None:
        msg = "Either arch or calibration must be specified"
        raise ValueError(msg)

    architecture = load_architecture(arch)
    load_calibration(architecture, calibration)

    config = Configuration()
    config.method = Method(method)
    config.initial_layout = InitialLayout(initial_layout)
    config.layering = Layering(layering)
    config.encoding = Encoding(encoding)
    config.commander_grouping = CommanderGrouping(commander_grouping)
    config.swap_reduction = SwapReduction(swap_reduction)
    config.swap_limit = swap_limit
    config.use_bdd = use_bdd
    config.include_WCNF = include_WCNF
    config.use_subsets = use_subsets
    config.subgraph = subgraph
    config.use_teleportation = use_teleportation
    config.teleportation_fake = teleportation_fake
    config.teleportation_seed = teleportation_seed
    config.pre_mapping_optimizations = pre_mapping_optimizations
    config.post_mapping_optimizations = post_mapping_optimizations
    config.add_measurements_to_mapped_circuit = add_measurements_to_mapped_circuit
    config.verbose = verbose
    config.debug = debug

    results = map(circ, architecture, config)

    circ = QuantumCircuit.from_qasm_str(results.mapped_circuit)
    layout = extract_initial_layout_from_qasm(results.mapped_circuit, circ.qregs)

    # qiskit-terra 0.22.0 introduced a breaking change in the `_layout` of the `QuantumCircuit` class.
    # To maintain backwards compatibility, the following `try... except` block is necessary.
    try:
        from qiskit.transpiler.layout import TranspileLayout

        circ._layout = TranspileLayout(initial_layout=layout, input_qubit_mapping=layout.get_virtual_bits())
    except ImportError:
        circ._layout = layout

    return circ, results
