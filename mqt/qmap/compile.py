#
# This file is part of MQT QMAP library which is released under the MIT license.
# See file README.md or go to http://iic.jku.at/eda/research/quantum_verification/ for more information.
#
from __future__ import annotations

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

from qiskit import QuantumCircuit, QuantumRegister
from qiskit.providers import Backend
from qiskit.providers.models import BackendProperties
from qiskit.transpiler import Layout
from qiskit.transpiler.target import Target


def extract_initial_layout_from_qasm(qasm: str, qregs: list[QuantumRegister]) -> Layout:
    """
    Extracts the initial layout resulting from compiling a circuit from a QASM file.
    :param qasm: QASM file
    :type qasm: str
    :param qregs: The quantum registers to apply the layout to.
    :type qregs: list[QuantumRegister]
    :return: layout to be used in Qiskit
    """
    for line in qasm.split("\n"):
        if line.startswith("// i "):
            # strip away initial part of line
            line = line[5:]
            # split line into tokens
            tokens = line.split(" ")
            # convert tokens to integers
            int_tokens = [int(token) for token in tokens]
            # create an empty layout
            layout = Layout().from_intlist(int_tokens, *qregs)
            return layout
    raise ValueError("No initial layout found in QASM file.")


def compile(
    circ: QuantumCircuit | str,
    arch: str | Arch | Architecture | Backend | None,
    calibration: str | BackendProperties | Target | None = None,
    method: str | Method = "heuristic",
    initial_layout: str | InitialLayout = "dynamic",
    layering: str | Layering = "individual_gates",
    use_teleportation: bool = False,
    teleportation_fake: bool = False,
    teleportation_seed: int = 0,
    encoding: str | Encoding = "naive",
    commander_grouping: str | CommanderGrouping = "halves",
    use_bdd: bool = False,
    swap_reduction: str | SwapReduction = "coupling_limit",
    swap_limit: int = 0,
    include_WCNF: bool = False,
    use_subsets: bool = True,
    subgraph: set[int] | None = None,
    pre_mapping_optimizations: bool = True,
    post_mapping_optimizations: bool = True,
    add_measurements_to_mapped_circuit: bool = True,
    verbose: bool = False,
) -> tuple[QuantumCircuit, MappingResults]:
    """Interface to the MQT QMAP tool for mapping quantum circuits

    :param circ: Qiskit QuantumCircuit object or path to circuit file
    :type circ: QuantumCircuit | str
    :param arch: Architecture to map to. Either a path to a file with architecture information, one of the available architectures (:py:mod:`mqt.qmap.Arch`), Architecture, or `qiskit.providers.backend` (if Qiskit is installed)
    :type arch: str | Arch | Architecture | Backend | None
    :param calibration: Path to file containing calibration information, `qiskit.providers.models.BackendProperties` object (if Qiskit is installed), or `qiskit.transpiler.target.Target` object (if Qiskit is installed)
    :type calibration: str | BackendProperties | Target | None
    :param method: Mapping technique to use (*heuristic* | exact)
    :type method: str | Method
    :param initial_layout: Strategy to use for determining initial layout in heuristic mapper (identity | static | *dynamic*)
    :type initial_layout: str | InitialLayout
    :param layering: Circuit layering strategy to use (*individual_gates* | disjoint_qubits | odd_qubits | qubit_triangle)
    :type layering: str | Layering
    :param encoding: Choose encoding for AMO and exactly one constraints (*naive* | commander | bimander)
    :type encoding: str | Encoding
    :param commander_grouping: Choose method of grouping in commander and bimander encoding (*halves* | fixed2 | fixed3 | logarithm)
    :type commander_grouping: str | CommanderGrouping
    :param use_bdd: Limit swaps per layer using BDDs, faster in some cases, but use with caution (default: False)
    :type use_bdd: bool
    :param swap_reduction: Choose method of limiting the search space (none | *coupling_limit* | custom | increasing)
    :type swap_reduction: str | SwapReduction
    :param swap_limit: Set a custom limit for max swaps per layer, for the increasing reduction strategy it sets the max swaps per layer
    :type swap_limit: int
    :param include_WCNF: Include WCNF file in the results (default: False)
    :type include_WCNF: bool
    :param use_subsets: Use qubit subsets, or consider all available physical qubits at once (default: True)
    :type use_subsets: bool
    :param subgraph: List of qubits to consider for mapping (in exact mapper), if None all qubits are considered
    :type subgraph: set[int] | None
    :param use_teleportation:  Use teleportation in addition to swaps
    :param teleportation_fake: Assign qubits as ancillary for teleportation in the initial placement but don't actually use them (used for comparisons)
    :param teleportation_seed: Fix a seed for the RNG in the initial ancilla placement (0 means the RNG will be seeded from /dev/urandom/ or similar)
    :param pre_mapping_optimizations: Run pre-mapping optimizations (default: True)
    :type pre_mapping_optimizations: bool
    :param post_mapping_optimizations: Run post-mapping optimizations (default: True)
    :type post_mapping_optimizations: bool
    :param add_measurements_to_mapped_circuit: Whether to add measurements at the end of the mapped circuit (default: True)
    :type add_measurements_to_mapped_circuit: bool
    :param verbose: Print more detailed information during the mapping process
    :type verbose: bool

    :return: Mapped circuit (as Qiskit `QuantumCircuit`) and results
    :rtype: tuple[QuantumCircuit, MappingResults]
    """

    if subgraph is None:
        subgraph = set()

    architecture = Architecture()
    if arch is None and calibration is None:
        raise ValueError("Either arch or calibration must be specified")

    if arch is not None:
        if isinstance(arch, str):
            try:
                architecture.load_coupling_map(Arch(arch))
            except ValueError:
                architecture.load_coupling_map(arch)
        elif isinstance(arch, Arch):
            architecture.load_coupling_map(arch)
        elif isinstance(arch, Architecture):
            architecture = arch
        elif isinstance(arch, Backend):
            from mqt.qmap.qiskit.backend import import_backend

            architecture = import_backend(arch)
        else:
            raise ValueError("No compatible type for architecture:", type(arch))

    if calibration is not None:
        if isinstance(calibration, str):
            architecture.load_properties(calibration)
        elif isinstance(calibration, BackendProperties):
            from mqt.qmap.qiskit.backend import import_backend_properties

            architecture.load_properties(import_backend_properties(calibration))
        elif isinstance(calibration, Target):
            from mqt.qmap.qiskit.backend import import_target

            architecture.load_properties(import_target(calibration))
        else:
            raise ValueError("No compatible type for calibration:", type(calibration))

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
