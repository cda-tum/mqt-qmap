#
# This file is part of MQT QMAP library which is released under the MIT license.
# See file README.md or go to http://iic.jku.at/eda/research/quantum_verification/ for more information.
#
import pickle
from pathlib import Path
from typing import Union, Optional, Set, List, Tuple
from mqt.qmap.pyqmap import map, Method, InitialLayout, Layering, Arch, Encoding, CommanderGrouping, SwapReduction, Configuration, MappingResults, Architecture

try:
    from qiskit import QuantumCircuit, QuantumRegister
    from qiskit.providers import Backend
    from qiskit.providers.models import BackendProperties
    from qiskit.transpiler.target import Target
    from qiskit.transpiler import Layout

    PossibleArchitectureTypes = Union[str, Arch, Architecture, Backend]
    PossibleCalibrationTypes = Union[str, BackendProperties, Target]
    CircuitReturnType = QuantumCircuit


    def extract_initial_layout_from_qasm(qasm: str, qregs: List[QuantumRegister]) -> Layout:
        """
        Extracts the initial layout resulting from compiling a circuit from a QASM file.
        :param qasm: QASM file
        :type qasm: str
        :param qregs: The quantum registers to apply the layout to.
        :type qregs: List[QuantumRegister]
        :return: layout to be used in Qiskit
        """
        for line in qasm.split("\n"):
            if line.startswith("// i "):
                # strip away initial part of line
                line = line[5:]
                # split line into tokens
                tokens = line.split(" ")
                # convert tokens to integers
                tokens = [int(token) for token in tokens]
                # create an empty layout
                layout = Layout().from_intlist(tokens, *qregs)
                return layout


except ModuleNotFoundError:
    PossibleArchitectureTypes = Union[str, Arch, Architecture]
    PossibleCalibrationTypes = str
    CircuitReturnType = str


def compile(circ, arch: Optional[PossibleArchitectureTypes],
            calibration: Optional[PossibleCalibrationTypes] = None,
            method: Union[str, Method] = "heuristic",
            initial_layout: Union[str, InitialLayout] = "dynamic",
            layering: Union[str, Layering] = "individual_gates",
            use_teleportation: bool = False,
            teleportation_fake: bool = False,
            teleportation_seed: int = 0,
            encoding: Union[str, Encoding] = "naive",
            commander_grouping: Union[str, CommanderGrouping] = "halves",
            use_bdd: bool = False,
            swap_reduction: Union[str, SwapReduction] = "coupling_limit",
            swap_limit: int = 0,
            include_WCNF: bool = False,
            use_subsets: bool = True,
            subgraph: Optional[Set[int]] = None,
            pre_mapping_optimizations: bool = True,
            post_mapping_optimizations: bool = True,
            verbose: bool = False
            ) -> Tuple[CircuitReturnType, MappingResults]:
    """Interface to the MQT QMAP tool for mapping quantum circuits

    :param circ: Path to circuit file, path to Qiskit QuantumCircuit pickle, or Qiskit QuantumCircuit object
    :param arch: Architecture to map to. Either a path to a file with architecture information, one of the available architectures (Arch), qmap.Architecture, or `qiskit.providers.backend` (if Qiskit is installed)
    :type arch: Optional[PossibleArchitectureTypes]
    :param calibration: Path to file containing calibration information, `qiskit.providers.models.BackendProperties` object (if Qiskit is installed), or `qiskit.transpiler.target.Target` object (if Qiskit is installed)
    :type calibration: Optional[PossibleCalibrationTypes]
    :param method: Mapping technique to use (*heuristic* | exact)
    :type method: Union[str, Method]
    :param initial_layout: Strategy to use for determining initial layout in heuristic mapper (identity | static | *dynamic*)
    :type initial_layout: Union[str, InitialLayout]
    :param layering: Circuit layering strategy to use (*individual_gates* | disjoint_qubits | odd_qubits | qubit_triangle)
    :type layering: Union[str, Layering]
    :param encoding - Choose encoding for AMO and exactly one constraints (*naive* | commander | bimander)
    :type encoding: Union[str, Encoding]
    :param commander_grouping - Choose method of grouping in commander and bimander encoding (*halves* | fixed2 | fixed3 | logarithm)
    :type commander_grouping: Union[str, CommanderGrouping]
    :param use_bdd: Limit swaps per layer using BDDs, faster in some cases, but use with caution (default: False)
    :type use_bdd: bool
    :param swap_reduction - Choose method of limiting the search space (none | *coupling_limit* | custom | increasing)
    :type swap_reduction: Union[str, SwapReduction]
    :param swap_limit - Set a custom limit for max swaps per layer, for the increasing reduction strategy it sets the max swaps per layer
    :type swap_limit: int
    :param include_WCNF: Include WCNF file in the results (default: False)
    :type include_WCNF: bool
    :param use_subsets: Use qubit subsets, or consider all available physical qubits at once (default: True)
    :type use_subsets: bool
    :param subgraph: List of qubits to consider for mapping (in exact mapper), if None all qubits are considered
    :type subgraph: Optional[List[int]]
    :param use_teleportation:  Use teleportation in addition to swaps
    :param teleportation_fake: Assign qubits as ancillary for teleportation in the initial placement but don't actually use them (used for comparisons)
    :param teleportation_seed: Fix a seed for the RNG in the initial ancilla placement (0 means the RNG will be seeded from /dev/urandom/ or similar)
    :param pre_mapping_optimizations: Run pre-mapping optimizations (default: True)
    :type pre_mapping_optimizations: bool
    :param post_mapping_optimizations: Run post-mapping optimizations (default: True)
    :type post_mapping_optimizations: bool
    :param verbose: Print more detailed information during the mapping process
    :type verbose: bool
    :return: Mapped circuit (either as Qiskit `QuantumCircuit`, if Qiskit is available, or `.qasm` string) and results
    :rtype: Tuple[CircuitReturnType, MappingResults]
    """

    if subgraph is None:
        subgraph = set()
    if type(circ) == str and Path(circ).suffix == '.pickle':
        circ = pickle.load(open(circ, "rb"))

    architecture = Architecture()
    if arch is None and calibration is None:
        raise ValueError("Either arch or calibration must be specified")

    if arch is not None:
        if type(arch) == str:
            try:
                architecture.load_coupling_map(Arch(arch))
            except ValueError:
                architecture.load_coupling_map(arch)
        elif type(arch) == Arch:
            architecture.load_coupling_map(arch)
        elif isinstance(arch, Architecture):
            architecture = arch
        else:
            try:
                from qiskit.providers.backend import Backend
                from mqt.qmap.qiskit.backend import import_backend
                if isinstance(arch, Backend):
                    architecture = import_backend(arch)
                else:
                    raise ValueError("No compatible type for architecture:", type(arch))
            except ModuleNotFoundError:
                raise ValueError("No compatible type for architecture:", type(arch))

    if calibration is not None:
        if type(calibration) == str:
            architecture.load_properties(calibration)
        else:
            try:
                from qiskit.providers.models import BackendProperties
                from qiskit.transpiler.target import Target
                from mqt.qmap.qiskit.backend import import_backend_properties, import_target

                if isinstance(calibration, BackendProperties):
                    architecture.load_properties(import_backend_properties(calibration))
                elif isinstance(calibration, Target):
                    architecture.load_properties(import_target(calibration))
                else:
                    raise ValueError("No compatible type for calibration:", type(calibration))
            except ModuleNotFoundError:
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
    config.verbose = verbose

    results = map(circ, architecture, config)

    try:
        from qiskit import QuantumCircuit

        circ = QuantumCircuit.from_qasm_str(results.mapped_circuit)
        layout = extract_initial_layout_from_qasm(results.mapped_circuit, circ.qregs)
        circ._layout = layout
    except ModuleNotFoundError:
        circ = results.mapped_circuit

    return circ, results
