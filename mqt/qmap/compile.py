#
# This file is part of MQT QMAP library which is released under the MIT license.
# See file README.md or go to http://iic.jku.at/eda/research/quantum_verification/ for more information.
#
import pickle
from pathlib import Path
from typing import Union
from mqt.qmap.pyqmap import map, Method, InitialLayout, Layering, Arch, Encoding, CommanderGrouping, SwapReduction, Configuration, MappingResults


def compile(circ, arch: Union[str, Arch],
            calibration: str = "",
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
            verbose: bool = False
            ) -> MappingResults:
    """Interface to the MQT QMAP tool for mapping quantum circuits

    :param circ: Path to first circuit file, path to Qiskit QuantumCircuit pickle, or Qiskit QuantumCircuit object
    :param arch: Path to architecture file or one of the available architectures (Arch)
    :type arch: Union[str, Arch]
    :param calibration: Path to file containing calibration information
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
    :type use_bdd: bool    :param swap_reduction - Choose method of limiting the search space (none | *coupling_limit* | custom | increasing)
    :type swap_reduction: Union[str, SwapReduction]
    :param swap_limit - Set a custom limit for max swaps per layer, for the increasing reduction strategy it sets the max swaps per layer
    :type swap_limit: int
    :param include_WCNF: Include WCNF file in the results (default: False)
    :type include_WCNF: bool
    :param use_subsets: Use qubit subsets, or consider all available physical qubits at once (default: True)
    :type use_subsets: bool
    :param use_teleportation:  Use teleportation in addition to swaps
    :param teleportation_fake: Assign qubits as ancillary for teleportation in the initial placement but don't actually use them (used for comparisons)
    :param teleportation_seed: Fix a seed for the RNG in the initial ancilla placement (0 means the RNG will be seeded from /dev/urandom/ or similar)
    :param verbose: Print more detailed information during the mapping process
    :type verbose: bool
    :return: Object containing all the results
    :rtype: MappingResults
    """

    if type(circ) == str and Path(circ).suffix == '.pickle':
        circ = pickle.load(open(circ, "rb"))

    config = Configuration()
    config.calibration = calibration
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
    config.use_teleportation = use_teleportation
    config.teleportation_fake = teleportation_fake
    config.teleportation_seed = teleportation_seed
    config.verbose = verbose

    return map(circ, arch, config)
