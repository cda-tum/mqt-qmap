#
# This file is part of MQT QMAP library which is released under the MIT license.
# See file README.md or go to http://iic.jku.at/eda/research/quantum_verification/ for more information.
#

from __future__ import annotations

from mqt.qmap.compile import extract_initial_layout_from_qasm
from mqt.qmap.pyqmap import (
    Arch,
    Architecture,
    OptimizationStrategy,
    SynthesisConfiguration,
    SynthesisResults,
    TargetMetric,
    synthesize,
)

from qiskit import QuantumCircuit
from qiskit.providers import Backend
from qiskit.providers.models import BackendProperties
from qiskit.quantum_info import Clifford, StabilizerTable
from qiskit.transpiler.target import Target

from .load_architecture import load_architecture
from .load_calibration import load_calibration


def optimize_clifford(
    description: str | QuantumCircuit,
    arch: str | Arch | Architecture | Backend | None = None,
    calibration: str | BackendProperties | Target | None = None,
    target: str | TargetMetric = "gates",
    strategy: str | OptimizationStrategy = "use_minimizer",
    choose_best: bool = False,
    initial_timestep: int = 10,
    fidelity_scaling: int = 1000,
    limit_finding_factor: float = 0.5,
    circuit_splitting_increase: float = 0.2,
    nthreads: int = 1,
    verbosity: int = 0,
) -> tuple[QuantumCircuit, SynthesisResults]:
    """
    Optimize a Clifford circuit using the Clifford synthesizer.

    :param description: Description of the Clifford functionality to be optimized. Either
     - a Qiskit :code:`QuantumCircuit`
     - a path to a qasm file
    :type description: str | StabilizerTable | Clifford | QuantumCircuit
    :param arch: Architecture to synthesize the circuit on, optional. Either a path to a file with architecture information, one of the available architectures (:py:mod:`mqt.qmap.Arch`), Architecture, or `qiskit.providers.backend`
    :type arch: str | Arch | Architecture | Backend | None
    :param calibration: Path to file containing calibration information, `qiskit.providers.models.BackendProperties` object (if Qiskit is installed), or `qiskit.transpiler.target.Target` object
    :type calibration: str | BackendProperties | Target | None
    :param target: Target metric to synthesize for. One of the available metrics (*gates* | depth | fidelity | two_qubit_gates)
    :type target: str | TargetMetric
    :param strategy: Optimization strategy to use. One of the available strategies (*use_minimizer* | minmax | start_low | start_high | split_iter)
    :type strategy: str | OptimizationStrategy
    :param choose_best: Whether to choose the fully connected subset from the architecture with the highest fidelity, or try all possible subsets, only relevant if architecture information is given
    :type choose_best: bool
    :param initial_timestep: Initial timesteps for the synthesis, lower limit for start_low and upper limit for start_high
    :type initial_timestep: int
    :param fidelity_scaling: Fidelity scaling factor for the synthesis *1000*, higher values are needed if the fidelities are very similar or high
    :type fidelity_scaling: int
    :param limit_finding_factor: Factor to multiply the initial guess for timesteps used in startLow (increase by 1 + factor) or startHigh (decrease by factor) *0.5*
    :type limit_finding_factor: float
    :param circuit_splitting_increase: Factor to multiply circuit splitting size by, calculation is as follows increase if a unsat instance is found is max(1, split * factor) *0.2*
    :type circuit_splitting_increase: float
    :param nthreads: Number of threads to use for the synthesis
    :type nthreads: int
    :param verbosity: Verbosity level of the debug output takes values from 0 (no output) to 5 (most output)
    :type verbosity: int

    :return: Synthesized circuit (as Qiskit `QuantumCircuit`) and results
    :rtype: tuple[QuantumCircuit, SynthesisResults]
    """
    architecture = load_architecture(arch)
    load_calibration(architecture, calibration)

    config = SynthesisConfiguration()
    config.target_metric = TargetMetric(target)
    config.optimization_strategy = OptimizationStrategy(strategy)
    config.choose_best = choose_best
    config.initial_timestep = initial_timestep
    config.fidelity_scaling = fidelity_scaling
    config.limit_finding_factor = limit_finding_factor
    config.circuit_splitting_increase = circuit_splitting_increase
    config.nthreads = nthreads
    config.verbosity = verbosity

    results = synthesize(description, architecture, config)

    return QuantumCircuit.from_qasm_str(results.result_circuit), results


def synthesize_clifford(
    description: str | StabilizerTable | Clifford,
    arch: str | Arch | Architecture | Backend | None = None,
    calibration: str | BackendProperties | Target | None = None,
    target: str | TargetMetric = "gates",
    strategy: str | OptimizationStrategy = "use_minimizer",
    choose_best: bool = False,
    initial_timestep: int = 10,
    fidelity_scaling: int = 1000,
    limit_finding_factor: float = 0.5,
    circuit_splitting_increase: float = 0.2,
    nthreads: int = 1,
    verbosity: int = 0,
) -> tuple[QuantumCircuit, SynthesisResults]:
    """
    Synthesize a Clifford circuit using the Clifford synthesizer.

    :param description: Description of the Clifford functionality to be synthesized. Either
     - a Qiskit :code`StabilizerTable` or :code:`Clifford` operator
     - a list of stabilizers (as exported by Qiskit) or as a semicolon separated binary matrix
    :type description: str | StabilizerTable | Clifford
    :param arch: Architecture to synthesize the circuit on, optional. Either a path to a file with architecture information, one of the available architectures (:py:mod:`mqt.qmap.Arch`), Architecture, or `qiskit.providers.backend`
    :type arch: str | Arch | Architecture | Backend | None
    :param calibration: Path to file containing calibration information, `qiskit.providers.models.BackendProperties` object (if Qiskit is installed), or `qiskit.transpiler.target.Target` object
    :type calibration: str | BackendProperties | Target | None
    :param target: Target metric to synthesize for. One of the available metrics (*gates* | depth | fidelity | two_qubit_gates)
    :type target: str | TargetMetric
    :param strategy: Optimization strategy to use. One of the available strategies (*use_minimizer* | minmax | start_low | start_high | split_iter)
    :type strategy: str | OptimizationStrategy
    :param choose_best: Whether to choose the fully connected subset from the architecture with the highest fidelity, or try all possible subsets, only relevant if architecture information is given
    :type choose_best: bool
    :param initial_timestep: Initial timesteps for the synthesis, lower limit for start_low and upper limit for start_high
    :type initial_timestep: int
    :param fidelity_scaling: Fidelity scaling factor for the synthesis *1000*, higher values are needed if the fidelities are very similar or high
    :type fidelity_scaling: int
    :param limit_finding_factor: Factor to multiply the initial guess for timesteps used in startLow (increase by 1 + factor) or startHigh (decrease by factor) *0.5*
    :type limit_finding_factor: float
    :param circuit_splitting_increase: Factor to multiply circuit splitting size by, calculation is as follows increase if a unsat instance is found is max(1, split * factor) *0.2*
    :type circuit_splitting_increase: float
    :param nthreads: Number of threads to use for the synthesis
    :type nthreads: int
    :param verbosity: Verbosity level of the debug output takes values from 0 (no output) to 5 (most output)
    :type verbosity: int

    :return: Synthesized circuit (as Qiskit `QuantumCircuit`) and results
    :rtype: tuple[QuantumCircuit, SynthesisResults]
    """
    architecture = load_architecture(arch)
    load_calibration(architecture, calibration)

    config = SynthesisConfiguration()
    config.target_metric = TargetMetric(target)
    config.optimization_strategy = OptimizationStrategy(strategy)
    config.choose_best = choose_best
    config.initial_timestep = initial_timestep
    config.fidelity_scaling = fidelity_scaling
    config.limit_finding_factor = limit_finding_factor
    config.circuit_splitting_increase = circuit_splitting_increase
    config.nthreads = nthreads
    config.verbosity = verbosity

    stabilizers: str = ""
    if isinstance(description, Clifford):
        stabilizers = str(description.stabilizer.to_labels())
    elif isinstance(description, StabilizerTable):
        stabilizers = str(description.to_labels())
    elif isinstance(description, str):
        stabilizers = description

    results = synthesize(stabilizers, architecture, config)

    circ = QuantumCircuit.from_qasm_str(results.result_circuit)
    layout = extract_initial_layout_from_qasm(results.result_circuit, circ.qregs)
    circ._layout = layout

    return circ, results
