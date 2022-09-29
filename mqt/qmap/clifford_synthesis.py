#
# This file is part of MQT QMAP library which is released under the MIT license.
# See file README.md or go to http://iic.jku.at/eda/research/quantum_verification/ for more information.
#

from __future__ import annotations

from mqt.qmap.load_architecture import load_architecture
from mqt.qmap.load_calibration import load_calibration
from mqt.qmap.pyqmap import (
    Arch,
    Architecture,
    OptimizationStrategy,
    SynthesisConfiguration,
    SynthesisResults,
    TargetMetric,
    optimize,
    synthesize,
)

from qiskit import QuantumCircuit
from qiskit.providers import Backend
from qiskit.providers.models import BackendProperties
from qiskit.transpiler.target import Target


def optimize_clifford(
    circ: QuantumCircuit | str,
    arch: str | Arch | Architecture | Backend | None = None,
    calibration: str | BackendProperties | Target | None = None,
    target: str | TargetMetric = "gates",
    strategy: str | OptimizationStrategy = "use_minimizer",
    choose_best: bool = False,
    initial_timestep: int = 10,
    nthreads: int = 1,
    verbosity: int = 0,
) -> tuple[QuantumCircuit, SynthesisResults]:
    """
    Optimize a circuit using the clifford synthesizer

    :param circ: Qiskit QuantumCircuit object or path to circuit file
    :type circ: QuantumCircuit | str
    :param arch: Architecture to synthesize the circuit on, optional. Either a path to a file with architecture information, one of the available architectures (:py:mod:`mqt.qmap.Arch`), Architecture, or `qiskit.providers.backend` (if Qiskit is installed)
    :type arch: str | Arch | Architecture | Backend | None
    :param calibration: Path to file containing calibration information, `qiskit.providers.models.BackendProperties` object (if Qiskit is installed), or `qiskit.transpiler.target.Target` object (if Qiskit is installed)
    :type calibration: str | BackendProperties | Target | None
    :param target: Target metric to optimize for. One of the available metrics (*gates* | depth | fidelity | two_qubit_gates)
    :type target: str | TargetMetric
    :param strategy: Optimization strategy to use. One of the available strategies (*use_minimizer* | minmax | start_low | start_high | split_iter)
    :type strategy: str | OptimizationStrategy
    :param choose_best: Whether to choose the fully connected subset from the architecture with the highest fidelity, or try all possible subsets, only relevant if architecture information is given
    :type choose_best: bool
    :param initial_timestep: Initial timesteps for the optimization, lower limit for start_low and upper limit for start_high
    :type initial_timestep: int
    :param nthreads: Number of threads to use for the optimization
    :type nthreads: int
    :param verbosity: Verbosity level of the debug output takes values from 0 (no output) to 5 (most output)
    :type verbosity: int

    :return: Optimized circuit (as Qiskit `QuantumCircuit`) and results
    :rtype: tuple[QuantumCircuit, SynthesisResults]
    """
    architecture = load_architecture(arch)

    architecture = load_calibration(calibration, architecture)

    config = SynthesisConfiguration()
    config.target_metric = TargetMetric(target)
    config.optimization_strategy = OptimizationStrategy(strategy)
    config.choose_best = choose_best
    config.initial_timestep = initial_timestep
    config.nthreads = nthreads
    config.verbosity = verbosity

    results = optimize(circ, architecture, config)

    return QuantumCircuit.from_qasm_str(results.resultCircuit), results


def synthesize_clifford(
    tableau: str,
    arch: str | Arch | Architecture | Backend | None = None,
    calibration: str | BackendProperties | Target | None = None,
    target: str | TargetMetric = "gates",
    strategy: str | OptimizationStrategy = "use_minimizer",
    choose_best: bool = False,
    initial_timestep: int = 10,
    nthreads: int = 1,
    verbosity: int = 0,
) -> tuple[QuantumCircuit, SynthesisResults]:
    """
    Synthesize a clifford circuit using the clifford synthesizer and a tableau input.

    :param tableau: String representation of the stabilizer tableau to be synthesized, either as semicolon separated binary matrix or output format of qiskit
    :type tableau: str
    :param arch: Architecture to synthesize the circuit on, optional. Either a path to a file with architecture information, one of the available architectures (:py:mod:`mqt.qmap.Arch`), Architecture, or `qiskit.providers.backend` (if Qiskit is installed)
    :type arch: str | Arch | Architecture | Backend | None
    :param calibration: Path to file containing calibration information, `qiskit.providers.models.BackendProperties` object (if Qiskit is installed), or `qiskit.transpiler.target.Target` object (if Qiskit is installed)
    :type calibration: str | BackendProperties | Target | None
    :param target: Target metric to synthesize for. One of the available metrics (*gates* | depth | fidelity | two_qubit_gates)
    :type target: str | TargetMetric
    :param strategy: Optimization strategy to use. One of the available strategies (*use_minimizer* | minmax | start_low | start_high | split_iter)
    :type strategy: str | OptimizationStrategy
    :param choose_best: Whether to choose the fully connected subset from the architecture with the highest fidelity, or try all possible subsets, only relevant if architecture information is given
    :type choose_best: bool
    :param initial_timestep: Initial timesteps for the synthesis, lower limit for start_low and upper limit for start_high
    :type initial_timestep: int
    :param nthreads: Number of threads to use for the synthesis
    :type nthreads: int
    :param verbosity: Verbosity level of the debug output takes values from 0 (no output) to 5 (most output)
    :type verbosity: int

    :return: Synthesized circuit (as Qiskit `QuantumCircuit`) and results
    :rtype: tuple[QuantumCircuit, SynthesisResults]
    """
    architecture = load_architecture(arch)

    architecture = load_calibration(calibration, architecture)

    config = SynthesisConfiguration()
    config.target_metric = TargetMetric(target)
    config.optimization_strategy = OptimizationStrategy(strategy)
    config.choose_best = choose_best
    config.initial_timestep = initial_timestep
    config.nthreads = nthreads
    config.verbosity = verbosity

    results = synthesize(tableau, architecture, config)

    return QuantumCircuit.from_qasm_str(results.resultCircuit), results
