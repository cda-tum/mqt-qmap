#
# This file is part of MQT QMAP library which is released under the MIT license.
# See file README.md or go to http://iic.jku.at/eda/research/quantum_verification/ for more information.
#

from __future__ import annotations

from mqt.qmap.pyqmap import (
    Arch,
    Architecture,
    CliffordOptResults,
    SynthesisStrategy,
    SynthesisTarget,
    optimize,
    synthesize,
)
from mqt.qmap.load_architecture import load_architecture
from mqt.qmap.load_calibration import load_calibration

from qiskit import QuantumCircuit
from qiskit.providers import Backend
from qiskit.providers.models import BackendProperties
from qiskit.transpiler.target import Target


def optimize_clifford(
        circ: QuantumCircuit | str,
        arch: str | Arch | Architecture | Backend | None = None,
        calibration: str | BackendProperties | Target | None = None,
        target: str | None = None,
        strategy: str | None = None,
) -> tuple[QuantumCircuit, CliffordOptResults]:
    """
    Optimize a circuit using the clifford synthesizer.

    """
    architecture = load_architecture(arch)

    architecture = load_calibration(calibration, architecture)

    strategy = SynthesisStrategy(strategy)
    target = SynthesisTarget(target)

    results = optimize(circ, architecture, strategy)

    return QuantumCircuit.from_qasm_str(results.resultCircuit), results


def synthesize_clifford(
        tableau: str,
        arch: str | Arch | Architecture | Backend | None = None,
        calibration: str | BackendProperties | Target | None = None,
        target: str | None = None,
        strategy: str | None = None,
) -> tuple[QuantumCircuit, CliffordOptResults]:
    """
    Synthesize a clifford circuit using the clifford synthesizer and a tableau input.

    """
    architecture = load_architecture(arch)

    architecture = load_calibration(calibration, architecture)

    strategy = SynthesisStrategy(strategy)
    target = SynthesisTarget(target)

    results = synthesize(tableau, architecture, strategy)

    return QuantumCircuit.from_qasm_str(results.resultCircuit), results
