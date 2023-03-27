"""Module for importing Qiskit backends."""

from __future__ import annotations

from typing import TYPE_CHECKING

from qiskit.providers import Backend, BackendV1, BackendV2

if TYPE_CHECKING:  # pragma: no cover
    from qiskit.providers.models import BackendProperties
    from qiskit.transpiler import Target

from mqt.qmap import Architecture


def import_backend(backend: Backend) -> Architecture:
    """Import a backend from qiskit.providers.Backend.

    Args:
        backend: The backend to import.

    Returns:
        The imported backend as an Architecture.

    """
    if isinstance(backend, BackendV1):
        return import_backend_v1(backend)
    if isinstance(backend, BackendV2):
        return import_backend_v2(backend)
    msg = f"Backend type {type(backend)} not supported."
    raise TypeError(msg)


def import_backend_properties(backend_properties: BackendProperties) -> Architecture.Properties:
    """Import backend properties from qiskit.providers.models.BackendProperties.

    Args:
        backend_properties: The backend properties to import.

    Returns:
        The imported backend properties as an Architecture.Properties object.
    """
    props = Architecture.Properties()
    props.name = backend_properties.backend_name
    props.num_qubits = len(backend_properties.qubits)
    for qubit in range(props.num_qubits):
        props.set_t1(qubit, backend_properties.t1(qubit))
        props.set_t2(qubit, backend_properties.t2(qubit))
        props.set_frequency(qubit, backend_properties.frequency(qubit))
        props.set_readout_error(qubit, backend_properties.readout_error(qubit))

    for gate in backend_properties.gates:
        if gate.gate == "reset":
            continue

        if len(gate.qubits) == 1:
            props.set_single_qubit_error(
                gate.qubits[0], gate.gate, backend_properties.gate_error(gate.gate, gate.qubits)
            )
        elif len(gate.qubits) == 2:
            props.set_two_qubit_error(
                gate.qubits[0], gate.qubits[1], backend_properties.gate_error(gate.gate, gate.qubits), gate.gate
            )
    return props


def import_backend_v1(backend: BackendV1) -> Architecture:
    """Import a backend from qiskit.providers.BackendV1."""
    arch = Architecture()
    arch.name = backend.name()
    arch.num_qubits = backend.configuration().n_qubits
    arch.coupling_map = {(a, b) for a, b in backend.configuration().coupling_map}
    arch.properties = import_backend_properties(backend.properties())

    return arch


def import_target(target: Target) -> Architecture.Properties:
    """Import a target from qiskit.transpiler.Target.

    Args:
        target: The target to import.

    Returns:
        The imported target as an Architecture.Properties object.
    """
    props = Architecture.Properties()
    props.num_qubits = len(target.qubit_properties)

    for i in range(props.num_qubits):
        qubit_props = target.qubit_properties[i]
        props.set_t1(i, qubit_props.t1)
        props.set_t2(i, qubit_props.t2)
        props.set_frequency(i, qubit_props.frequency)

    for instruction, qargs in target.instructions:
        if instruction.name == "reset" or instruction.name == "delay":
            continue

        instruction_props = target[instruction.name][qargs]
        if instruction.name == "measure":
            props.set_readout_error(qargs[0], instruction_props.error)
        elif len(qargs) == 1:
            props.set_single_qubit_error(qargs[0], instruction.name, instruction_props.error)
        elif len(qargs) == 2:
            props.set_two_qubit_error(qargs[0], qargs[1], instruction_props.error, instruction.name)

    return props


def import_backend_v2(backend: BackendV2) -> Architecture:
    """Import a backend from qiskit.providers.BackendV2."""
    arch = Architecture()
    arch.name = backend.name
    arch.num_qubits = backend.num_qubits
    arch.coupling_map = set(backend.coupling_map.get_edges())
    arch.properties = import_target(backend.target)
    arch.properties.name = backend.name

    return arch
