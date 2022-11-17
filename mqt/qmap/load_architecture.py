from __future__ import annotations

from mqt.qmap.pyqmap import Arch, Architecture

from qiskit.providers import Backend


def load_architecture(arch: str | Arch | Architecture | Backend | None = None) -> Architecture:
    """
    Load an architecture from a string, Arch, Architecture, or Backend. If None is passed, no architecture is loaded.
    :param arch: Architecture to map to. Either a path to a file with architecture information, one of the available architectures (:py:mod:`mqt.qmap.Arch`), Architecture, or `qiskit.providers.backend`
    :type arch: str | Arch | Architecture | Backend | None

    :return: Architecture
    :rtype: Architecture
    """
    architecture = Architecture()

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

    return architecture
