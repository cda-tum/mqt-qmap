#
# This file is part of MQT QMAP library which is released under the MIT license.
# See file README.md or go to http://iic.jku.at/eda/research/quantum_verification/ for more information.
#

from __future__ import annotations

from mqt.qmap.pyqmap import Architecture, Target

from qiskit.providers.models import BackendProperties


def load_calibration(
    calibration: str | BackendProperties | Target | None = None, architecture: Architecture = None
) -> Architecture:
    """
    Load a calibration from a string, BackendProperties, or Target.
    :param calibration: Path to file containing calibration information, `qiskit.providers.models.BackendProperties` object (if Qiskit is installed), or `qiskit.transpiler.target.Target` object (if Qiskit is installed)
    :type calibration: str | BackendProperties | Target | None
    :param architecture: Architecture to load the calibration for
    :type architecture: Architecture

    :return: Architecture
    :rtype: Architecture
    """
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
    return architecture
