from __future__ import annotations

from mqt.qmap.pyqmap import Architecture

from qiskit.providers.models import BackendProperties
from qiskit.transpiler.target import Target


def load_calibration(architecture: Architecture, calibration: str | BackendProperties | Target | None = None) -> None:
    """
    Load a calibration from a string, BackendProperties, or Target.
    :param architecture: Architecture to load the calibration for
    :type architecture: Architecture
    :param calibration: Path to file containing calibration information, `qiskit.providers.models.BackendProperties` object, or `qiskit.transpiler.target.Target` object
    :type calibration: str | BackendProperties | Target | None
    """
    if calibration is None:
        return

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
