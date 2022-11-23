Results
=======

This class captures additional results from the :func:`~mqt.qmap.synthesize_clifford` and the :func:`~mqt.qmap.optimize_clifford` method.

    .. currentmodule:: mqt.qmap
    .. autoclass:: SynthesisResults

In addition to the synthesized/optimized circuit

    .. autoattribute:: SynthesisResults.result_circuit

it also includes runtime information of the mapping procedure

    .. autoattribute:: SynthesisResults.total_seconds

and information about the synthesized circuit

    .. autoattribute:: SynthesisResults.single_qubit_gates
    .. autoattribute:: SynthesisResults.two_qubit_gates
    .. autoattribute:: SynthesisResults.depth
    .. autoattribute:: SynthesisResults.fidelity
    .. autoattribute:: SynthesisResults.qubits

it also includes information about the synthesis configuration

    .. autoattribute:: SynthesisResults.choose_best
    .. autoattribute:: SynthesisResults.initial_timestep
    .. autoattribute:: SynthesisResults.verbosity
    .. autoattribute:: SynthesisResults.strategy
    .. autoattribute:: SynthesisResults.target

In addition, the class provides methods to export to other formats.

    .. automethod:: SynthesisResults.json
