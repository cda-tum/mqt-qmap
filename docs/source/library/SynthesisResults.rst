Results
=======

This class captures additional results from the :func:`~mqt.qmap.synthesize_clifford` and the :func:`~mqt.qmap.optimize_clifford` method.

    .. currentmodule:: mqt.qmap
    .. autoclass:: SynthesisResults

In addition to the synthesized/optimized circuit

    .. autoattribute:: SynthesisResults.circuit

it also includes runtime information of the mapping procedure and calls to the SAT solver.

    .. autoattribute:: SynthesisResults.runtime
    .. autoattribute:: SynthesisResults.solver_calls

and information about the synthesized circuit

    .. autoattribute:: SynthesisResults.gates
    .. autoattribute:: SynthesisResults.single_qubit_gates
    .. autoattribute:: SynthesisResults.two_qubit_gates
    .. autoattribute:: SynthesisResults.depth

it also includes the tableau of the resulting synthesized circuit

    .. autoattribute:: SynthesisResults.tableau

in addition it provides methods to check if the synthesis was successful or unsuccessful

    .. automethod:: SynthesisResults.sat
    .. automethod:: SynthesisResults.unsat

In addition, the class provides methods to export to other formats.

    .. automethod:: SynthesisResults.json
