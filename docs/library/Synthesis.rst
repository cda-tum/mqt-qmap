Quantum Circuit Synthesis
=========================

The main functionality for synthesizing quantum circuits is provided by the :func:`~mqt.qmap.synthesize_clifford` and the :func:`~mqt.qmap.optimize_clifford` methods.



    .. currentmodule:: mqt.qmap
    .. autofunction:: synthesize_clifford
    .. autofunction:: optimize_clifford



Synthesis Configuration
#######################

The following classes provide a more explicit way of initializing the respective parameters of the :code:`synthesize_clifford` and the :code:`optimize_clifford` method. Instead of setting parameters via :code:`str` members, they can be set with :code:`Enumeration.member` where :code:`Enumeration` is the enumeration class.

    .. autoclass:: TargetMetric

    .. autoclass:: Verbosity

The :code:`SynthesisConfiguration` class is used to provide several parameters to the :code:`synthesize_clifford` and the :code:`optimize_clifford` method. The following table lists the parameters and their default values.

    .. autoclass:: SynthesisConfiguration
        :members:
