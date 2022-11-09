Quantum Circuit Synthesis
===========================

The main functionality for synthesizing quantum circuits is the :code:`synthesize_clifford` and the :code:`optimize_clifford` method.



    .. currentmodule:: mqt.qmap
    .. automethod:: mqt.qmap::synthesize_clifford



Compile Setting Enumerations
############################

The following classes provide a more explicit way of initializing the respective parameters of the :code:`synthesize_clifford` and the :code:`optimize_clifford` method. Instead of setting parameters via :code:`str` members, they can be set with :code:`Enumeration.member` where :code:`Enumeration` is the enumeration class.

    .. currentmodule:: mqt.qmap
    .. autoclass:: TargetMetric

    .. autoclass:: OptimizationStrategy
