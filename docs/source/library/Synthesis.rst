Quantum Circuit Synthesis
===========================

The main functionality for synthesizing quantum circuits is provided by the :func:`~mqt.qmap.synthesize_clifford` and the :func:`~mqt.qmap.optimize_clifford` methods.



    .. currentmodule:: mqt.qmap
    .. automethod:: mqt.qmap::synthesize_clifford
    .. automethod:: mqt.qmap::optimize_clifford



Synthesis Configuration
#######################

The following classes provide a more explicit way of initializing the respective parameters of the :code:`synthesize_clifford` and the :code:`optimize_clifford` method. Instead of setting parameters via :code:`str` members, they can be set with :code:`Enumeration.member` where :code:`Enumeration` is the enumeration class.

    .. currentmodule:: mqt.qmap
    .. autoclass:: TargetMetric


    .. currentmodule:: mqt.qmap
    .. autoclass:: SynthesisConfiguration

The :code:`SynthesisConfiguration` class is used to provide several parameters to the :code:`synthesize_clifford` and the :code:`optimize_clifford` method. The following table lists the parameters and their default values.

    .. autoattribute:: SynthesisConfiguration.initial_timestep_limit
    .. autoattribute:: SynthesisConfiguration.use_maxsat
    .. autoattribute:: SynthesisConfiguration.target_metric
    .. autoattribute:: SynthesisConfiguration.use_symmetry_breaking
    .. autoattribute:: SynthesisConfiguration.n_threads
    .. autoattribute:: SynthesisConfiguration.minimize_gates_after_depth_optimization
    .. autoattribute:: SynthesisConfiguration.try_higher_gate_limit_for_two_qubit_gate_optimization
    .. autoattribute:: SynthesisConfiguration.gate_limit_factor
    .. autoattribute:: SynthesisConfiguration.minimize_gates_after_two_qubit_gate_optimization