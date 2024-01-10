Quantum Circuit Compilation
===========================

The main functionality for mapping quantum circuits is the :code:`compile` method.

    .. note::
        This method assumes that the input quantum circuit is already decomposed into elementary gates for the specified architecture.

    .. currentmodule:: mqt.qmap
    .. autofunction:: compile



Compile Setting Enumerations
############################

The following classes provide a more explicit way of initializing the respective parameters of the :code:`compile` method. Instead of setting parameters via :code:`str` members, they can be set with :code:`Enumeration.member` where :code:`Enumeration` is the enumeration class.

    .. currentmodule:: mqt.qmap
    .. autoclass:: Method
        
    .. autoclass:: Heuristic
        
    .. autoclass:: LookaheadHeuristic

    .. autoclass:: Layering

    .. autoclass:: InitialLayout
        
    .. autoclass:: EarlyTermination

    .. autoclass:: Encoding

    .. autoclass:: CommanderGrouping

    .. autoclass:: SwapReduction
