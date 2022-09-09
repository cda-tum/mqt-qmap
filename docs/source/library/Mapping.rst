Quantum Circuit Compilation
===========================

The main functionality for mapping quantum circuits is the :code:`compile` method.

    .. currentmodule:: mqt.qmap
    .. automethod:: mqt.qmap::compile

Some of the minor settings of the :code:`compile` method are given here.


Compile Setting Enumerations
############################

The following classes provide a more explicit way of initializing the respective parameters of the :code:`compile` method. Instead of setting parameters via :code:`str` members, they can be set with :code:`Enumeration.member` where :code:`Enumeration` is the enumeration class.

    .. currentmodule:: mqt.qmap
    .. autoclass:: Method

    .. autoclass:: Layering

    .. autoclass:: InitialLayout

    .. autoclass:: Encoding

    .. autoclass:: CommanderGrouping

    .. autoclass:: SwapReduction
