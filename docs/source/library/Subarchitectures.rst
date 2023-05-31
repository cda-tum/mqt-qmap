Optimal Subarchitectures
========================

To compute (near-)optimal subarchitectures of quantum computing architectures with restricted connectivity as described in :cite:labelpar:`peham2023OptimalSubarchitectures` the :code:`SubarchitectureOrder` class is provided. This class has functionality to compute the quasi-order that allows for fast computation of optimal subarchitectures.

Note that the initial construction of the ordering might take a while for larger architectures.

    .. currentmodule:: mqt.qmap
    .. autoclass:: SubarchitectureOrder
       :members:

QMAP also provides precomputed subarchitecture libraries. The available libraries are available via:

    .. autoattribute:: subarchitectures.precomputed_backends

Convenience methods are provided to import these precomputed orderings:

    .. automethod:: subarchitectures.ibm_guadalupe_subarchitectures
    .. automethod:: subarchitectures.rigetti_16_subarchitectures
