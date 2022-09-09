Architecture
============

The *architecture* describes the properties of the target device when mapping a quantum circuit.

Constructing an Architecture
############################

An architecture can be constructed in various ways. The most convenient is to simply define the coupling map and pass it to the constructor of the Architecture class.

    .. code-block:: python3

        from mqt import qmap

        coupling_map = set([(0, 1), (1, 2), (2, 3), (2, 4)])
        num_qubits = 5

        arch = qmap.Architecture(num_qubits, coupling_map)

The architecture above is a directional architecture, meaning that control and target between connected qubits cannot be arbitrarily assigned. In the above example, a CNOT with qubit 0 as control and qubit 1 as target is possible, while the reverse CNOT with qubit 1 as control and qubit 0 as target is invalid.

A bidirectional architecture where any qubit can be either control or target of a CNOT gate can be defined by adding both directions for each edge.

    .. code-block:: python3

        from mqt import qmap

        coupling_map = set([(0, 1), (1, 0), (1, 2), (2, 1), (2, 3), (3, 2), (2, 4), (4, 2)])
        num_qubits = 5

        arch = qmap.Architecture(num_qubits, coupling_map)


In addition to the number of qubits and the coupling map, QMAP's Architecture class also holds further device information such as readout errors (see below for details). Note that QMAP does not consider this information during mapping at this point. Noise-aware mapping is planned for future releases of QMAP.

    .. currentmodule:: mqt.qmap

    .. automethod:: Architecture.__init__

    .. autoclass:: Architecture
        :undoc-members:
        :members:


mqt.qmap.Arch
#############

For convenience, this module provides several architectures.

    .. automodule:: mqt.qmap.Arch
        :undoc-members:
        :members:
