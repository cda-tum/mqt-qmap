Architecture
============

The *architecture* describes the properties of the target device when mapping a quantum circuit.

An Architecture can be constructed in various ways:

- From a Qiskit :code:`Backend` instance such as those defined under :code:`qiskit.providers.fake_provider` **(recommended)**,
- From a *coupling map* defined as a set of tuples of qubits.
- From the path to a file containing the number of qubits and a line-by-line enumeration of the qubit connections.


Constructing an Architecture
############################

From a Qiskit Backend
^^^^^^^^^^^^^^^^^^^^^

Architectures can be imported directly from Qiskit backends. Assuming :code:`backend` is a Qiskit backend, converting it to an Architecture is done as follows:

    .. code-block:: python3

        from mqt.qmap.qiskit.backend import import_backend

        architecture = import_backend(backend)

From a Coupling Map
^^^^^^^^^^^^^^^^^^^

A coupling map is defined as a set of tuples of qubits. The qubits are given as integers ranging from :math:`0` to :math:`num\_qubits - 1`.

    .. code-block:: python3

        from mqt import qmap

        coupling_map = {(0, 1), (1, 2), (2, 3)}
        num_qubits = 4

        arch = qmap.Architecture(num_qubits, coupling_map)

The architecture above is a directional architecture, meaning that control and target between connected qubits cannot be arbitrarily assigned. In the above example, a CNOT with qubit 0 as control and qubit 1 as target is possible, while the reverse CNOT with qubit 1 as control and qubit 0 as target is invalid.

A bidirectional architecture where any qubit can be either control or target of a CNOT gate can be defined by adding both directions for each edge.

    .. code-block:: python3

        from mqt import qmap

        coupling_map = {(0, 1), (1, 0), (1, 2), (2, 1), (2, 3), (3, 2)}
        num_qubits = 4

        arch = qmap.Architecture(num_qubits, coupling_map)


In addition to the number of qubits and the coupling map, QMAP's Architecture class also holds further device information such as gate and readout errors (see below for details). Note that QMAP does not consider this information during mapping at this point. Noise-aware mapping is planned for future releases of QMAP.

From a File
^^^^^^^^^^^

An architecture can be imported from a file. Assuming :code:`path` is a string variable containing the path to an architecture file, this is done as follows:

    .. code-block:: python3

        arch = Architecture()
        arch.load_coupling_map(path)

The first line of an architecture file is the number of qubits of the architecture. The remaining lines are tuples defining the connections of the coupling map. More explicitly, the file format is defined by the following Backus-Naur-Form.

    .. code-block:: console

        <coupling_map> ::= <integer>"\n"(<qubit>" "<qubit>"\n")*
        <qubit>        ::= 0|1|2| ... |nqubits-2|nqubits-1

Here the first integer defines the number of qubits of the architecture.

Full API of the Architecture class
##################################
    .. currentmodule:: mqt.qmap
    .. autoclass:: Architecture
        :special-members: __init__
        :undoc-members:
        :members:


mqt.qmap.Arch
#############

For convenience, this module provides several pre-defined architectures.

    .. automodule:: mqt.qmap.Arch
        :undoc-members:
        :members:
