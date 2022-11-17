Synthesis of Clifford Circuits
================================

Executing quantum circuits on a quantum computer requires compilation to representations that conform to all restrictions imposed by the device.
Due to device's limited coherence times and gate fidelity, the compilation process has to be optimized as much as possible.
To this end, an algorithm's description first has to be *synthesized* using the device's gate library.
In addition, circuits have to be *mapped* to the target quantum device to satisfy its connectivity constraints.
Even though Clifford circuits form a finite subgroup of all quantum circuits -- one that is not even universal for quantum computing -- the search space for these problems grows exponentially with respect to the number of considered qubits.

The *Clifford synthesis approach* in QMAP can be used to produce optimal Clifford circuits based on the methods proposed in :cite:labelpar:`schneider2023satEncodingOptimalClifford`.
To this end, it encodes the underlying task as a satisfiability (SAT) problem and solves it using the `SMT solver Z3 <https://github.com/Z3Prover/z3>`_ in conjunction with a binary search scheme.

The following gives a brief overview on Clifford circuits and how QMAP can be used for their synthesis.

Clifford Circuits
#################

Clifford circuits, i.e., circuits generated from the set :math:`\{H, S, \mathit{CNOT}\}`, form an important subclass of quantum circuits.
This is due to several factors

- According to the Gottesman-Knill theorem, they can be simulated in polynomial time and space on classical computers using the *stabilizer* formalism.
- They can be used to describe several quantum phenomena such as superposition, entanglement, superdense coding, and teleportation.
- Many error correcting codes rely on them.

Quantum states that can be obtained from the all-zero basis state :math:`|0\dots 0\rangle` by applying Clifford operations are called stabilizer states.
The name originates from the fact that such a state is uniquely and efficiently described by the set of operators that generate the group of its stabilizers.
Specifically, any *n*-qubit stabilizer state can be described by a set of *n* Pauli strings :math:`\pm P_{i,0}P_{i,1}P_{i,2}\dots P_{i,n-1}`, with :math:`P_{i,j}\in\{I, X, Y, Z\}` and :math:`i, j\in 0,\dots, n-1`.
Hence,two bits per qubit are needed to identify the Pauli operator, as well as one additional bit for the phase, which leads to a total of :math:`n(2n+1)` bits needed to uniquely describe a particular stabilizer state.

The stabilizer representation of a quantum state is conveniently described by a *tableau*:

.. math::

    \begin{bmatrix}
        x_{0,0}   & \cdots & x_{0,n-1}   & z_{0,0}    & \cdots & z_{0,n-1}   & r_0    \\
        \vdots    & \ddots &  \vdots         & \vdots     & \ddots &    \vdots         & \vdots \\
        x_{n-1,0} & \cdots & x_{n-1,n-1} & z_{n-1,0} & \cdots & z_{n-1,n-1} & r_{n-1}  \\
    \end{bmatrix}

Here, the binary variables :math:`x_{ij}` and :math:`z_{ij}` specify whether the Pauli term :math:`P_{i,j}` is :math:`X` or :math`Z`, respectively.
Since :math:`Y = iXZ`, setting :math:`x_{ij} = z_{ij} = 1` corresponds to :math:`P_{i,j}=Y`.
Finally, :math:`r_i` describes whether the generator has a negative phase.

Consider the following quantum circuit

    .. code-block:: python3

        from qiskit import QuantumCircuit

        qc = QuantumCircuit(2)
        qc.h(0)
        qc.cx(0, 1)
        qc.h(0)
        qc.h(1)

        print(qc.draw(fold=-1))

Then, the corresponding stabilizer tableau is

    .. code-block:: console

        0 0 | 1 1 | 0
        1 1 | 0 0 | 0

which corresponds to the stabilizers

    .. code-block:: python3

        stabilizers = ["+XX", "+ZZ"]

Using QMAP for Optimal Synthesis
################################

*QMAP* can be used in a multitude of ways to efficiently synthesize Clifford circuits:

- Starting from an initial (Clifford) circuit :code:`qc`, an optimal realization of that circuit's functionality can be determined as follows

    .. code-block:: python3

        from qiskit import QuantumCircuit
        from mqt import qmap

        qc = QuantumCircuit(2)
        qc.h(0)
        qc.cx(0, 1)
        qc.h(0)
        qc.h(1)

        qc_opt, results = qmap.optimize_clifford(qc)

        print(qc_opt.draw(fold=-1))

    .. code-block:: console

        TODO!

- Starting from a functional description, e.g., a list of stabilizers, an optimal realization of that functionality can be determined as follows

    .. code-block:: python3

        from qiskit.quantum_info import StabilizerTable
        from mqt import qmap

        stabilizers = ["+XX", "+ZZ"]
        table = StabilizerTable.from_labels(stabilizers)
        qc_synth, results = qmap.synthesize_clifford(table)

        print(qc_synth.draw(fold=-1))

    .. code-block:: console

        TODO!

The synthesis method offers lots of configuration options to fine-tune the synthesis procedure, change the cost metric, or to perform architecture-aware synthesis. Further details on that will be included in a future update.
