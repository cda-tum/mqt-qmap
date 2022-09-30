Quantum Circuit Synthesis
=======================

Many quantum algorithms are not available in circuit representation, that is directly executable on a quantum computer.
They rather are represented in a *functional representation*, which first needs to be translated into machine executable instructions.
The functional representation considered here is called *tableau*, which is a representation of a limited set of quantum functionality, called the *clifford group*.
This translation is commonly called *synthesis*.

Consider the following circuit.

    .. code-block:: python3

        from mqt import qmap
        from qiskit import QuantumCircuit

        qc = QuantumCircuit(2)
        qc.h(0)
        qc.cx(0, 1)
        qc.h(0)
        qc.h(1)

        print(qc.draw(fold=-1))

It has a corresponding functional representation, the *tableau* looks like the following:

    .. code-block:: console

        0 0 | 1 1 | 0
        1 1 | 0 0 | 0

with the corresponding *generators*:

    .. code-block:: python3

        tableau = "Stabilizer = ['+XX', '+ZZ']"


Keeping the number of used gates as small as possible is key for ensuring the successful execution of a quantum circuit. Finding an optimal circuit even for a finite group such as the Clifford Group grows with complexity 2^O(n^2).
In recent years the focus has shifted from gate optimal synthesis to optimizing depth or expected fidelity of a circuit. The complexity for these goals is of equal magnitude to gate optimal synthesis.

Circuit Synthesis
#################

The *clifford circuit synthesis* implemented in *QMAP* synthesizes quantum circuits using either minimal *gates*, *depth* or maximal expected *fidelity*. The latter is only available if the corresponding architecture information is given.
To this end, it encodes the mapping task as a MaxSAT problem and subsequently solves it using the `SMT solver Z3 <https://github.com/Z3Prover/z3>`_. Due to the complexity of 2^O(n^2) the mapping task, this approach is only scalable up to 26 qubits for most application.

Using the circuit synthesis is as simple as:

    .. code-block:: python3

        qc_mapped, results = qmap.synthesize_clifford(tableau, target="gates")

        print(qc_mapped.draw(fold=-1))

    .. code-block:: console



By default, the :code:`synthesize_clifford` method synthesizes a circuit with the minimum number of gates.
There are are also options for optimizing depth instead using :code:`target='depth'`.

To maximize expected fidelity if architecture/fidelity information is given :code:`target='fidelity'`:

    .. code-block:: python3

        arch = "ibmq_london.arch"
        calibration = "ibmq_london.csv"

        qc_mapped, results = qmap.synthesize_clifford(
            tableau, arch, calibration, target="fidelity"
        )

        print(qc_mapped.draw(fold=-1))


Circuit Optimization
####################

The *clifford circuit optimization* implemented in *QMAP* optimizes quantum circuits using either minimal *gates*, *depth* or maximal expected *fidelity*. The latter is only available if the corresponding architecture information is given.
Similar to synthesis, it encodes the mapping task as a MaxSAT problem and subsequently solves it using the `SMT solver Z3 <https://github.com/Z3Prover/z3>`_. Due to the complexity of 2^O(n^2) the mapping task, this approach is only scalable up to 26 qubits for most application.

Using the circuit optimization is as simple as:

    .. code-block:: python3

        qc_mapped, results = qmap.optimize_clifford(qc, target="gates")

        print(qc_mapped.draw(fold=-1))

    .. code-block:: console


By default, the :code:`optimize_clifford` method synthesizes a circuit with the minimum number of gates.
There are are also options for optimizing depth instead using :code:`target='depth'`.

To maximize expected fidelity if architecture/fidelity information is given :code:`target='fidelity'`:

    .. code-block:: python3

        arch = "ibmq_london.arch"
        calibration = "ibmq_london.csv"

        qc_mapped, results = qmap.synthesize_clifford(
            tableau, arch, calibration, target="fidelity"
        )

        print(qc_mapped.draw(fold=-1))

If architecture information is given, both the synthesis and optimization of a clifford circuit produce an already mapped circuit, in what is called *architecture aware synthesis*.
