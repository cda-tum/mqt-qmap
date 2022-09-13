Quantum Circuit Mapping
=======================

Many quantum computing architectures limit the pairs of qubits that two-qubit operations can be applied to.
This is commonly described by a device's *coupling map*.
To execute a generic quantum circuit (with arbitrary interactions between its qubits) on such an architecture, the circuit needs to be *mapped*.
This involves *qubit allocation*, where logical qubits are assigned to physical qubits in an *initial layout*, and *routing*, where the original circuit is augmented with SWAP gates such that it adheres to the target device's *coupling map*.

Consider the following circuit.

    .. code-block:: python3

        from mqt import qmap
        from qiskit import QuantumCircuit

        qc = QuantumCircuit(4)
        qc.h(0)
        qc.cx(0, 1)
        qc.cx(0, 2)
        qc.cx(0, 3)

        qc.barrier()

        qc.t(0)
        qc.t(1)
        qc.t(2)
        qc.t(3)

        qc.barrier()

        qc.cx(0, 3)
        qc.cx(0, 2)
        qc.cx(0, 1)

        qc.measure_all()

        print(qc.draw(fold=-1))

    .. code-block:: console

               ┌───┐                ░ ┌───┐ ░                 ░ ┌─┐
          q_0: ┤ H ├──■────■────■───░─┤ T ├─░───■────■────■───░─┤M├─────────
               └───┘┌─┴─┐  │    │   ░ ├───┤ ░   │    │  ┌─┴─┐ ░ └╥┘┌─┐
          q_1: ─────┤ X ├──┼────┼───░─┤ T ├─░───┼────┼──┤ X ├─░──╫─┤M├──────
                    └───┘┌─┴─┐  │   ░ ├───┤ ░   │  ┌─┴─┐└───┘ ░  ║ └╥┘┌─┐
          q_2: ──────────┤ X ├──┼───░─┤ T ├─░───┼──┤ X ├──────░──╫──╫─┤M├───
                         └───┘┌─┴─┐ ░ ├───┤ ░ ┌─┴─┐└───┘      ░  ║  ║ └╥┘┌─┐
          q_3: ───────────────┤ X ├─░─┤ T ├─░─┤ X ├───────────░──╫──╫──╫─┤M├
                              └───┘ ░ └───┘ ░ └───┘           ░  ║  ║  ║ └╥┘
       meas: 4/══════════════════════════════════════════════════╩══╩══╩══╩═
                                                           0  1  2  3


Now assume this circuit shall be mapped to a 4-qubit architecture defined by the following coupling map:

.. image:: /images/linear_arch.svg
   :width: 60%
   :alt: 4-qubit Architecture
   :align: center

|

In *QMAP* this architecture can be manually defined.

    .. code-block:: python3

       arch = qmap.Architecture(
           4,
           {
               (0, 1),
               (1, 0),
               (1, 2),
               (2, 1),
               (2, 3),
               (3, 2),
           },
       )

The quantum circuit :code:`qc` can not be run directly on this architecture since it contains gates that act on qubits not connected on the device architecture.
Naively inserting SWAP gates that permute the logical-to-physical qubit mapping on the fly may yield the following compiled circuit.

    .. code-block:: console

                      ┌───┐                      ░ ┌───┐ ░                       ░ ┌─┐
          q_0 -> q_0: ┤ H ├──■───X───────────────░─┤ T ├─░───────────────X───■───░─┤M├─────────
                      └───┘┌─┴─┐ │               ░ ├───┤ ░               │ ┌─┴─┐ ░ └╥┘┌─┐
          q_1 -> q_1: ─────┤ X ├─X───■───X───────░─┤ T ├─░───────X───■───X─┤ X ├─░──╫─┤M├──────
                           └───┘   ┌─┴─┐ │       ░ ├───┤ ░       │ ┌─┴─┐   └───┘ ░  ║ └╥┘┌─┐
          q_2 -> q_2: ─────────────┤ X ├─X───■───░─┤ T ├─░───■───X─┤ X ├─────────░──╫──╫─┤M├───
                                   └───┘   ┌─┴─┐ ░ ├───┤ ░ ┌─┴─┐   └───┘         ░  ║  ║ └╥┘┌─┐
          q_3 -> q_3: ─────────────────────┤ X ├─░─┤ T ├─░─┤ X ├─────────────────░──╫──╫──╫─┤M├
                                           └───┘ ░ └───┘ ░ └───┘                 ░  ║  ║  ║ └╥┘
              meas: 4/══════════════════════════════════════════════════════════════╩══╩══╩══╩═
                                                                                    0  1  2  3

Over the course of the mapping, *four* SWAP gates have been introduce to satisfy the connectivity constraints of the device's architecture.
Since every additional gate increases the probability of errors, this is a very costly overhead for such a small circuit.

Keeping the number of additionally introduced gates as small as possible is key for ensuring the successful execution of the quantum circuit. Finding an optimal mapping for a quantum circuit is an NP-hard problem :cite:labelpar:`boteaComplexityQuantumCircuit2018`.
*QMAP* offers two dedicated techniques for tackling that problem:
- An *exact* mapping approach (based on :cite:labelpar:`willeMappingQuantumCircuits2019`, :cite:labelpar:`burgholzer2022limitingSearchSpace`) that guarantees (gate-optimal) solutions and is typically suitable for up to 8 qubits.
- A *heuristic* mapping approach (based on :cite:labelpar:`zulehnerEfficientMethodologyMapping2019`, :cite:labelpar:`hillmichExlpoitingQuantumTeleportation2021`) that allows to determine efficient mapping solutions in a scalable fashion for up to hundreds of qubits.

Exact Mapping
#############

The *exact mapper* implemented in *QMAP* maps quantum circuits using the *minimal* number of SWAP gates.
To this end, it encodes the mapping task as a MaxSAT problem and subsequently solves it using the `SMT solver Z3 <https://github.com/Z3Prover/z3>`_. Due to the NP-hardness of the mapping task, this approach is only scalable up to roughly eight qubits in most scenarios.

    .. note::
        On directional architectures, it can be significantly cheaper to surround a CNOT gate with four Hadamard operations (effectively exchanging its control and target) instead of adding a SWAP gate. For these architectures, QMAP minimizes the number of additional SWAP *and* H gates.

Using the exact mapper is as simple as:

    .. code-block:: python3

        qc_mapped, res = qmap.compile(qc, arch, method="exact")

        print(qc_mapped.draw(fold=-1))

        print("Additional gates: %d" % res.json()["statistics"]["additional_gates"])
        print("Runtime:          %f" % res.json()["statistics"]["mapping_time"])

    .. code-block:: console

                                       ┌───┐┌───┐┌───┐     ┌─┐
          q_3 -> 0 ────────────────────┤ X ├┤ T ├┤ X ├─────┤M├───────────────────
                                  ┌───┐└─┬─┘├───┤└─┬─┘┌───┐└╥┘          ┌─┐
          q_2 -> 1 ────────────■──┤ X ├──■──┤ T ├──■──┤ X ├─╫───■───────┤M├──────
                   ┌───┐     ┌─┴─┐└─┬─┘┌───┐└───┘     └─┬─┘ ║ ┌─┴─┐     └╥┘┌─┐
          q_0 -> 2 ┤ H ├──■──┤ X ├──■──┤ T ├────────────■───╫─┤ X ├──■───╫─┤M├───
                   └───┘┌─┴─┐├───┤     └───┘                ║ └───┘┌─┴─┐ ║ └╥┘┌─┐
          q_1 -> 3 ─────┤ X ├┤ T ├──────────────────────────╫──────┤ X ├─╫──╫─┤M├
                        └───┘└───┘                          ║      └───┘ ║  ║ └╥┘
              c: 4/═════════════════════════════════════════╩════════════╩══╩══╩═
                                                            3            2  0  1
          Additional gates: 2
          Runtime:          0.015844

By default, the :code:`compile` method optimizes sequences of CNOTs preceded/followed by SWAP gates in a post-processing step.
As a result, *no* SWAP gates are needed at all for mapping the circuit. Instead, only two additional CNOT operations are necessary to map the circuit.

Heuristic Mapping
#################

The *heuristic mapper* implemented in *QMAP* uses A\*-search to efficiently traverse the immense search space of the mapping problem.
It effectively trades optimality for runtime.
This allows to reliably determine suitable mappings for circuits with up to hundreds of qubits.
Using the heuristic mapper works completely analogous to the exact mapper.

    .. code-block:: python3

        qc_mapped, res = qmap.compile(qc, arch, method="heuristic")

        print(qc_mapped.draw(fold=-1))

        print("Additional gates: %d" % res.json()["statistics"]["additional_gates"])
        print("Runtime:          %f" % res.json()["statistics"]["mapping_time"])

    .. code-block:: console

	         ┌───┐┌───┐     ┌───┐                              ┌───┐   ┌─┐
	q_0 -> 0 ┤ H ├┤ X ├──■──┤ T ├──────────────────────────────┤ X ├───┤M├
	         └───┘└─┬─┘┌─┴─┐├───┤     ┌───┐               ┌───┐└─┬─┘┌─┐└╥┘
	q_1 -> 1 ───────■──┤ X ├┤ X ├──■──┤ T ├────────────■──┤ X ├──■──┤M├─╫─
	                   └───┘└─┬─┘┌─┴─┐└───┘┌───┐     ┌─┴─┐└─┬─┘ ┌─┐ └╥┘ ║
	q_2 -> 2 ─────────────────■──┤ X ├──■──┤ T ├──■──┤ X ├──■───┤M├──╫──╫─
	                             └───┘┌─┴─┐├───┤┌─┴─┐└┬─┬┘      └╥┘  ║  ║
	q_3 -> 3 ─────────────────────────┤ X ├┤ T ├┤ X ├─┤M├────────╫───╫──╫─
	                                  └───┘└───┘└───┘ └╥┘        ║   ║  ║
	    c: 4/══════════════════════════════════════════╩═════════╩═══╩══╩═
	                                                   3         2   0  1
	Additional gates: 3
        Runtime:          0.000065

While this solution is not optimal anymore it only requires one more gate and even for such a small example the heuristic mapper is orders of magnitudes faster than the exact mapper.
