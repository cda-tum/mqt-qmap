Quantum Circuit Mapping
=======================

Quantum computing architectures often provide only limited connectivity between physical qubits, meaning that two-qubit gates can not be executed with arbitrary qubits. Quantum algorithms are usually developed and implemented without device-specific restrictions in mind. Before running a quantum circuit on a quantum computer the circuit needs to be *mapped*. This involves *qubit allocation* where logical qubits are assigned to physical qubits in an *initial layout* and *routing* where the original circuit is amended with SWAP gates such that it adheres to the target device's *coupling map*.

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


Now assume this circuit shall be mapped to a 4-qubit architecture with the following coupling map:

.. image:: /images/linear_arch_4-1.png
   :width: 80%
   :alt: 4-qubit Architecture
   :align: center


In *QMAP* this architecture can be manually defined.

    .. code-block:: python3

        arch = qmap.Architecture(4, set([(0, 1), (1, 2), (2, 3), (3, 2), (2, 1), (1, 0)]))

The quantum circuit :code:`qc` can not be run directly on this architecture. Naive mapping yields the following possible compiled circuit.

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

This mapping uses 4 SWAPs. As every SWAP decomposes into 3 CNOT gates, this is a very costly overhead for such a small circuit.
Minimizing the number of additional gates introduced by mapping a quantum circuit is precisely the job of *QMAP*.

Exact Mapping
#############

The *exact mapper* implemented in *QMAP* maps quantum circuits using the *minimal* number of SWAP gates. For directional architectures it furthermore tries to minimize the number of additional Hadamard gates. Using the exact mapper is as simple as:

    .. code-block:: python3

        qc_mapped, res = qmap.compile(qc, arch, method="exact")

        print(qc_mapped.draw(fold=-1))

        print("Additional gates: %d" % res.json()["statistics"]["additional_gates"])
	print("Additional SWAPs: %d" % res.json()['statistics']['swaps'])
        print("Runtime:          %f" % res.json()['statistics']['mapping_time'])

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
	Additional SWAPs: 2
        Runtime:          0.015844

Not only does this solution yield less SWAPs, by default the :code:`compile` method optimizes SWAP gates in a post-processing step.

Heuristic Mapping
#################

Finding an optimal mapping for a quantum circuit is an NP-hard problem. *QMAP* employs the `SMT solver Z3 <https://github.com/Z3Prover/z3>`_ to encode and solve the problem of finding an optimal mapping as a Boolean formula. This, of course, brings scalablity issues along with it.

The *heuristic mapper* implemented in *QMAP* provides an alternative approach to quantum circuit mapping that trades off quality of the mapping result against runtime. Using the heuristic mapper works completely analogous to the exact mapper.

    .. code-block:: python3

        qc_mapped, res = qmap.compile(qc, arch, method="heuristic")

        print(qc_mapped.draw(fold=-1))

        print("Additional gates: %d" % res.json()["statistics"]["additional_gates"])
        print("Additional SWAPs: %d" % res.json()["statistics"]["swaps"])
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
	Additional SWAPs: 3
        Runtime:          0.000065

While this solution is not optimal anymore it only requires one more gate and even for such a small example the heuristic mapper is orders of magnitudes faster than the exact mapper.
