from qiskit import QuantumCircuit
from qiskit.test.mock.backends import FakeLondon

from mqt import qmap

if __name__ == '__main__':
    qc = QuantumCircuit(3)
    qc.h(0)
    qc.cx(0, 1)
    qc.cx(1, 2)
    qc.measure_all()
    print(qc.draw(fold=-1))

    # compile the circuit
    results = qmap.compile(qc, arch=FakeLondon())

    # get the mapped circuit
    qc_mapped = QuantumCircuit.from_qasm_str(results.mapped_circuit)
    print(qc_mapped.draw(fold=-1))
