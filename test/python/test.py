from mqt import qmap

from qiskit import QuantumCircuit
from qiskit.providers.fake_provider import FakeLondon

if __name__ == "__main__":
    qc = QuantumCircuit(3)
    qc.h(0)
    qc.cx(0, 1)
    qc.cx(1, 2)
    qc.measure_all()
    print(qc.draw(fold=-1))

    # compile the circuit
    qc_mapped, results = qmap.compile(qc, arch=FakeLondon())
    print(qc_mapped.draw(fold=-1))
