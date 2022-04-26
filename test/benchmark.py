from typing import Tuple

from mqt import qmap
from qiskit import QuantumCircuit
from os import walk


def parse_filename(filename: str) -> Tuple[str, int]:
    # strip trailing .qasm
    filename = filename.replace('.qasm', '')

    # split on '_'
    name, _, _, _, qubits = filename.split('_')

    return name, int(qubits)


if __name__ == '__main__':
    benchmark_location = "./benchmarks/"
    original_benchmark_dir = benchmark_location + "original/"

    mapped_benchmark_dir = benchmark_location + "mapped/"
    heuristic_mapped_benchmark_dir = mapped_benchmark_dir + "heuristic/"
    exact_mapped_benchmark_dir = mapped_benchmark_dir + "exact/"

    benchmark_categories = ["small", "medium", "large"]
    benchmarks = {}

    for category in benchmark_categories:
        _, _, benchmarks[category] = next(walk(original_benchmark_dir + category), (None, None, []))
        benchmarks[category] = sorted(benchmarks[category])

    for benchmark in benchmarks["small"]:
        name, qubits = parse_filename(benchmark)

        print('Starting benchmark: {} with {} qubits ... '.format(name, qubits), end='')
        qc = QuantumCircuit.from_qasm_file(original_benchmark_dir + "small/" + benchmark)
        qc.name = name

        result = qmap.compile(circ=qc, arch="IBMQ_Ehningen")
        print('finished!')

        qc_mapped = QuantumCircuit.from_qasm_str(result.mapped_circuit)
