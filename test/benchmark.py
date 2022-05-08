from typing import Tuple, List, Optional

from mqt import qmap
from qiskit import QuantumCircuit
from os import walk


def parse_filename(filename: str) -> Tuple[str, int]:
    # strip trailing .qasm
    filename = filename.replace('.qasm', '')

    # split on '_'
    name, _, _, _, qubits = filename.split('_')

    return name, int(qubits)


def heuristic_mapping(benchmark_location: str, mapped_circuit_location: str, benchmark_categories: Optional[List[str]] = None) -> None:
    if benchmark_categories is None:
        benchmark_categories = ["small", "medium", "large"]

    with open(benchmark_location + "heuristic_mapping_results.csv", 'w+') as heuristic_mapping_results:
        for category in benchmark_categories:
            _, _, benchmarks = next(walk(original_benchmark_dir + category), (None, None, []))
            benchmarks = sorted(benchmarks)

            for benchmark in benchmarks:
                name, qubits = parse_filename(benchmark)

                print('Starting benchmark: {} with {} qubits ... '.format(name, qubits), end='')
                qc = QuantumCircuit.from_qasm_file(original_benchmark_dir + category + "/" + benchmark)
                qc.name = name

                result = qmap.compile(circ=qc, arch="IBMQ_Ehningen")
                print('finished! Required {} swaps.'.format(result.output.swaps))
                heuristic_mapping_results.write(result.csv() + '\n')

                # write resulting circuit to file
                with open(mapped_circuit_location + category + "/" + benchmark, 'w+') as f:
                    f.write(result.mapped_circuit)


def exact_mapping(benchmark_location: str, mapped_circuit_location: str, benchmark_categories: Optional[List[str]] = None) -> None:
    if benchmark_categories is None:
        benchmark_categories = ["small"]

    with open(benchmark_location + "exact_mapping_results.csv", 'w+') as exact_mapping_results:
        for category in benchmark_categories:
            _, _, benchmarks = next(walk(original_benchmark_dir + category), (None, None, []))
            benchmarks = sorted(benchmarks)

            for benchmark in benchmarks:
                name, qubits = parse_filename(benchmark)
                print('Starting benchmark: {} with {} qubits ... '.format(name, qubits), end='')

                # orchestrate the exact mapper
                qc = QuantumCircuit.from_qasm_file(original_benchmark_dir + category + "/" + benchmark)
                qc.name = name

                subgraph_file = "subgraphs/ibmq_ehningen.txt"
                subgraphs = qmap.load_subgraphs_from_file(subgraph_file, qubits)
                best_result = None
                for subgraph in subgraphs:
                    result = qmap.compile(circ=qc, arch="IBMQ_Ehningen", method="exact", subgraph=subgraph)
                    if best_result is None or result.output.swaps < best_result.output.swaps:
                        best_result = result

                if best_result is not None:
                    print(' finished! Required {} swaps.'.format(best_result.output.swaps))

                    exact_mapping_results.write(best_result.csv() + '\n')

                    # write resulting circuit to file
                    with open(mapped_circuit_location + category + "/" + benchmark, 'w+') as f:
                        f.write(best_result.mapped_circuit)


if __name__ == '__main__':
    benchmark_location = "./benchmarks/"
    original_benchmark_dir = benchmark_location + "original/"

    mapped_benchmark_dir = benchmark_location + "mapped/"
    heuristic_mapped_benchmark_dir = mapped_benchmark_dir + "heuristic/"
    exact_mapped_benchmark_dir = mapped_benchmark_dir + "exact/"

    print("Starting heuristic mapping ...")
    heuristic_mapping(benchmark_location, heuristic_mapped_benchmark_dir)
    print("Starting exact mapping ...")
    exact_mapping(benchmark_location, exact_mapped_benchmark_dir)
