import os
import re
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


def normalize_subgraph(filename: str):
    with open(filename, 'r') as f:
        lines = f.readlines()
        qubits = int(lines[0])
        edges = lines[1:]
        edges = [(int(a), int(b)) for a, b in (edge.rstrip('\n').split(' ') for edge in edges)]
        # create an empty mapping between logical and physical qubits
        logical_to_physical_mapping = {}
        for i in range(qubits):
            logical_to_physical_mapping[i] = None

        # collect which physical qubits are used in the architecture
        used_physical_qubits = set()
        for a, b in edges:
            used_physical_qubits.add(a)
            used_physical_qubits.add(b)

        # create an empty mapping between physical and logical qubits
        physical_to_logical_mapping = {}
        for used_physical_qubit in used_physical_qubits:
            physical_to_logical_mapping[used_physical_qubit] = None

        current_qubit = 0
        for a, b in edges:
            if physical_to_logical_mapping[a] is None:
                physical_to_logical_mapping[a] = current_qubit
                logical_to_physical_mapping[current_qubit] = a
                current_qubit += 1

            if physical_to_logical_mapping[b] is None:
                physical_to_logical_mapping[b] = current_qubit
                logical_to_physical_mapping[current_qubit] = b
                current_qubit += 1

        normalized_edges = [(physical_to_logical_mapping[a], physical_to_logical_mapping[b]) for a, b in edges]
        normalized_filename = "IBMQ_Ehningen.arch"
        with open(normalized_filename, 'w+') as normalized_f:
            normalized_f.write(str(qubits) + '\n')
            for a, b in normalized_edges:
                normalized_f.write(str(a) + ' ' + str(b) + '\n')
        return logical_to_physical_mapping


def remap_qasm(device_qubits: int, logical_to_physical_mapping: dict, qasm: str):
    qasm_lines = qasm.split('\n')
    new_qasm = ''
    qubit_register_match = re.compile(r'q\[(\d+)]')

    for line in qasm_lines:
        # deal with initial and output permutation info
        if line.startswith('// i') or line.startswith('// o'):
            prefix = line[:5]
            line = line[5:]
            mapping = [logical_to_physical_mapping[int(i)] for i in line.split(' ')]
            new_qasm += prefix + ' '.join(str(i) for i in mapping) + '\n'
        elif line.startswith('qreg'):
            new_qasm += 'qreg q[' + str(device_qubits) + '];\n'
        elif line.startswith('creg') or line.startswith('include') or line.startswith('OPENQASM'):
            new_qasm += line + '\n'
        else:
            # try to match a gate
            new_qasm += re.sub(qubit_register_match, lambda x: 'q[' + str(logical_to_physical_mapping[int(x.group(1))]) + ']', line) + '\n'

    return new_qasm


def exact_mapping(device_qubits: int, benchmark_location: str, mapped_circuit_location: str, benchmark_categories: Optional[List[str]] = None) -> None:
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
                    result = qmap.compile(circ=qc, arch="IBMQ_Ehningen", method="exact", subgraph=subgraph, verbose=True)
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
    exact_mapping(27, benchmark_location, exact_mapped_benchmark_dir)
