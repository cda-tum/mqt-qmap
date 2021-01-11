[![CI](https://github.com/iic-jku/qmap/workflows/CI/badge.svg)](https://github.com/iic-jku/qmap/actions?query=workflow%3A%22CI%22)
[![codecov](https://codecov.io/gh/iic-jku/qmap/branch/master/graph/badge.svg?token=TSFLDIO7HX)](https://codecov.io/gh/iic-jku/qmap)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![toolset: JKQ](https://img.shields.io/badge/toolset-JKQ-blue)](https://github.com/iic-jku/jkq)

# QMAP - A JKQ tool for Quantum Circuit Mapping written in C++
A [JKQ](https://github.com/iic-jku/jkq) tool for quantum circuit mapping by the [Institute for Integrated Circuits](https://iic.jku.at/eda/) at the [Johannes Kepler University Linz](https://jku.at) based on methods proposed in [[1]](https://iic.jku.at/files/eda/2018_tcad_efficient_mapping_of_quantum_circuits_to_ibm_qx_architectures.pdf), [[2]](https://iic.jku.at/files/eda/2019_dac_mapping_quantum_circuits_ibm_architectures_using_minimal_number_swap_h_gates.pdf).

[[1]](https://iic.jku.at/files/eda/2018_tcad_efficient_mapping_of_quantum_circuits_to_ibm_qx_architectures.pdf)
A. Zulehner, A. Paler, and R. Wille. **"An Efficient Methodology for Mapping Quantum Circuits to the IBM QX Architectures"**.
IEEE Transactions on Computer Aided Design of Integrated Circuits and Systems (TCAD), 2018.

[[2]](https://iic.jku.at/files/eda/2019_dac_mapping_quantum_circuits_ibm_architectures_using_minimal_number_swap_h_gates.pdf)
R. Wille, L. Burgholzer, and A. Zulehner. **"Mapping Quantum Circuits to IBM QX Architectures
Using the Minimal Number of SWAP and H Operations"**. In Design Automation Conference (DAC), 2019.

The tool can be used for mapping quantum circuits in any of the following formats:
* `Real` (e.g. from [RevLib](http://revlib.org)),
* `OpenQASM` (e.g. used by IBM's [Qiskit](https://github.com/Qiskit/qiskit)),
* `TFC` (e.g. from [Reversible Logic Synthesis Benchmarks Page](http://webhome.cs.uvic.ca/~dmaslov/mach-read.html))
* `QC` (e.g. from [Feynman](https://github.com/meamy/feynman))

to any given architecture, e.g., the IBM Q London architecture, which is specified by the coupling map
```
5
0 1
1 0
1 2
2 1
1 3
3 1
3 4
4 3
```
with the following available methods:
- **Heuristic Mapper**:  Heuristic solution based on A* search. For details see [[1]](https://iic.jku.at/files/eda/2018_tcad_efficient_mapping_of_quantum_circuits_to_ibm_qx_architectures.pdf).
- **Exact Mapper**: Exact solution utilizing the SMT Solver Z3. For details see [[2]](https://iic.jku.at/files/eda/2019_dac_mapping_quantum_circuits_ibm_architectures_using_minimal_number_swap_h_gates.pdf).

For more information, please visit [iic.jku.at/eda/research/ibm_qx_mapping/](https://iic.jku.at/eda/research/ibm_qx_mapping/).

If you have any questions, feel free to contact us via [iic-quantum@jku.at](mailto:iic-quantum@jku.at) or by creating an issue on [GitHub](https://github.com/iic-jku/qmap/issues).

## Usage

```commandline
$ ./qmap_heuristic --in grover_2.qasm --out grover_2m.qasm --arch ibmq_london.arch --ps
{
	"circuit": {
		"name": "grover_2",
		"qubits": 3,
		"gates": 30,
	},
	"mapped_circuit": {
		"name": "grover_2m",
		"qubits": 5,
		"gates": 33,
	},
	"statistics": {
		"mapping_time": 0.001638,
		"additional_gates": 3,
		"method": "heuristic",
		"arch": "ibmq_london"
	}
}
```
The heuristic mapping tool `qmap_heuristic` also offers the `--initiallayout` option, which allows to choose one of the following strategies for choosing an initial layout:
- `identity`: map logical qubit q_i to physical qubit Q_i,
- `static`: determine fixed initial layout statically at the start of mapping,
- `dynamic` (*default*): determine initial layout on demand during the mapping.

The exact mapping tool `qmap_exact` also offers the `--layering` option, which allows to choose one of the following strategies for partitioning the circuit:
- `individual` (*default*): consider each gate separately,
- `disjoint`: consider gates acting on disjoint qubits as a layer,
- `odd`: group pairs of gates. (Note that this strategy was only tested for IBM QX4 and may not work on different architectures)
- `triangle`: add gates to a layer, as long as no more than three qubits are involved. (Note that this strategy only works if the architecture's coupling map contains a triangle, e.g. IBM QX4)

### System Requirements
Building (and running) is continuously tested under Linux and MacOS using the [latest available system versions for GitHub Actions](https://github.com/actions/virtual-environments).
However, the implementation should be compatible with any current C++ compiler supporting C++14 and a minimum CMake version of 3.10.

`boost/program_options >= 1.50` is required for building the the commandline applications of the mapping tool.

In order to build the exact mapping tool, the SMT Solver [Z3 >= 4.8.3](https://github.com/Z3Prover/z3) has to be installed and on the path.

### Build and Run

For building the commandline applications the following commands should be used
```
$ cmake -DCMAKE_BUILD_TYPE=Release -S . -B build
$ cmake --build build --config Release
```
Windows users need to configure CMake by calling
```
$ cmake -G "Visual Studio 15 2017" -A x64 -DCMAKE_BUILD_TYPE=Release -S . -B build
```
instead.

Afterwards `./build/qmap_heuristic` and `./build/qmap_exact` (requires [Z3 >= 4.8.3](https://github.com/Z3Prover/z3)) are available.

## Reference

If you use our tool for your research, we will be thankful if you refer to it by citing one or both of the following publications.

If you used the heuristic mapping, please cite
```bibtex
@article{zulehner2019efficient,
  title={An efficient methodology for mapping quantum circuits to the {IBM} {QX} architectures},
  author={Zulehner, Alwin and Paler, Alexandru and Wille, Robert},
  journal={{IEEE} Transactions on Computer-Aided Design of Integrated Circuits and Systems},
  volume={38},
  number={7},
  pages={1226--1236},
  year={2019}
}
```

If you used the exact mapping, please cite
```bibtex
@inproceedings{wille2019mapping,
    title={Mapping Quantum Circuits to {IBM QX} Architectures Using the Minimal Number of {SWAP} and {H} Operations},
    author={Wille, Robert and Burgholzer, Lukas and Zulehner, Alwin},
    booktitle={Design Automation Conference},
    year={2019}
}
````
