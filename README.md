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
* `QuantumCircuit` object from IBM's [Qiskit](https://github.com/Qiskit/qiskit) (only through the JKQ QMAP Python bindings)
* `OpenQASM` (e.g. used by IBM's [Qiskit](https://github.com/Qiskit/qiskit)),
* `Real` (e.g. from [RevLib](http://revlib.org)),
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

Note that, at the moment, circuits to be mapped are assumed to be already decomposed into elementary gates supported by the targeted device. More specifically, circuits must not contain gates acting on more than two qubits.

For more information, please visit [iic.jku.at/eda/research/ibm_qx_mapping/](https://iic.jku.at/eda/research/ibm_qx_mapping/).

If you have any questions, feel free to contact us via [iic-quantum@jku.at](mailto:iic-quantum@jku.at) or by creating an issue on [GitHub](https://github.com/iic-jku/qmap/issues).

## Usage

JKQ QMAP is mainly developed as a C++ library with a [commandline interface](#command-line-executable). However, using it in Python is as easy as
```bash
pip install jkq.qmap
```
and then in Python
```python
from jkq import qmap
qmap.compile(...)
```
where the `compile` function is defined as follows:
```python
"""
Interface to the JKQ QMAP tool for mapping quantum circuits

Params:
    circ – Path to first circuit file, path to Qiskit QuantumCircuit pickle, or Qiskit QuantumCircuit object
    arch – Path to architecture file or one of the available architectures (Arch)
    calibration – Path to file containing calibration information
    method – Mapping technique to use (*heuristic* | exact)
    initial_layout – Strategy to use for determining initial layout (only relevant for heuristic mapper)
    layering – Circuit layering strategy to use (*individual_gates* | disjoint_qubits | odd_qubits | qubit_triangle)
    save_mapped_circuit – Include .qasm string of the mapped circuit in result
    csv – Create CSV string for result
    statistics – Print statistics
    verbose – Print more detailed information during the mapping process
Returns:
    JSON object containing results
"""
def compile(circ, arch: Union[str, Arch],
            calibration = "",
            method: Method = Method.heuristic,
            initial_layout: InitialLayoutStrategy = InitialLayoutStrategy.dynamic,
            layering: LayeringStrategy = LayeringStrategy.individual_gates,
            save_mapped_circuit: bool = False,
            csv: bool = False,
            statistics: bool = False,
            verbose: bool = False
            ) -> object
```

Note that in order for the bindings to work the SMT Solver [Z3 >= 4.8.3](https://github.com/Z3Prover/z3) has to be installed on the system and the dynamic linker has to be able to find the library.
This can be accomplished in a multitude of ways:
- Under Ubuntu 20.04 and newer: `sudo apt-get install z3`
- Under macOS: `brew install z3`
- Alternatively: `pip install z3-solver` and then append the corresponding path to the library path (`LD_LIBRARY_PATH` under Linux, `DYLD_LIBRARY_PATH` under macOS), e.g. via
    ```bash
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(python -c "import z3; print(z3.__path__[0]+'/lib')")
    ```
- Download pre-built binaries from https://github.com/Z3Prover/z3/releases and copy the files to the respective system directories
- Build Z3 from source and install it to the system
### Command-line Executable
JKQ QMAP also provides two **standalone executables** with command-line interface called `qmap_heuristic` and `qmap_exact`.
They provide the same options as the Python module as flags (e.g., `--ps` for printing statistics). Per default, this produces JSON formatted output.

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

Both, the exact and the heuristic mapping tool also offer the `--layering` option, which allows to choose one of the following strategies for partitioning the circuit:
- `individual` (*default*): consider each gate separately,
- `disjoint`: consider gates acting on disjoint qubits as a layer,
- `odd`: group pairs of gates. (Note that this strategy was only tested for IBM QX4 with the exact mapping tool and may not work on different architectures)
- `triangle`: add gates to a layer, as long as no more than three qubits are involved. (Note that this strategy only works if the architecture's coupling map contains a triangle, e.g. IBM QX4, and was only tested using the exact mapping tool)

### System Requirements
Building (and running) is continuously tested under Linux, MacOS, and Windows using the [latest available system versions for GitHub Actions](https://github.com/actions/virtual-environments).
However, the implementation should be compatible with any current C++ compiler supporting C++14 and a minimum CMake version of 3.10.

`boost/program_options >= 1.50` is required for building the the commandline applications of the mapping tool.

In order to build the exact mapping tool, the SMT Solver [Z3 >= 4.8.3](https://github.com/Z3Prover/z3) has to be installed and on the path.

### Library Organisation
Internally the JKQ QMAP library works in the following way
- Import input file into a `qc::QuantumComputation` object
    ```c++
    qc::QuantumComputation qc{};
    std::string circ = "<PATH_TO_CIRCUIT_FILE>";
    qc.import(circ);
    ```
- Import architecture file into a `Architecture` object
    ```c++
    Architecture arch{};
    std::string cm = "<PATH_TO_ARCH_FILE>";
    arch.loadCouplingMap(cm);
    ```
- (Optional) Import calibration file into `arch` object
    ```c++
    std::string cal = "<PATH_TO_CAL_FILE>";
    arch.loadCalibrationData(cal);
    ```
- Depending on `Method`, instantiate a `HeuristicMapper` or `ExactMapper` object with the circuit and the architecture
    ```c++
    HeuristicMapper mapper(qc, arch);
    ```
  or
    ```c++ 
    ExactMapper mapper(qc, arch);
    ```
- Set configuration options, e.g.,
    ```c++
    MappingSettings ms{};
    ms.layeringStrategy = LayeringStrategy::DisjointQubits;
    ```
- Perform the actual mapping
    ```c++
    mapper.map(ms);
    ```
- Dump the mapped circuit
  ```c++
  mapper.dumpResult("<PATH_TO_OUTPUT_FILE>");
  ```
- Print the results (include statistics by setting `printStatistics=true`)
    ```c++
    mapper.printResult(std::cout, printStatistics);
    ```

### Configure, Build, and Install

In order to build the library execute the following in the project's main directory
1) Configure CMake
    ```commandline
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    ```
   Windows users using Visual Studio and the MSVC compiler may try
   ```commandline
   cmake -S . -B build -G "Visual Studio 16 2019" -A x64 -DCMAKE_BUILD_TYPE=Release
   ```
   Older CMake versions not supporting the above syntax (< 3.13) may be used with
   ```commandline
   mkdir build && cd build
   cmake .. -DCMAKE_BUILD_TYPE=Release
   ```
2) Build the respective target.
    ```commandline
   cmake --build ./build --config Release --target <target>
   ```
   The following CMake targets are available
    - `qmap_heuristic`: The heuristic mapper commandline executable (only available if Boost is found)
    - `qmap_exact`: The exact mapper commandline executable (only available if Boost and Z3 is found)
    - `qmap_heuristic_lib`: The standalone heuristic mapper library
    - `qmap_exact_lib`: The standalone exact mapper library (only available if Z3 is found)
    - `qmap_heuristic_test`: Unit tests for heuristic maper using GoogleTest
    - `qmap_exact_test`: Unit tests for exact maper using GoogleTest (only available if Z3 is found)

3) Optional: The QMAP library and tool may be installed on the system by executing

    ```commandline
    cmake --build ./build --config Release --target install
    ```

   It can then also be included in other projects using the following CMake snippet

    ```cmake
    find_package(qmap)
    target_link_libraries(${TARGET_NAME} PRIVATE JKQ::qmap)
    ```

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
