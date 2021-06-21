[![PyPI](https://img.shields.io/pypi/v/jkq.qmap?logo=pypi&style=plastic)](https://pypi.org/project/jkq.qmap/)
[![GitHub Workflow Status](https://img.shields.io/github/workflow/status/iic-jku/qmap/CI?logo=github&style=plastic)](https://github.com/iic-jku/qmap/actions?query=workflow%3A%22CI%22)
[![Codecov branch](https://img.shields.io/codecov/c/github/iic-jku/qmap/master?label=codecov&logo=codecov&style=plastic)](https://codecov.io/gh/iic-jku/qmap)
![GitHub](https://img.shields.io/github/license/iic-jku/qmap?style=plastic)
[![toolset: JKQ](https://img.shields.io/static/v1?label=toolset&message=JKQ&color=blue&style=plastic)](https://github.com/iic-jku/jkq)
[![arXiv](https://img.shields.io/static/v1?label=arXiv&message=1712.04722&color=inactive&style=plastic)](https://arxiv.org/abs/1712.04722)
[![arXiv](https://img.shields.io/static/v1?label=arXiv&message=2009.02376&color=inactive&style=plastic)](https://arxiv.org/abs/2009.02376)
[![arXiv](https://img.shields.io/static/v1?label=arXiv&message=2011.07314&color=inactive&style=plastic)](https://arxiv.org/abs/2011.07314)

# QMAP - A JKQ tool for Quantum Circuit Mapping written in C++
A [JKQ](https://github.com/iic-jku/jkq) tool for quantum circuit mapping by the [Institute for Integrated Circuits](https://iic.jku.at/eda/) at the [Johannes Kepler University Linz](https://jku.at) based on methods proposed in [[1]](https://iic.jku.at/files/eda/2018_tcad_efficient_mapping_of_quantum_circuits_to_ibm_qx_architectures.pdf), [[2]](https://iic.jku.at/files/eda/2019_dac_mapping_quantum_circuits_ibm_architectures_using_minimal_number_swap_h_gates.pdf), [[3]](https://iic.jku.at/files/eda/2021_aspdac_exploiting_teleportation_in_quantum_circuit_mappping.pdf).

[[1]](https://iic.jku.at/files/eda/2018_tcad_efficient_mapping_of_quantum_circuits_to_ibm_qx_architectures.pdf)
A. Zulehner, A. Paler, and R. Wille. An Efficient Methodology for Mapping Quantum Circuits to the IBM QX Architectures.
*IEEE Transactions on Computer Aided Design of Integrated Circuits and Systems (TCAD)*, 2018.

[[2]](https://iic.jku.at/files/eda/2019_dac_mapping_quantum_circuits_ibm_architectures_using_minimal_number_swap_h_gates.pdf)
R. Wille, L. Burgholzer, and A. Zulehner. Mapping Quantum Circuits to IBM QX Architectures
Using the Minimal Number of SWAP and H Operations. In *Design Automation Conference (DAC)*, 2019.

[[3]](https://iic.jku.at/files/eda/2021_aspdac_exploiting_teleportation_in_quantum_circuit_mappping.pdf)
S. Hillmich, A. Zulehner, and R. Wille. Exploiting Quantum Teleportation in Quantum Circuit Mapping. 
*In Asia and South Pacific Design Automation Conference (ASP-DAC)*, 2021.

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
- **Heuristic Mapper**:  Heuristic solution based on A* search. For details see [[1]](https://iic.jku.at/files/eda/2018_tcad_efficient_mapping_of_quantum_circuits_to_ibm_qx_architectures.pdf) and [[3]](https://iic.jku.at/files/eda/2021_aspdac_exploiting_teleportation_in_quantum_circuit_mappping.pdf).
- **Exact Mapper**: Exact solution utilizing the SMT Solver Z3. For details see [[2]](https://iic.jku.at/files/eda/2019_dac_mapping_quantum_circuits_ibm_architectures_using_minimal_number_swap_h_gates.pdf).

Note that, at the moment, circuits to be mapped are assumed to be already decomposed into elementary gates supported by the targeted device. 
More specifically, circuits must not contain gates acting on more than two qubits.

For more information, please visit [iic.jku.at/eda/research/ibm_qx_mapping/](https://iic.jku.at/eda/research/ibm_qx_mapping/).

If you have any questions, feel free to contact us via [iic-quantum@jku.at](mailto:iic-quantum@jku.at) or by creating an issue on [GitHub](https://github.com/iic-jku/qmap/issues).

## Usage

JKQ QMAP is mainly developed as a C++ library with a [commandline interface](#command-line-executable). 
However, using it in Python is as easy as
```bash
pip install jkq.qmap
```
and then in Python
```python
from jkq import qmap
qmap.compile(circ, arch, ...)
```
where `circ` is either a Qiskit `QuantumCircuit` object or the path to an input file (supporting various formats, such as `.qasm`, `.real`,...)
and `arch` is either one of the pre-defined architectures (see below) or the path to a file containing the number of qubits and a line-by-line enumeration of the qubit connections.

Architectures that are available per default (under `qmap.Arch.<...>`) include (corresponding files are available in `extern/architectures/`):
- `IBM_QX4` (5 qubit, directed bow tie layout)
- `IBM_QX5` (16 qubit, directed ladder layout)
- `IBMQ_Yorktown` (5 qubit, undirected bow tie layout)
- `IBMQ_London` (5 qubit, undirected T-shape layout)
- `IBMQ_Bogota` (5 qubit, undirected linear chain layout)
- `IBMQ_Tokyo` (20 qubit, undirected brick-like layout)

Whether the heuristic (*default*) or the exact mapper is used can be controlled by passing `method=qmap.Method.heuristic` or `method=qmap.Method.exact` to the `compile` function.

There are several options that can be passed to the `compile` function:
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
    use_teleportation - Use teleportation in addition to swaps
    teleportation_fake - Assign qubits as ancillary for teleportation but don't actually use them
    teleportationSeed - Fix a seed for the initial ancilla placement (0 means no fixed seed)
    save_mapped_circuit – Include .qasm string of the mapped circuit in result
    csv – Create CSV string for result
    statistics – Print statistics
    verbose – Print more detailed information during the mapping process
Returns:
    JSON object containing results
"""
def compile(circ, arch: Union[str, Arch],
            calibration: str = "",
            method: Method = Method.heuristic,
            initial_layout: InitialLayoutStrategy = InitialLayoutStrategy.dynamic,
            layering: LayeringStrategy = LayeringStrategy.individual_gates,
            use_teleportation: bool = False,
            teleportation_fake: bool = False,
            teleportation_seed: int = 0,
            save_mapped_circuit: bool = True,
            csv: bool = False,
            statistics: bool = False,
            verbose: bool = False
            ) -> Dict[str, Any]:
    ...
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
$ ./qmap_heuristic --in grover_2.qasm --out grover_2m.qasm --arch extern/architectures/ibm_london_5qubit.arch --ps
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
- `dynamic` (*default*): determine initial layout on demand during the mapping (this is the only one compatible with teleportation).

Both, the exact and the heuristic mapping tool also offer the `--layering` option, which allows to choose one of the following strategies for partitioning the circuit:
- `individual` (*default*): consider each gate separately,
- `disjoint`: consider gates acting on disjoint qubits as a layer,
- `odd`: group pairs of gates. (Note that this strategy was only tested for IBM QX4 with the exact mapping tool and may not work on different architectures)
- `triangle`: add gates to a layer, as long as no more than three qubits are involved. (Note that this strategy only works if the architecture's coupling map contains a triangle, e.g. IBM QX4, and was only tested using the exact mapping tool)

### System Requirements
Building (and running) is continuously tested under Linux, MacOS, and Windows using the [latest available system versions for GitHub Actions](https://github.com/actions/virtual-environments).
However, the implementation should be compatible with any current C++ compiler supporting C++17 and a minimum CMake version of 3.14.

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

To start off, clone this repository using
```shell
$ git clone --recurse-submodules -j8 https://github.com/iic-jku/qmap 
```
Note the `--recurse-submodules` flag. It is required to also clone all the required submodules. 
If you happen to forget passing the flag on your initial clone, you can initialize all the submodules by executing `git submodule update --init --recursive` in the main project directory.

Our projects use CMake as the main build configuration tool. Building a project using CMake is a two-stage process. First, CMake needs to be *configured* by calling
```shell 
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```
This tells CMake to search the current directory `.` (passed via `-S`) for a *CMakeLists.txt* file and process it into a directory `build` (passed via `-B`).
The flag `-DCMAKE_BUILD_TYPE=Release` tells CMake to configure a *Release* build (as opposed to, e.g., a *Debug* build).

After configuring with CMake, the project can be built by calling
```shell
$ cmake --build build --config Release
```
This tries to build the project in the `build` directory (passed via `--build`).
Some operating systems and developer environments explicitly require a configuration to be set, which is why the `--config` flag is also passed to the build command. The flag `--parallel <NUMBER_OF_THREADS>` may be added to trigger a parallel build.

Building the project this way generates
- the heuristic library `libqmap_heuristic_lib.a` (Unix) / `qmap_heuristic_lib.lib` (Windows) in the `build/src` directory
- the heuristic mapper commandline executable `qmap_heuristic` in the `build/apps` directory (only available if Boost is found)
- a test executable `qmap_heuristic_test` containing a small set of unit tests for the heuristic mapper in the `build/test` directory
- the exact library `libqmap_exact_lib.a` (Unix) / `qmap_exact_lib.lib` (Windows) in the `build/src` directory (only available if Z3 is found)
- the exact mapper commandline executable `qmap_exact` in the `build/apps` directory (only available if Boost and Z3 is found)
- a test executable `qmap_exact_test` containing a small set of unit tests for the exact mapper in the `build/test` directory (only available if Z3 is found)

### Extending the Python Bindings

To extend the Python bindings you can locally install the package in edit mode, so that changes in the Python code are instantly available.
The following example assumes you have a [virtual environment](https://docs.python.org/3/library/venv.html) set up and activated.

```commandline
(venv) $ pip install cmake
(venv) $ pip install --editable .
```

If you change parts of the C++ code, you have to run the second line to make the changes visible in Python.

## Reference

If you use our tool for your research, we will be thankful if you refer to it by citing the appropriate publications.

For the [heuristic mapping](https://dblp.org/rec/journals/tcad/ZulehnerPW19.html?view=bibtex), please cite
```bibtex
@article{DBLP:journals/tcad/ZulehnerPW19,
  author    = {Alwin Zulehner and Alexandru Paler and Robert Wille},
  title     = {An Efficient Methodology for Mapping Quantum Circuits to the {IBM QX} Architectures},
  journal   = {{IEEE} Transactions on Computer-Aided Design of Integrated Circuits and Systems},
  volume    = {38},
  number    = {7},
  pages     = {1226--1236},
  year      = {2019}
}
```

For the [teleportation in the heuristic mapping](https://dblp.org/rec/conf/aspdac/HillmichZW21.html?view=bibtex), please cite
```bibtex
@inproceedings{DBLP:conf/aspdac/HillmichZW21,
  author    = {Stefan Hillmich and Alwin Zulehner and Robert Wille},
  title     = {Exploiting Quantum Teleportation in Quantum Circuit Mapping},
  booktitle = {Asia and South Pacific Design Automation Conference},
  pages     = {792--797},
  publisher = {{ACM}},
  year      = {2021}
}
```

For the [exact mapping](https://dblp.org/rec/conf/dac/WilleBZ19.html?view=bibtex), please cite
```bibtex
@inproceedings{DBLP:conf/dac/WilleBZ19,
  author    = {Robert Wille and Lukas Burgholzer and Alwin Zulehner},
  title     = {Mapping Quantum Circuits to {IBM QX} Architectures Using the Minimal Number of {SWAP} and {H} Operations},
  booktitle = {Design Automation Conference},
  publisher = {{ACM}},
  year      = {2019}
}
````
