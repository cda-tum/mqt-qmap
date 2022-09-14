[![PyPI](https://img.shields.io/pypi/v/mqt.qmap?logo=pypi&style=flat-square)](https://pypi.org/project/mqt.qmap/)
![OS](https://img.shields.io/badge/os-linux%20%7C%20macos%20%7C%20windows-blue?style=flat-square)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square)](https://opensource.org/licenses/MIT)
[![CI](https://img.shields.io/github/workflow/status/cda-tum/qmap/C++?style=flat-square&logo=github&label=c%2B%2B)](https://github.com/cda-tum/qmap/actions/workflows/ci.yml)
[![Python-CI](https://img.shields.io/github/workflow/status/cda-tum/qmap/Python?style=flat-square&logo=github&label=python)](https://github.com/cda-tum/qmap/actions/workflows/python-ci.yml)
[![Bindings](https://img.shields.io/github/workflow/status/cda-tum/qmap/Python%20Packaging?style=flat-square&logo=github&label=packaging)](https://github.com/cda-tum/qmap/actions/workflows/deploy.yml)
[![codecov](https://img.shields.io/codecov/c/github/cda-tum/qmap?style=flat-square&logo=codecov)](https://codecov.io/gh/cda-tum/qmap)

# MQT QMAP - A tool for Quantum Circuit Mapping written in C++

A tool for quantum circuit mapping developed by the [Chair for Design Automation](https://www.cda.cit.tum.de/) at the [Technical University of Munich](https://www.tum.de/) based on methods proposed
in [[1]](https://www.cda.cit.tum.de/files/eda/2018_tcad_efficient_mapping_of_quantum_circuits_to_ibm_qx_architectures.pdf)
, [[2]](https://www.cda.cit.tum.de/files/eda/2019_dac_mapping_quantum_circuits_ibm_architectures_using_minimal_number_swap_h_gates.pdf)
, [[3]](https://www.cda.cit.tum.de/files/eda/2021_aspdac_exploiting_teleportation_in_quantum_circuit_mappping.pdf)
, [[4]](https://www.cda.cit.tum.de/files/eda/2022_aspdac_limiting_search_space_optimal_quantum_circuit_mapping.pdf).

[[1]](https://www.cda.cit.tum.de/files/eda/2018_tcad_efficient_mapping_of_quantum_circuits_to_ibm_qx_architectures.pdf)
A. Zulehner, A. Paler, and R. Wille. An Efficient Methodology for Mapping Quantum Circuits to the IBM QX Architectures.
_IEEE Transactions on Computer Aided Design of Integrated Circuits and Systems (TCAD)_, 2018.

[[2]](https://www.cda.cit.tum.de/files/eda/2019_dac_mapping_quantum_circuits_ibm_architectures_using_minimal_number_swap_h_gates.pdf)
R. Wille, L. Burgholzer, and A. Zulehner. Mapping Quantum Circuits to IBM QX Architectures Using the Minimal Number of SWAP and H Operations. In _Design Automation Conference (DAC)_, 2019.

[[3]](https://www.cda.cit.tum.de/files/eda/2021_aspdac_exploiting_teleportation_in_quantum_circuit_mappping.pdf)
S. Hillmich, A. Zulehner, and R. Wille. Exploiting Quantum Teleportation in Quantum Circuit Mapping.
_In Asia and South Pacific Design Automation Conference (ASP-DAC)_, 2021.

[[4]](https://www.cda.cit.tum.de/files/eda/2022_aspdac_limiting_search_space_optimal_quantum_circuit_mapping.pdf)
L. Burgholzer, S. Schneider, and R. Wille. Limiting the Search Space in Optimal Quantum Circuit Mapping.
_In Asia and South Pacific Design Automation Conference (ASP-DAC)_, 2022.

QMAP is part of the Munich Quantum Toolkit (MQT; formerly known as JKQ and developed by the [Institute for Integrated Circuits](https://iic.jku.at/eda/) at the [Johannes Kepler University Linz](https://jku.at)).
It builds upon [our quantum functionality representation (QFR)](https://github.com/cda-tum/qfr) and can be used for mapping quantum circuits in any of the following formats:

- `QuantumCircuit` object from IBM's [Qiskit](https://github.com/Qiskit/qiskit) (only through the MQT QMAP Python bindings)
- `OpenQASM` (e.g. used by IBM's [Qiskit](https://github.com/Qiskit/qiskit)),
- `Real` (e.g. from [RevLib](http://revlib.org)),
- `TFC` (e.g. from [Reversible Logic Synthesis Benchmarks Page](http://webhome.cs.uvic.ca/~dmaslov/mach-read.html))
- `QC` (e.g. from [Feynman](https://github.com/meamy/feynman))

to any given architecture with the following available methods:

- **Heuristic Mapper**: Heuristic solution based on A\* search. For details see [[1]](https://www.cda.cit.tum.de/files/eda/2018_tcad_efficient_mapping_of_quantum_circuits_to_ibm_qx_architectures.pdf)
  and [[3]](https://www.cda.cit.tum.de/files/eda/2021_aspdac_exploiting_teleportation_in_quantum_circuit_mappping.pdf).
- **Exact Mapper**: Exact solution utilizing the SMT Solver Z3. For details see [[2]](https://www.cda.cit.tum.de/files/eda/2019_dac_mapping_quantum_circuits_ibm_architectures_using_minimal_number_swap_h_gates.pdf)
  and [[4]](https://www.cda.cit.tum.de/files/eda/2022_aspdac_limiting_search_space_optimal_quantum_circuit_mapping.pdf).

Note that, at the moment, circuits to be mapped are assumed to be already decomposed into elementary gates supported by the targeted device. More specifically, circuits must not contain gates acting on more than two qubits.

For more information, please visit [cda.cit.tum.de/research/ibm_qx_mapping/](https://www.cda.cit.tum.de/research/ibm_qx_mapping/).

If you have any questions, feel free to contact us via [quantum.cda@xcit.tum.de](mailto:quantum.cda@xcit.tum.de) or by creating an issue on [GitHub](https://github.com/cda-tum/qmap/issues).

## Usage

MQT QMAP is developed as a C++ library with an easy to use Python interface.

- In order to make the library as easy to use as possible (without compilation), we provide pre-built wheels for most common platforms (64-bit Linux, MacOS, Windows). These can be installed using
  ```bash
  pip install mqt.qmap
  ```
  However, in order to get the best performance out of QMAP, it is recommended to build it locally from the source distribution (see [system requirements](#system-requirements)) via
  ```bash
  pip install  mqt.qmap --no-binary mqt.qmap
  ```
  This enables platform specific compiler optimizations that cannot be enabled on portable wheels.
- Once installed, start using it in Python:

  ```python
  from mqt import qmap

  circ_mapped, results = qmap.compile(circ, arch)
  ```

where `circ` is either a Qiskit `QuantumCircuit` object or the path to an input file (in any of the formats listed above)
and `arch` is either

- a Qiskit `Backend` instance such as those defined under `qiskit.providers.fake_provider` **(recommended)**,
- one of the pre-defined architectures (see below), or
- the path to a file containing the number of qubits and a line-by-line enumeration of the qubit connections.

Architectures that are available per default (either as strings or under `qmap.Arch.<...>`) include:

- `IBM_QX4` (5 qubit, directed bow tie layout)
- `IBM_QX5` (16 qubit, directed ladder layout)
- `IBMQ_Yorktown` (5 qubit, undirected bow tie layout)
- `IBMQ_London` (5 qubit, undirected T-shape layout)
- `IBMQ_Bogota` (5 qubit, undirected linear chain layout)
- `IBMQ_Casablanca` (7 qubit, undirected H-shape layout)
- `IBMQ_Tokyo` (20 qubit, undirected brick-like layout)
- `Rigetti_Agave` (8 qubit, undirected ring layout)
- `Rigetti_Aspen` (16 qubit, undirected dumbbell layout)

Whether the heuristic (_default_) or the exact mapper is used can be controlled by passing `method="heuristic"` or `method="exact"` to the `compile` function.

There are several configuration options that can be passed to the `compile` function:

- The heuristic mapper offers the `initial_layout` option, which allows to choose one of the following strategies for
  choosing an initial layout:

  - `identity`: map logical qubit q_i to physical qubit Q_i,
  - `static`: determine fixed initial layout statically at the start of mapping,
  - `dynamic` (_default_): determine initial layout on demand during the mapping (this is the only one compatible with
    teleportation).

- Both, the exact and the heuristic mapper also offer the `layering` option, which allows to choose one of the following
  strategies for partitioning the circuit:

  - `individual_gates` (_default_): consider each gate separately,
  - `disjoint_qubits`: consider gates acting on disjoint qubits as a layer,
  - `odd_gates`: group pairs of gates. (Note that this strategy was only tested for IBM QX4 with the exact mapping tool
    and may not work on different architectures)
  - `qubit_triangle`: add gates to a layer, as long as no more than three qubits are involved. (Note that this strategy
    only works if the architecture's coupling map contains a triangle, e.g. IBM QX4, and was only tested using the exact
    mapping tool)

- The exact mapper offers the `encoding` option, which allows to choose a different encoding for at-most-one and
  exactly-one constraints:

  - `naive` (_default_): use naive encoding for constraints
  - `commander`: use commander encoding for at-most-one and exactly-one constraints
  - `bimander`: use bimander encoding for at-most-one and commander for exactly-one constraints

  As the commander encoding can use different strategies to group the variables, there are
  different `commander_grouping` options:

  - `halves` (_default_): each group contains half of the total variables
  - `logarithm`: each group contains at most log2 of the total variables
  - `fixed2`: each group contains exactly two variables
  - `fixed3`: each group contains exactly three variables

- Per default, the exact mapper searches for a suitable mapping by considering every possible (connected) subset of
  qubits instead of the whole architecture at once. This can be disabled by setting `use_subsets=False`.

- The exact mapper also offers the `swap_reduction` option to enable limiting the number of swaps considered per layer (
  as proposed
  in [[4]](https://www.cda.cit.tum.de/files/eda/2022_aspdac_limiting_search_space_optimal_quantum_circuit_mapping.pdf)
  , which offers the following options:

  - `none`: consider whole search space
  - `coupling_limit` (_default_): calculate the max swaps per layer based on the longest path of current choice of
    qubits, or if `use_subsets` is disabled considers the whole architecture
  - `increasing`: start with 0 swaps and geometrically increase the number of swaps per layer
  - `custom`: set a custom limit, needs the `swap_limit` option to set the limit

  Using the `use_bdd` option, the mapping utilizes BDDs instead of simply removing the permutations from the core
  routine. This option is not generally advised, as it is more resource intensive in most cases, but is something to try
  in
  cases of timeout.

### Command-line Executable

QMAP also provides two **standalone executables** with command-line interface called `qmap_heuristic` and `qmap_exact`. They provide the same options as the Python module as flags. Per default, this produces JSON formatted output. In
general, we recommend to use the Python approach described above, as these commandline executables might not be maintained in the future.

### System Requirements

Building (and running) is continuously tested under Linux, MacOS, and Windows using the [latest available system versions for GitHub Actions](https://github.com/actions/virtual-environments). However, the implementation should be compatible
with any current C++ compiler supporting C++17 and a minimum CMake version of 3.18.

`boost/program_options >= 1.50` is required for building the commandline applications of the mapping tool.

In order to build the exact mapping tool and for the Python bindings to work, the SMT Solver [Z3 >= 4.8.15](https://github.com/Z3Prover/z3) has to be installed and the dynamic linker has to be able to find the library. This can be
accomplished in a multitude of ways:

- Under macOS: `brew install z3`
- Alternatively: `pip install z3-solver` in the environment used for development
- Download pre-built binaries from https://github.com/Z3Prover/z3/releases and copy the files to the respective system directories
- Build Z3 from source and install it to the system

### Library Organisation

Internally the MQT QMAP library works in the following way

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
  arch.loadProperties(cal);
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
  Configuration config{};
  config.layeringStrategy = Layering::DisjointQubits;
  ```
- Perform the actual mapping
  ```c++
  mapper.map(config);
  ```
- Dump the mapped circuit
  ```c++
  mapper.dumpResult("<PATH_TO_OUTPUT_FILE>");
  ```
- Print the results
  ```c++
  mapper.printResult(std::cout);
  ```

### Configure, Build, and Install

To start off, clone this repository using

```shell
git clone --recurse-submodules -j8 https://github.com/cda-tum/qmap
```

Note the `--recurse-submodules` flag. It is required to also clone all the required submodules.
If you happen to forget passing the flag on your initial clone, you can initialize all the submodules by
executing `git submodule update --init --recursive` in the main project directory.

Our projects use CMake as the main build configuration tool. Building a project using CMake is a two-stage process.
First, CMake needs to be _configured_ by calling

```shell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

This tells CMake to search the current directory `.` (passed via `-S`) for a _CMakeLists.txt_ file and process it into a
directory `build` (passed via `-B`).
The flag `-DCMAKE_BUILD_TYPE=Release` tells CMake to configure a _Release_ build (as opposed to, e.g., a _Debug_ build).

After configuring with CMake, the project can be built by calling

```shell
cmake --build build --config Release
```

This tries to build the project in the `build` directory (passed via `--build`).
Some operating systems and developer environments explicitly require a configuration to be set, which is why
the `--config` flag is also passed to the build command. The flag `--parallel <NUMBER_OF_THREADS>` may be added to
trigger a parallel build.

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
```

For the [search space limitation in the exact mapping](https://www.cda.cit.tum.de/files/eda/2022_aspdac_limiting_search_space_optimal_quantum_circuit_mapping.pdf), please cite

```bibtex
@inproceedings{burgholzer2022limitingSearchSpace,
  author    = {Lukas Burgholzer and Sarah Schneider and Robert Wille},
  title     = {Limiting the Search Space in Optimal Quantum Circuit Mapping},
  booktitle = {Asia and South Pacific Design Automation Conference},
  year      = {2022}
}
```
