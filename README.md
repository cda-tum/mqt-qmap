# A tool to map quantum circuits to the IBM QX architecture (developed in C++)
Copyright (c) 2017 by Alwin Zulehner (alwin.zulehner@jku.at), Stefan Hillmich (stefan.hillmich@jku.at), Alexandru Paler (alexandru.paler@jku.at), and Robert Wille (robert.wille@jku.at)  
Johannes Kepler University Linz, Austria  
http://www.jku.at/iic/eda/ibm_qx_mapping

The software is intellectual property of the above mentioned authors. You can freely redistribute this software for non-commercial purposes as long as it includes a reference to its origin (e.g. by referring to the corresponding paper cited below).

Use at your own risk!
In no event shall the authors be liable for any damages whatsoever (including without limitation damages for loss of business profits, business interruption, loss of business information, or any other pecuniary loss) arising from the use of or inability to use the software, even if the authors have been advised of the possibility of such damages.

If you have any questions feel free to contact us using alwin.zulehner@jku.at, stefan.hillmich@jku.at, alexandru.paler@jku.at, and robert.wille@jku.at

## Overview

In March 2017, IBM launched the project IBM Q with the goal to provide access to quantum computers for a broad audience. This allowed users to conduct quantum experiments on a 5-qubit and, since June 2017, also on a 16-qubit quantum computer (called IBM QX2 and IBM QX3, respectively). Later, further architectures have been added. In order to use these, the desired quantum functionality (e.g. provided in terms of a quantum circuit) has to properly be mapped so that the underlying physical constraints are satisfied – a complex task. This demands for solutions to automatically and efficiently conduct this mapping process. Here, we propose such an approach which satisfies all constraints given by the architecture and, at the same time, aims to keep the overhead in terms of additionally required quantum gates minimal. The proposed approach is generic and can easily be configured for future architectures. Experimental evaluations show that the proposed approach clearly outperforms IBM’s own mapping solution. In fact, for many quantum circuits, the proposed approach determines a mapping to the IBM architecture within less than five minutes (and within a fraction of a second in most cases), while IBM’s solution suffers from long runtimes and runs into a timeout of 1 hour in several cases. As an additional benefit, the proposed approach yields mapped circuits with smaller costs (i.e. fewer additional gates are required).

## Usage

### System Requirements

The package has been tested under Linux (Ubuntu 17.04, 64-bit) and should be compatible with any current version of g++/cmake. No additional packages are required.

### Build and Run

To build the quantum simulator type:

```commandline
mkdir build
cd build 
cmake ..
cmake --build .
cd ..
./build/ibm_qx_mapping "examples/4gt11_84.qasm" "mapped.qasm"
```


The mapping of a quantum circuit can be conducted as follows:

`./build/ibm_qx_mapping <input_file> <output_file>` maps the circuit `<input_file>` (given in the OpenQASM 2.0 format) to the IBM QX5 quantum processor.
Note that this implementation contains only one certain aspect of the mapping procedure, namely satisfying the architectural constraints.
Therefore, it assumes that the circuit is already decomposed into elementary operations.
The resulting circuit is written to `<output_file>` and can then be executed on the IBM QX architecture.
Our implementation does not perform post mapping optimization (as done e.g. in IBM's Python SDK).
However, they can be easily conducted by passing the resulting circuit to IBM's SDK. 

## Reference

If you use out mapping algorithm for your research, we would be thankful if you referred to it by citing the following publication: 

```
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

## A tool to map quantum circuits to the IBM QX architecture using the minimal number of SWAP and H operations (developed in C++)
Developers: Robert Wille, Lukas Burgholzer and Alwin Zulehner

For more information, please visit http://iic.jku.at/eda/research/ibm_qx_mapping/

If you have any questions, feel free to contact us via iic_quantum@jku.at or by creating an issue on GitHub.

## Overview

The recent progress in the physical realization of quantum computers (the first publicly available ones—IBM’s QX architectures—have been launched in 2017) has motivated research on automatic methods that aid users in running quantum circuits on them. Here, certain physical constraints given by the architectures which restrict the allowed interactions of the involved qubits have to be satisfied.
 
Thus far, this has been addressed by inserting SWAP and H operations. However, it remains unknown whether existing methods add a minimum number of SWAP and H operations or, if not, how far they are away from that minimum—an NP-complete problem. 
 
This tool addresses this by formulating the mapping task as a symbolic optimization problem that is solved using the Z3 theorem solver. While only applicable for circuits with a small number of qubits, experimental evaluations show that the number of operations added by IBM’s heuristic solution exceeds the lower bound obtained with the proposed approach by more than 100% on average.

## Usage

### System Requirements 

The package has been tested under Linux (Ubuntu 18.10, 64-bit) and should be compatible with any current version of g++/cmake.
Additionally the following packages are required:
* The Z3 Theorem Solver: https://github.com/Z3Prover/z3
* qelib.inc Gate Library: included with this package
  
### Build and Run 

To build the mapping tool type

    mkdir build
    cd build 
    cmake ..
    make -j4
    cd ..
    
It is important to include the `qelib1.inc` file in the build directory in order to allow the parser to find the gate declarations.

    cp ./qelib.inc ./build/

Executing 

    ./build/minimal_ibm_qx_mapping
     
 runs the mapping task for the provided test circuit
mapped to IBM's QX4; showing the results for all available strategies. 

To run your own computations for mapping to IBM's QX4, either execute 

    ./build/minimal_ibm_qx_mapping <input file>

 or specify the filename in ``main.cpp``.

The strategy to use can be configured in ``main.cpp``.

To run the mapping task for a different architecture, specify the architectures' coupling map, as well as the correspondingly required number of SWAP gates for each permutation in ``main.cpp``.

## Reference

If you use our mapping algorithm for your research, we would be thankful if you referred to it
by citing the following publication: 
````
@inproceedings{wille2019mapping,
    title={Mapping Quantum Circuits to {IBM QX} Architectures Using the Minimal Number of {SWAP} and {H} Operations},
    author={Wille, Robert and Burgholzer, Lukas and Zulehner, Alwin},
    booktitle={Design Automation Conference},
    year={2019}
}
````

