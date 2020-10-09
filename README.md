[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![toolset: JKQ](https://img.shields.io/badge/toolset-JKQ-blue)](https://github.com/iic-jku/jkq)

# JKQ QMAP - A tool for mapping quantum circuits to quantum architectures (developed in C++)
A tool for quantum circuit simulation by the [Institute for Integrated Circuits](https://iic.jku.at/eda/) at the [Johannes Kepler University Linz](https://jku.at) 
and a part of the [JKQ toolset](https://github.com/iic-jku/jkq).

For more information, please visit [iic.jku.at/eda/research/ibm_qx_mapping/](https://iic.jku.at/eda/research/ibm_qx_mapping/).

If you have any questions, feel free to contact us via [iic-quantum@jku.at](mailto:iic-quantum@jku.at) or by creating an [issue](https://github.com/iic-jku/qmap/issues) on GitHub.

## Usage

```commandline
$ ./heuristic_app --in grover_2.qasm --out grover_2m.qasm --arch ibmq_london.cm --ps
{
	"circuit": {
		"name": "grover_2",
		"n_qubits": 3,
		"n_gates": 30,
	},
	"mapped_circuit": {
		"name": "grover_2m",
		"n_qubits": 5,
		"n_gates": 33,
	},
	"statistics": {
		"mapping_time": 0.001638,
		"additional_gates": 3,
		"method": "heuristic",
		"arch": "ibmq_london.cm"
	}
}
```

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

Afterwards `./build/exact_app` and `./build/heuristic_app` are available.

### System Requirements
Building (and running) should work under Linux, MacOS, and Windows with any current C++ compiler supporting C++14 and a minimum CMake version of 3.10.

## Reference

If you use our tool for your research, we will be thankful if you refer to it by citing the one or both of following publications.

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
