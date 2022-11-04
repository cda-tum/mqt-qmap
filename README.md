[![PyPI](https://img.shields.io/pypi/v/mqt.qmap?logo=pypi&style=flat-square)](https://pypi.org/project/mqt.qmap/)
![OS](https://img.shields.io/badge/os-linux%20%7C%20macos%20%7C%20windows-blue?style=flat-square)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square)](https://opensource.org/licenses/MIT)
[![CI](https://img.shields.io/github/workflow/status/cda-tum/qmap/C++?style=flat-square&logo=github&label=c%2B%2B)](https://github.com/cda-tum/qmap/actions/workflows/ci.yml)
[![Python-CI](https://img.shields.io/github/workflow/status/cda-tum/qmap/Python?style=flat-square&logo=github&label=python)](https://github.com/cda-tum/qmap/actions/workflows/python-ci.yml)
[![Bindings](https://img.shields.io/github/workflow/status/cda-tum/qmap/Python%20Packaging?style=flat-square&logo=github&label=packaging)](https://github.com/cda-tum/qmap/actions/workflows/deploy.yml)
[![Documentation](https://img.shields.io/readthedocs/mqtqmap?logo=readthedocs&style=flat-square)](https://mqtqmap.readthedocs.io/en/latest/)
[![codecov](https://img.shields.io/codecov/c/github/cda-tum/qmap?style=flat-square&logo=codecov)](https://codecov.io/gh/cda-tum/qmap)

# MQT QMAP - A tool for Quantum Circuit Compilation

A tool for quantum circuit compilation developed as part of the _Munich Quantum Toolkit_ (_MQT_)[^1] by the [Chair for Design Automation](https://www.cda.cit.tum.de/) at the [Technical University of Munich](https://www.tum.de/).
It builds upon [our quantum functionality representation (QFR)](https://github.com/cda-tum/qfr).
.

<p align="center">
  <a href="https://mqtqmap.readthedocs.io/en/latest/">
  <img width=30% src="https://img.shields.io/badge/documentation-blue?style=for-the-badge&logo=read%20the%20docs" alt="Documentation" />
  </a>
</p>

If you have any questions, feel free to contact us via [quantum.cda@xcit.tum.de](mailto:quantum.cda@xcit.tum.de) or by creating an issue on [GitHub](https://github.com/cda-tum/qmap/issues).

## Getting Started

QMAP is available via [PyPI](https://pypi.org/project/mqt.qmap/) for Linux, macOS, and Windows.

```console
(venv) $ pip install mqt.qmap
```

Compiling a given quantum circuit to a certain device is as easy as

```python3
from mqt import qmap
from qiskit import QuantumCircuit
from qiskit.providers.fake_provider import FakeLondon

circ = QuantumCircuit(3)
circ.h(0)
circ.cx(0, 1)
circ.cx(0, 2)

circ_mapped, results = qmap.compile(circ, arch=FakeLondon())
```

**Detailed documentation on all available methods, options, and input formats is available at [ReadTheDocs](https://mqtqmap.readthedocs.io/en/latest/).**

## System Requirements and Building

The implementation is compatible with any C++17 compiler and a minimum CMake version of 3.19.
Please refer to the [documentation](https://mqtqmap.readthedocs.io/en/latest/) on how to build the project.

Building (and running) is continuously tested under Linux, macOS, and Windows using the [latest available system versions for GitHub Actions](https://github.com/actions/virtual-environments).

## References

QMAP has been developed based on methods proposed in the following papers:

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

[[5]](https://arxiv.org/pdf/2210.09321.pdf)
T. Peham, L. Burgholzer, and R. Wille. On Optimal Subarchitectures for Quantum Circuit Mapping.
_arXiv:2210.09321_, 2022.

[^1]: The Munich Quantum Toolkit was formerly known under the acronym _JKQ_ and developed by the [Institute for Integrated Circuits](https://iic.jku.at/eda/) at the [Johannes Kepler University Linz](https://jku.at)).
