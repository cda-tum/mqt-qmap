---
file_format: mystnb
kernelspec:
  name: python3
mystnb:
  number_source_lines: true
---

```{code-cell} ipython3
:tags: [remove-cell]
%config InlineBackend.figure_formats = ['svg']
```

<style>.widget-subarea{display:none;} /*hide widgets as they do not work with sphinx*/</style>

# Zoned Neutral Atom Compiler

Successful quantum computation requires advanced software, especially compilers that optimize quantum algorithms for
hardware execution.
Zoned neutral atom architectures execute operations in designated spatially separated zones.
The zones facilitate higher coherence times and overall fidelity of the quantum computation.
However, the zones also require the rearrangement of atoms during the quantum computation.
This tool provides a routing-aware placement method to improve the efficiency of atom rearrangements during quantum
computation.

## Example: GHZ State on Neutral Atom Architecture

In this example, we will demonstrate how to use the zoned neutral atom compiler to generate a sequence of
target-specific instructions for a quantum circuit.
The circuit prepares a GHZ state of 8 qubits using the native gates of the zoned neutral atom architecture, i.e., global
ry-gates, local rz-gates, and controlled z-gates.

```{note}
Even though other one-qubit gates may not be supported by the hardware, the compiler can handle arbitrary one-qubit gates.
It will translate them to generic u3 gates and include them in the output.
```

```{code-cell} ipython3
from qiskit import QuantumCircuit
from numpy import pi

def global_ry(theta, num_qubits):
    """:returns: a global ry gate"""
    qc = QuantumCircuit(num_qubits)
    qc.ry(theta, range(num_qubits))
    return qc.to_gate(label = f"Ry({theta:.3f})")


qc = QuantumCircuit(8)
qc.append(global_ry(-pi/4, 8), range(8))
qc.z(range(8))
qc.append(global_ry(pi/4, 8), range(8))
qc.cz(0, 4)
qc.append(global_ry(-pi/4, 8), range(8))
qc.z(4)
qc.append(global_ry(pi/4, 8), range(8))
qc.cz(0, 2)
qc.cz(4, 6)
qc.append(global_ry(-pi/4, 8), range(8))
qc.z([2, 6])
qc.append(global_ry(pi/4, 8), range(8))
qc.cz(0, 1)
qc.cz(2, 3)
qc.cz(4, 5)
qc.cz(6, 7)
qc.append(global_ry(-pi/4, 8), range(8))
qc.z([1,3,5,7])
qc.append(global_ry(pi/4, 8), range(8))

qc.draw(output="mpl")
```

On the considered architecture, the one-qubit gates, i.e., the global ry and local rz-gates can be executed everywhere.
However, the controlled z-gates can only be executed between nearby atoms in the so-called entanglement zone.
This entanglement zone is spatially separated from the storage zone, where all atoms not involved in a cz-gate are
stored.
To find an optimized sequence of target-specific instructions, we use the zoned neutral atom compiler.
This compiler requires first a specification of the architecture.

```{code-cell} ipython3
from json import loads as parse_json
from mqt.qmap.na.azac import AZACArchitecture

arch = AZACArchitecture(parse_json("""{
  "name": "Architecture with one entanglement and one storage zone",
  "operation_duration": {"rydberg": 0.36, "1qGate": 52, "atom_transfer": 15},
  "operation_fidelity": {"two_qubit_gate": 0.995, "single_qubit_gate": 0.9997, "atom_transfer": 0.999},
  "qubit_spec": {"T": 1.5e6},
  "storage_zones": [{
    "zone_id": 0,
    "slms": [{"id": 0, "site_separation": [3, 3], "r": 20, "c": 100, "location": [0, 0]}],
    "offset": [0, 0],
    "dimension": [297, 57]
  }],
  "entanglement_zones": [{
    "zone_id": 0,
    "slms": [
      {"id": 1, "site_separation": [12, 10], "r": 7, "c": 20, "location": [35, 67]},
      {"id": 2, "site_separation": [12, 10], "r": 7, "c": 20, "location": [37, 67]}
    ],
    "offset": [35, 67],
    "dimension": [230, 60]
  }],
  "aods": [{"id": 0, "site_separation": 2, "r": 100, "c": 100}],
  "arch_range": [[0, 0], [297, 162]],
  "rydberg_range": [[[30, 62], [270, 132]]]
}"""))
```

Furthermore, the different stages of the compiler can be configured with a set of parameters.

```{code-cell} ipython3
from mqt.qmap.na.azac import AZACompiler

compiler = AZACompiler(arch, parse_json("""{
  "code_generator": {
    "parking_offset": 1,
    "warn_unsupported_gates": true
  },
  "a_star_placer" : {
    "use_window" : true,
    "window_min_width" : 8,
    "window_ratio" : 1.5,
    "window_share" : 0.6,
    "deepening_factor" : 0.6,
    "deepening_value" : 0.2,
    "lookahead_factor": 0.2,
    "reuse_level": 5.0,
    "max_nodes": 50000000
  }
}"""))
```

Now, the created compiler can be used to compile the circuit from above.
The output is in the `.naviz` format that can be read by the `NAViz` tool
at [github.com/cda-tum/mqt-naviz](https://github.com/cda-tum/mqt-naviz).
This tool allows visualizing the resulting quantum computation.

```{code-cell} ipython3
from mqt.core import load
circ = load(qc)
code = compiler.compile(circ)
print(code)
```
