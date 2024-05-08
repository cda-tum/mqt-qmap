# Hybrid Neutral Atom Mapper

This file gives some "big picture" information about the Hybrid Neutral Atom Mapper (NAMAP) project.
The complete mapping process has the following structure:

1. Neutral Atom Mapper
   1. Neutral Atom Architecture: Reads a json file containing hardware/architecture information and provides methods to access this information.
   2. Neutral Atom Layer: Manages the creation of the front and lookahead layer, given a circuit.
   3. Mapping Loop: The gates are assigned to gate-based or shuttling-based mapping. Then suitable swaps or shuttling moves are found and applied.
2. AOD Scheduler: The abstract move operations are bundled into parallelize move groups and then converted into native AOD sequences (activate, move, deactivate).
3. Neutral Atom Scheduler: Aligns gates in a as-soon-as-possible fashion, taking into account the blocking constraint of multi-qubit gates.

## Neutral Atom Mapper

### Neutral Atom Architecture

The architecture has two main tasks:

- Provide information about the hardware, such as the gate times etc.
- Precompute the connectivity graph regarding the coordinates
  It represents the possible SLM traps as coordinates, either by a row-first running index or by a 2D coordinate tuple.

### Neutral Atom Layer

The layer is the interface between the mapper and the quantum circuit. It provides the following functionality:

- Start the layer creation from a certain point in the circuit and retrieve the gates in the layer.
- Provide gates that have been executed. These are removed from the layer, and it is updated accordingly.
  This allows the layer to be used in a loop, where the gates are executed one after another, without any knowledge of the mapping criteria.

### Mapping Loop

#### Gate-based mapping

The gate-based mapping has the following structure:

- For all multi-qubit gates:
  - find a suitable position (close by using BFS) and store the necessary swaps as "exact" swaps. This means, the qubits need to be swapped exactly to the position.
  - The swaps are weighted (larger weight for more qubits and swaps that would "finish" the gate)
- For two-qubit gates the swaps are stored as closeBy swaps, as they only need to be next to each other.
- All possible swaps are evaluated based on the cost function (distance reduction \* weight) and the best one is applied.

#### Shuttling-based mapping

The shuttling-based mapping has the following structure:

- For all multi-qubit gates:
  - find a suitable position
  - Compute move combinations (sequences of moves) that would bring the qubits to the position
  - The moves are evaluated based on the cost function (distance reduction + parallelism) and the best one is applied.

#### Result

The mapping output is a circuit including SWAP gates and shuttling move operations. This circuit can be dumped into a qasm-similar file.

### AOD Scheduler

To convert the abstract move operations into native AOD sequences, the following steps are performed:

- The move operations are bundled into move groups. A move group is a set of move operations that can be executed in parallel.
- For both, the activation of the move group and the deactivation, an AOD activation helper is used. This helper provides some functionality to facilitate the management of the AOD moves in form of AOD activations.
  An activations consists of an initial position, the small offset move needed and the actual move delta of the operation, both for the x and y direction.
- These activations are computed depending on if the loading/unloading can be done in parallel or not (indicated by ActivationMergeType)
- The activations are then converted into AOD sequences, which can then be dumped into a file.

### Neutral Atom Scheduler

The neutral atom scheduler aligns the gates in a as-soon-as-possible fashion, taking into account the blocking constraint of multi-qubit gates.
It also computes the total necessary execution time of the circuit, based on the hardware information of the architecture.
Currently, gates and shuttling moves are treated equally, and can be executed in parallel.
