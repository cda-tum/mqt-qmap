//
// Created by Ludwig Schmid on 19.10.23.
//

#pragma once

#include "operations/AodOperation.hpp"
#include "utils.hpp"

namespace qc {
// A CoordIndex corresponds to node in the SLM grid, where an atom can be placed
// (or not).
using CoordIndex   = std::uint32_t;
using CoordIndices = std::vector<CoordIndex>;
// A HwQubit corresponds to an atom in the neutral atom architecture. It can be
// used as qubit or not and occupies a certain position in the architecture.
using HwQubit  = uint32_t;
using HwQubits = std::set<HwQubit>;
using HwPositions =
    std::vector<HwQubits>; // TODO make consistent use between set and vector
// A Qubit corresponds to a qubit in the quantum circuit. It can be mapped to a
// hardware qubit.
using Qubits = std::set<Qubit>;

// Swaps are between hardware Qubits (one of them can be unmapped).
using Swap          = std::pair<HwQubit, HwQubit>;
using Swaps         = std::vector<Swap>;
using WeightedSwap  = std::pair<Swap, fp>;
using WeightedSwaps = std::vector<WeightedSwap>;
// Moves are between coordinates (the first is occupied, the second is not).
using AtomMove = std::pair<CoordIndex, CoordIndex>;

// Used to represent operations
using OpPointer = std::unique_ptr<AodOperation>;
using GateList  = std::vector<const Operation*>;

} // namespace qc
