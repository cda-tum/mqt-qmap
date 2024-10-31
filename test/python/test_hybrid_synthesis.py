"""Test the hybrid Neutral Atom synthesis mapping."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from qiskit import QuantumCircuit

from mqt.qmap import HybridSynthesisMapper, NeutralAtomHybridArchitecture

arch_dir = Path(__file__).parent.parent / "hybridmap" / "architectures"
circuit_dir = Path(__file__).parent.parent / "hybridmap" / "circuits"

qc1 = QuantumCircuit(3)
qc1.h(0)
qc1.cx(0, 1)
qc1.cx(1, 2)

qc2 = QuantumCircuit(3)
qc2.cx(0, 2)
qc2.cx(1, 2)


@pytest.mark.parametrize(
    "arch_filename",
    [
        "rubidium.json",
        "rubidium_hybrid.json",
        "rubidium_shuttling.json",
    ],
)
def test_hybrid_synthesis(arch_filename: str) -> None:
    """Test the hybrid Neutral Atom synthesis mapper evaluation of different circuits."""
    arch = NeutralAtomHybridArchitecture(str(arch_dir / arch_filename))

    synthesis_mapper = HybridSynthesisMapper(arch)
    synthesis_mapper.init_mapping(3)
    best_circuit = synthesis_mapper.evaluate_synthesis_steps([qc1, qc2], True)

    assert best_circuit is not None


@pytest.mark.parametrize(
    "arch_filename",
    [
        "rubidium.json",
        "rubidium_hybrid.json",
        "rubidium_shuttling.json",
    ],
)
def test_hybrid_synthesis_input_output(arch_filename: str) -> None:
    """Test printing and saving the produced circuits."""
    arch = NeutralAtomHybridArchitecture(str(arch_dir / arch_filename))
    synthesis_mapper = HybridSynthesisMapper(arch)
    synthesis_mapper.init_mapping(3)

    synthesis_mapper.append_with_mapping(qc1)
    synthesis_mapper.append_without_mapping(qc2)

    qasm = synthesis_mapper.get_mapped_qc()
    assert qasm is not None

    filename_mapped = Path(__file__).parent / f"{arch_filename}_mapped.qasm"
    synthesis_mapper.save_mapped_qc(str(filename_mapped))

    synthesis_mapper.convert_to_aod()
    qasm_aod = synthesis_mapper.get_mapped_qc_aod()
    assert qasm_aod is not None

    filename_mapped_aod = Path(__file__).parent / f"{arch_filename}_mapped_aod.qasm"
    synthesis_mapper.save_mapped_qc_aod(str(filename_mapped_aod))

    qasm_synth = synthesis_mapper.get_synthesized_qc()
    assert qasm_synth is not None

    filename_synth = Path(__file__).parent / f"{arch_filename}_synthesized.qasm"
    synthesis_mapper.save_synthesized_qc(str(filename_synth))


def test_adjacency_matrix() -> None:
    """Test the adjacency matrix of the hybrid Neutral Atom synthesis mapper."""
    arch = NeutralAtomHybridArchitecture(str(arch_dir / "rubidium.json"))
    synthesis_mapper = HybridSynthesisMapper(arch)
    circ_size = 3
    synthesis_mapper.init_mapping(circ_size)
    synthesis_mapper.append_with_mapping(qc1)
    adj_mat = np.array(synthesis_mapper.get_circuit_adjacency_matrix())
    assert adj_mat is not None
    assert adj_mat.shape == (circ_size, circ_size)
    for i in range(circ_size):
        for j in range(circ_size):
            assert adj_mat[i, j] == adj_mat[j, i]


def help_create_arch(arch_filename: str) -> NeutralAtomHybridArchitecture:
    """Helper function to create a hybrid Neutral Atom architecture."""
    return NeutralAtomHybridArchitecture(str(arch_dir / arch_filename))


def help_create_mapper(arch_filename: str) -> HybridSynthesisMapper:
    """Helper function to create a hybrid synthesis mapper."""
    arch = help_create_arch(arch_filename)
    synthesis_mapper = HybridSynthesisMapper(arch)
    synthesis_mapper.init_mapping(3)
    return synthesis_mapper


def test_keep_alive() -> None:
    """Test the keep alive functionality of the hybrid Neutral Atom synthesis mapper."""
    synthesis_mapper = help_create_mapper("rubidium.json")
    synthesis_mapper.append_with_mapping(qc1)
    _ = synthesis_mapper.get_circuit_adjacency_matrix()
