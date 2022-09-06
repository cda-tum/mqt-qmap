"""
Functionality for computing good subarchitectures for quantum circuit mapping.

The function provided in this module... TODO
"""

from typing import NewType, Union
from itertools import combinations
from collections import defaultdict
import rustworkx as rx
import matplotlib.pyplot as plt

Architecture = NewType('Architecture',
                       Union[rx.PyGraph, list[tuple[int, int]]])


class SubarchitectureOrder:
    """
    Class representing partial order of Subarchitectures.

    Attributes
    ----------
    arch : Architecture
        quantum computing architecture whose subarchitectures should be ordered
    order : dict[int, list[int]]
        ordering of variables according to subarchitecture order
    sgs: list[list[Architecture]]
        subarchitectures of arch. sgs[i][j]: the j-th subarchitecture of size i

    Methods
    -------
    optimal_candidate(nqubits)
        return optimal candidate for mapping quantum circuits of a given size
    plot_subarchitectures()
        plot subarchitectures on grid
    covering(nqubits, size)
        return a covering for nqubit circuits limited by size
    """

    def __init__(self, arch: Union[Architecture, str]):
        """
        Initialize the partial order.

        If an architecture is given, the order will be computed for this
        specific architecture.
        If a str is given instead, the ordering will be loaded from the
        subarchitecture library of that name.
        """
        if type(arch) is str:
            self.__load_library(arch)
            return

        self.arch = arch if type(arch) is rx.PyGraph else rx.PyGraph(arch)
        self.subarch_order = defaultdict(lambda: [])
        self.path_order = defaultdict(lambda: [])

        self.__isomorphisms = defaultdict(lambda: {})

        self.__compute_subgraphs()
        self.__compute_subarch_order()
        self.__compute_path_order()
        self.__transitive_closure()

    def plot_subarchitectures(self) -> plt.Figure:
        """Draw subarchitectures on a Grid."""
        pass

    def optimal_candidate(self, nqubits: int) -> Architecture:
        """Return optimal subarchitecture candidate."""
        pass

    def covering(nqubits: int, size: int) -> list[Architecture]:
        """
        Return covering for nqubit circuits.

        The size of the covering is limited by size.
        Note that a smaller covering might be found.
        """
        pass

    def store_library(self, lib_name: str):
        """Store ordering."""
        pass

    def __load_library(self, lib_name_str):
        pass

    def __compute_subgraphs(self):
        self.sgs = [[] for i in range(self.arch.num_nodes()+1)]

        for i in range(1, self.arch.num_nodes()+1):
            node_combinations = combinations(range(self.arch.num_nodes()), i)
            for sg in (self.arch.subgraph(selected_nodes)
                       for selected_nodes in node_combinations):
                if rx.is_connected(sg):
                    new_class = True
                    for g in self.sgs[i]:
                        if rx.is_isomorphic(sg, g):
                            new_class = False
                            break
                    if new_class:
                        self.sgs[i].append(sg)

        return

    def __compute_subarch_order(self):
        for n, sgs_n in enumerate(self.sgs[:-1]):
            for i, sg in enumerate(sgs_n):
                for j, parent_sg in enumerate(self.sgs[n+1]):
                    matcher = rx.graph_vf2_mapping(parent_sg, sg, subgraph=True)
                    for iso in matcher:
                        self.subarch_order[i].append(j)
                        iso_rev = {}
                        for key, val in iso.items():
                            iso_rev[val] = key
                        self.__isomorphisms[(n, i)][(n+1, j)] = iso_rev
                     #   self.__isomorphisms[(n+1,j)][(n, i)] = iso
                        # TODO: might need isos
                        break  # One isomorphism suffices
        self.__complete_isos()

    def __complete_isos(self):
        for n in reversed(range(len(self.sgs[:-2]))):
            for i in range(len(self.sgs[n])):
                for j in self.order[n][i]:
                    self.__combine_isos(n, i, j)

    def __combine_isos(self, nqubits, i, j):
        """Combine all isomorphisms from sgs[n][i] with those frome sgs[n+1][j]."""
        first = self.__isomorphisms[(n, i)][(n+1, j)]
        for row, k, second in self.__isomorphisms[(n+1,j)].items():
            combined = {}
            for src, img in second.items():
                combined[src] = first[val]
            self.__isomorphisms[(n, i)][(row, k)] = combined

    def __compute_path_order(self):
        for n, sgs_n in enumerate(reversed(self.sgs[:-1])):
            for i, sgs_par in enumerate(self.sgs[n+1]):


    def __transitive_closure(self):
        pass
