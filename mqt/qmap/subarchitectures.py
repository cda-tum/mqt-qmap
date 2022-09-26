"""Functionality for computing good subarchitectures for quantum circuit mapping."""

from __future__ import annotations

import sys

if sys.version_info < (3, 10, 0):
    import importlib_resources as resources
else:
    from importlib import resources

import pathlib
import pickle
from itertools import combinations
from typing import Dict, List, NewType, Tuple, Union, Set

import networkx as nx
import retworkx as rx

Subarchitecture = NewType("Subarchitecture", Union[rx.PyGraph, List[Tuple[int, int]]])
PartialOrder = NewType("PartialOrder", Dict[Tuple[int, int], Set[Tuple[int, int]]])

precomputed_backends = ["rigetti_16", "ibm_guadalupe_16", "sycamore_23"]


class SubarchitectureOrder:
    """
    Class representing partial order of Subarchitectures.

    Attributes
    ----------
    arch : Architecture
        quantum computing architecture whose subarchitectures should be ordered
    sgs: list[list[Architecture]]
        subarchitectures of arch. sgs[i][j]: the j-th subarchitecture of size i
    subarch_order : dict[tuple[int, int], list[tuple[int, int]]]
        ordering of variables according to subarchitecture order.
        tuples (n, i) are indices into self.sgs
    desirable_subarchitectures: dict[tuple[int, int], list[tuple[int, int]]]
        mapping of indices (n, i) to its desirable subarchitectures
        tuples (n, i) are indices into self.sgs

    Methods
    -------
    optimal_candidates(nqubits)
        return optimal candidate for mapping quantum circuits of a given size
    covering(nqubits, size)
        return a covering for nqubit circuits limited by size
    store_library(lib_name)
        serialize this object to avoid recomputing the ordering in the future
    """

    def __init__(self, arch: Subarchitecture | str | pathlib.path):
        """
        Initialize the partial order.

        If an architecture is given, the order will be computed for this
        specific architecture.
        If a str or a path is given instead, the ordering will be loaded from the
        subarchitecture library of that name.
        """
        if type(arch) is str:
            if arch in precomputed_backends:
                ref = resources.files("mqt.qmap") / "libs" / (arch+".pickle")
                with resources.as_file(ref) as path:
                    self.__load_library(path)
            else:
                self.__load_library(path)
            return

        if isinstance(arch, pathlib.Path):
            self.__load_library(arch)
            return

        print(type(arch))

        if isinstance(arch, rx.PyGraph):
            self.arch = arch
        else:
            num_nodes = max(map(lambda edge: max(edge), arch))
            self.arch = rx.PyGraph()
            self.arch.add_nodes_from(list(range(num_nodes+1)))
            self.arch.add_edges_from_no_data(list(map(lambda edge: tuple(edge), arch)))

        self.subarch_order = {}
        self.desirable_subarchitectures = {}
        self.__isomorphisms = {}

        self.__compute_subarchs()
        self.__compute_subarch_order()
        self.__compute_desirable_subarchitectures()

    def optimal_candidates(self, nqubits: int) -> list[Subarchitecture]:
        """Return optimal subarchitecture candidate.

        nqubits : int
            size of circuit for which the optimal candidate should be given.
        """
        if nqubits <= 0 or nqubits > self.arch.num_nodes():
            raise ValueError(
                "Number of qubits must not be smaller or equal 0 or larger then number of physical qubits of architecture."
            )

        if nqubits == self.arch.num_nodes():
            return [self.arch]

        cands = self.__cand(nqubits)
        trans_ord = self.__transitive_closure(self.subarch_order)
        ref_ord = self.__reflexive_closure(trans_ord)

        opt_cands = set(ref_ord[next(iter(cands))])
        for cand in cands:
            opt_cands = opt_cands.intersection(set(ref_ord[cand]))

        ordered_cands = list(opt_cands)
        ordered_cands.sort()
        for cand in ordered_cands:
            opt_cands = opt_cands.difference(trans_ord[cand])

        return [self.sgs[n][i] for (n, i) in opt_cands]

    def covering(self, nqubits: int, size: int) -> list[Subarchitecture]:
        """
        Return covering for nqubit circuits.

        The size of the covering is limited by size.
        Note that a smaller covering might be found.
        """
        cov = self.__cand(nqubits)
        po_trans = self.__transitive_closure(self.subarch_order)
        ref_trans_po = self.__reflexive_closure(po_trans)
        queue = list({el for cand in cov for el in ref_trans_po[cand]})
        queue.sort(reverse=True)

        po_inv = self.__inverse_relation(po_trans)

        while len(cov) > size:
            d = queue.pop()
            cov_d = cov.intersection(po_inv[d])
            if len(cov_d) > 1:
                cov = cov.difference(cov_d)
                cov.add(d)

        return [self.sgs[n][i] for n, i in cov]

    def store_library(self, lib_name: Union[str, pathlib.Path]):
        """Store ordering."""
        if type(lib_name) is str:
            with open(lib_name+".pickle", "wb") as f:
                pickle.dump(self, file=f)
        else:
            pickle.dump(self, file=lib_name)

    def __load_library(self, lib_name: Union[str, pathlib.Path]) -> SubarchitectureOrder:
        temp = None
        if type(lib_name) is str:
            lib_name += ".pickle"
        with open(lib_name, "rb") as f:
            temp = pickle.load(f)

        self.__dict__.update(temp.__dict__)

        # self.__isomorphisms = temp.__isomorphisms
        # self.arch = temp.arch
        # self.subarch_order = temp.subarch_order
        # self.desirable_subarchitectures = temp.desirable_subarchitectures
        # self.sgs = temp.sgs

    def __compute_subarchs(self) -> None:
        self.sgs = [[] for i in range(self.arch.num_nodes() + 1)]

        for i in range(1, self.arch.num_nodes() + 1):
            node_combinations = combinations(range(self.arch.num_nodes()), i)
            for sg in (self.arch.subgraph(selected_nodes) for selected_nodes in node_combinations):
                if rx.is_connected(sg):
                    new_class = True
                    for g in self.sgs[i]:
                        if rx.is_isomorphic(g, sg):
                            new_class = False
                            break
                    if new_class:
                        self.sgs[i].append(sg)
        # init orders
        for n in range(self.arch.num_nodes()+1):
            for i in range(len(self.sgs[n])):
                self.subarch_order[(n,i)] = set()
                self.desirable_subarchitectures[(n,i)] = set()
                self.__isomorphisms[(n,i)] = {}

    def __compute_subarch_order(self) -> None:
        for n, sgs_n in enumerate(self.sgs[:-1]):
            for i, sg in enumerate(sgs_n):
                for j, parent_sg in enumerate(self.sgs[n + 1]):
                    matcher = rx.graph_vf2_mapping(parent_sg, sg, subgraph=True)
                    for iso in matcher:
                        self.subarch_order[(n, i)].add((n + 1, j))
                        iso_rev = {}
                        for key, val in iso.items():
                            iso_rev[val] = key
                        self.__isomorphisms[(n, i)][(n + 1, j)] = iso_rev
                        break  # One isomorphism suffices

    def __complete_isos(self) -> None:
        for n in reversed(range(1, len(self.sgs[:-1]))):
            for i in range(len(self.sgs[n])):
                for _, i_prime in self.subarch_order[(n, i)]:
                    self.__combine_iso_with_parent(n, i, i_prime)

    def __combine_iso_with_parent(self, n, i, j) -> None:
        """Combine all isomorphisms from sgs[n][i] with those from sgs[n+1][j]."""
        first = self.__isomorphisms[(n, i)][(n + 1, j)]
        for (row, k), second in self.__isomorphisms[(n + 1, j)].items():
            self.__isomorphisms[(n, i)][(row, k)] = SubarchitectureOrder.__combine_isos(first, second)

    @staticmethod
    def __combine_isos(first: dict[int, int], second: dict[int, int]) -> dict[int, int]:
        combined = {}
        for src, img in first.items():
            combined[src] = second[img]
        return combined

    def __transitive_closure(self, po: PartialOrder) -> PartialOrder:
        po_trans = {}
        po_trans[self.arch.num_nodes(), 0] = set()

        for n in reversed(range(1, len(self.sgs[:-1]))):
            for i in range(len(self.sgs[n])):
                new_rel = set(po[(n, i)])
                po_trans[(n, i)] = new_rel.copy()
                for n_prime, i_prime in po_trans[(n, i)]:
                    new_rel = new_rel.union(po_trans[(n_prime, i_prime)])
                po_trans[(n, i)] = new_rel

        return po_trans

    def __reflexive_closure(self, po: PartialOrder) -> PartialOrder:
        po_ref = dict()
        for k, v in po.items():
            v_copy = v.copy()
            v_copy.add(k)
            po_ref[k] = v_copy
        return po_ref

    def __inverse_relation(self, po: PartialOrder) -> PartialOrder:
        po_inv = dict()
        for n in range(self.arch.num_nodes() + 1):
            for i in range(len(self.sgs[n])):
                po_inv[(n,i)] = set()
        for k, v in po.items():
            for e in v:
                po_inv[e].add(k)
        return po_inv

    def __path_order_less(self, n, i, n_prime, i_prime) -> bool:
        lhs = self.sgs[n][i]
        rhs = self.sgs[n_prime][i_prime]
        iso = self.__isomorphisms[(n, i)][(n_prime, i_prime)]
        for v in range(lhs.num_nodes()):
            for w in range(lhs.num_nodes()):
                if v is w:
                    continue
                if (
                    rx.dijkstra_shortest_path_lengths(lhs, v, lambda x: 1, goal=w)[w]
                    > rx.dijkstra_shortest_path_lengths(rhs, iso[v], lambda x: 1, goal=iso[w])[iso[w]]
                ):
                    return True
        return False

    def __compute_desirable_subarchitectures(self):
        self.__complete_isos()
        for n in reversed(range(1, len(self.sgs[:-1]))):
            for i in range(len(self.sgs[n])):
                val = self.__isomorphisms[(n, i)]
                # for (n, i), val in self.__isomorphisms.items():
                for n_prime, i_prime in val.keys():
                    if self.__path_order_less(n, i, n_prime, i_prime):
                        self.desirable_subarchitectures[(n, i)].add((n_prime, i_prime))
                des = list(self.desirable_subarchitectures[(n, i)])
                des.sort()
                new_des = []
                for j, (n_prime, i_prime) in enumerate(reversed(des)):
                    j = len(des) - j - 1
                    if not any(
                        map(
                            lambda k, n_prime=n_prime, i_prime=i_prime: (n_prime, i_prime) in self.subarch_order[k],
                            des[:j],
                        )
                    ):
                        new_des.append((n_prime, i_prime))

                self.desirable_subarchitectures[(n, i)] = new_des
                if len(self.desirable_subarchitectures[(n, i)]) == 0:
                    self.desirable_subarchitectures[(n, i)].append((n, i))
        self.desirable_subarchitectures[self.arch.num_nodes(), 0] = [(self.arch.num_nodes(), 0)]

    def __cand(self, nqubits: int) -> set[Subarchitecture]:
        return {des for (n, i), desirables in self.desirable_subarchitectures.items() if n == nqubits for des in desirables}


def ibm_guadalupe_subarchitectures() -> SubarchitectureOrder:
    """Load the precomputed ibm guadalupe subarchitectures."""
    ref = resources.files("mqt.qmap") / "libs" / "ibm_guadalupe_16.pickle"
    with resources.as_file(ref) as path:
        return SubarchitectureOrder(path)


def rigetti_16_subarchitectures() -> SubarchitectureOrder:
    """Load the precomputed rigetti subarchitectures."""
    ref = resources.files("mqt.qmap") / "libs" / "rigetti_16.pickle"
    with resources.as_file(ref) as path:
        return SubarchitectureOrder(path)
