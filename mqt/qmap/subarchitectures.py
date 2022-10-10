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
from typing import Dict, NewType, Set, Tuple

import retworkx as rx
from mqt.qmap import Architecture

from qiskit.providers import Backend

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

    # def __init__(self, arch: rx.PyGraph | list[tuple[int, int]] | str | pathlib.Path | Backend | Architecture):
    #     """
    #     Initialize the partial order.

    #     If an architecture is given, the order will be computed for this
    #     specific architecture.
    #     If a str or a path is given instead, the ordering will be loaded from the
    #     subarchitecture library of that name.
    #     """
    #     if type(arch) is str:
    #         if arch in precomputed_backends:
    #             ref = resources.files("mqt.qmap") / "libs" / (arch + ".pickle")
    #             with resources.as_file(ref) as path:
    #                 self.__load_library(path)
    #         else:
    #             self.__load_library(arch)
    #         return

    #     if isinstance(arch, pathlib.Path):
    #         self.__load_library(arch)
    #         return

    #     print(type(arch))

    #     if isinstance(arch, rx.PyGraph):
    #         self.arch = arch
    #     elif isinstance(arch, list) or isinstance(arch, Backend) or isinstance(arch, Architecture):
    #         if isinstance(arch, Backend):
    #             arch = {(a, b) for a, b in arch.configuration().coupling_map}
    #         elif isinstance(arch, Architecture):
    #             arch = arch.coupling_map
    #         num_nodes = max(max(int(u), int(v)) for u, v in arch)
    #         self.arch = rx.PyGraph()
    #         self.arch.add_nodes_from(list(range(num_nodes + 1)))
    #         self.arch.add_edges_from_no_data([tuple(edge) for edge in arch])

    #     self.subarch_order: PartialOrder = PartialOrder({})
    #     self.desirable_subarchitectures: PartialOrder = PartialOrder({})
    #     self.__isomorphisms: dict[tuple[int, int], dict[tuple[int, int], dict[int, int]]] = {}

    #     self.__compute_subarchs()
    #     self.__compute_subarch_order()
    #     self.__compute_desirable_subarchitectures()
    #     return

    @classmethod
    def from_retworkx_graph(cls, graph: rx.PyGraph) -> SubarchitectureOrder:
        """Construct SubarchitectureOrder from retworkx graph."""
        so = SubarchitectureOrder()
        so.arch = graph
        so.subarch_order: PartialOrder = PartialOrder({})
        so.desirable_subarchitectures: PartialOrder = PartialOrder({})
        so.__isomorphisms: dict[tuple[int, int], dict[tuple[int, int], dict[int, int]]] = {}

        so.__compute_subarchs()
        so.__compute_subarch_order()
        so.__compute_desirable_subarchitectures()
        return so

    @classmethod
    def from_coupling_map(cls, coupling_map: set[tuple[int, int]] | list[tuple[int, int]]) -> SubarchitectureOrder:
        """Construct SubarchitectureOrder from coupling map defined as set of tuples of connected qubits."""
        num_nodes = max(max(int(u), int(v)) for u, v in coupling_map)
        graph = rx.PyGraph()
        graph.add_nodes_from(list(range(num_nodes + 1)))
        graph.add_edges_from_no_data([tuple(edge) for edge in coupling_map])

        return cls.from_retworkx_graph(graph)

    @classmethod
    def from_backend(cls, backend: Backend) -> SubarchitectureOrder:
        """Construct SubarchitectureOrder from coupling map defined by qiskit backend."""
        coupling_map = {(a, b) for a, b in backend.configuration().coupling_map}
        return cls.from_coupling_map(coupling_map)

    @classmethod
    def from_qmap_architecture(cls, arch: Architecture) -> SubarchitectureOrder:
        """Construct SubarchitectureOrder from qmap Architecture object."""
        return cls.from_coupling_map(arch.coupling_map)

    @classmethod
    def from_library(cls, path: pathlib.Path) -> SubarchitectureOrder:
        """Construct SubarchitectureOrder from stored library."""
        temp = None
        with pathlib.Path(path).open("rb") as f:
            temp = pickle.load(f)

        so = SubarchitectureOrder()
        so.__dict__.update(temp.__dict__)

        return so

    @classmethod
    def from_string(cls, path: str) -> SubarchitectureOrder:
        """Construct SubarchitectureOrder from library name."""
        if path in precomputed_backends:
            ref = resources.files("mqt.qmap") / "libs" / (path + ".pickle")
            with resources.as_file(ref) as lib_path:
                return cls.from_library(lib_path)

    def optimal_candidates(self, nqubits: int) -> list[rx.PyGraph]:
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

    def covering(self, nqubits: int, size: int) -> list[rx.PyGraph]:
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

    def store_library(self, lib_name: str | pathlib.Path) -> None:
        """Store ordering."""
        if type(lib_name) is str:
            with pathlib.Path(lib_name + ".pickle").open("wb") as f:
                pickle.dump(self, file=f)
        else:
            with pathlib.Path(lib_name).open("wb") as f:
                pickle.dump(self, file=f)

    def __compute_subarchs(self) -> None:
        self.sgs: list[list[rx.PyGraph]] = [[] for i in range(self.arch.num_nodes() + 1)]

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
        for n in range(self.arch.num_nodes() + 1):
            for i in range(len(self.sgs[n])):
                self.subarch_order[(n, i)] = set()
                self.desirable_subarchitectures[(n, i)] = set()
                self.__isomorphisms[(n, i)] = {}

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

    def __combine_iso_with_parent(self, n: int, i: int, j: int) -> None:
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
        po_trans: PartialOrder = PartialOrder({})
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
        po_ref = PartialOrder({})
        for k, v in po.items():
            v_copy = v.copy()
            v_copy.add(k)
            po_ref[k] = v_copy
        return po_ref

    def __inverse_relation(self, po: PartialOrder) -> PartialOrder:
        po_inv = PartialOrder({})
        for n in range(self.arch.num_nodes() + 1):
            for i in range(len(self.sgs[n])):
                po_inv[(n, i)] = set()
        for k, v in po.items():
            for e in v:
                po_inv[e].add(k)
        return po_inv

    def __path_order_less(self, n: int, i: int, n_prime: int, i_prime: int) -> bool:
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

    def __compute_desirable_subarchitectures(self) -> None:
        self.__complete_isos()
        for n in reversed(range(1, len(self.sgs[:-1]))):
            for i in range(len(self.sgs[n])):
                val = self.__isomorphisms[(n, i)]
                for n_prime, i_prime in val.keys():
                    if self.__path_order_less(n, i, n_prime, i_prime):
                        self.desirable_subarchitectures[(n, i)].add((n_prime, i_prime))
                des = list(self.desirable_subarchitectures[(n, i)])
                des.sort()
                new_des: set[tuple[int, int]] = set()
                for j, (n_prime, i_prime) in enumerate(reversed(des)):
                    j = len(des) - j - 1
                    if not any([(n_prime, i_prime) in self.subarch_order[k] for k in des[:j]]):
                        new_des.add((n_prime, i_prime))

                self.desirable_subarchitectures[(n, i)] = new_des
                if len(self.desirable_subarchitectures[(n, i)]) == 0:
                    self.desirable_subarchitectures[(n, i)].add((n, i))
        self.desirable_subarchitectures[self.arch.num_nodes(), 0] = {(self.arch.num_nodes(), 0)}

    def __cand(self, nqubits: int) -> set[tuple[int, int]]:
        return {
            des for (n, i), desirables in self.desirable_subarchitectures.items() if n == nqubits for des in desirables
        }


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
