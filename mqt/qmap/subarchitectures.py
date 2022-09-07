"""
Functionality for computing good subarchitectures for quantum circuit mapping.

The function provided in this module... TODO
"""

from collections import defaultdict
from itertools import combinations
from typing import NewType, Union

import networkx as nx
import retworkx as rx
import pickle

Subarchitecture = NewType('Subarchitecture', Union[rx.PyGraph, list[tuple[int, int]]])
PartialOrder = NewType('PartialOrder', defaultdict[tuple[int, int], tuple[int,int]])


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

    def __init__(self, arch: Union[Subarchitecture, str]):
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

        if isinstance(arch, rx.PyGraph):
            self.arch = arch
        else:
            self.arch = rx.networkx_converter(nx.from_edgelist(arch))

        self.subarch_order = defaultdict(lambda: [])
        self.desirable_subarchitectures = defaultdict(lambda: [])

        self.__isomorphisms = defaultdict(lambda: {})

        self.__compute_subgraphs()
        self.__compute_subarch_order()
        self.__compute_desirable_subarchitectures()
        self.subarch_order = dict(self.subarch_order)
        self.desirable_subarchitectures = dict(self.desirable_subarchitectures)
        self.__isomorphisms = dict(self.__isomorphisms)

    def optimal_candidates(self, nqubits: int) -> list[Subarchitecture]:
        """Return optimal subarchitecture candidate."""
        if nqubits <= 0 or nqubits > self.arch.num_nodes():
            raise ValueError('Number of qubits must not be smaller or equal 0 or larger then number of physical qubits of architecture.')

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
            opt_cands = opt_cands.difference(set(trans_ord[cand]))

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
        queue = list(set(el for cand in cov for el in ref_trans_po[cand]))
        queue.sort()
        queue.reverse()

        po_inv = self.__inverse_relation(po_trans)

        while len(cov) > size:
            d = queue.pop()
            print(d)
            cov_d = cov.intersection(set(po_inv[d]))
            print(cov_d)
            if len(cov_d) > 1:
                cov = cov.difference(cov_d)
                cov.add(d)
            print(cov)
            print()

        return [self.sgs[n][i] for n, i in cov]

    def store_library(self, lib_name: str):
        """Store ordering."""
        with open(lib_name, 'wb') as f:
            pickle.dump(self, file=f)

    def __load_library(self, lib_name_str):
        with open(lib_name_str, 'rb') as f:
            temp = pickle.load(f)
            self.__isomorphisms = temp.__isomorphisms
            self.arch = temp.arch
            self.subarch_order = temp.subarch_order
            self.desirable_subarchitectures = temp.desirable_subarchitectures
            self.sgs = temp.sgs

    def __compute_subgraphs(self) -> None:
        self.sgs = [[] for i in range(self.arch.num_nodes()+1)]

        for i in range(1, self.arch.num_nodes()+1):
            node_combinations = combinations(range(self.arch.num_nodes()), i)
            for sg in (self.arch.subgraph(selected_nodes)
                       for selected_nodes in node_combinations):
                if rx.is_connected(sg):
                    new_class = True
                    for g in self.sgs[i]:
                        if rx.is_subgraph_isomorphic(sg, g):
                            new_class = False
                            break
                    if new_class:
                        self.sgs[i].append(sg)

    def __compute_subarch_order(self) -> None:
        for n, sgs_n in enumerate(self.sgs[:-1]):
            for i, sg in enumerate(sgs_n):
                for j, parent_sg in enumerate(self.sgs[n+1]):
                    matcher = rx.graph_vf2_mapping(parent_sg, sg, subgraph=True)
                    for iso in matcher:
                        self.subarch_order[(n, i)].append((n+1,j))
                        iso_rev = {}
                        for key, val in iso.items():
                            iso_rev[val] = key
                        self.__isomorphisms[(n, i)][(n+1, j)] = iso_rev
                        break  # One isomorphism suffices

    def __complete_isos(self) -> None:
        for n in reversed(range(1, len(self.sgs[:-1]))):
            for i in range(len(self.sgs[n])):
                for _, i_prime in self.subarch_order[(n, i)]:
                    self.__combine_iso_with_parent(n, i, i_prime)

    def __combine_iso_with_parent(self, n, i, j) -> None:
        """Combine all isomorphisms from sgs[n][i] with those frome sgs[n+1][j]."""
        first = self.__isomorphisms[(n, i)][(n+1, j)]
        for (row, k), second in self.__isomorphisms[(n+1,j)].items():
            self.__isomorphisms[(n, i)][(row, k)] = SubarchitectureOrder.__combine_isos(first, second)

    @staticmethod
    def __combine_isos(first: dict[int,int], second: dict[int, int]) -> dict[int, int]:
        combined = {}
        for src, img in first.items():
            combined[src] = second[img]
        return combined

    def __transitive_closure(self, po: PartialOrder) -> PartialOrder:
        po_trans = defaultdict(lambda: [])
        for n in reversed(range(1, len(self.sgs[:-1]))):
            for i in range(len(self.sgs[n])):
                new_rel = po[(n,i)].copy()
                po_trans[(n,i)] = new_rel.copy()
                for n_prime, i_prime in po_trans[(n, i)]:
                    new_rel += po_trans[(n_prime, i_prime)]
                po_trans[(n,i)] = new_rel

        return po_trans

    def __reflexive_closure(self, po: PartialOrder) -> PartialOrder:
        po_ref = defaultdict(lambda: [])
        for k, v in po.items():
            v_copy = v.copy()
            v_copy.append(k)
            po_ref[k] = v_copy
        return po_ref

    def __inverse_relation(self, po: PartialOrder) -> PartialOrder:
        po_inv = defaultdict(lambda: [])
        for k, v in po.items():
            for e in v:
                po_inv[e].append(k)
        return po_inv

    def __path_order_less(self, n, i, n_prime, i_prime) -> bool:
        lhs = self.sgs[n][i]
        rhs = self.sgs[n_prime][i_prime]
        iso = self.__isomorphisms[(n,i)][(n_prime,i_prime)]
        for v in range(lhs.num_nodes()):
            for w in range(lhs.num_nodes()):
                if v is w:
                    continue
                if rx.dijkstra_shortest_path_lengths(lhs, v, lambda x: 1, goal=w)[w] > rx.dijkstra_shortest_path_lengths(rhs, iso[v], lambda x: 1, goal=iso[w])[iso[w]]:
                    return True
        return False

    def __compute_desirable_subarchitectures(self):
        self.__complete_isos()
        for n in reversed(range(1, len(self.sgs[:-1]))):
            for i in range(len(self.sgs[n])):
                val = self.__isomorphisms[(n,i)]
        # for (n, i), val in self.__isomorphisms.items():
                for n_prime, i_prime in val.keys():
                    if self.__path_order_less(n, i, n_prime, i_prime):
                        self.desirable_subarchitectures[(n,i)].append((n_prime, i_prime))
                des = self.desirable_subarchitectures[(n, i)].copy()
                des.sort()
                new_des = []
                for j, (n_prime, i_prime) in (enumerate(reversed(des))):
                    j = len(des) - j - 1
                    if not any(map(lambda k, n_prime=n_prime, i_prime=i_prime: (n_prime, i_prime) in self.subarch_order[k], des[:j])):
                        new_des.append((n_prime, i_prime))

                self.desirable_subarchitectures[(n, i)] = new_des
                if len(self.desirable_subarchitectures[(n,i)]) == 0:
                    self.desirable_subarchitectures[(n,i)].append((n,i))
        self.desirable_subarchitectures[self.arch.num_nodes(),0] = [(self.arch.num_nodes(),0)]

    def __cand(self, nqubits: int) -> set[Subarchitecture]:
        all_desirables = [desirables for (n,i), desirables in self.desirable_subarchitectures.items() if n == nqubits]
        return set(des for desirables in all_desirables for des in desirables)
