from pathlib import Path

import pytest
import retworkx as rx
from mqt.qmap import Architecture
from mqt.qmap.subarchitectures import (
    SubarchitectureOrder,
    ibm_guadalupe_subarchitectures,
    rigetti_16_subarchitectures,
)


@pytest.fixture
def ibm_guadalupe() -> SubarchitectureOrder:
    return SubarchitectureOrder.from_coupling_map(
        [
            [0, 1],
            [1, 2],
            [2, 3],
            [3, 5],
            [1, 4],
            [5, 8],
            [4, 7],
            [6, 7],
            [8, 9],
            [7, 10],
            [8, 11],
            [10, 12],
            [12, 15],
            [12, 13],
            [13, 14],
            [11, 14],
        ]
    )


@pytest.fixture
def rigetti16() -> SubarchitectureOrder:
    return SubarchitectureOrder.from_coupling_map(
        [
            [0, 1],
            [1, 2],
            [2, 3],
            [3, 4],
            [4, 5],
            [5, 6],
            [6, 7],
            [7, 8],
            [8, 9],
            [9, 10],
            [10, 11],
            [11, 12],
            [12, 13],
            [13, 14],
            [14, 15],
            [0, 15],
            [3, 12],
            [4, 11],
        ]
    )


@pytest.fixture
def rigetti16_opt() -> rx.PyGraph:
    cm = [
        [0, 1],
        [1, 2],
        [2, 3],
        [3, 4],
        [4, 5],
        [5, 6],
        [6, 7],
        [7, 8],
        [8, 9],
        [9, 10],
        [10, 11],
        [11, 12],
        [12, 13],
        [3, 12],
        [4, 11],
    ]

    num_nodes = max(max(int(u), int(v)) for u, v in cm)
    graph = rx.PyGraph()
    graph.add_nodes_from(list(range(num_nodes + 1)))
    graph.add_edges_from_no_data([tuple(edge) for edge in cm])
    return graph


@pytest.fixture
def singleton_graph() -> rx.PyGraph:
    g = rx.PyGraph()
    g.add_node(0)
    return g


def test_singleton_graph(singleton_graph: rx.PyGraph) -> None:
    """Verify that singleton graph has trivial ordering."""

    order = SubarchitectureOrder.from_retworkx_graph(singleton_graph)

    assert len(order.sgs) == 2
    assert len(order.sgs[0]) == 0
    assert len(order.sgs[1]) == 1
    assert rx.is_isomorphic(order.optimal_candidates(1)[0], singleton_graph)


def test_two_node_graph(singleton_graph: rx.PyGraph) -> None:
    """Verify ordering for graph with two nodes and one edge."""

    order = SubarchitectureOrder.from_coupling_map([[0, 1]])
    assert len(order.sgs) == 3
    assert len(order.sgs[0]) == 0
    assert len(order.sgs[1]) == 1
    assert len(order.sgs[2]) == 1
    assert rx.is_isomorphic(order.optimal_candidates(2)[0], order.sgs[2][0])
    assert rx.is_isomorphic(order.optimal_candidates(1)[0], singleton_graph)


def test_ibm_guadalupe_opt(ibm_guadalupe: SubarchitectureOrder) -> None:
    opt_cand_9 = ibm_guadalupe.optimal_candidates(9)
    assert len(opt_cand_9) == 2
    assert opt_cand_9[0].num_nodes() == 15
    assert opt_cand_9[1].num_nodes() == 15
    assert not rx.is_isomorphic(opt_cand_9[0], opt_cand_9[1])


def test_ibm_guadalupe_cov(ibm_guadalupe: SubarchitectureOrder) -> None:
    cov = ibm_guadalupe.covering(9, 2)
    assert 1 <= len(cov) <= 2

    for sg in ibm_guadalupe.sgs[9]:
        covered = False
        for co in cov:
            if rx.is_subgraph_isomorphic(co, sg):
                covered = True
                break
        assert covered


def test_rigetti16_opt(rigetti16: SubarchitectureOrder, rigetti16_opt: rx.PyGraph) -> None:
    opt = rigetti16.optimal_candidates(10)
    assert len(opt) == 1

    opt_cand = opt[0]
    assert rx.is_isomorphic(opt_cand, rigetti16_opt)


def test_rigetti16_opt_library(rigetti16_opt: rx.PyGraph) -> None:
    opt = rigetti_16_subarchitectures().optimal_candidates(10)
    assert len(opt) == 1

    opt_cand = opt[0]
    assert rx.is_isomorphic(opt_cand, rigetti16_opt)


def test_rigetti16_opt_library_from_str(rigetti16_opt: rx.PyGraph) -> None:
    opt = SubarchitectureOrder.from_string("rigetti_16").optimal_candidates(10)
    assert len(opt) == 1

    opt_cand = opt[0]
    assert rx.is_isomorphic(opt_cand, rigetti16_opt)


def test_ibm_guadalupe_library() -> None:
    opt_cand_9 = ibm_guadalupe_subarchitectures().optimal_candidates(9)
    assert len(opt_cand_9) == 2
    assert opt_cand_9[0].num_nodes() == 15
    assert opt_cand_9[1].num_nodes() == 15
    assert not rx.is_isomorphic(opt_cand_9[0], opt_cand_9[1])


def test_store_subarch(ibm_guadalupe: SubarchitectureOrder) -> None:
    ibm_guadalupe.store_library("tmp")

    p = Path("tmp.pickle")

    loaded_tmp = SubarchitectureOrder.from_library(p)

    if p.exists():
        p.unlink()

    opt_origin = ibm_guadalupe.optimal_candidates(8)
    opt_loaded = loaded_tmp.optimal_candidates(8)

    assert len(opt_origin) == len(opt_loaded)
    for opt_cand_orig, opt_cand_load in zip(opt_origin, opt_loaded):
        assert rx.is_isomorphic(opt_cand_load, opt_cand_orig)


def test_subarchitecture_from_qmap_arch() -> None:
    cm = {(0, 1), (1, 0), (1, 2), (2, 1)}
    arch = Architecture(cm)
    so_arch = SubarchitectureOrder.from_qmap_architecture(arch)
    so_cm = SubarchitectureOrder.from_coupling_map(cm)

    for sg_arch, sg_cm in zip(so_arch, so_cm):
        rx.is_isomorphic(sg_arch, sg_cm)


def test_invalid_opt_cand_arg(ibm_guadalupe: SubarchitectureOrder) -> None:
    with pytest.raises(ValueError):
        ibm_guadalupe.optimal_candidates(100)
