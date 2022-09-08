import retworkx as rx
from mqt.qmap.subarchitectures import SubarchitectureOrder


def test_singleton_graph():
    """Verify that singleton graph has trivial ordering."""
    g = rx.PyGraph()
    g.add_node(0)

    order = SubarchitectureOrder(g)

    assert len(order.sgs) == 2
    assert len(order.sgs[0]) == 0
    assert len(order.sgs[1]) == 1
    assert rx.is_isomorphic(order.optimal_candidates(1)[0], g)


def test_two_node_graph():
    """Verify ordering for graph with two nodes and one edge."""

    order = SubarchitectureOrder([[0, 1]])
    assert len(order.sgs) == 3
    assert len(order.sgs[0]) == 0
    assert len(order.sgs[1]) == 1
    assert len(order.sgs[2]) == 1
    assert rx.is_isomorphic(order.optimal_candidates(2)[0], order.sgs[2][0])

    g = rx.PyGraph()
    g.add_node(0)  # g is a singleton graph.

    assert rx.is_isomorphic(order.optimal_candidates(1)[0], g)


def test_ibm_guadalupe_opt():
    so = SubarchitectureOrder(
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

    opt_cand_9 = so.optimal_candidates(9)
    assert len(opt_cand_9) == 2
    assert opt_cand_9[0].num_nodes() == 15
    assert opt_cand_9[1].num_nodes() == 15
    assert not rx.is_isomorphic(opt_cand_9[0], opt_cand_9[1])


# def test_ibm_guadalupe_cov():
#     so = SubarchitectureOrder(
#         [
#             [0, 1],
#             [1, 2],
#             [2, 3],
#             [3, 5],
#             [1, 4],
#             [5, 8],
#             [4, 7],
#             [6, 7],
#             [8, 9],
#             [7, 10],
#             [8, 11],
#             [10, 12],
#             [12, 15],
#             [12, 13],
#             [13, 14],
#             [11, 14],
#         ]
#     )

#     cov = so.covering(9, 2)
#     assert len(cov) == 2
#     assert cov[0].num_nodes() == 14
#     assert cov[1].num_nodes() == 11
