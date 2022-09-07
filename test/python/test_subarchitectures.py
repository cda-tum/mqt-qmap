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
