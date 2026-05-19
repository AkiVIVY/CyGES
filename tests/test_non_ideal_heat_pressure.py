"""非理想换热边：组内层号与压力偏移。"""

import pytest

from core.closed_cycle_layer import Node, SimplifiedEdge, SimplifiedTopology
from core.non_ideal_closed_cycle_layer import (
    NonIdealClosedCycleLayer,
    build_directed_groups,
)


def _heat_chain_topology() -> tuple[SimplifiedTopology, dict[int, Node]]:
    """三节点换热链 0→1→2，理想 P 分别为 100、200、300 kPa。"""
    edges = {
        "SH1": SimplifiedEdge("heat", 0, 1, ("H1",), (), 1.0),
        "SH2": SimplifiedEdge("heat", 1, 2, ("H2",), (), 1.0),
    }
    topo = SimplifiedTopology(
        kept_nodes=frozenset({0, 1, 2}),
        simplified_edges=tuple(sorted(edges.items())),
        merged_into=(),
    )
    nodes = {
        0: Node(0, 300.0, 100.0, 1000.0, 3.0),
        1: Node(1, 350.0, 200.0, 1100.0, 3.1),
        2: Node(2, 400.0, 300.0, 1200.0, 3.2),
    }
    return topo, nodes


def test_heat_group_upstream_layer():
    topo, _ = _heat_chain_topology()
    groups = build_directed_groups(topo, "heat")
    assert len(groups) == 1
    g = groups[0]
    assert g.depth_dict() == {0: 0, 1: 1, 2: 2}
    assert g.upstream_special_nodes == {0}
    assert g.max_depth == 2


class _PSStubSolver:
    """仅实现 ``PS`` 闪蒸：返回 ``H = P + 1000·S``，``T = 300 + P·0.1``，``P,S`` 原样回写。"""

    fluid = "He"

    def __init__(self) -> None:
        self.calls: list[tuple[str, float, float]] = []

    def state(self, pair: str, x: float, y: float):
        self.calls.append((pair, x, y))
        if pair != "PS":
            raise NotImplementedError(pair)
        p, s = x, y
        return {"T": 300.0 + p * 0.1, "P": p, "H": p + 1000.0 * s, "S": s}


def test_apply_heat_pressure_offsets():
    topo, ideal = _heat_chain_topology()
    solver = _PSStubSolver()
    layer = type("L", (), {"simplified": topo, "nodes": ideal, "properties": solver})()
    ni = NonIdealClosedCycleLayer.from_closed_cycle_layer(layer)
    eta = 0.9
    ni.apply_heat_pressure_offsets(heat_efficiency=eta)

    assert ni.nodes is not None
    layers = ni.heat_groups[0].depth_dict()
    for v in (0, 1, 2):
        p_ideal = ideal[v].P
        expected_p = p_ideal * (eta ** layers[v])
        assert ni.nodes[v].P == pytest.approx(expected_p)
        assert ideal[v].P == p_ideal
        assert ni.nodes[v].S == pytest.approx(ideal[v].S)
        assert ni.nodes[v].H == pytest.approx(expected_p + 1000.0 * ideal[v].S)

    assert ni.nodes[0].P == pytest.approx(100.0)
    assert ni.nodes[1].P == pytest.approx(200.0 * 0.9)
    assert ni.nodes[2].P == pytest.approx(300.0 * 0.9**2)
    assert ni.heat_efficiency == eta
    assert all(c[0] == "PS" for c in solver.calls)
    assert len(solver.calls) == len(ideal)
