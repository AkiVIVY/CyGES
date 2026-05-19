"""非理想机械边：等熵效率焓修正与 HP 物性闭合。"""

import random

import pytest

from core.closed_cycle_layer import ClosedCycleLayer, ClosedCycleTPInput, Node, SimplifiedEdge, SimplifiedTopology
from core.fluid_property_solver import ThermoStateTPHS
from core.non_ideal_closed_cycle_layer import NonIdealClosedCycleLayer


class MockSolver:
    """实现 ``HP`` 与 ``PS`` 闪蒸：

    - ``PS(P,S) → H = P + 1000·S``，``T = 300 + P·0.1``；
    - ``HP(H,P) → S = 3.0 + H·0.001``，``T = 300 + H·0.01``。

    两套关系互为反函数仅在 ``H = P + 1000·(3 + H·0.001) = P + 3000 + H`` 时同时成立（无解），
    因此用于测试时只断言「先 PS 得 H1、再 HP 写回」这一调用顺序与显式数值，不要求物理一致。
    """

    fluid = "He"

    def __init__(self) -> None:
        self.calls: list[tuple[str, float, float]] = []

    def state(self, pair: str, x: float, y: float) -> ThermoStateTPHS:
        self.calls.append((pair, x, y))
        if pair == "HP":
            h, p = x, y
            return ThermoStateTPHS(T=300.0 + h * 0.01, P=p, H=h, S=3.0 + h * 0.001)
        if pair == "PS":
            p, s = x, y
            return ThermoStateTPHS(T=300.0 + p * 0.1, P=p, H=p + 1000.0 * s, S=s)
        raise NotImplementedError(pair)


def _layer_stub(simplified: SimplifiedTopology, nodes: dict[int, Node], solver: MockSolver):
    return type("L", (), {"simplified": simplified, "nodes": nodes, "properties": solver})()


def test_mechanical_anchor_prefers_primary():
    """一级节点为锚点；从锚点 PS 得 H1，压缩公式得 H2，再 HP 写回。"""
    edges = {"SM1": SimplifiedEdge("mechanical", 0, 1, ("M1",), (), 1.0)}
    topo = SimplifiedTopology(
        kept_nodes=frozenset({0, 1}),
        simplified_edges=(("SM1", edges["SM1"]),),
        merged_into=(),
    )
    ideal = {
        0: Node(0, 300.0, 100.0, 1000.0, 3.0, parent=None),
        1: Node(1, 320.0, 200.0, 1200.0, 3.1, parent=0),
    }
    solver = MockSolver()
    ni = NonIdealClosedCycleLayer.from_closed_cycle_layer(_layer_stub(topo, ideal, solver))
    eta = 0.85
    ni.apply_mechanical_isentropic_offsets(mechanical_efficiency=eta)

    assert ni.nodes is not None
    assert ni.nodes[0].S == pytest.approx(3.0)
    h1 = 200.0 + 1000.0 * 3.0
    expected_h = (h1 - 1000.0) / eta + 1000.0
    assert ni.nodes[1].H == pytest.approx(expected_h)
    assert ("PS", 200.0, 3.0) in solver.calls
    assert solver.calls[-1] == ("HP", expected_h, 200.0)


def test_compression_enthalpy_formula():
    """单条压缩边：H1=PS(P_head,S_tail)，H2=(H1-H_tail)/η_is + H_tail。"""
    edges = {"SM1": SimplifiedEdge("mechanical", 10, 11, ("M1",), (), 1.0)}
    topo = SimplifiedTopology(
        kept_nodes=frozenset({10, 11}),
        simplified_edges=(("SM1", edges["SM1"]),),
        merged_into=(),
    )
    ideal = {
        10: Node(10, 300.0, 100.0, 500.0, 2.0, parent=None),
        11: Node(11, 350.0, 250.0, 800.0, 2.5),
    }
    solver = MockSolver()
    ni = NonIdealClosedCycleLayer.from_closed_cycle_layer(_layer_stub(topo, ideal, solver))
    eta = 0.8
    ni.apply_mechanical_isentropic_offsets(mechanical_efficiency=eta)

    h1 = 250.0 + 1000.0 * 2.0
    expected = (h1 - 500.0) / eta + 500.0
    assert ni.nodes[11].H == pytest.approx(expected)


def test_expansion_enthalpy_formula():
    """单条膨胀边：H1=PS(P_head,S_tail)，H2=(H1-H_tail)*η_is + H_tail。"""
    edges = {"SM1": SimplifiedEdge("mechanical", 20, 21, ("M1",), (), 1.0)}
    topo = SimplifiedTopology(
        kept_nodes=frozenset({20, 21}),
        simplified_edges=(("SM1", edges["SM1"]),),
        merged_into=(),
    )
    ideal = {
        20: Node(20, 400.0, 300.0, 900.0, 3.0, parent=None),
        21: Node(21, 350.0, 150.0, 700.0, 3.2),
    }
    solver = MockSolver()
    ni = NonIdealClosedCycleLayer.from_closed_cycle_layer(_layer_stub(topo, ideal, solver))
    eta = 0.9
    ni.apply_mechanical_isentropic_offsets(mechanical_efficiency=eta)

    h1 = 150.0 + 1000.0 * 3.0
    expected = (h1 - 900.0) * eta + 900.0
    assert ni.nodes[21].H == pytest.approx(expected)


def test_equal_pressure_hp_only():
    """等压边：等熵假设下 H1≡H_known，故 H2=H_known，仅 HP 闭合。"""
    edges = {"SM1": SimplifiedEdge("mechanical", 30, 31, ("M1",), (), 1.0)}
    topo = SimplifiedTopology(
        kept_nodes=frozenset({30, 31}),
        simplified_edges=(("SM1", edges["SM1"]),),
        merged_into=(),
    )
    ideal = {
        30: Node(30, 300.0, 150.0, 600.0, 2.8, parent=None),
        31: Node(31, 305.0, 150.0, 650.0, 2.9),
    }
    solver = MockSolver()
    ni = NonIdealClosedCycleLayer.from_closed_cycle_layer(_layer_stub(topo, ideal, solver))
    ni.apply_mechanical_isentropic_offsets(mechanical_efficiency=0.85)

    assert ni.nodes is not None
    assert solver.calls[-1] == ("HP", 600.0, 150.0)
    assert not any(c[0] == "PS" for c in solver.calls)
    assert ni.nodes[31].H == pytest.approx(600.0)


def test_branch_walk_updates_upstream_branch():
    """A→D 与 E→D 形成「侧枝」：基准 A 沿无向支路依次更新 D 再更新 E。"""
    edges = {
        "SM1": SimplifiedEdge("mechanical", 0, 3, ("M1",), (), 1.0),
        "SM2": SimplifiedEdge("mechanical", 4, 3, ("M2",), (), 1.0),
    }
    topo = SimplifiedTopology(
        kept_nodes=frozenset({0, 3, 4}),
        simplified_edges=tuple(sorted(edges.items())),
        merged_into=(),
    )
    ideal = {
        0: Node(0, 300.0, 100.0, 1000.0, 3.0, parent=None),
        3: Node(3, 320.0, 200.0, 1200.0, 3.1),
        4: Node(4, 310.0, 120.0, 1100.0, 3.05),
    }
    solver = MockSolver()
    ni = NonIdealClosedCycleLayer.from_closed_cycle_layer(_layer_stub(topo, ideal, solver))
    eta = 0.85
    ni.apply_mechanical_isentropic_offsets(mechanical_efficiency=eta)

    assert ni.nodes is not None
    assert ni.nodes[0].S == pytest.approx(3.0)

    h1_d = 200.0 + 1000.0 * 3.0
    expected_h_d = (h1_d - 1000.0) / eta + 1000.0
    assert ni.nodes[3].H == pytest.approx(expected_h_d)
    expected_s_d = 3.0 + expected_h_d * 0.001

    h1_e = 120.0 + 1000.0 * expected_s_d
    expected_h_e = (h1_e - expected_h_d) / eta + expected_h_d
    assert ni.nodes[4].H == pytest.approx(expected_h_e)


def test_ideal_nodes_unchanged():
    ideal = {
        0: Node(0, 300.0, 100.0, 1000.0, 3.0, parent=None),
        1: Node(1, 320.0, 200.0, 1100.0, 3.1),
    }
    edges = {"SM1": SimplifiedEdge("mechanical", 0, 1, ("M1",), (), 1.0)}
    topo = SimplifiedTopology(
        kept_nodes=frozenset({0, 1}),
        simplified_edges=(("SM1", edges["SM1"]),),
        merged_into=(),
    )
    solver = MockSolver()
    layer = _layer_stub(topo, ideal, solver)
    ni = NonIdealClosedCycleLayer.from_closed_cycle_layer(layer)
    ni.apply_mechanical_isentropic_offsets(mechanical_efficiency=0.85)

    assert layer.nodes[1].H == pytest.approx(1100.0)
    assert ni.ideal_nodes[1].H == pytest.approx(1100.0)


def test_forked_mechanical_group_entropy_increases_along_sm6():
    """
    He 分叉组 {SM5:22→4, SM6:23→22, SM7:23→24}：一级节点 4 在侧枝末端，
    基准须为 upstream_special 23，沿 SM6 流动方向熵不减。
    """
    inp = ClosedCycleTPInput(
        fluid="He",
        t_min=100.0,
        t_max=900.0,
        p_min=1000.0,
        p_max=9000.0,
        t_quantiles=(0.3, 0.7),
        p_quantiles=(0.3, 0.7),
        max_mass_flow=10.0,
    )
    layer = ClosedCycleLayer(inp)
    rng = random.Random(42)
    idxs = list(range(len(layer.subcycles)))
    rng.shuffle(idxs)
    mf = float(inp.max_mass_flow)
    for j in idxs[: min(8, max(3, len(layer.subcycles) // 2))]:
        layer.subcycle_mass_flows[j] = round(rng.uniform(-0.45 * mf, 0.45 * mf), 3)
    layer.commit_subcycle_mass_flows_to_topology()
    ni = layer.ensure_non_ideal()
    g = next(g for g in ni.mechanical_groups if 23 in g.nodes())
    assert g.upstream_special_nodes == {23}
    assert min(g.upstream_special_nodes) == 23
    ni.apply_heat_pressure_offsets()
    ni.apply_mechanical_isentropic_offsets()
    assert ni.nodes is not None
    se = ni.simplified.edges_dict()["SM6"]
    assert se.tail == 23 and se.head == 22
    assert ni.nodes[se.head].S >= ni.nodes[se.tail].S - 1e-9
