"""闭式循环 TP 拓扑：一级网格 + 等熵二级节点。"""

from pathlib import Path

import pytest

from core.closed_cycle_layer import ClosedCycleLayer, ClosedCycleTPInput, build_axis
from core.fluid_property_solver import CoolPropFluidPropertySolver


def test_axis_bounds_and_quantiles():
    assert build_axis(300.0, 400.0, (0.5,)) == [300.0, 350.0, 400.0]
    assert build_axis(100.0, 500.0, (0.25, 0.5, 0.75)) == [100.0, 200.0, 300.0, 400.0, 500.0]


def test_grid_and_isentropic_secondary_water_via_layer():
    inp = ClosedCycleTPInput(
        fluid="Water",
        t_min=400.0,
        t_max=500.0,
        p_min=200.0,
        p_max=2000.0,
        t_quantiles=(0.25, 0.75),
        p_quantiles=(0.25, 0.75),
        max_mass_flow=10.0,
    )
    layer = ClosedCycleLayer(inp)
    layer.analyze_topology()
    primary = [n for n in layer.nodes.values() if n.parent is None]
    secondary = [n for n in layer.nodes.values() if n.parent is not None]
    by_index = {n.index: n for n in primary}
    assert len(primary) >= 1
    for n in primary:
        st = layer.properties.state("TP", n.T, n.P)
        assert st["T"] == pytest.approx(n.T, rel=1e-12)
        assert st["P"] == pytest.approx(n.P, rel=1e-12)
        assert n.H == pytest.approx(st["H"], rel=1e-9)
        assert n.S == pytest.approx(st["S"], rel=1e-9)
    for snode in secondary:
        parent = by_index[snode.parent]
        assert snode.S == pytest.approx(parent.S, rel=1e-6)
        assert inp.t_min <= snode.T <= inp.t_max

    for e in layer.mechanical_edges + layer.heat_edges:
        assert e.tail in layer.nodes and e.head in layer.nodes
        assert e.mass_flow is None


def test_analyze_accepts_injected_solver():
    inp = ClosedCycleTPInput(
        fluid="Water",
        t_min=400.0,
        t_max=420.0,
        p_min=200.0,
        p_max=400.0,
    )
    solver = CoolPropFluidPropertySolver("Water")
    layer = ClosedCycleLayer(inp, properties=solver)
    layer.analyze_topology()
    assert len([n for n in layer.nodes.values() if n.parent is None]) >= 1


def test_helium_tp_topology_ts_plot():
    """He 工质宽温压范围拓扑，并输出 T–S 图（节点 + 机械边 + 换热边）。"""
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    inp = ClosedCycleTPInput(
        fluid="He",
        t_min=100.0,
        t_max=900.0,
        p_min=1000.0,
        p_max=9000.0,
        t_quantiles=(0.5,),
        p_quantiles=(0.5,),
    )
    layer = ClosedCycleLayer(inp)
    layer.analyze_topology()
    primary = [n for n in layer.nodes.values() if n.parent is None]
    secondary = [n for n in layer.nodes.values() if n.parent is not None]
    assert len(primary) >= 1

    fig, ax = plt.subplots(figsize=(8, 6))
    z_edge = 2
    z_pt = 3

    # 机械边（等熵离散链）：绿线
    for k, e in enumerate(layer.mechanical_edges):
        nt, nh = layer.nodes[e.tail], layer.nodes[e.head]
        ax.plot(
            [nt.S, nh.S],
            [nt.T, nh.T],
            color="tab:green",
            linewidth=1.2,
            alpha=0.85,
            zorder=z_edge,
            label="Mechanical edges" if k == 0 else None,
        )

    # 换热边（等压、温度相邻）：红虚线
    for k, e in enumerate(layer.heat_edges):
        nt, nh = layer.nodes[e.tail], layer.nodes[e.head]
        ax.plot(
            [nt.S, nh.S],
            [nt.T, nh.T],
            color="tab:red",
            linewidth=1.0,
            linestyle="--",
            alpha=0.85,
            zorder=z_edge,
            label="Heat edges" if k == 0 else None,
        )

    ax.scatter(
        [n.S for n in primary],
        [n.T for n in primary],
        c="tab:blue",
        s=36,
        label="Level-1 (TP grid)",
        zorder=z_pt,
    )
    if secondary:
        ax.scatter(
            [n.S for n in secondary],
            [n.T for n in secondary],
            c="tab:orange",
            s=40,
            marker="x",
            label="Level-2 (isentropic)",
            zorder=z_pt,
        )
    ax.set_xlabel("S [kJ/(kg·K)]")
    ax.set_ylabel("T [K]")
    ax.set_title(
        "He closed-cycle TP topology (T in [100,900] K, P in [1,9] MPa, quantile 0.5)"
    )
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = Path(__file__).resolve().parent / "ts_topology_he.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    assert out.is_file()


def test_fluid_solver_hp_hs_ps_match_tp_reference():
    """HP/HS/PS 与 TP 参考态在数值上应一致（过冷水单相区）。"""
    sol = CoolPropFluidPropertySolver("Water")
    ref = sol.state("TP", 450.0, 500.0)
    from_hp = sol.state("HP", ref["H"], ref["P"])
    from_hs = sol.state("HS", ref["H"], ref["S"])
    from_ps = sol.state("PS", ref["P"], ref["S"])
    for st in (from_hp, from_hs, from_ps):
        assert st["T"] == pytest.approx(ref["T"], rel=1e-6)
        assert st["P"] == pytest.approx(ref["P"], rel=1e-6)
        assert st["H"] == pytest.approx(ref["H"], rel=1e-6)
        assert st["S"] == pytest.approx(ref["S"], rel=1e-6)
