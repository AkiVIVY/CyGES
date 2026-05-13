"""闭式循环 TP 拓扑：氦气工况单测，输出含节点、边与子循环的综合图。"""

from pathlib import Path

import config
import pytest

from core.closed_cycle_layer import ClosedCycleLayer, ClosedCycleTPInput


def test_subcycle_mass_flow_defaults_zero_when_max_mass_flow_none():
    """未给 max_mass_flow 时子循环初值系数与 0 相乘，列表与 SubCycle 均为 0。"""
    inp = ClosedCycleTPInput(
        fluid="He",
        t_min=100.0,
        t_max=900.0,
        p_min=1000.0,
        p_max=9000.0,
        t_quantiles=(0.3, 0.7),
        p_quantiles=(0.3, 0.7),
        max_mass_flow=None,
    )
    layer = ClosedCycleLayer(inp)
    assert len(layer.subcycles) >= 1
    assert all(q == 0.0 for q in layer.subcycle_mass_flows)
    assert all(sc.mass_flow == 0.0 for sc in layer.subcycles)
    for ek in layer.subcycles[0].edges:
        assert layer.edges[ek].mass_flow == pytest.approx(0.0)


def test_subcycle_mass_flow_step_fraction_defaults_from_config():
    """未传入 subcycle_mass_flow_step_fraction 时与 config 默认一致。"""
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
    assert inp.subcycle_mass_flow_step_fraction == config.SUBCYCLE_MASS_FLOW_STEP_FRACTION_DEFAULT


def test_auto_analyze_false_then_manual_analyze():
    """``auto_analyze=False`` 时构造后拓扑为空，显式 ``analyze_topology()`` 后填充。"""
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
    layer = ClosedCycleLayer(inp, auto_analyze=False)
    assert len(layer.nodes) == 0
    assert len(layer.edges) == 0
    assert len(layer.subcycles) == 0
    layer.analyze_topology()
    assert len(layer.nodes) >= 1
    assert len(layer.subcycles) >= 1


def test_non_ideal_cleared_on_analyze():
    """``ensure_non_ideal`` 挂载后，``analyze_topology`` 清空 ``non_ideal``；再次 ensure 为新实例。"""
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
    assert layer.non_ideal is None
    ni = layer.ensure_non_ideal()
    assert layer.non_ideal is ni
    layer.analyze_topology()
    assert layer.non_ideal is None
    ni2 = layer.ensure_non_ideal()
    assert ni2 is not ni


def test_baseline_snapshot_edge_mass_flow_detached():
    """基准边记录为值拷贝；父层 ``Edge.mass_flow`` 原地修改不改变已拍快照。"""
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
    ni = layer.ensure_non_ideal()
    rec0 = ni.baseline.edge_records[0]
    k0 = rec0.key
    snap_mf = rec0.mass_flow
    layer.edges[k0].mass_flow = 99999.0 if snap_mf != 99999.0 else 88888.0
    rec_again = next(r for r in ni.baseline.edge_records if r.key == k0)
    assert rec_again.mass_flow == snap_mf


def test_non_ideal_cleared_after_commit():
    """``commit_subcycle_mass_flows_to_topology`` 后 ``non_ideal`` 被清空。"""
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
    layer.ensure_non_ideal()
    assert layer.non_ideal is not None
    layer.commit_subcycle_mass_flows_to_topology()
    assert layer.non_ideal is None


def test_commit_subcycle_mass_flows_len_mismatch_raises():
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
    layer.subcycle_mass_flows = [1.0]
    with pytest.raises(ValueError, match="不一致"):
        layer.commit_subcycle_mass_flows_to_topology()


def test_helium_topology_overview_plot():
    """
    He 宽温压拓扑：断言存在子循环；输出双子图（T–S 与 P–S），含一级/二级节点、
    机械/换热边（箭头表示 ``ṁ`` 方向，标注 ``|ṁ|``；``ṁ<0`` 时箭头与 tail→head 相反）、最小子循环多边形及编号。
    """
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import math

    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, Polygon

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
    primary = [n for n in layer.nodes.values() if n.parent is None]
    secondary = [n for n in layer.nodes.values() if n.parent is not None]
    assert len(primary) >= 1
    assert len(layer.subcycles) >= 1
    n_sc = len(layer.subcycles)
    assert len(layer.subcycle_mass_flows) == n_sc
    default_q = config.SUBCYCLE_INITIAL_MASS_FLOW_FRACTION_OF_MAX * inp.max_mass_flow
    assert all(q == default_q for q in layer.subcycle_mass_flows)
    assert all(sc.mass_flow == default_q for sc in layer.subcycles)
    # 量化步长 = 1% * max_mass_flow = 0.1 kg/s
    layer.subcycle_mass_flows[0] = 1.27
    layer.commit_subcycle_mass_flows_to_topology()
    assert layer.subcycle_mass_flows[0] == pytest.approx(1.3)
    assert layer.subcycles[0].mass_flow == pytest.approx(1.3)

    assert any(e.mass_flow is not None for e in layer.edges.values())

    if n_sc >= 2:
        layer.subcycle_mass_flows[1] = -float(default_q)
        layer.commit_subcycle_mass_flows_to_topology()

    fig, (ax_ts, ax_ps) = plt.subplots(1, 2, figsize=(14, 6.5))
    z_edge = 2
    z_poly = 3
    z_pt = 4

    def draw_edges(ax, xkey: str, ykey: str) -> None:
        """边：无 ``mass_flow`` 时画线段；有时画箭头（负流量时箭头与 tail→head 相反），标注 ``|ṁ|`` 为正。"""

        def inset(ax0: float, ay0: float, bx0: float, by0: float, shrink_f: float) -> tuple[tuple[float, float], tuple[float, float]]:
            dx, dy = bx0 - ax0, by0 - ay0
            L = math.hypot(dx, dy)
            if L < 1e-15:
                return (ax0, ay0), (bx0, by0)
            ux, uy = dx / L, dy / L
            s = shrink_f * L
            return (ax0 + ux * s, ay0 + uy * s), (bx0 - ux * s, by0 - uy * s)

        shrink_fr = 0.08
        for e in layer.edges.values():
            nt, nh = layer.nodes[e.tail], layer.nodes[e.head]
            xt, yt = getattr(nt, xkey), getattr(nt, ykey)
            xh, yh = getattr(nh, xkey), getattr(nh, ykey)
            mf = e.mass_flow
            color = "tab:green" if e.kind == "mechanical" else "tab:red"
            ls = "-" if e.kind == "mechanical" else "--"

            if mf is None:
                ax.plot(
                    (xt, xh),
                    (yt, yh),
                    color=color,
                    linestyle=ls,
                    linewidth=1.1 if e.kind == "mechanical" else 1.0,
                    alpha=0.75,
                    zorder=z_edge,
                )
                continue

            mag = abs(mf)
            if mf >= 0.0:
                (xa, ya), (xb, yb) = inset(xt, yt, xh, yh, shrink_fr)
            else:
                (xa, ya), (xb, yb) = inset(xh, yh, xt, yt, shrink_fr)

            arr = FancyArrowPatch(
                (xa, ya),
                (xb, yb),
                arrowstyle="-|>",
                mutation_scale=11,
                linewidth=1.25,
                color=color,
                linestyle=ls,
                alpha=0.88,
                zorder=z_edge,
            )
            ax.add_patch(arr)
            mx, my = (xa + xb) / 2.0, (ya + yb) / 2.0
            ax.annotate(
                f"{mag:.4g}",
                (mx, my),
                xytext=(0, 3),
                textcoords="offset points",
                fontsize=6.5,
                color="0.2",
                ha="center",
                va="bottom",
                zorder=z_edge + 1,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.72),
            )

    def draw_nodes(ax, xkey: str, ykey: str) -> None:
        ax.scatter(
            [getattr(n, xkey) for n in primary],
            [getattr(n, ykey) for n in primary],
            c="tab:blue",
            s=38,
            label="Level-1 (TP grid)",
            zorder=z_pt,
            edgecolors="white",
            linewidths=0.45,
        )
        if secondary:
            ax.scatter(
                [getattr(n, xkey) for n in secondary],
                [getattr(n, ykey) for n in secondary],
                c="tab:orange",
                s=42,
                marker="x",
                label="Level-2 (isentropic)",
                zorder=z_pt,
            )

    def draw_subcycles(ax, xkey: str, ykey: str) -> None:
        for k, sc in enumerate(layer.subcycles):
            corners = [layer.nodes[i] for i in sc.nodes]
            xs = [getattr(c, xkey) for c in corners]
            ys = [getattr(c, ykey) for c in corners]
            color = f"C{k % 10}"
            poly = Polygon(
                list(zip(xs, ys)),
                closed=True,
                facecolor=color,
                edgecolor=color,
                linewidth=2.0,
                alpha=0.2,
                zorder=z_poly,
            )
            ax.add_patch(poly)
            cx = sum(xs) / 4.0
            cy = sum(ys) / 4.0
            ax.annotate(
                str(k + 1),
                (cx, cy),
                ha="center",
                va="center",
                fontsize=9,
                fontweight="bold",
                color="0.2",
                zorder=z_poly + 1,
            )

    for ax, xk, yk, xl, yl, title in (
        (
            ax_ts,
            "S",
            "T",
            "S [kJ/(kg·K)]",
            "T [K]",
            "T–S",
        ),
        (
            ax_ps,
            "S",
            "P",
            "S [kJ/(kg·K)]",
            "P [kPa]",
            "P–S",
        ),
    ):
        draw_edges(ax, xk, yk)
        draw_subcycles(ax, xk, yk)
        draw_nodes(ax, xk, yk)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=8)

    green = plt.Line2D([0], [0], color="tab:green", lw=1.5, label="Mechanical (arrow = ṁ dir, label = |ṁ|)")
    red = plt.Line2D([0], [0], color="tab:red", lw=1.5, ls="--", label="Heat (arrow = ṁ dir, label = |ṁ|)")
    fig.legend(
        handles=[green, red],
        loc="upper center",
        ncol=2,
        fontsize=8,
        frameon=True,
        bbox_to_anchor=(0.5, 1.02),
    )
    fig.suptitle(
        f"He closed-cycle topology (T=[100,900] K, P=[1,9] MPa); "
        f"{len(layer.subcycles)} minimal subcycle(s)",
        fontsize=11,
        y=1.06,
    )
    fig.tight_layout()
    out = Path(__file__).resolve().parent / "ts_topology_he.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    assert out.is_file()