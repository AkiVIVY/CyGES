"""理想闭式循环 He 工况绘图。"""

from pathlib import Path

import config
import pytest

from core.closed_cycle_layer import ClosedCycleLayer, ClosedCycleTPInput


def test_helium_topology_overview_plot():
    """
    He 宽温压理想拓扑：输出双子图（T–S 与 P–S），含一级/二级节点、
    机械/换热边（箭头表示 ``ṁ`` 方向，标注 ``|ṁ|``）、最小子循环多边形及编号。
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
    default_q = config.SUBCYCLE_INITIAL_MASS_FLOW_FRACTION_OF_MAX * inp.max_mass_flow
    layer.subcycle_mass_flows[0] = 1.27
    layer.commit_subcycle_mass_flows_to_topology()
    if n_sc >= 2:
        layer.subcycle_mass_flows[1] = -float(default_q)
        layer.commit_subcycle_mass_flows_to_topology()

    fig, (ax_ts, ax_ps) = plt.subplots(1, 2, figsize=(14, 6.5))
    z_edge = 2
    z_poly = 3
    z_pt = 4

    def draw_edges(ax, xkey: str, ykey: str) -> None:
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
        (ax_ts, "S", "T", "S [kJ/(kg·K)]", "T [K]", "T–S"),
        (ax_ps, "S", "P", "S [kJ/(kg·K)]", "P [kPa]", "P–S"),
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
