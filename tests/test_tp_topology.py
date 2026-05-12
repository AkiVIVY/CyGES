"""闭式循环 TP 拓扑：氦气工况单测，输出含节点、边与子循环的综合图。"""

from pathlib import Path

import pytest

from core.closed_cycle_layer import ClosedCycleLayer, ClosedCycleTPInput


def test_helium_topology_overview_plot():
    """
    He 宽温压拓扑：断言存在子循环；输出双子图（T–S 与 P–S），含一级/二级节点、
    机械/换热边、最小子循环多边形及编号。
    """
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    inp = ClosedCycleTPInput(
        fluid="He",
        t_min=100.0,
        t_max=900.0,
        p_min=1000.0,
        p_max=9000.0,
        t_quantiles=(0.3,  0.7),
        p_quantiles=(0.3,  0.7),
    )
    layer = ClosedCycleLayer(inp)
    layer.analyze_topology()
    primary = [n for n in layer.nodes.values() if n.parent is None]
    secondary = [n for n in layer.nodes.values() if n.parent is not None]
    assert len(primary) >= 1
    assert len(layer.subcycles) >= 1

    fig, (ax_ts, ax_ps) = plt.subplots(1, 2, figsize=(14, 6.5))
    z_edge = 2
    z_poly = 3
    z_pt = 4

    def draw_edges(ax, xkey: str, ykey: str) -> None:
        for e in layer.edges.values():
            nt, nh = layer.nodes[e.tail], layer.nodes[e.head]
            xs = (getattr(nt, xkey), getattr(nh, xkey))
            ys = (getattr(nt, ykey), getattr(nh, ykey))
            if e.kind == "mechanical":
                ax.plot(
                    xs,
                    ys,
                    color="tab:green",
                    linestyle="-",
                    linewidth=1.1,
                    alpha=0.75,
                    zorder=z_edge,
                )
            else:
                ax.plot(
                    xs,
                    ys,
                    color="tab:red",
                    linestyle="--",
                    linewidth=1.0,
                    alpha=0.75,
                    zorder=z_edge,
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

    green = plt.Line2D([0], [0], color="tab:green", lw=1.5, label="Mechanical edges")
    red = plt.Line2D([0], [0], color="tab:red", lw=1.5, ls="--", label="Heat edges")
    fig.legend(
        handles=[green, red],
        loc="upper center",
        ncol=2,
        fontsize=9,
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