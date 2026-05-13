"""非理想精简拓扑：He 工况绘图用例（理想层生成后对若干子循环赋随机流量，再 commit 并 ensure_non_ideal 观察简化效果）。"""

from pathlib import Path
import random

import pytest

from core.closed_cycle_layer import ClosedCycleLayer, ClosedCycleTPInput
from core.non_ideal_closed_cycle_layer import MERGED_ISOLATED_NODE_EDGE_KEY


def test_helium_non_ideal_simplified_topology_plot():
    """
    He 理想拓扑 ``ClosedCycleLayer`` 默认分析后，对若干子循环 ``subcycle_mass_flows`` 写入随机值，
    ``commit_subcycle_mass_flows_to_topology()`` 写回边流量，再 ``ensure_non_ideal()`` 得到 ``simplified``。

    输出双子图（T–S 与 P–S）：淡色 baseline 边；保留节点 / 链上合并点 / 孤立占位合并点分色；
    精简边粗线，非零 ``mass_flow`` 时箭头沿 tail→head；中点仅标数值流量（无单位前缀），多段合并时另起一行标 ``(×n)``。
    """
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch

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
    n_sc = len(layer.subcycles)
    assert n_sc >= 1

    rng = random.Random(42)
    mf = float(inp.max_mass_flow)
    idxs = list(range(n_sc))
    rng.shuffle(idxs)
    n_pick = min(8, max(3, n_sc // 2))
    for j in idxs[:n_pick]:
        layer.subcycle_mass_flows[j] = round(rng.uniform(-0.45 * mf, 0.45 * mf), 3)
    layer.commit_subcycle_mass_flows_to_topology()

    ni = layer.ensure_non_ideal()
    simp = ni.simplified

    merged_keys = set(simp.merged_dict().keys())
    merged_isolated = {
        i for i, k in simp.merged_dict().items() if k == MERGED_ISOLATED_NODE_EDGE_KEY
    }
    merged_chain = merged_keys - merged_isolated
    kept_keys = simp.kept_nodes
    primary = [n for n in layer.nodes.values() if n.parent is None]
    secondary = [n for n in layer.nodes.values() if n.parent is not None]

    fig, (ax_ts, ax_ps) = plt.subplots(1, 2, figsize=(14, 6.5))
    z_base_edge = 1
    z_base_pt = 2
    z_simp_edge = 3
    z_simp_pt = 4

    def draw_baseline_edges(ax, xkey: str, ykey: str) -> None:
        for _, e in layer.edges.items():
            nt, nh = layer.nodes[e.tail], layer.nodes[e.head]
            xt, yt = getattr(nt, xkey), getattr(nt, ykey)
            xh, yh = getattr(nh, xkey), getattr(nh, ykey)
            color = "tab:green" if e.kind == "mechanical" else "tab:red"
            ls = "-" if e.kind == "mechanical" else "--"
            ax.plot(
                [xt, xh],
                [yt, yh],
                color=color,
                linestyle=ls,
                linewidth=0.7,
                alpha=0.22,
                zorder=z_base_edge,
            )

    def draw_baseline_nodes(ax, xkey: str, ykey: str) -> None:
        kept_primary = [n for n in primary if n.index in kept_keys]
        kept_secondary = [n for n in secondary if n.index in kept_keys]
        mch_pri = [n for n in primary if n.index in merged_chain]
        mch_sec = [n for n in secondary if n.index in merged_chain]
        iso_pri = [n for n in primary if n.index in merged_isolated]
        iso_sec = [n for n in secondary if n.index in merged_isolated]

        if kept_primary:
            ax.scatter(
                [getattr(n, xkey) for n in kept_primary],
                [getattr(n, ykey) for n in kept_primary],
                c="tab:blue",
                s=34,
                label="Level-1 kept",
                zorder=z_base_pt,
                edgecolors="white",
                linewidths=0.4,
                alpha=0.9,
            )
        if kept_secondary:
            ax.scatter(
                [getattr(n, xkey) for n in kept_secondary],
                [getattr(n, ykey) for n in kept_secondary],
                c="tab:orange",
                s=38,
                marker="x",
                label="Level-2 kept",
                zorder=z_base_pt,
                alpha=0.9,
            )
        mch_all = mch_pri + mch_sec
        if mch_all:
            ax.scatter(
                [getattr(n, xkey) for n in mch_all],
                [getattr(n, ykey) for n in mch_all],
                c="0.35",
                s=70,
                marker="X",
                linewidths=1.6,
                label=f"Merged chain ({len(mch_all)})",
                zorder=z_simp_pt,
            )
        iso_all = iso_pri + iso_sec
        if iso_all:
            ax.scatter(
                [getattr(n, xkey) for n in iso_all],
                [getattr(n, ykey) for n in iso_all],
                c="tab:purple",
                s=55,
                marker="s",
                linewidths=0.8,
                label=f"Isolated placeholder ({len(iso_all)})",
                zorder=z_simp_pt,
                alpha=0.85,
            )

    def draw_simplified_edges(ax, xkey: str, ykey: str) -> None:
        for sk, se in simp.simplified_edges:
            nt, nh = layer.nodes[se.tail], layer.nodes[se.head]
            xt, yt = getattr(nt, xkey), getattr(nt, ykey)
            xh, yh = getattr(nh, xkey), getattr(nh, ykey)
            color = "darkgreen" if se.kind == "mechanical" else "darkred"
            ls = "-" if se.kind == "mechanical" else "--"
            mf_edge = se.mass_flow
            if mf_edge is not None and mf_edge > 0:
                arr = FancyArrowPatch(
                    (xt, yt),
                    (xh, yh),
                    arrowstyle="-|>",
                    mutation_scale=10,
                    linewidth=1.6,
                    color=color,
                    linestyle=ls,
                    alpha=0.92,
                    zorder=z_simp_edge,
                )
                ax.add_patch(arr)
            else:
                ax.plot(
                    [xt, xh],
                    [yt, yh],
                    color=color,
                    linestyle=ls,
                    linewidth=1.6,
                    alpha=0.92,
                    zorder=z_simp_edge,
                )
            n_const = len(se.constituent_edges)
            mx = (xt + xh) / 2.0
            my = (yt + yh) / 2.0
            if mf_edge is not None:
                mf_line = f"{mf_edge:.3f}"
            else:
                mf_line = "—"
            label_lines = [mf_line]
            if n_const > 1:
                label_lines.append(f"{sk} (×{n_const})")
            ax.annotate(
                "\n".join(label_lines),
                (mx, my),
                fontsize=6.5,
                color="0.15",
                ha="center",
                va="center",
                bbox=dict(
                    boxstyle="round,pad=0.15",
                    facecolor="white",
                    edgecolor="none",
                    alpha=0.82,
                ),
                zorder=z_simp_edge + 1,
            )

    for ax, xk, yk, xl, yl, title in (
        (ax_ts, "S", "T", "S [kJ/(kg·K)]", "T [K]", "T–S"),
        (ax_ps, "S", "P", "S [kJ/(kg·K)]", "P [kPa]", "P–S"),
    ):
        draw_baseline_edges(ax, xk, yk)
        draw_baseline_nodes(ax, xk, yk)
        draw_simplified_edges(ax, xk, yk)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=8)

    green = plt.Line2D(
        [0], [0], color="darkgreen", lw=1.6,
        label="Simplified Mechanical (arrow = flow dir)",
    )
    red = plt.Line2D(
        [0], [0], color="darkred", lw=1.6, ls="--",
        label="Simplified Heat (arrow = flow dir)",
    )
    fig.legend(
        handles=[green, red],
        loc="upper center",
        ncol=2,
        fontsize=8,
        frameon=True,
        bbox_to_anchor=(0.5, 1.02),
    )
    fig.suptitle(
        f"He: random subcycle mass flow on {n_pick}/{n_sc} cells → commit → simplified; "
        f"{len(simp.simplified_edges)} simp. edges, {len(merged_chain)} chain-merged, "
        f"{len(merged_isolated)} isolated, {len(layer.nodes)} baseline nodes",
        fontsize=10.5,
        y=1.06,
    )
    fig.tight_layout()

    out = Path(__file__).resolve().parent / "simplified_topology_he.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    assert out.is_file()
