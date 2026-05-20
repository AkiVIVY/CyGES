"""He 精简拓扑绘图：理想 vs 非理想（换热 σ + 机械 η_is）节点与有向边。"""

from pathlib import Path
import random
from typing import Mapping

import pytest

from core.closed_cycle_layer import ClosedCycleLayer, ClosedCycleTPInput, Node, SimplifiedTopology
from core.non_ideal_closed_cycle_layer import NonIdealClosedCycleLayer
from core.non_ideal_node_offsets_legacy import (
    apply_heat_pressure_offsets,
    apply_mechanical_isentropic_offsets,
)

FOCUS_NODE_INDICES = frozenset({9, 14, 35, 44})


def _helium_tp_input() -> ClosedCycleTPInput:
    return ClosedCycleTPInput(
        fluid="He",
        t_min=100.0,
        t_max=900.0,
        p_min=1000.0,
        p_max=9000.0,
        t_quantiles=(0.3, 0.7),
        p_quantiles=(0.3, 0.7),
        max_mass_flow=10.0,
    )


def _assign_random_subcycle_flows(layer: ClosedCycleLayer, seed: int = 42) -> int:
    n_sc = len(layer.subcycles)
    rng = random.Random(seed)
    mf = float(layer.input.max_mass_flow)
    idxs = list(range(n_sc))
    rng.shuffle(idxs)
    n_pick = min(8, max(3, n_sc // 2))
    for j in idxs[:n_pick]:
        layer.subcycle_mass_flows[j] = round(rng.uniform(-0.45 * mf, 0.45 * mf), 3)
    layer.commit_subcycle_mass_flows_to_topology()
    return n_pick


def _build_offsets_case() -> tuple[ClosedCycleLayer, NonIdealClosedCycleLayer]:
    layer = ClosedCycleLayer(_helium_tp_input())
    _assign_random_subcycle_flows(layer)
    ni = layer.ensure_non_ideal()
    apply_heat_pressure_offsets(ni)
    apply_mechanical_isentropic_offsets(ni)
    assert ni.nodes is not None
    return layer, ni


def _draw_simplified_topology(
    ax,
    *,
    nodes: Mapping[int, Node],
    simp: SimplifiedTopology,
    xkey: str,
    ykey: str,
    xlabel: str,
    ylabel: str,
    title: str,
    node_indices: frozenset[int] | None = None,
    zoom_to_nodes: bool = False,
) -> None:
    from matplotlib.patches import FancyArrowPatch

    indices = (
        sorted(node_indices & simp.kept_nodes)
        if node_indices is not None
        else sorted(simp.kept_nodes)
    )

    xs: list[float] = []
    ys: list[float] = []
    for i in indices:
        n = nodes[i]
        x, y = getattr(n, xkey), getattr(n, ykey)
        xs.append(x)
        ys.append(y)
        ax.scatter(
            x,
            y,
            c="0.25",
            s=36 if node_indices is not None else 28,
            zorder=2,
            edgecolors="white",
            linewidths=0.3,
        )
        ax.annotate(
            str(i),
            (x, y),
            fontsize=5,
            color="0.1",
            ha="center",
            va="bottom",
            xytext=(0, 2),
            textcoords="offset points",
            zorder=3,
        )

    for _, se in simp.simplified_edges:
        if node_indices is not None and (
            se.tail not in node_indices or se.head not in node_indices
        ):
            continue
        nt, nh = nodes[se.tail], nodes[se.head]
        xt, yt = getattr(nt, xkey), getattr(nt, ykey)
        xh, yh = getattr(nh, xkey), getattr(nh, ykey)
        color = "tab:green" if se.kind == "mechanical" else "tab:red"
        ls = "-" if se.kind == "mechanical" else "--"
        mf = se.mass_flow
        if mf is not None and mf > 0:
            ax.add_patch(
                FancyArrowPatch(
                    (xt, yt),
                    (xh, yh),
                    arrowstyle="-|>",
                    mutation_scale=9,
                    linewidth=1.4,
                    color=color,
                    linestyle=ls,
                    zorder=1,
                )
            )
        else:
            ax.plot([xt, xh], [yt, yh], color=color, linestyle=ls, linewidth=1.4, zorder=1)

    if zoom_to_nodes and xs:
        pad_x = (max(xs) - min(xs)) * 0.12 + 1e-6
        pad_y = (max(ys) - min(ys)) * 0.12 + 1e-6
        ax.set_xlim(min(xs) - pad_x, max(xs) + pad_x)
        ax.set_ylim(min(ys) - pad_y, max(ys) + pad_y)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.25)


def test_helium_ideal_and_non_ideal_offsets_topology_plot():
    """
    He 工况：随机子循环流量 → commit → ensure_non_ideal → 换热 σ 与机械 η_is。

    输出 ``tests/non_ideal_offsets_he.png``：2×2 子图（上行理想 ``layer.nodes``，下行非理想 ``ni.nodes``；
    左列 T–S，右列 P–S）。仅保留节点 index 与有向精简边（无流量数字）。
    """
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    layer, ni = _build_offsets_case()
    simp = ni.simplified
    ideal_nodes = layer.nodes

    fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex="col")

    panels = (
        (axes[0, 0], ideal_nodes, "S", "T", "T–S", "ideal"),
        (axes[0, 1], ideal_nodes, "S", "P", "P–S", "ideal"),
        (axes[1, 0], ni.nodes, "S", "T", "T–S", "non-ideal"),
        (axes[1, 1], ni.nodes, "S", "P", "P–S", "non-ideal"),
    )
    for ax, node_map, xk, yk, plane, row_label in panels:
        _draw_simplified_topology(
            ax,
            nodes=node_map,
            simp=simp,
            xkey=xk,
            ykey=yk,
            xlabel="S [kJ/(kg·K)]",
            ylabel=f"{'T [K]' if yk == 'T' else 'P [kPa]'}",
            title=f"{plane} ({row_label})",
        )

    fig.tight_layout()
    out = Path(__file__).resolve().parent / "non_ideal_offsets_he.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    assert out.is_file()


def test_helium_focus_nodes_9_14_35_44_offsets_plot():
    """
    与全图相同工况，额外输出节点 {9, 14, 35, 44} 的局部 2×2 图（端点均在子集内的精简边）。

    输出 ``tests/non_ideal_offsets_nodes_9_14_35_44.png``。
    """
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    layer, ni = _build_offsets_case()
    simp = ni.simplified
    focus = FOCUS_NODE_INDICES & simp.kept_nodes
    assert focus, "聚焦节点应至少有一个在精简 kept_nodes 中"

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    panels = (
        (axes[0, 0], layer.nodes, "S", "T", "T–S", "ideal"),
        (axes[0, 1], layer.nodes, "S", "P", "P–S", "ideal"),
        (axes[1, 0], ni.nodes, "S", "T", "T–S", "non-ideal"),
        (axes[1, 1], ni.nodes, "S", "P", "P–S", "non-ideal"),
    )
    for ax, node_map, xk, yk, plane, row_label in panels:
        _draw_simplified_topology(
            ax,
            nodes=node_map,
            simp=simp,
            xkey=xk,
            ykey=yk,
            xlabel="S [kJ/(kg·K)]",
            ylabel=f"{'T [K]' if yk == 'T' else 'P [kPa]'}",
            title=f"{plane} ({row_label})",
            node_indices=FOCUS_NODE_INDICES,
            zoom_to_nodes=True,
        )

    fig.suptitle("nodes 9, 14, 35, 44", fontsize=9, y=1.02)
    fig.tight_layout()
    out = Path(__file__).resolve().parent / "non_ideal_offsets_nodes_9_14_35_44.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    assert out.is_file()
