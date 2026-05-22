"""非理想偏置两工况对比图（2×2：理想 vs 非理想，T–S / P–S）。

pytest（仓库根、``PYTHONPATH=.``）::

    pytest tests/test_non_ideal_two_cases.py -q

或直接运行::

    python tests/test_non_ideal_two_cases.py

输出 PNG（不纳入 git）：

- ``tests/non_ideal_offsets_case_a_he_wide.png``
- ``tests/non_ideal_offsets_case_b_he_mid_fine.png``
"""

from __future__ import annotations

import random
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import pytest

TESTS_DIR = Path(__file__).resolve().parent
ROOT = TESTS_DIR.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.closed_cycle_layer import ClosedCycleLayer, ClosedCycleTPInput, Node, SimplifiedTopology
from core.non_ideal_bias import NonIdealClosedCycleLayer, apply_combined_offsets


@dataclass(frozen=True)
class PlotCase:
    """一组工况：输入参数、随机流量种子、输出文件名。"""

    name: str
    inp: ClosedCycleTPInput
    flow_seed: int
    outfile: str


CASES: tuple[PlotCase, ...] = (
    PlotCase(
        name="Case A: He wide T/P",
        inp=ClosedCycleTPInput(
            fluid="He",
            t_min=100.0,
            t_max=900.0,
            p_min=1000.0,
            p_max=9000.0,
            t_quantiles=(0.3, 0.7),
            p_quantiles=(0.3, 0.7),
            mass_flow_max=10.0,
        ),
        flow_seed=42,
        outfile="non_ideal_offsets_case_a_he_wide.png",
    ),
    PlotCase(
        name="Case B: He mid T/P (t_q=0.5, p_q=0.3/0.7)",
        inp=ClosedCycleTPInput(
            fluid="He",
            t_min=250.0,
            t_max=750.0,
            p_min=2500.0,
            p_max=7500.0,
            t_quantiles=(0.5,),
            p_quantiles=(0.3, 0.7),
            mass_flow_max=6.0,
        ),
        flow_seed=2026,
        outfile="non_ideal_offsets_case_b_he_mid_fine.png",
    ),
)


def _assign_random_subcycle_flows(layer: ClosedCycleLayer, seed: int) -> int:
    n_sc = len(layer.subcycles)
    rng = random.Random(seed)
    mf_min = float(layer.input.mass_flow_min)
    mf_max = float(layer.input.mass_flow_max)
    idxs = list(range(n_sc))
    rng.shuffle(idxs)
    n_pick = min(8, max(3, n_sc // 2))
    for j in idxs[:n_pick]:
        layer.subcycle_mass_flows[j] = round(rng.uniform(mf_min, mf_max), 3)
    layer.commit_subcycle_mass_flows_to_topology()
    return n_pick


def _build_offsets_case(inp: ClosedCycleTPInput, flow_seed: int) -> tuple[ClosedCycleLayer, NonIdealClosedCycleLayer]:
    layer = ClosedCycleLayer(inp)
    _assign_random_subcycle_flows(layer, seed=flow_seed)
    ni = layer.ensure_non_ideal()
    apply_combined_offsets(ni)
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
) -> None:
    from matplotlib.patches import FancyArrowPatch

    indices = sorted(simp.kept_nodes)
    for i in indices:
        n = nodes[i]
        x, y = getattr(n, xkey), getattr(n, ykey)
        ax.scatter(
            x,
            y,
            c="0.25",
            s=28,
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

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.25)


def _plot_case(case: PlotCase, out_dir: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    layer, ni = _build_offsets_case(case.inp, case.flow_seed)
    simp = ni.simplified
    n_kept = len(simp.kept_nodes)
    n_sc = len(layer.subcycles)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex="col")
    fig.suptitle(
        f"{case.name}\n"
        f"T=[{case.inp.t_min},{case.inp.t_max}] K  P=[{case.inp.p_min},{case.inp.p_max}] kPa  "
        f"ṁ∈[{case.inp.mass_flow_min},{case.inp.mass_flow_max}]  kept={n_kept}  subcycles={n_sc}  flow_seed={case.flow_seed}",
        fontsize=9,
        y=1.02,
    )

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
            ylabel="T [K]" if yk == "T" else "P [kPa]",
            title=f"{plane} ({row_label})",
        )

    fig.tight_layout()
    out = out_dir / case.outfile
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def test_non_ideal_offsets_case_a_he_wide() -> None:
    """Case A：宽温压 He；输出 ``non_ideal_offsets_case_a_he_wide.png``。"""
    pytest.importorskip("matplotlib")
    out = _plot_case(CASES[0], TESTS_DIR)
    assert out.is_file()


def test_non_ideal_offsets_case_b_he_mid_fine() -> None:
    """Case B：中等温压；输出 ``non_ideal_offsets_case_b_he_mid_fine.png``。"""
    pytest.importorskip("matplotlib")
    out = _plot_case(CASES[1], TESTS_DIR)
    assert out.is_file()


def main() -> None:
    for case in CASES:
        path = _plot_case(case, TESTS_DIR)
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
