"""非理想单工况（Case A）性能可视化：

1) 机械过程功率柱状图（压缩为负、膨胀为正），按 |功率| 升序；
2) 四类过程汇总柱状图（压缩负、膨胀正、吸热负、放热正），柱顶显示 1 位小数；
3) T-Q 图：按温度节点离散后汇总区间能量，仅绘制吸热线与放热线两条连续折线。
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from core import (
    ClosedCycleLayer,
    ClosedCycleTPInput,
    ProcessCategory,
    ProcessRecord,
    apply_combined_offsets,
    compute_cycle_performance,
    resolve_performance_context,
)

TESTS_DIR = Path(__file__).resolve().parent


def _case_a_input() -> ClosedCycleTPInput:
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


def _assign_random_subcycle_flows(layer: ClosedCycleLayer, seed: int) -> None:
    n_sc = len(layer.subcycles)
    rng = random.Random(seed)
    mf = float(layer.input.max_mass_flow)
    idxs = list(range(n_sc))
    rng.shuffle(idxs)
    n_pick = min(8, max(3, n_sc // 2))
    for j in idxs[:n_pick]:
        layer.subcycle_mass_flows[j] = round(rng.uniform(0.05 * mf, 0.45 * mf), 3)
    layer.commit_subcycle_mass_flows_to_topology()


def _signed_power_by_rule(category: ProcessCategory, power_rate: float) -> float:
    if category == ProcessCategory.COMPRESSION:
        return -abs(power_rate)
    if category == ProcessCategory.EXPANSION:
        return abs(power_rate)
    if category == ProcessCategory.HEAT_ABSORPTION:
        return -abs(power_rate)
    if category == ProcessCategory.HEAT_REJECTION:
        return abs(power_rate)
    return power_rate


def _plot_case_a_visualization(out_dir: Path) -> Path:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    layer = ClosedCycleLayer(_case_a_input())
    _assign_random_subcycle_flows(layer, seed=42)
    ni = layer.ensure_non_ideal()
    apply_combined_offsets(ni)
    report = compute_cycle_performance(resolve_performance_context(layer))

    by_edge = dict(report.by_edge)

    mech_rows: list[tuple[str, float, float]] = []
    for edge_key, rec in sorted(by_edge.items()):
        if rec.kind != "mechanical" or rec.power_rate is None:
            continue
        signed = _signed_power_by_rule(rec.category, rec.power_rate)
        if rec.category in (ProcessCategory.COMPRESSION, ProcessCategory.EXPANSION):
            mech_rows.append((edge_key, signed, abs(rec.power_rate)))

    mech_rows.sort(key=lambda x: x[2])
    mech_labels = [r[0] for r in mech_rows]
    mech_values = [r[1] for r in mech_rows]

    cat_order = [
        ProcessCategory.COMPRESSION,
        ProcessCategory.EXPANSION,
        ProcessCategory.HEAT_ABSORPTION,
        ProcessCategory.HEAT_REJECTION,
    ]
    cat_names = ["compression", "expansion", "heat_absorption", "heat_rejection"]
    cat_sums: list[float] = []
    for cat in cat_order:
        s = 0.0
        for _, rec in by_edge.items():
            if rec.category != cat or rec.power_rate is None:
                continue
            s += _signed_power_by_rule(cat, rec.power_rate)
        cat_sums.append(s)

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5))
    fig.suptitle("Case A (non-ideal): process performance visualization", fontsize=10)

    ax0 = axes[0]
    if mech_labels:
        colors = ["tab:blue" if v >= 0 else "tab:orange" for v in mech_values]
        ax0.bar(range(len(mech_labels)), mech_values, color=colors, alpha=0.9)
        ax0.set_xticks(range(len(mech_labels)))
        ax0.set_xticklabels(mech_labels, rotation=60, ha="right", fontsize=8)
    ax0.axhline(0.0, color="0.2", linewidth=1.0)
    ax0.set_title("Mechanical edge power (sorted by |Ẇ| ascending)")
    ax0.set_ylabel("Power [kW]")
    ax0.grid(True, axis="y", alpha=0.25)

    ax1 = axes[1]
    colors1 = ["tab:orange", "tab:blue", "tab:red", "tab:green"]
    bars = ax1.bar(cat_names, cat_sums, color=colors1, alpha=0.9)
    ax1.axhline(0.0, color="0.2", linewidth=1.0)
    ax1.set_title("Category sums")
    ax1.set_ylabel("Power / Heat rate [kW]")
    ax1.tick_params(axis="x", rotation=25)
    ax1.grid(True, axis="y", alpha=0.25)
    ymax = max((abs(v) for v in cat_sums), default=1.0)
    pad = 0.06 * ymax + 1e-9
    for rect, val in zip(bars, cat_sums):
        h = rect.get_height()
        y_text = h + pad if h >= 0 else h - pad
        va = "bottom" if h >= 0 else "top"
        ax1.annotate(
            f"{val:.1f}",
            xy=(rect.get_x() + rect.get_width() / 2.0, y_text),
            ha="center",
            va=va,
            fontsize=9,
        )

    ax2 = axes[2]
    for curve in report.heat_tq_curves:
        color = "tab:red" if curve.category == ProcessCategory.HEAT_ABSORPTION else "tab:green"
        label = "absorption" if curve.category == ProcessCategory.HEAT_ABSORPTION else "rejection"
        ax2.plot(curve.q_points, curve.t_points, marker="o", color=color, linewidth=2.0, label=label)
    ax2.set_title("Heat-process T-Q (two merged polylines)")
    ax2.set_xlabel("Cumulative |Q| [kW] (per line)")
    ax2.set_ylabel("T [K]")
    ax2.grid(True, alpha=0.25)
    ax2.legend(loc="best", fontsize=8)

    fig.tight_layout()
    out = out_dir / "cycle_performance_case_a_non_ideal.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def test_cycle_performance_visualization_case_a_non_ideal() -> None:
    """输出 Case A 非理想性能统计图。"""
    out = _plot_case_a_visualization(TESTS_DIR)
    assert out.is_file()
