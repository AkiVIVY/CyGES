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


def _interp_pressure_by_temperature(rec: ProcessRecord, T: float) -> float:
    """在线性温度参数下插值压力，供区间端点 TP 查焓使用。"""
    t0 = float(rec.tail_state.T)
    t1 = float(rec.head_state.T)
    p0 = float(rec.tail_state.P)
    p1 = float(rec.head_state.P)
    if abs(t1 - t0) <= 1e-12:
        return p0
    alpha = (T - t0) / (t1 - t0)
    return p0 + alpha * (p1 - p0)


def _build_tq_polyline_by_temperature_nodes(
    layer: ClosedCycleLayer,
    heat_records: list[tuple[str, ProcessRecord]],
    category: ProcessCategory,
) -> tuple[list[float], list[float]]:
    """按温度节点离散并填充区间能量（TP 查焓），生成单条连续 T-Q 折线。

    对每条换热过程，在每个重叠温区的两个端点 `(T_a, T_b)` 通过 TP 查焓并计算：

    ``dq = |mass_flow| * |h(T_b, P_b) - h(T_a, P_a)|``

    各区间累计后，得到单条连续折线；Q 统一按绝对值，故曲线位于 ``Q>=0`` 右半轴。
    """
    selected = [(ek, rec) for ek, rec in heat_records if rec.category == category]
    if not selected:
        return [0.0], [0.0]

    temp_nodes = sorted(
        {
            float(rec.tail_state.T)
            for _, rec in selected
        }
        | {
            float(rec.head_state.T)
            for _, rec in selected
        }
    )
    if len(temp_nodes) == 1:
        return [0.0], [temp_nodes[0]]

    interval_q = [0.0] * (len(temp_nodes) - 1)
    for _, rec in selected:
        if rec.mass_flow is None:
            continue
        m_abs = abs(float(rec.mass_flow))
        t0 = float(rec.tail_state.T)
        t1 = float(rec.head_state.T)
        lo = min(t0, t1)
        hi = max(t0, t1)
        width = hi - lo
        if width <= 1e-12:
            continue
        for i in range(len(temp_nodes) - 1):
            a = temp_nodes[i]
            b = temp_nodes[i + 1]
            overlap = max(0.0, min(hi, b) - max(lo, a))
            if overlap <= 0.0:
                continue
            Ta = max(lo, a)
            Tb = min(hi, b)
            Pa = _interp_pressure_by_temperature(rec, Ta)
            Pb = _interp_pressure_by_temperature(rec, Tb)
            h_a = float(layer.properties.state("TP", Ta, Pa)["H"])
            h_b = float(layer.properties.state("TP", Tb, Pb)["H"])
            dq = m_abs * abs(h_b - h_a)
            interval_q[i] += dq

    q_points = [0.0]
    for dq in interval_q:
        q_points.append(q_points[-1] + dq)
    return q_points, temp_nodes


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

    heat_records: list[tuple[str, ProcessRecord]] = []
    for edge_key, rec in sorted(by_edge.items()):
        if rec.kind == "heat" and rec.power_rate is not None:
            heat_records.append((edge_key, rec))

    q_abs, t_abs = _build_tq_polyline_by_temperature_nodes(
        layer, heat_records, ProcessCategory.HEAT_ABSORPTION
    )
    q_rej, t_rej = _build_tq_polyline_by_temperature_nodes(
        layer, heat_records, ProcessCategory.HEAT_REJECTION
    )

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
    ax2.plot(q_abs, t_abs, marker="o", color="tab:red", linewidth=2.0, label="absorption")
    ax2.plot(q_rej, t_rej, marker="o", color="tab:green", linewidth=2.0, label="rejection")
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
