"""精简拓扑过程性能统计：理想/非理想 source 与四类过程归类。

绘图（不纳入 git）::

    pytest tests/test_cycle_performance.py::test_ideal_and_non_ideal_cycle_totals_same_figure -q

输出 ``tests/ideal_vs_non_ideal_cycle_performance.png``：3×2 子图，上排理想、下排非理想；
每排为机械过程柱图、四类合计柱图、换热 T-Q 折线（同输入与同子循环流量）。
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
from core.cycle_performance import CyclePerformanceReport
from core.non_ideal_bias import NonIdealClosedCycleLayer

TESTS_DIR = Path(__file__).resolve().parent


def _he_wide_input() -> ClosedCycleTPInput:
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


def _categories(report: CyclePerformanceReport) -> set[ProcessCategory]:
    return {cat for cat, _ in report.by_category}


def _edge_dict(report: CyclePerformanceReport) -> dict:
    return dict(report.by_edge)


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
    selected = [(ek, rec) for ek, rec in heat_records if rec.category == category]
    if not selected:
        return [0.0], [0.0]

    temp_nodes = sorted(
        {float(rec.tail_state.T) for _, rec in selected}
        | {float(rec.head_state.T) for _, rec in selected}
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
        if hi - lo <= 1e-12:
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
            interval_q[i] += m_abs * abs(h_b - h_a)

    q_points = [0.0]
    for dq in interval_q:
        q_points.append(q_points[-1] + dq)
    return q_points, temp_nodes


def _draw_performance_row(
    axes,
    *,
    layer: ClosedCycleLayer,
    report: CyclePerformanceReport,
    row_label: str,
) -> None:
    by_edge = dict(report.by_edge)
    ax_mech, ax_cat, ax_tq = axes

    mech_rows: list[tuple[str, float, float]] = []
    for edge_key, rec in sorted(by_edge.items()):
        if rec.kind != "mechanical" or rec.power_rate is None:
            continue
        if rec.category not in (ProcessCategory.COMPRESSION, ProcessCategory.EXPANSION):
            continue
        signed = _signed_power_by_rule(rec.category, rec.power_rate)
        mech_rows.append((edge_key, signed, abs(rec.power_rate)))
    mech_rows.sort(key=lambda x: x[2])

    if mech_rows:
        labels = [r[0] for r in mech_rows]
        values = [r[1] for r in mech_rows]
        colors = ["tab:blue" if v >= 0 else "tab:orange" for v in values]
        ax_mech.bar(range(len(labels)), values, color=colors, alpha=0.9)
        ax_mech.set_xticks(range(len(labels)))
        ax_mech.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
    ax_mech.axhline(0.0, color="0.2", linewidth=1.0)
    ax_mech.set_title(f"{row_label}: mechanical (|Ẇ|↑)")
    ax_mech.set_ylabel("Power [kW]")
    ax_mech.grid(True, axis="y", alpha=0.25)

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

    colors1 = ["tab:orange", "tab:blue", "tab:red", "tab:green"]
    bars = ax_cat.bar(cat_names, cat_sums, color=colors1, alpha=0.9)
    ax_cat.axhline(0.0, color="0.2", linewidth=1.0)
    ax_cat.set_title(f"{row_label}: category sums")
    ax_cat.set_ylabel("Power / Heat [kW]")
    ax_cat.tick_params(axis="x", rotation=25, labelsize=8)
    ax_cat.grid(True, axis="y", alpha=0.25)
    ymax = max((abs(v) for v in cat_sums), default=1.0)
    pad = 0.06 * ymax + 1e-9
    for rect, val in zip(bars, cat_sums):
        h = rect.get_height()
        y_text = h + pad if h >= 0 else h - pad
        va = "bottom" if h >= 0 else "top"
        ax_cat.annotate(
            f"{val:.1f}",
            xy=(rect.get_x() + rect.get_width() / 2.0, y_text),
            ha="center",
            va=va,
            fontsize=8,
        )

    heat_records = [
        (ek, rec)
        for ek, rec in sorted(by_edge.items())
        if rec.kind == "heat" and rec.power_rate is not None
    ]
    q_abs, t_abs = _build_tq_polyline_by_temperature_nodes(
        layer, heat_records, ProcessCategory.HEAT_ABSORPTION
    )
    q_rej, t_rej = _build_tq_polyline_by_temperature_nodes(
        layer, heat_records, ProcessCategory.HEAT_REJECTION
    )
    ax_tq.plot(q_abs, t_abs, marker="o", color="tab:red", linewidth=1.8, label="absorption")
    ax_tq.plot(q_rej, t_rej, marker="o", color="tab:green", linewidth=1.8, label="rejection")
    ax_tq.set_title(f"{row_label}: heat T-Q")
    ax_tq.set_xlabel("Cumulative |Q| [kW]")
    ax_tq.set_ylabel("T [K]")
    ax_tq.grid(True, alpha=0.25)
    ax_tq.legend(loc="best", fontsize=7)


def test_ideal_performance_report_categories():
    layer = ClosedCycleLayer(_he_wide_input())
    _assign_random_subcycle_flows(layer, seed=42)
    assert layer.simplified is not None
    assert len(layer.simplified.simplified_edges) >= 1

    report = layer.performance_report()
    assert report.source == "ideal"
    assert len(report.by_edge) == len(layer.simplified.simplified_edges)
    assert report.nodes
    assert all(idx in layer.simplified.kept_nodes for idx, _ in report.nodes)

    cats = _categories(report)
    assert ProcessCategory.COMPRESSION in cats
    assert ProcessCategory.EXPANSION in cats
    assert ProcessCategory.HEAT_ABSORPTION in cats
    assert ProcessCategory.HEAT_REJECTION in cats

    for _, rec in report.by_edge:
        assert rec.tail_state.index == rec.tail
        assert rec.head_state.index == rec.head
        if rec.mass_flow is not None:
            assert rec.power_rate == pytest.approx(rec.mass_flow * rec.delta_H)


def test_non_ideal_performance_uses_offset_nodes():
    layer = ClosedCycleLayer(_he_wide_input())
    _assign_random_subcycle_flows(layer, seed=42)
    ni = layer.ensure_non_ideal()
    apply_combined_offsets(ni)

    ctx = resolve_performance_context(layer)
    assert ctx.source == "non_ideal"

    report = compute_cycle_performance(ctx)
    assert report.source == "non_ideal"

    cats = _categories(report)
    assert ProcessCategory.COMPRESSION in cats
    assert ProcessCategory.EXPANSION in cats
    assert ProcessCategory.HEAT_ABSORPTION in cats
    assert ProcessCategory.HEAT_REJECTION in cats

    fresh_ni = NonIdealClosedCycleLayer.from_closed_cycle_layer(layer)
    ideal_report = compute_cycle_performance(resolve_performance_context(layer, non_ideal=fresh_ni))
    assert ideal_report.source == "ideal"

    by_edge = _edge_dict(report)
    ideal_by_edge = _edge_dict(ideal_report)
    mech_keys = [k for k, e in layer.simplified.edges_dict().items() if e.kind == "mechanical"]
    assert mech_keys
    key = mech_keys[0]
    assert by_edge[key].delta_H != ideal_by_edge[key].delta_H


def test_ideal_and_non_ideal_cycle_totals_same_figure() -> None:
    """3×2 性能对比图：上排理想、下排非理想（机械柱 / 四类合计 / T-Q）。"""
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    layer = ClosedCycleLayer(_he_wide_input())
    _assign_random_subcycle_flows(layer, seed=42)

    snap = NonIdealClosedCycleLayer.from_closed_cycle_layer(layer)
    ideal_report = compute_cycle_performance(
        resolve_performance_context(layer, non_ideal=snap)
    )
    assert ideal_report.source == "ideal"

    ni = layer.ensure_non_ideal()
    apply_combined_offsets(ni)
    ni_report = compute_cycle_performance(resolve_performance_context(layer))
    assert ni_report.source == "non_ideal"

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(
        "Ideal vs non-ideal cycle performance (He wide, same subcycle flows, seed=42)",
        fontsize=11,
        y=0.98,
    )

    _draw_performance_row(axes[0], layer=layer, report=ideal_report, row_label="ideal")
    _draw_performance_row(axes[1], layer=layer, report=ni_report, row_label="non-ideal")

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = TESTS_DIR / "ideal_vs_non_ideal_cycle_performance.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    assert out.is_file()
