"""夹点分析测试：理想/非理想工况下计算 pinch，绘制平移前/后对比图。

绘图输出 ``tests/pinch_analysis_ideal_vs_non_ideal.png``（2×2 子图）。
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from core import (
    ClosedCycleLayer,
    ClosedCycleTPInput,
    HeatTQCurve,
    PinchAnalysisResult,
    ProcessCategory,
    apply_combined_offsets,
    build_heat_tq_curves,
    compute_cycle_performance,
    compute_pinch,
    resolve_performance_context,
)
from core.non_ideal_bias import NonIdealClosedCycleLayer
from core.postprocess import _interp_T_at_Q

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
        mass_flow_max=10.0,
    )


def _assign_random_subcycle_flows(layer: ClosedCycleLayer, seed: int) -> None:
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


def _run_pinch_on_report(report_curves: tuple[HeatTQCurve, ...], delta_T_min: float) -> PinchAnalysisResult:
    rej = next((c for c in report_curves if c.category == ProcessCategory.HEAT_REJECTION), None)
    abs_ = next((c for c in report_curves if c.category == ProcessCategory.HEAT_ABSORPTION), None)
    if rej is None or abs_ is None:
        raise RuntimeError("报告中缺少放热或吸热 T-Q 曲线")
    return compute_pinch(rej, abs_, delta_T_min=delta_T_min)


def _draw_pinch_row(
    axes,
    *,
    row_label: str,
    rej_curve: HeatTQCurve,
    abs_curve: HeatTQCurve,
    pinch: PinchAnalysisResult,
) -> None:
    ax_before, ax_after = axes

    # ── 左图：平移前 ──
    ax_before.plot(rej_curve.q_points, rej_curve.t_points, "s-", color="tab:red", linewidth=1.8, label="rejection")
    ax_before.plot(abs_curve.q_points, abs_curve.t_points, "o-", color="tab:blue", linewidth=1.8, label="absorption")
    ax_before.set_title(f"{row_label}: before pinch")
    ax_before.set_xlabel("Q [kW]")
    ax_before.set_ylabel("T [K]")
    ax_before.grid(True, alpha=0.25)
    ax_before.legend(loc="best", fontsize=7)

    # ── 右图：平移后 ──
    ax_after.plot(
        pinch.rejection.q_points, pinch.rejection.t_points,
        "s-", color="tab:red", linewidth=1.8, label="rejection",
    )
    ax_after.plot(
        pinch.absorption_shifted.q_points, pinch.absorption_shifted.t_points,
        "o-", color="tab:blue", linewidth=1.8, label="absorption (shifted)",
    )

    # 重叠区填充
    if pinch.overlap_rejection is not None and pinch.overlap_absorption is not None:
        x_q = pinch.overlap_rejection.q_points
        y_abs_shifted = tuple(_interp_T_at_Q(pinch.absorption_shifted, q) for q in x_q)
        ax_after.fill_between(
            x_q, y_abs_shifted, pinch.overlap_rejection.t_points,
            alpha=0.15, color="gray", label="overlap",
        )

    # 夹点标记
    ax_after.axvline(x=pinch.pinch_Q_hot, color="black", linestyle="--", linewidth=1.0, alpha=0.6)
    ax_after.plot(
        pinch.pinch_Q_hot, pinch.pinch_T_hot,
        "D", color="darkred", markersize=8, zorder=5,
    )
    ax_after.plot(
        pinch.pinch_Q_hot, pinch.pinch_T_cold,
        "D", color="darkblue", markersize=8, zorder=5,
    )
    ax_after.annotate(
        f"pinch\nΔT={pinch.delta_T_min:.0f} K",
        xy=(pinch.pinch_Q_hot, pinch.pinch_T_hot),
        xytext=(10, 15), textcoords="offset points",
        fontsize=8, color="black",
        arrowprops=dict(arrowstyle="->", color="0.4", lw=0.8),
    )

    ax_after.set_title(f"{row_label}: after pinch (shift={pinch.delta_Q:.2f} kW)")
    ax_after.set_xlabel("Q [kW]")
    ax_after.set_ylabel("T [K]")
    ax_after.grid(True, alpha=0.25)
    ax_after.legend(loc="best", fontsize=7)


def test_pinch_analysis_ideal_and_non_ideal() -> None:
    """2×2 夹点分析图：上排理想、下排非理想；左列平移前、右列平移后。"""
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    delta_T_min = 20.0  # 目标夹点温差 [K]

    # ── 构建理想层 ──
    layer = ClosedCycleLayer(_he_wide_input())
    _assign_random_subcycle_flows(layer, seed=42)
    enthalpy_fn = lambda f, T, P: layer.properties.state("TP", T, P)["H"]

    # 理想
    snap = NonIdealClosedCycleLayer.from_closed_cycle_layer(layer)
    ideal_ctx = resolve_performance_context(layer, non_ideal=snap)
    ideal_report = compute_cycle_performance(ideal_ctx)
    ideal_tq = build_heat_tq_curves(ideal_report, enthalpy_fn)
    ideal_pinch = _run_pinch_on_report(ideal_tq, delta_T_min)

    # 非理想
    ni = layer.ensure_non_ideal()
    apply_combined_offsets(ni)
    ni_report = compute_cycle_performance(resolve_performance_context(layer))
    ni_tq = build_heat_tq_curves(ni_report, enthalpy_fn)
    ni_pinch = _run_pinch_on_report(ni_tq, delta_T_min)

    # ── 绘图 ──
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        f"Pinch analysis (ΔT_min = {delta_T_min:.0f} K) — ideal vs non-ideal (seed=42)",
        fontsize=11, y=0.98,
    )

    # 理想行：放热原曲线 & 吸热原曲线 → pinch
    _draw_pinch_row(
        axes[0],
        row_label="ideal",
        rej_curve=ideal_pinch.rejection,
        abs_curve=next(c for c in ideal_tq if c.category == ProcessCategory.HEAT_ABSORPTION),
        pinch=ideal_pinch,
    )

    # 非理想行
    _draw_pinch_row(
        axes[1],
        row_label="non-ideal",
        rej_curve=ni_pinch.rejection,
        abs_curve=next(c for c in ni_tq if c.category == ProcessCategory.HEAT_ABSORPTION),
        pinch=ni_pinch,
    )

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = TESTS_DIR / "pinch_analysis_ideal_vs_non_ideal.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    assert out.is_file()

    # ── 基本断言 ──
    for pinch_res, label in [(ideal_pinch, "ideal"), (ni_pinch, "non-ideal")]:
        assert pinch_res.delta_T_min == pytest.approx(delta_T_min)
        # 平移后最小温差应 >= delta_T_min（容差内）
        assert pinch_res.pinch_T_hot - pinch_res.pinch_T_cold == pytest.approx(delta_T_min, abs=1e-6)
        assert pinch_res.overlap_rejection is not None, f"{label}: 应有重叠区"
        assert pinch_res.overlap_absorption is not None, f"{label}: 应有重叠区"
