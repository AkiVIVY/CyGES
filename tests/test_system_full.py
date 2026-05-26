"""全管线可视化测试：热源/冷源 TP-PS → 循环 TP-PS → 循环夹点 → 循环性能 → 系统夹点。

输出 5 张 PNG（``tests/`` 下）。
"""

from __future__ import annotations

from pathlib import Path

import pytest

from core import (
    ClosedCycleLayer,
    ClosedCycleTPInput,
    CycleConfig,
    ExternalSourceInput,
    HeatTQCurve,
    PinchResult,
    ProcessCategory,
    ProcessRecord,
    PropertyRegistry,
    SystemInput,
    SystemPipeline,
    SystemResult,
    build_heat_tq_curves,
)

TESTS_DIR = Path(__file__).resolve().parent


# ──────────────────────────────────────────────────────────────
# 输入构造
# ──────────────────────────────────────────────────────────────


def _make_system_input(use_non_ideal: bool) -> SystemInput:
    """热源 Air + 冷源 H₂ + He 循环。"""

    hot_src = ExternalSourceInput(
        fluid="Air", mass_flow=100.0, T_in=1250.0, P_in=200.0, T_out=500.0, P_out=180.0,
    )
    cold_src = ExternalSourceInput(
        fluid="Hydrogen", mass_flow=3.1, T_in=20.0, P_in=5000.0, T_out=900.0, P_out=4500.0,
    )

    cycle_input = ClosedCycleTPInput(
        fluid="He",
        t_min=40.0,
        t_max=1000.0,
        p_min=2000.0,
        p_max=20000.0,
        t_quantiles=(0.5,),
        p_quantiles=(0.5,),
        subcycle_mass_flow_initial=20.0,
    )

    cycle_cfg = CycleConfig(
        input=cycle_input,
        use_non_ideal=use_non_ideal,
        delta_T_min=20.0,
        heat_method="pinch",
    )

    return SystemInput(
        heat_sources=(hot_src,),
        cold_sources=(cold_src,),
        cycles=(cycle_cfg,),
        delta_T_min=20.0,
        heat_method="pinch",
    )


# ──────────────────────────────────────────────────────────────
# 图 1：冷热源 TP / PS — 纵坐标 T/P，横坐标 P/S
# ──────────────────────────────────────────────────────────────


def _draw_source_ts_ps(props: PropertyRegistry, out_path: Path) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sys_inp = _make_system_input(use_non_ideal=False)

    sources: list[tuple[str, ExternalSourceInput, str]] = [
        ("heat", sys_inp.heat_sources[0], "tab:red"),
        ("cold", sys_inp.cold_sources[0], "tab:blue"),
    ]

    fig, (ax_ts, ax_ps) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("External sources: T\u2013S & P\u2013S", fontsize=11)

    for label, src, color in sources:
        in_st = props(src.fluid, "TP", src.T_in, src.P_in)
        out_st = props(src.fluid, "TP", src.T_out, src.P_out)

        # TS: 横 S / 纵 T
        ax_ts.plot([in_st["S"], out_st["S"]], [src.T_in, src.T_out], "o-", color=color, linewidth=2, label=f"{label} ({src.fluid})")
        ax_ts.annotate("in", (in_st["S"], src.T_in), textcoords="offset points", xytext=(5, 5), fontsize=8, color=color)
        ax_ts.annotate("out", (out_st["S"], src.T_out), textcoords="offset points", xytext=(5, 5), fontsize=8, color=color)

        # PS: 横 S / 纵 P
        ax_ps.plot([in_st["S"], out_st["S"]], [src.P_in, src.P_out], "o-", color=color, linewidth=2, label=f"{label} ({src.fluid})")
        ax_ps.annotate("in", (in_st["S"], src.P_in), textcoords="offset points", xytext=(5, 5), fontsize=8, color=color)
        ax_ps.annotate("out", (out_st["S"], src.P_out), textcoords="offset points", xytext=(5, 5), fontsize=8, color=color)

    ax_ts.set_xlabel("S [kJ/(kg·K)]"); ax_ts.set_ylabel("T [K]"); ax_ts.grid(True, alpha=0.25); ax_ts.legend(fontsize=8)
    ax_ps.set_xlabel("S [kJ/(kg·K)]"); ax_ps.set_ylabel("P [kPa]"); ax_ps.grid(True, alpha=0.25); ax_ps.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ──────────────────────────────────────────────────────────────
# 图 2：循环 TP / PS — 标注过程名
# ──────────────────────────────────────────────────────────────


def _draw_cycle_ts_ps(ideal_report, ni_report, out_path: Path) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle("Cycle topology: T\u2013S & P\u2013S", fontsize=11)

    for col, (report, label) in enumerate([(ideal_report, "ideal"), (ni_report, "non-ideal")]):
        ax_ts = axes[0, col]
        ax_ps = axes[1, col]

        # kept_nodes 散点
        for _, snap in report.nodes:
            ax_ts.scatter(snap.S, snap.T, s=10, color="0.3", zorder=2)
            ax_ps.scatter(snap.S, snap.P, s=10, color="0.3", zorder=2)

        # 边连线 + 箭头 + 标注
        for _, rec in report.by_edge:
            c = "tab:orange" if rec.kind == "mechanical" else "tab:green"
            # TS
            ax_ts.plot(
                [rec.tail_state.S, rec.head_state.S],
                [rec.tail_state.T, rec.head_state.T],
                color=c, linewidth=1.0, alpha=0.7,
            )
            # 方向箭头：从 tail 40% 处画到 tail 60% 处（沿流向前进 20% 边长）
            ts_fr = 0.38
            ts_to = 0.62
            ax_ts.annotate("", xytext=(
                rec.tail_state.S + ts_fr * (rec.head_state.S - rec.tail_state.S),
                rec.tail_state.T + ts_fr * (rec.head_state.T - rec.tail_state.T),
            ), xy=(
                rec.tail_state.S + ts_to * (rec.head_state.S - rec.tail_state.S),
                rec.tail_state.T + ts_to * (rec.head_state.T - rec.tail_state.T),
            ), arrowprops=dict(arrowstyle="->", color="0.3", lw=1.2, alpha=0.7))
            mx_s = (rec.tail_state.S + rec.head_state.S) / 2
            mx_t = (rec.tail_state.T + rec.head_state.T) / 2
            ax_ts.annotate(rec.edge_key, (mx_s, mx_t), textcoords="offset points", xytext=(3, 3), fontsize=5, color="0.4")

            # PS
            ax_ps.plot(
                [rec.tail_state.S, rec.head_state.S],
                [rec.tail_state.P, rec.head_state.P],
                color=c, linewidth=1.0, alpha=0.7,
            )
            ax_ps.annotate("", xytext=(
                rec.tail_state.S + ts_fr * (rec.head_state.S - rec.tail_state.S),
                rec.tail_state.P + ts_fr * (rec.head_state.P - rec.tail_state.P),
            ), xy=(
                rec.tail_state.S + ts_to * (rec.head_state.S - rec.tail_state.S),
                rec.tail_state.P + ts_to * (rec.head_state.P - rec.tail_state.P),
            ), arrowprops=dict(arrowstyle="->", color="0.3", lw=1.2, alpha=0.7))
            my_s = (rec.tail_state.S + rec.head_state.S) / 2
            my_p = (rec.tail_state.P + rec.head_state.P) / 2
            ax_ps.annotate(rec.edge_key, (my_s, my_p), textcoords="offset points", xytext=(3, 3), fontsize=5, color="0.4")

        ax_ts.set_title(f"{label} T\u2013S"); ax_ts.set_xlabel("S [kJ/(kg·K)]"); ax_ts.set_ylabel("T [K]"); ax_ts.grid(True, alpha=0.25)
        ax_ps.set_title(f"{label} P\u2013S"); ax_ps.set_xlabel("S [kJ/(kg·K)]"); ax_ps.set_ylabel("P [kPa]"); ax_ps.grid(True, alpha=0.25)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ──────────────────────────────────────────────────────────────
# 图 3/5 共用一行：T-Q（左）+ 柱状图（右），柱横坐标 = 过程名，纵坐标 = T
# ──────────────────────────────────────────────────────────────


def _draw_pinch_row(
    ax_tq, ax_bar,
    all_curves: list[HeatTQCurve],
    label: str,
    *,
    pinch: PinchResult | None = None,
) -> None:
    """一行：左 T-Q 曲线，右过程组成柱状图（横坐标=过程名，纵坐标=T）。

    如果 ``pinch`` 非空，在 T-Q 图上标注夹点位置和匹配/额外区域。
    """
    import matplotlib.pyplot as plt

    colors_tq = {ProcessCategory.HEAT_ABSORPTION: "tab:red", ProcessCategory.HEAT_REJECTION: "tab:green"}
    labels_tq = {ProcessCategory.HEAT_ABSORPTION: "absorption", ProcessCategory.HEAT_REJECTION: "rejection"}

    # ── 左：T-Q 曲线（吸热曲线平移 delta_Q 对齐放热）──
    for curve in all_curves:
        ls = "-"
        q = curve.q_points
        # 额外曲线用虚线；吸热曲线在绘图中平移 delta_Q
        if pinch is not None:
            if (curve.category == ProcessCategory.HEAT_ABSORPTION and pinch.extra_absorption is curve) or \
               (curve.category == ProcessCategory.HEAT_REJECTION and pinch.extra_rejection is curve):
                ls = "--"
            if curve.category == ProcessCategory.HEAT_ABSORPTION:
                q = tuple(qp + pinch.delta_Q for qp in q)
        ax_tq.plot(q, curve.t_points, marker="o", color=colors_tq.get(curve.category, "gray"),
                   linewidth=1.8, linestyle=ls, markersize=4,
                   label=labels_tq.get(curve.category, str(curve.category)))
    ax_tq.set_title(f"{label} T-Q"); ax_tq.set_xlabel("Q [kW]"); ax_tq.set_ylabel("T [K]")
    ax_tq.grid(True, alpha=0.25); ax_tq.legend(fontsize=7)

    # ── 夹点标记 ──
    if pinch is not None:
        pq = pinch.pinch_Q_hot
        ax_tq.axvline(x=pq, color="black", linestyle=":", linewidth=1.0, alpha=0.5)
        ax_tq.plot(pq, pinch.pinch_T_hot, "D", color="darkred", markersize=8, zorder=5, label="pinch hot")
        ax_tq.plot(pq, pinch.pinch_T_cold, "D", color="darkblue", markersize=8, zorder=5, label="pinch cold")

    # ── 右：过程组成柱状图（横坐标=过程名，纵坐标=T 范围）──
    edge_keys: list[str] = []
    edge_rects: list[tuple[float, float, str, ProcessCategory]] = []

    for curve in all_curves:
        for seg in curve.segments:
            for rec, _dq in seg.contributions:
                k = rec.edge_key
                if k not in edge_keys:
                    edge_keys.append(k)
                edge_rects.append((seg.T_from, seg.T_to, k, rec.category))

    if not edge_rects:
        return

    seen: set[str] = set()
    ordered_keys: list[str] = []
    for _, _, k, _ in sorted(edge_rects, key=lambda x: x[0]):
        if k not in seen:
            seen.add(k)
            ordered_keys.append(k)

    x_positions = range(len(ordered_keys))
    x_map = {k: i for i, k in enumerate(ordered_keys)}
    bar_width = 0.6

    for t1, t2, k, cat in edge_rects:
        x = x_map[k]
        c = "tab:red" if cat == ProcessCategory.HEAT_ABSORPTION else "tab:green"
        ax_bar.bar(x, t2 - t1, width=bar_width, bottom=t1, color=c, alpha=0.7)

    ax_bar.set_xticks(list(x_positions))
    ax_bar.set_xticklabels(ordered_keys, rotation=60, ha="right", fontsize=6)
    ax_bar.set_ylabel("T [K]")
    ax_bar.set_title(f"{label} composition")
    ax_bar.grid(True, axis="y", alpha=0.25)
    ax_bar.set_ylim(ax_tq.get_ylim())


# ──────────────────────────────────────────────────────────────
# 图 4：循环性能
# ──────────────────────────────────────────────────────────────


def _draw_performance(ideal_report, ni_report, out_path: Path) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle("Cycle performance", fontsize=11)

    for col, (report, label) in enumerate([(ideal_report, "ideal"), (ni_report, "non-ideal")]):
        ax_mech = axes[0, col]
        ax_cat = axes[1, col]

        # 机械功率柱
        mech_data: list[tuple[str, float, float]] = []
        for _, rec in report.by_edge:
            if rec.kind != "mechanical" or rec.power_rate is None:
                continue
            signed = abs(rec.power_rate) if rec.category == ProcessCategory.COMPRESSION else -abs(rec.power_rate)
            mech_data.append((rec.edge_key, signed, abs(rec.power_rate)))
        mech_data.sort(key=lambda x: x[2])

        if mech_data:
            labels = [r[0] for r in mech_data]
            values = [r[1] for r in mech_data]
            colors = ["tab:blue" if v >= 0 else "tab:orange" for v in values]
            ax_mech.bar(range(len(labels)), values, color=colors, alpha=0.9)
            ax_mech.set_xticks(range(len(labels)))
            ax_mech.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
        ax_mech.axhline(0, color="0.2", linewidth=1)
        ax_mech.set_title(f"{label} mechanical (|\u1e44|↑)"); ax_mech.set_ylabel("Power [kW]"); ax_mech.grid(True, axis="y", alpha=0.25)

        # 类别合计柱
        cat_sums = {c: 0.0 for c in (ProcessCategory.COMPRESSION, ProcessCategory.EXPANSION, ProcessCategory.HEAT_ABSORPTION, ProcessCategory.HEAT_REJECTION)}
        for _, rec in report.by_edge:
            if rec.power_rate is None:
                continue
            if rec.category == ProcessCategory.COMPRESSION:
                cat_sums[rec.category] += abs(rec.power_rate)
            elif rec.category == ProcessCategory.EXPANSION:
                cat_sums[rec.category] -= abs(rec.power_rate)
            elif rec.category == ProcessCategory.HEAT_ABSORPTION:
                cat_sums[rec.category] += abs(rec.power_rate)
            elif rec.category == ProcessCategory.HEAT_REJECTION:
                cat_sums[rec.category] -= abs(rec.power_rate)

        cat_names = ["compression", "expansion", "absorption", "rejection"]
        cat_vals = [cat_sums[c] for c in (ProcessCategory.COMPRESSION, ProcessCategory.EXPANSION, ProcessCategory.HEAT_ABSORPTION, ProcessCategory.HEAT_REJECTION)]
        ax_cat.bar(cat_names, cat_vals, color=["tab:orange", "tab:blue", "tab:red", "tab:green"], alpha=0.9)
        ax_cat.axhline(0, color="0.2", linewidth=1)
        ax_cat.set_title(f"{label} category sums"); ax_cat.set_ylabel("Power / Heat [kW]")
        ax_cat.tick_params(axis="x", rotation=25, labelsize=8); ax_cat.grid(True, axis="y", alpha=0.25)

        net_mech = cat_vals[0] + cat_vals[1]
        net_heat = cat_vals[2] + cat_vals[3]
        ymax = max(abs(v) for v in cat_vals) if cat_vals else 1
        ax_cat.text(3.5, ymax * 0.9, f"\u03a3 mech={net_mech:.1f}\n\u03a3 heat={net_heat:.1f}", fontsize=8, ha="center", bbox=dict(boxstyle="round", alpha=0.2))

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ──────────────────────────────────────────────────────────────
# 主测试
# ──────────────────────────────────────────────────────────────


def test_system_full_visualization() -> None:
    props = PropertyRegistry()

    # ── 理想 ──
    ideal_sys = _make_system_input(use_non_ideal=False)
    ideal_result = SystemPipeline(ideal_sys).run(props)
    assert ideal_result.cycle_reports
    assert ideal_result.system_pinch is not None

    # ── 非理想 ──
    ni_sys = _make_system_input(use_non_ideal=True)
    ni_result = SystemPipeline(ni_sys).run(props)
    assert ni_result.cycle_reports
    assert ni_result.system_pinch is not None

    # ── 图 1：冷热源 TP/PS ──
    _draw_source_ts_ps(props, TESTS_DIR / "full_01_source_ts_ps.png")

    # ── 图 2：循环 TP/PS ──
    ideal_report = ideal_result.cycle_reports[0]
    ni_report = ni_result.cycle_reports[0]
    _draw_cycle_ts_ps(ideal_report, ni_report, TESTS_DIR / "full_02_cycle_ts_ps.png")

    # ── 图 3：循环夹点 T-Q + 过程组成（2 行：理想/非理想，每行左 T-Q 右柱状图）──
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig3, axes3 = plt.subplots(2, 2, figsize=(16, 10))
    fig3.suptitle("Cycle pinch: T-Q & segment composition", fontsize=11)

    for row_idx, (result, label) in enumerate([(ideal_result, "ideal"), (ni_result, "non-ideal")]):
        cycle_pinch = result.cycle_pinch
        if cycle_pinch is not None:
            curves: list[HeatTQCurve] = [cycle_pinch.matched_absorption, cycle_pinch.matched_rejection]
            if cycle_pinch.extra_absorption is not None:
                curves.append(cycle_pinch.extra_absorption)
            if cycle_pinch.extra_rejection is not None:
                curves.append(cycle_pinch.extra_rejection)
            curves = [c for c in curves if not (c.q_points == (0.0,) and c.t_points == (0.0,))]
            _draw_pinch_row(axes3[row_idx, 0], axes3[row_idx, 1], curves, label, pinch=cycle_pinch)
        else:
            report = result.cycle_reports[0]
            tq = list(build_heat_tq_curves(report, props))
            _draw_pinch_row(axes3[row_idx, 0], axes3[row_idx, 1], tq, label)

    fig3.tight_layout()
    fig3.savefig(TESTS_DIR / "full_03_cycle_pinch_composition.png", dpi=150)
    plt.close(fig3)

    # ── 图 4：循环性能 ──
    _draw_performance(ideal_report, ni_report, TESTS_DIR / "full_04_cycle_performance.png")

    # ── 图 5：系统夹点 T-Q + 过程组成（2 行：理想/非理想）──
    fig5, axes5 = plt.subplots(2, 2, figsize=(16, 10))
    fig5.suptitle("System pinch: T-Q & segment composition", fontsize=11)

    for row_idx, (result, label) in enumerate([(ideal_result, "ideal"), (ni_result, "non-ideal")]):
        pinch = result.system_pinch
        assert pinch is not None
        all_curves_sys: list[HeatTQCurve] = [pinch.matched_absorption, pinch.matched_rejection]
        if pinch.extra_absorption is not None:
            all_curves_sys.append(pinch.extra_absorption)
        if pinch.extra_rejection is not None:
            all_curves_sys.append(pinch.extra_rejection)
        all_curves_sys = [c for c in all_curves_sys if not (c.q_points == (0.0,) and c.t_points == (0.0,))]
        _draw_pinch_row(axes5[row_idx, 0], axes5[row_idx, 1], all_curves_sys, label, pinch=pinch)

    fig5.tight_layout()
    fig5.savefig(TESTS_DIR / "full_05_system_pinch_composition.png", dpi=150)
    plt.close(fig5)
test_system_full_visualization()
