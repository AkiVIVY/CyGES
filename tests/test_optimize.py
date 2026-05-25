"""差分进化优化测试：Air/H₂/He 系统，非理想，min_max_utility 闭环目标。

每代打印最优值；绘制残差下降曲线。
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from core import (
    ClosedCycleLayer,
    ClosedCycleTPInput,
    CycleConfig,
    ExternalSourceInput,
    HeatTQCurve,
    ProcessCategory,
    PropertyRegistry,
    SystemInput,
    build_heat_tq_curves,
)
from optimize import Optimizer

TESTS_DIR = Path(__file__).resolve().parent


def _make_system_input() -> SystemInput:
    """Air/H₂/He 系统，非理想。"""
    hot = ExternalSourceInput(
        fluid="Air", mass_flow=100.0,
        T_in=1250.0, P_in=200.0, T_out=500.0, P_out=180.0,
    )
    cold = ExternalSourceInput(
        fluid="Hydrogen", mass_flow=3.5,
        T_in=20.0, P_in=5000.0, T_out=900.0, P_out=4500.0,
    )
    cycle_input = ClosedCycleTPInput(
        fluid="He", t_min=40.0, t_max=1000.0, p_min=2000.0, p_max=20000.0,
        t_quantiles=(0.5,), p_quantiles=(0.5,),
        mass_flow_min=-10.0, mass_flow_max=50.0,
        subcycle_mass_flow_initial=20.0,
    )
    return SystemInput(
        heat_sources=(hot,), cold_sources=(cold,),
        cycles=(CycleConfig(
            input=cycle_input, use_non_ideal=True,
            delta_T_min=20.0, heat_method="pinch",
        ),),
        delta_T_min=20.0, heat_method="pinch",
    )


def test_optimize_min_max_utility() -> None:
    """差分进化优化，每代打印进度并绘制残差下降曲线。"""
    sys_inp = _make_system_input()
    props = PropertyRegistry()
    opt = Optimizer(base_input=sys_inp, props=props, objective="min_max_utility")

    print(f"\n维度: {len(opt.bounds)} (t_max + t_q + p_q + {opt._max_sc}×mf)")
    print("CMA-ES, popsize=30, maxiter=300, seed=42\n")

    history: list[tuple[int, float]] = []  # (gen, best_val)

    def _on_gen(gen: int, best_x: list[float], best_val: float, n_evals: int) -> None:
        history.append((gen, best_val))
        t_max, t_min, tq, pq = best_x[0], best_x[1], best_x[2], best_x[3]
        mf_str = ", ".join(f"{best_x[4+i]:.1f}" for i in range(6) if abs(best_x[4+i]) > 1e-6)
        sys.stdout.write(
            f"  gen {gen:3d} | obj={best_val:.5f} | t_max={t_max:.0f} t_min={t_min:.0f} t_q={tq:.4f} p_q={pq:.4f} "
            f"| mf=[{mf_str}] | evals={n_evals}\n"
        )
        sys.stdout.flush()

    print("  开始优化...\n")
    result = opt.run(
        method="cma", maxiter=300, seed=42, early_stop=30, sigma0=0.3,
        callback=_on_gen,
    )

    sp = result.system_result.system_pinch
    assert sp is not None, "系统夹点不应为空"

    report = result.system_result.cycle_reports[0]
    ct = report.cycle_totals

    print("\n" + "=" * 55)
    print("最优解")
    print("=" * 55)
    print(f"  t_max / t_min     = {result.x_opt[0]:.1f} / {result.x_opt[1]:.1f} K")
    print(f"  t_q / p_q         = {result.x_opt[2]:.4f} / {result.x_opt[3]:.4f}")
    mfs = [f"{result.x_opt[4+i]:.2f}" for i in range(min(8, len(result.x_opt)-4))]
    print(f"  mass_flows (前8)   = [{', '.join(mfs)}]")
    print(f"  objective          = {result.objective:.5f}")
    print(f"  hot_utility_demand = {sp.hot_utility_demand:.1f} kW ({sp.hot_utility_ratio*100:.2f}%)")
    print(f"  cold_utility_demand= {sp.cold_utility_demand:.1f} kW ({sp.cold_utility_ratio*100:.2f}%)")
    print(f"  net_mechanical     = {ct.net_mechanical_power:.1f} kW")
    print(f"  net_heat           = {ct.net_heat_rate:.1f} kW")
    print(f"  n_evaluations      = {result.n_evaluations}")

    assert result.n_evaluations > 0

    # ── 绘制残差下降曲线 ──
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    gens, vals = zip(*history)
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(gens, vals, "o-", markersize=2, linewidth=1.0, color="tab:blue")
    ax.set_xlabel("Generation"); ax.set_ylabel("Objective (min-max utility ratio)")
    ax.set_title("DE Convergence — Air/H₂/He (non-ideal, seed=42)")
    ax.grid(True, alpha=0.25)
    ax.set_ylim(bottom=0)

    # 标注初始和最终值
    ax.annotate(f"{vals[0]:.4f}", (gens[0], vals[0]), textcoords="offset points", xytext=(5, 8), fontsize=8, color="0.3")
    ax.annotate(f"{vals[-1]:.4f}", (gens[-1], vals[-1]), textcoords="offset points", xytext=(5, -12), fontsize=8, color="0.3")

    out = TESTS_DIR / "optimize_convergence.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    assert out.is_file()

    # ── 用最优参数绘图（复用 test_system_full 的 _draw_pinch_row）──
    from tests.test_system_full import _draw_pinch_row

    best_result = result.system_result
    best_report = best_result.cycle_reports[0]
    best_pinch = best_result.cycle_pinch or best_result.system_pinch
    assert best_pinch is not None

    tq_curves = list(build_heat_tq_curves(best_report, props))

    # 图 A: 循环 TS/PS（单张）
    fig_a, (ax_ts, ax_ps) = plt.subplots(1, 2, figsize=(12, 5))
    fig_a.suptitle("Optimized cycle: T–S & P–S", fontsize=11)
    for _, snap in best_report.nodes:
        ax_ts.scatter(snap.S, snap.T, s=10, color="0.3", zorder=2)
        ax_ps.scatter(snap.S, snap.P, s=10, color="0.3", zorder=2)
    for _, rec in best_report.by_edge:
        c = "tab:orange" if rec.kind == "mechanical" else "tab:green"
        ax_ts.plot([rec.tail_state.S, rec.head_state.S], [rec.tail_state.T, rec.head_state.T], color=c, lw=1.0, alpha=0.7)
        ax_ps.plot([rec.tail_state.S, rec.head_state.S], [rec.tail_state.P, rec.head_state.P], color=c, lw=1.0, alpha=0.7)
        ts_fr, ts_to = 0.38, 0.62
        for ax, y1, y2 in [(ax_ts, rec.tail_state.T, rec.head_state.T), (ax_ps, rec.tail_state.P, rec.head_state.P)]:
            ax.annotate("", xytext=(rec.tail_state.S + ts_fr*(rec.head_state.S-rec.tail_state.S), y1 + ts_fr*(y2-y1)),
                        xy=(rec.tail_state.S + ts_to*(rec.head_state.S-rec.tail_state.S), y1 + ts_to*(y2-y1)),
                        arrowprops=dict(arrowstyle="->", color="0.3", lw=1.2, alpha=0.7))
    ax_ts.set_xlabel("S"); ax_ts.set_ylabel("T [K]"); ax_ts.grid(True, alpha=0.25)
    ax_ps.set_xlabel("S"); ax_ps.set_ylabel("P [kPa]"); ax_ps.grid(True, alpha=0.25)
    fig_a.tight_layout(); fig_a.savefig(TESTS_DIR / "opt_cycle_ts_ps.png", dpi=150); plt.close(fig_a)

    # 图 B: 循环夹点 T-Q + 组成
    curves_for_pinch: list[HeatTQCurve] = []
    if best_pinch.matched_absorption.q_points != (0.0,) or best_pinch.matched_absorption.t_points != (0.0,):
        curves_for_pinch.append(best_pinch.matched_absorption)
    if best_pinch.matched_rejection.q_points != (0.0,) or best_pinch.matched_rejection.t_points != (0.0,):
        curves_for_pinch.append(best_pinch.matched_rejection)
    if best_pinch.extra_absorption is not None:
        curves_for_pinch.append(best_pinch.extra_absorption)
    if best_pinch.extra_rejection is not None:
        curves_for_pinch.append(best_pinch.extra_rejection)
    fig_b, (ax_tq_b, ax_bar_b) = plt.subplots(1, 2, figsize=(14, 6))
    _draw_pinch_row(ax_tq_b, ax_bar_b, curves_for_pinch, "opt", pinch=best_pinch)
    fig_b.savefig(TESTS_DIR / "opt_pinch_composition.png", dpi=150); plt.close(fig_b)

    # 图 C: 循环性能
    fig_c, (ax_mech, ax_cat) = plt.subplots(1, 2, figsize=(12, 5))
    fig_c.suptitle("Optimized cycle performance", fontsize=11)
    mech_data = []; cat_sums = {pc: 0.0 for pc in (ProcessCategory.COMPRESSION, ProcessCategory.EXPANSION, ProcessCategory.HEAT_ABSORPTION, ProcessCategory.HEAT_REJECTION)}
    for _, rec in best_report.by_edge:
        if rec.power_rate is None: continue
        cat_sums[rec.category] += rec.power_rate
        if rec.kind == "mechanical":
            signed = abs(rec.power_rate) if rec.category == ProcessCategory.COMPRESSION else -abs(rec.power_rate)
            mech_data.append((rec.edge_key, signed, abs(rec.power_rate)))
    mech_data.sort(key=lambda x: x[2])
    if mech_data:
        labels, vals_m = zip(*[(r[0], r[1]) for r in mech_data])
        colors = ["tab:blue" if v >= 0 else "tab:orange" for v in vals_m]
        ax_mech.bar(range(len(labels)), vals_m, color=colors, alpha=0.9)
        ax_mech.set_xticks(range(len(labels))); ax_mech.set_xticklabels(labels, rotation=60, ha="right", fontsize=6)
    ax_mech.axhline(0, color="0.2"); ax_mech.set_title("mechanical"); ax_mech.set_ylabel("kW"); ax_mech.grid(True, axis="y", alpha=0.25)
    cn = ["compression", "expansion", "absorption", "rejection"]
    cv = [cat_sums[pc] for pc in (ProcessCategory.COMPRESSION, ProcessCategory.EXPANSION, ProcessCategory.HEAT_ABSORPTION, ProcessCategory.HEAT_REJECTION)]
    ax_cat.bar(cn, cv, color=["tab:orange","tab:blue","tab:red","tab:green"], alpha=0.9)
    ax_cat.axhline(0, color="0.2"); ax_cat.set_title("category sums"); ax_cat.grid(True, axis="y", alpha=0.25)
    ax_cat.text(3.5, max(abs(v) for v in cv)*0.9, f"Σ mech={cv[0]+cv[1]:.0f}\nΣ heat={cv[2]+cv[3]:.0f}", fontsize=8, ha="center", bbox=dict(boxstyle="round", alpha=0.2))
    fig_c.tight_layout(); fig_c.savefig(TESTS_DIR / "opt_performance.png", dpi=150); plt.close(fig_c)

    # 图 D: 系统夹点
    sys_pinch = best_result.system_pinch
    if sys_pinch is not None:
        sys_curves: list[HeatTQCurve] = [sys_pinch.matched_absorption, sys_pinch.matched_rejection]
        if sys_pinch.extra_absorption is not None:
            sys_curves.append(sys_pinch.extra_absorption)
        if sys_pinch.extra_rejection is not None:
            sys_curves.append(sys_pinch.extra_rejection)
        sys_curves = [c for c in sys_curves if not (c.q_points == (0.0,) and c.t_points == (0.0,))]
        fig_d, (ax_tq_d, ax_bar_d) = plt.subplots(1, 2, figsize=(14, 6))
        _draw_pinch_row(ax_tq_d, ax_bar_d, sys_curves, "opt system", pinch=sys_pinch)
        fig_d.savefig(TESTS_DIR / "opt_system_pinch.png", dpi=150); plt.close(fig_d)

    # 图 E: 子循环网格 + 流量（TS/PS，非简化拓扑）
    t_q_opt = (result.x_opt[2],)
    p_q_opt = (result.x_opt[3],)
    t_max_opt, t_min_opt = result.x_opt[0], result.x_opt[1]
    tp_opt = ClosedCycleTPInput(
        fluid="He", t_min=t_min_opt, t_max=t_max_opt, p_min=2000.0, p_max=20000.0,
        t_quantiles=t_q_opt, p_quantiles=p_q_opt,
        mass_flow_min=-10.0, mass_flow_max=50.0,
    )
    layer_opt = ClosedCycleLayer(tp_opt)
    mf_opt = [float(result.x_opt[4+i]) for i in range(min(len(layer_opt.subcycles), len(result.x_opt)-4))]
    layer_opt.subcycle_mass_flows = mf_opt
    layer_opt.commit_subcycle_mass_flows_to_topology()

    fig_e, (ax_ts_e, ax_ps_e) = plt.subplots(1, 2, figsize=(14, 6))
    fig_e.suptitle("Optimized subcycle grid & flows", fontsize=11)

    # 散点: 所有节点
    for n in layer_opt.nodes.values():
        ax_ts_e.scatter(n.S, n.T, s=8, color="0.4", zorder=1)
        ax_ps_e.scatter(n.S, n.P, s=8, color="0.4", zorder=1)

    # 子循环多边形 + 流量标注
    colors_sc = plt.cm.tab10.colors
    for i, sc in enumerate(layer_opt.subcycles):
        c = colors_sc[i % len(colors_sc)]
        q = sc.mass_flow if sc.mass_flow is not None else 0.0
        n0, n1, n2, n3 = [layer_opt.nodes[idx] for idx in sc.nodes]
        for ax, y_get in [(ax_ts_e, lambda n: n.T), (ax_ps_e, lambda n: n.P)]:
            xs = [n0.S, n1.S, n2.S, n3.S, n0.S]
            ys = [y_get(n0), y_get(n1), y_get(n2), y_get(n3), y_get(n0)]
            ax.plot(xs, ys, "-", color=c, linewidth=1.2, alpha=0.5)
            cx = (n0.S + n2.S) / 2
            cy = (y_get(n0) + y_get(n2)) / 2
            ax.annotate(f"SC{i}\n{q:.1f}", (cx, cy), ha="center", va="center", fontsize=6, color=c, fontweight="bold")

    ax_ts_e.set_xlabel("S [kJ/(kg·K)]"); ax_ts_e.set_ylabel("T [K]"); ax_ts_e.grid(True, alpha=0.25)
    ax_ps_e.set_xlabel("S [kJ/(kg·K)]"); ax_ps_e.set_ylabel("P [kPa]"); ax_ps_e.grid(True, alpha=0.25)
    fig_e.tight_layout(); fig_e.savefig(TESTS_DIR / "opt_subcycle_grid.png", dpi=150); plt.close(fig_e)

    print(f"  最优解图已保存: opt_cycle_ts_ps.png, opt_pinch_composition.png, opt_performance.png, opt_system_pinch.png, opt_subcycle_grid.png")

test_optimize_min_max_utility()