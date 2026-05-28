"""hx_unmatched 优化 + HX 匹配测试: Air/H₂/He 系统, CMA-ES 200 代, skip_pinch。

HX 匹配已集成入目标函数 (dT_min=10K)；优化后打印匹配明细并输出 3 张 PNG。
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import pytest

from core import (
    ClosedCycleLayer,
    ClosedCycleTPInput,
    CycleConfig,
    ExternalSourceInput,
    ProcessCategory,
    ProcessRecord,
    PropertyRegistry,
    SystemInput,
)
from core.heat_exchanger import match_heat_exchanger_groups, HXMatchResult
from optimize import Optimizer

TESTS_DIR = Path(__file__).resolve().parent


# ============================================================
# 输入构造（与 test_optimize 一致）
# ============================================================


def _make_system_input() -> SystemInput:
    hot = ExternalSourceInput(
        fluid="Air", mass_flow=100.0,
        T_in=1250.0, P_in=200.0, T_out=500.0, P_out=180.0,
    )
    cold = ExternalSourceInput(
        fluid="Hydrogen", mass_flow=4.3,
        T_in=20.0, P_in=5000.0, T_out=900.0, P_out=4500.0,
    )
    cycle_input = ClosedCycleTPInput(
        fluid="He", t_min=40.0, t_max=1000.0, p_min=2000.0, p_max=10000.0,
        t_quantiles=(0.33, 0.67), p_quantiles=(0.5,),
        subcycle_mass_flow_initial=20.0,
    )
    return SystemInput(
        heat_sources=(hot,), cold_sources=(cold,),
        cycles=(CycleConfig(
            input=cycle_input, use_non_ideal=True,
            delta_T_min=20.0, heat_method=None,
        ),),
        delta_T_min=20.0, heat_method="system_pinch",
    )


# ============================================================
# 图 1: 理想骨架 TP/PS（子循环多边形 + 流量）
# ============================================================


def _draw_ideal_grid(layer: ClosedCycleLayer, out_path: Path) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax_ts, ax_ps) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("Ideal Topology — Subcycle Grid & Flows", fontsize=12, fontweight="bold")

    colors_sc = plt.cm.tab10.colors
    for n in layer.nodes.values():
        ax_ts.scatter(n.S, n.T, s=10, color="0.4", zorder=1)
        ax_ps.scatter(n.S, n.P, s=10, color="0.4", zorder=1)

    for i, sc in enumerate(layer.subcycles):
        c = colors_sc[i % len(colors_sc)]
        q = sc.mass_flow if sc.mass_flow is not None else 0.0
        n0, n1, n2, n3 = [layer.nodes[idx] for idx in sc.nodes]
        for ax, y_fn in [(ax_ts, lambda n: n.T), (ax_ps, lambda n: n.P)]:
            xs = [n0.S, n1.S, n2.S, n3.S, n0.S]
            ys = [y_fn(n0), y_fn(n1), y_fn(n2), y_fn(n3), y_fn(n0)]
            ax.plot(xs, ys, "-", color=c, lw=1.2, alpha=0.5)
            cx = (n0.S + n2.S) / 2
            cy = (y_fn(n0) + y_fn(n2)) / 2
            ax.annotate(f"SC{i}\n{q:.1f}", (cx, cy), ha="center", va="center",
                        fontsize=6, color=c, fontweight="bold")

    ax_ts.set_xlabel("S [kJ/(kg·K)]"); ax_ts.set_ylabel("T [K]")
    ax_ts.grid(True, alpha=0.25); ax_ts.set_title("T–S", fontsize=10)
    ax_ps.set_xlabel("S [kJ/(kg·K)]"); ax_ps.set_ylabel("P [kPa]")
    ax_ps.grid(True, alpha=0.25); ax_ps.set_title("P–S", fontsize=10)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ============================================================
# 图 2: 非理想 TS + PS（kept_nodes + simplified_edges）
# ============================================================


def _draw_non_ideal(report, out_path: Path) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax_ts, ax_ps) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("Non-Ideal Cycle — T–S & P–S", fontsize=12, fontweight="bold")

    for _, snap in report.nodes:
        ax_ts.scatter(snap.S, snap.T, s=12, color="0.3", zorder=2)
        ax_ps.scatter(snap.S, snap.P, s=12, color="0.3", zorder=2)

    for _, rec in report.by_edge:
        c = "tab:orange" if rec.kind == "mechanical" else "tab:green"
        ax_ts.plot([rec.tail_state.S, rec.head_state.S],
                    [rec.tail_state.T, rec.head_state.T],
                    color=c, lw=1.2, alpha=0.7)
        ax_ps.plot([rec.tail_state.S, rec.head_state.S],
                    [rec.tail_state.P, rec.head_state.P],
                    color=c, lw=1.2, alpha=0.7)
        fr, to = 0.38, 0.62
        ax_ts.annotate("", xytext=(
            rec.tail_state.S + fr * (rec.head_state.S - rec.tail_state.S),
            rec.tail_state.T + fr * (rec.head_state.T - rec.tail_state.T),
        ), xy=(
            rec.tail_state.S + to * (rec.head_state.S - rec.tail_state.S),
            rec.tail_state.T + to * (rec.head_state.T - rec.tail_state.T),
        ), arrowprops=dict(arrowstyle="->", color="0.3", lw=1.0, alpha=0.6))
        ax_ps.annotate("", xytext=(
            rec.tail_state.S + fr * (rec.head_state.S - rec.tail_state.S),
            rec.tail_state.P + fr * (rec.head_state.P - rec.tail_state.P),
        ), xy=(
            rec.tail_state.S + to * (rec.head_state.S - rec.tail_state.S),
            rec.tail_state.P + to * (rec.head_state.P - rec.tail_state.P),
        ), arrowprops=dict(arrowstyle="->", color="0.3", lw=1.0, alpha=0.6))
        mx_s = (rec.tail_state.S + rec.head_state.S) / 2
        mx_t = (rec.tail_state.T + rec.head_state.T) / 2
        mx_p = (rec.tail_state.P + rec.head_state.P) / 2
        ax_ts.annotate(rec.edge_key, (mx_s, mx_t), fontsize=5, color="0.4",
                       ha="center", xytext=(0, 4), textcoords="offset points")
        ax_ps.annotate(rec.edge_key, (mx_s, mx_p), fontsize=5, color="0.4",
                       ha="center", xytext=(0, 4), textcoords="offset points")

    ax_ts.set_xlabel("S [kJ/(kg·K)]"); ax_ts.set_ylabel("T [K]")
    ax_ts.grid(True, alpha=0.25); ax_ts.set_title("T–S", fontsize=10)
    ax_ps.set_xlabel("S [kJ/(kg·K)]"); ax_ps.set_ylabel("P [kPa]")
    ax_ps.grid(True, alpha=0.25); ax_ps.set_title("P–S", fontsize=10)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ============================================================
# 图 3: HX 匹配 T-Q 图（概览 + 每组详情）
# ============================================================


def _draw_hx_match(
    result: HXMatchResult,
    hot_recs, cold_recs,
    out_path: Path,
) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    n_units = max(len(result.units), 1)
    n_rows = (n_units + 2) // 3
    if n_rows == 0:
        n_rows = 1

    fig = plt.figure(figsize=(max(14, 5 * min(n_units, 3)), 3.0 + 3.5 * n_rows))
    fig.suptitle(
        f"HX Chain Matching | matched={result.total_matched:.0f}kW  "
        f"unmatched={result.total_unmatched:.0f}kW  {result.num_units} units",
        fontsize=12, fontweight="bold",
    )

    gs = fig.add_gridspec(n_rows + 1, 1, height_ratios=[1] + [2] * n_rows,
                          hspace=0.45, top=0.92)

    # ── 概览（所有记录沿 T_high 降序）──
    ax_ov = fig.add_subplot(gs[0])

    from core.heat_exchanger import _normalize_records
    recs, _, _ = _normalize_records(hot_recs, cold_recs)
    sorted_recs = sorted(recs, key=lambda r: (r.T_high, 0 if r.is_hot else 1), reverse=True)

    unit_colors = ["#FFD0D0", "#D0D0FF", "#D0FFD0", "#FFFFD0",
                   "#FFD0FF", "#D0FFFF", "#FFE8D0", "#E8D0FF"]
    unit_keys: list[set[str]] = []
    for u in result.units:
        keys: set[str] = set()
        for r in u.hot_records:
            keys.add(r.edge_key)
        for r in u.cold_records:
            keys.add(r.edge_key)
        unit_keys.append(keys)

    ua_keys: set[str] = set()
    for r in result.unassigned_hots:
        ua_keys.add(r.edge_key)
    for r in result.unassigned_colds:
        ua_keys.add(r.edge_key)

    cum_q = 0.0
    for r in sorted_recs:
        color = "tab:red" if r.is_hot else "tab:blue"
        label = r.record.edge_key
        ax_ov.plot([cum_q, cum_q + r.Q], [r.T_high, r.T_low],
                    color=color, lw=2.2, marker="o", ms=4)
        mid_q = cum_q + r.Q / 2
        mid_t = (r.T_high + r.T_low) / 2

        in_unit = -1
        for ui, keys in enumerate(unit_keys):
            if label in keys:
                in_unit = ui
                break
        if in_unit >= 0:
            bg = unit_colors[in_unit % len(unit_colors)]
        elif label in ua_keys:
            bg = "#EEEEEE"
        else:
            bg = "white"
        ax_ov.annotate(f"{label} {r.Q:.0f}kW", (mid_q, mid_t),
                        fontsize=5.5, ha="center", color=color,
                        bbox=dict(boxstyle="round,pad=0.1", facecolor=bg, alpha=0.85))
        cum_q += r.Q

    ax_ov.set_xlabel("Cumulative Q (T_high ↓) [kW]")
    ax_ov.set_ylabel("T [K]")
    ax_ov.grid(True, alpha=0.2)
    ax_ov.legend(handles=[
        mpatches.Patch(color="tab:red", label="hot"),
        mpatches.Patch(color="tab:blue", label="cold"),
    ], fontsize=8, loc="upper left")
    ax_ov.set_title("Overview — all records (sorted by T_high)", fontsize=9, loc="left")

    # ── 每组详情 ──
    for ui, unit in enumerate(result.units):
        row = 1 + ui // 3
        col = ui % 3
        n_sub = min(n_units - (ui // 3) * 3, 3)
        sub_gs = gs[row].subgridspec(1, n_sub, wspace=0.35)
        ax = fig.add_subplot(sub_gs[0, col])

        h_sorted = sorted(unit.hot_records,
                          key=lambda rr: max(rr.tail_state.T, rr.head_state.T), reverse=True)
        c_sorted = sorted(unit.cold_records,
                          key=lambda rr: max(rr.tail_state.T, rr.head_state.T), reverse=True)

        x = 0.0
        for hr in h_sorted:
            Th = max(hr.tail_state.T, hr.head_state.T)
            Tl = min(hr.tail_state.T, hr.head_state.T)
            Qh = abs(hr.power_rate) if hr.power_rate else 0.0
            ax.plot([x, x + Qh], [Th, Tl], "o-", color="tab:red", lw=2.2, ms=3)
            ax.annotate(f"{hr.edge_key}\n({Qh:.0f}kW)", (x + Qh / 2, (Th + Tl) / 2),
                        fontsize=5, color="darkred", ha="center",
                        xytext=(0, 6), textcoords="offset points")
            x += Qh

        x = 0.0
        for cr in c_sorted:
            Tc_l = min(cr.tail_state.T, cr.head_state.T)
            Tc_h = max(cr.tail_state.T, cr.head_state.T)
            Qc = abs(cr.power_rate) if cr.power_rate else 0.0
            # 逆流: 冷出口(左,高温) → 冷入口(右,低温)
            ax.plot([x, x + Qc], [Tc_h, Tc_l], "s--", color="tab:blue", lw=1.8, ms=3)
            ax.annotate(f"{cr.edge_key}\n({Qc:.0f}kW)", (x + Qc / 2, (Tc_h + Tc_l) / 2),
                        fontsize=5, color="darkblue", ha="center",
                        xytext=(0, -10), textcoords="offset points")
            x += Qc

        ax.set_title(f"U{ui + 1}  matched={unit.matched_heat:.0f}kW  "
                     f"resid={unit.residual:.1f}kW  pinch={unit.internal_pinch:.0f}K",
                     fontsize=7.5, fontweight="bold")
        ax.set_xlabel("Q [kW]"); ax.set_ylabel("T [K]")
        ax.grid(True, alpha=0.2)

    fig.subplots_adjust(left=0.05, right=0.97, bottom=0.05)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ============================================================
# 辅助: 收集换热 ProcessRecord
# ============================================================


def _collect_heat_records(system_result, props):
    """从 SystemResult 收集所有换热过程 (HOT=REJECTION, COLD=ABSORPTION)。"""
    hots: list[ProcessRecord] = []
    colds: list[ProcessRecord] = []

    from core.fluid_property_solver import PropertyRegistry as PR

    for rec in system_result.heat_source_records:
        hots.append(rec)
    for rec in system_result.cold_source_records:
        colds.append(rec)

    for report in system_result.cycle_reports:
        for _, rec in report.by_edge:
            if rec.kind != "heat" or rec.power_rate is None:
                continue
            if rec.category == ProcessCategory.HEAT_REJECTION:
                hots.append(rec)
            elif rec.category == ProcessCategory.HEAT_ABSORPTION:
                colds.append(rec)

    return hots, colds


# ============================================================
# 主测试
# ============================================================


def test_optimize_heat_balance_hx() -> None:
    sys_inp = _make_system_input()
    props = PropertyRegistry()

    hx_dT = 10.0

    run_dir = TESTS_DIR / "run_hx_opt"
    run_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 60)
    print("hx_unmatched 优化 (HX matching in obj, skip_pinch=True, CMA-ES 200代)")
    print("=" * 60)
    opt = Optimizer(
        base_input=sys_inp, props=props, objective="hx_unmatched",
        mf_step_fraction=0.01, quantile_step=0.01,
        basis_encoding=True, basis_s=4, basis_p=4,
        mf_bounds=(0.0, 50.0), quantile_merge_ratio=0.0,
        n_workers=6, skip_pinch=True, hx_dT_min=hx_dT, hx_max_group_size=3,
    )

    dim = len(opt.bounds)
    print(f"维度: {dim} (t_max + t_min + {opt._n_t_q}t_q + {opt._n_p_q}p_q + {opt._mf_dim}×basis)")
    print(f"CMA-ES maxiter=200 seed=42, {opt._n_workers} workers, "
          f"dT_min={opt._hx_dT_min}K, skip_pinch={opt._skip_pinch}\n")

    best_val_log: list[float] = []

    def _on_gen(gen: int, restart: int, best_x: list[float],
                best_val: float, n_evals: int) -> None:
        best_val_log.append(best_val)
        idx = 0
        t_max, t_min = best_x[idx], best_x[idx + 1]; idx += 2
        t_q_s = ", ".join(f"{best_x[idx + i]:.3f}" for i in range(opt._n_t_q))
        idx += opt._n_t_q
        p_q_s = ", ".join(f"{best_x[idx + i]:.3f}" for i in range(opt._n_p_q))
        idx += opt._n_p_q
        w_mean = sum(best_x[idx:]) / max(1, opt._mf_dim)
        sys.stdout.write(
            f"  gen {gen:3d} | r{restart} | obj={best_val:.5f} | "
            f"t_max={t_max:.0f} t_min={t_min:.0f} | "
            f"t_q=[{t_q_s}] p_q=[{p_q_s}] | w_mean={w_mean:.1f} | evals={n_evals}\n"
        )
        sys.stdout.flush()

    print("  开始优化...\n")
    t0 = time.time()
    opt_result = opt.run(
        method="cma", maxiter=200, seed=42, early_stop=30, sigma0=0.3,
        callback=_on_gen,
    )
    elapsed = time.time() - t0

    raw = opt_result.system_result

    # ── 单次系统评估 + HX 匹配计时 ──
    _hots_tmp, _colds_tmp = _collect_heat_records(raw, props)
    _t0 = time.time()
    for _ in range(100):
        match_heat_exchanger_groups(_hots_tmp, _colds_tmp, dT_min=hx_dT)
    _t_hx_ms = (time.time() - _t0) / 100 * 1000
    print(f"\n  计时: {len(_hots_tmp)+len(_colds_tmp)}条记录, 100次HX匹配 = {_t_hx_ms*100:.0f}ms → 单次 {_t_hx_ms:.2f}ms")
    print(f"  优化总耗时: {elapsed:.1f}s, 评估{opt_result.n_evaluations}次 → 每评 {elapsed/opt_result.n_evaluations*1000:.0f}ms")

    report = raw.cycle_reports[0]
    ct = report.cycle_totals

    print("\n" + "=" * 55)
    print("最优解")
    print("=" * 55)
    idx = 0
    t_max_opt = opt_result.x_opt[idx]; idx += 1
    t_min_opt = opt_result.x_opt[idx]; idx += 1
    t_q_vals = [f"{opt_result.x_opt[idx + i]:.4f}" for i in range(opt._n_t_q)]
    idx += opt._n_t_q
    p_q_vals = [f"{opt_result.x_opt[idx + i]:.4f}" for i in range(opt._n_p_q)]
    idx += opt._n_p_q
    print(f"  t_max / t_min  = {t_max_opt:.1f} / {t_min_opt:.1f} K")
    print(f"  t_q = [{', '.join(t_q_vals)}]")
    print(f"  p_q = [{', '.join(p_q_vals)}]")
    print(f"  objective       = {opt_result.objective:.5f}  (hx_unmatched/total_Q)")
    print(f"  net_mechanical  = {ct.net_mechanical_power:.1f} kW")
    print(f"  net_heat        = {ct.net_heat_rate:.1f} kW")
    print(f"  n_evaluations   = {opt_result.n_evaluations}")
    print(f"  runtime         = {elapsed:.1f}s")

    # 收集换热过程 + HX 匹配（优化时已内部运行，此处重复仅用于打印/绘图）
    hots, colds = _collect_heat_records(raw, props)
    Q_h = sum(abs(float(r.power_rate)) for r in hots if r.power_rate)
    Q_c = sum(abs(float(r.power_rate)) for r in colds if r.power_rate)
    print(f"\n  总放热 = {Q_h:.0f}kW  总吸热 = {Q_c:.0f}kW  "
          f"不平衡 = {abs(Q_h - Q_c):.0f}kW")
    print(f"  hot 记录 ({len(hots)}): {[r.edge_key for r in hots]}")
    print(f"  cold 记录 ({len(colds)}): {[r.edge_key for r in colds]}")

    # ── HX 匹配（计时）──
    t_hx0 = time.time()
    hx_result = match_heat_exchanger_groups(hots, colds, dT_min=hx_dT, max_group_size=3)
    t_hx = (time.time() - t_hx0) * 1000
    total_N = len(hots) + len(colds)
    print(f"\n  ── HX 匹配 (dT_min={hx_dT}K) ──")
    print(f"  记录数: {total_N} (热{len(hots)} + 冷{len(colds)})")
    print(f"  单次 HX 匹配耗时: {t_hx:.1f} ms")
    print(f"  优化总耗时: {elapsed:.1f}s")
    print(f"  每轮评估 ≈ 系统计算 + {t_hx:.1f}ms HX匹配")
    print(f"  换热器: {hx_result.num_units} 组")
    print(f"  matched   = {hx_result.total_matched:.0f} kW")
    print(f"  unmatched = {hx_result.total_unmatched:.0f} kW")
    print(f"  unmatched/total = {hx_result.total_unmatched/(Q_h+Q_c):.4f}")
    for ui, unit in enumerate(hx_result.units):
        h_labels = ", ".join(r.edge_key for r in unit.hot_records)
        c_labels = ", ".join(r.edge_key for r in unit.cold_records)
        print(f"    U{ui + 1}: hot=[{h_labels}]  cold=[{c_labels}]  "
              f"matched={unit.matched_heat:.0f}kW  resid={unit.residual:.1f}kW  "
              f"pinch={unit.internal_pinch:.0f}K")
    if hx_result.unassigned_hots:
        print(f"    未分配热: {[r.edge_key for r in hx_result.unassigned_hots]}"
              f"  ({sum(abs(float(r.power_rate)) for r in hx_result.unassigned_hots if r.power_rate):.0f}kW)")
    if hx_result.unassigned_colds:
        print(f"    未分配冷: {[r.edge_key for r in hx_result.unassigned_colds]}"
              f"  ({sum(abs(float(r.power_rate)) for r in hx_result.unassigned_colds if r.power_rate):.0f}kW)")

    # ── 构建最优参数的理想层 + 非理想层（用于绘图）──
    from optimize.solver import _round_and_dedup
    raw_tq = tuple(opt_result.x_opt[2 + i] for i in range(opt._n_t_q))
    raw_pq = tuple(opt_result.x_opt[2 + opt._n_t_q + i] for i in range(opt._n_p_q))
    tq_dedup = _round_and_dedup(raw_tq, opt._qstep, opt._qmerge)
    pq_dedup = _round_and_dedup(raw_pq, opt._qstep, opt._qmerge)

    tp_opt = ClosedCycleTPInput(
        fluid="He", t_min=t_min_opt, t_max=t_max_opt,
        p_min=opt._base_tp.p_min, p_max=opt._base_tp.p_max,
        t_quantiles=tq_dedup, p_quantiles=pq_dedup,
    )
    layer_ideal = ClosedCycleLayer(tp_opt)
    weights = list(opt_result.x_opt[4:4 + opt._mf_dim])
    mf_opt = opt._decode_basis_flows(layer_ideal, weights)
    layer_ideal.subcycle_mass_flows = mf_opt
    layer_ideal.commit_subcycle_mass_flows_to_topology()
    print(f"\n  子循环数: {len(layer_ideal.subcycles)}, "
          f"流量 (kg/s): {[f'{f:.2f}' for f in mf_opt]}")

    # 非理想层
    ni = layer_ideal.ensure_non_ideal()
    ni.apply_offsets()
    ni_report = layer_ideal.performance_report()

    # ── 图 1: 理想骨架 ──
    _draw_ideal_grid(layer_ideal, run_dir / "01_ideal_grid.png")
    print(f"  [1/3] 理想骨架: {run_dir / '01_ideal_grid.png'}")

    # ── 图 2: 非理想 TS+PS ──
    _draw_non_ideal(ni_report, run_dir / "02_non_ideal_ts_ps.png")
    print(f"  [2/3] 非理想 TS+PS: {run_dir / '02_non_ideal_ts_ps.png'}")

    # ── 图 3: HX 匹配 ──
    _draw_hx_match(hx_result, hots, colds, run_dir / "03_hx_match.png")
    print(f"  [3/4] HX 匹配: {run_dir / '03_hx_match.png'}")

    # ── 图 4: 收敛曲线 ──
    if best_val_log:
        matplotlib = pytest.importorskip("matplotlib")
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig4, ax4 = plt.subplots(figsize=(8, 4))
        ax4.plot(range(1, len(best_val_log) + 1), best_val_log, "o-", ms=2, lw=1, color="tab:blue")
        ax4.set_xlabel("Generation (callback)")
        ax4.set_ylabel("Objective (hx_unmatched)")
        ax4.set_title("Convergence — 2P quantile, 4x4 basis, CMA-ES 500")
        ax4.grid(True, alpha=0.25)
        ax4.set_ylim(bottom=0)
        if len(best_val_log) > 1:
            ax4.annotate(f"{best_val_log[0]:.4f}", (1, best_val_log[0]),
                         fontsize=7, xytext=(5, 8), textcoords="offset points")
            ax4.annotate(f"{best_val_log[-1]:.4f}", (len(best_val_log), best_val_log[-1]),
                         fontsize=7, xytext=(5, -12), textcoords="offset points")
        conv_path = run_dir / "04_convergence.png"
        fig4.tight_layout()
        fig4.savefig(conv_path, dpi=150)
        plt.close(fig4)
        print(f"  [4/4] 收敛曲线: {conv_path}")

    assert opt_result.n_evaluations > 0
    assert hx_result.total_unmatched >= 0.0


if __name__ == "__main__":
    test_optimize_heat_balance_hx()
