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
        t_quantiles=(0.5,), p_quantiles=(0.5,),
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
    """占位测试：验证 Optimizer 可正常创建（具体优化参数待讨论后配置）。"""
    sys_inp = _make_system_input()
    props = PropertyRegistry()
    opt = Optimizer(
        base_input=sys_inp, props=props, objective="hx_unmatched",
        mf_step_fraction=0.01, quantile_step=0.01,
        basis_encoding=True, basis_s=4, basis_p=4,
        mf_bounds=(0.0, 50.0), quantile_merge_ratio=0.0,
        n_workers=1, skip_pinch=True, hx_dT_min=10.0, hx_max_group_size=3,
        use_interp_he=True,
    )
    assert len(opt.bounds) > 0
    assert opt._n_t_q == 1
    assert opt._n_p_q == 1
    print(f"  Optimizer OK: dim={len(opt.bounds)}, objective={opt._objective_name}")


if __name__ == "__main__":
    test_optimize_heat_balance_hx()
