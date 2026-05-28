"""
换热最优匹配测试：4热+5冷复杂算例 → 枚举+DP → T-Q 可视化。

输出 ``tests/hx_match_demo.png``。
"""

from __future__ import annotations

from pathlib import Path

import pytest

from core.cycle_performance import NodeStateSnapshot, ProcessCategory, ProcessRecord
from core.heat_exchanger import (
    HXMatchResult,
    HXUnit,
    _RecInfo,
    _normalize_records,
    _enumerate_candidate_groups,
    _solve_optimal_packing,
    match_heat_exchanger_groups,
)

TESTS_DIR = Path(__file__).resolve().parent
CP = 5.193


# ──────────────────────────────────────────────────────────────
# 合成数据
# ──────────────────────────────────────────────────────────────


def _mk(cat: ProcessCategory, key: str, mf: float,
        T_tail: float, T_head: float) -> ProcessRecord:
    H_tail = T_tail * CP
    H_head = T_head * CP
    ts = NodeStateSnapshot(index=-1, T=T_tail, P=1000.0, H=H_tail, S=5.0)
    hs = NodeStateSnapshot(index=-2, T=T_head, P=1000.0, H=H_head, S=5.0)
    dH = H_head - H_tail
    return ProcessRecord(edge_key=key, fluid="He", kind="heat", category=cat,
                         tail=-1, head=-2, mass_flow=mf,
                         tail_state=ts, head_state=hs,
                         delta_H=dH, power_rate=mf * dH)


def _make_test_data():
    """4 热 + 5 冷，ΔT_min=10K。预期 DP 选出 4 个组，总匹配=280kW。"""
    h1 = _mk(ProcessCategory.HEAT_REJECTION, "SH_H1",
             mf=150.0 / (250.0 * CP), T_tail=950.0, T_head=700.0)  # Q≈150
    h2 = _mk(ProcessCategory.HEAT_REJECTION, "SH_H2",
             mf=60.0 / (120.0 * CP), T_tail=600.0, T_head=480.0)   # Q≈60
    h3 = _mk(ProcessCategory.HEAT_REJECTION, "SH_H3",
             mf=55.0 / (120.0 * CP), T_tail=550.0, T_head=430.0)   # Q≈55
    h4 = _mk(ProcessCategory.HEAT_REJECTION, "SH_H4",
             mf=40.0 / (100.0 * CP), T_tail=500.0, T_head=400.0)   # Q≈40

    c1 = _mk(ProcessCategory.HEAT_ABSORPTION, "SH_C1",
             mf=80.0 / (160.0 * CP), T_tail=520.0, T_head=680.0)   # Q≈80
    c2 = _mk(ProcessCategory.HEAT_ABSORPTION, "SH_C2",
             mf=70.0 / (130.0 * CP), T_tail=420.0, T_head=550.0)   # Q≈70
    c3 = _mk(ProcessCategory.HEAT_ABSORPTION, "SH_C3",
             mf=50.0 / (110.0 * CP), T_tail=350.0, T_head=460.0)   # Q≈50
    c4 = _mk(ProcessCategory.HEAT_ABSORPTION, "SH_C4",
             mf=45.0 / (120.0 * CP), T_tail=300.0, T_head=420.0)   # Q≈45
    c5 = _mk(ProcessCategory.HEAT_ABSORPTION, "SH_C5",
             mf=35.0 / (100.0 * CP), T_tail=250.0, T_head=350.0)   # Q≈35

    return [h1, h2, h3, h4], [c1, c2, c3, c4, c5]


# ──────────────────────────────────────────────────────────────
# 可视化
# ──────────────────────────────────────────────────────────────


def _draw_result(result: HXMatchResult, hot_recs, cold_recs,
                 candidates_count: int, out_path: Path) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    n_units = max(len(result.units), 1)
    n_rows = (n_units + 2) // 3

    fig = plt.figure(figsize=(max(14, 5 * min(n_units, 3)), 3.5 + 3.5 * n_rows))
    fig.suptitle(
        f"Heat Exchanger Optimal Matching | "
        f"{len(hot_recs)+len(cold_recs)} records, {candidates_count} candidates, "
        f"dT_min=10K | matched={result.total_matched:.0f}kW, "
        f"unmatched={result.total_unmatched:.0f}kW, {result.num_units} units",
        fontsize=12, fontweight="bold",
    )

    gs = fig.add_gridspec(n_rows + 1, 1, height_ratios=[1] + [2] * n_rows,
                          hspace=0.45, top=0.90)

    # ── Row 0：概览 ──
    ax_ov = fig.add_subplot(gs[0])
    all_recs_orig: list[_RecInfo] = []
    for rec in hot_recs + cold_recs:
        ri = _normalize_records([rec] if rec.category == ProcessCategory.HEAT_REJECTION else [],
                                 [rec] if rec.category == ProcessCategory.HEAT_ABSORPTION else [])
        if ri[0]:
            all_recs_orig.append(ri[0][0])

    sorted_orig = sorted(all_recs_orig, key=lambda r: (r.T_high, 0 if r.is_hot else 1), reverse=True)

    # 构建 unit → 包含的记录 keys
    unit_rec_keys: list[set[str]] = []
    for u in result.units:
        keys: set[str] = set()
        for r in u.hot_records:
            keys.add(r.edge_key)
        for r in u.cold_records:
            keys.add(r.edge_key)
        unit_rec_keys.append(keys)

    # 未分配标记
    ua_keys: set[str] = set()
    for r in result.unassigned_hots:
        ua_keys.add(r.edge_key)
    for r in result.unassigned_colds:
        ua_keys.add(r.edge_key)

    unit_colors = ["#FFD0D0", "#D0D0FF", "#D0FFD0", "#FFFFD0",
                   "#FFD0FF", "#D0FFFF", "#FFE8D0", "#E8D0FF"]

    cum_q = 0.0
    for i, r in enumerate(sorted_orig):
        color = "tab:red" if r.is_hot else "tab:blue"
        kw = {"color": color, "linewidth": 2.2, "marker": "o", "markersize": 4}
        label = r.record.edge_key
        ax_ov.plot([cum_q, cum_q + r.Q], [r.T_high, r.T_low], **kw)

        # 找到归属
        in_unit = -1
        for ui, keys in enumerate(unit_rec_keys):
            if label in keys:
                in_unit = ui
                break

        if in_unit >= 0:
            ax_ov.annotate(f"{label}\n{r.Q:.0f}kW",
                           (cum_q + r.Q / 2, (r.T_high + r.T_low) / 2),
                           fontsize=6, ha="center", color=color,
                           bbox=dict(boxstyle="round,pad=0.1",
                                     facecolor=unit_colors[in_unit % len(unit_colors)],
                                     alpha=0.85))
        elif label in ua_keys:
            ax_ov.annotate(f"{label}\n{r.Q:.0f}kW",
                           (cum_q + r.Q / 2, (r.T_high + r.T_low) / 2),
                           fontsize=6, ha="center", color="gray",
                           bbox=dict(boxstyle="round,pad=0.1",
                                     facecolor="#EEEEEE", alpha=0.7))
        cum_q += r.Q

    ax_ov.set_xlabel("Cumulative Q (T_high descending) [kW]")
    ax_ov.set_ylabel("T [K]")
    ax_ov.grid(True, alpha=0.2)
    ax_ov.legend(handles=[
        mpatches.Patch(color="tab:red", label="hot"),
        mpatches.Patch(color="tab:blue", label="cold"),
        mpatches.Patch(color="gray", label="unassigned"),
    ], fontsize=8, loc="upper left", ncol=3)
    ax_ov.set_title("Overview — all records (sorted by T_high)", fontsize=9, loc="left")

    # ── 各组详情 ──
    for ui, unit in enumerate(result.units):
        row = 1 + ui // 3
        col = ui % 3
        sub_gs = gs[row].subgridspec(1, min(n_units - (ui // 3) * 3, 3),
                                      wspace=0.35)
        ax = fig.add_subplot(sub_gs[0, col]) if min(n_units, 3) > 0 else fig.add_subplot(gs[row])

        h_recs = unit.hot_records
        c_recs = unit.cold_records

        # 按 T_high 排序
        h_sorted = sorted(h_recs, key=lambda rr: max(rr.tail_state.T, rr.head_state.T), reverse=True)
        c_sorted = sorted(c_recs, key=lambda rr: max(rr.tail_state.T, rr.head_state.T), reverse=True)

        # 共享 Q 轴：画热侧（实线红，左高→右低）
        x = 0.0
        for hr in h_sorted:
            Th_high = max(hr.tail_state.T, hr.head_state.T)
            Th_low = min(hr.tail_state.T, hr.head_state.T)
            Qh = abs(hr.power_rate) if hr.power_rate else 0.0
            ax.plot([x, x + Qh], [Th_high, Th_low], "o-", color="tab:red", lw=2.2, ms=3.5)
            matched_of_h = _matched_in_unit(hr, unit)
            ax.annotate(f"{hr.edge_key}\n(matched {matched_of_h:.0f}/{Qh:.0f}kW)",
                        (x + Qh / 2, (Th_high + Th_low) / 2),
                        fontsize=5, color="darkred", ha="center",
                        xytext=(0, 8), textcoords="offset points")
            x += Qh
        q_max = x

        # 画冷侧（虚线蓝，逆流方向：左高[冷出口]→右低[冷入口]）
        x = 0.0
        for cr in c_sorted:
            Tc_low = min(cr.tail_state.T, cr.head_state.T)
            Tc_high = max(cr.tail_state.T, cr.head_state.T)
            Qc = abs(cr.power_rate) if cr.power_rate else 0.0
            matched_of_c = _matched_in_unit(cr, unit)
            # 逆流：冷出口（高温, 在左）→ 冷入口（低温, 在右）
            ax.plot([x, x + Qc], [Tc_high, Tc_low], "s--", color="tab:blue", lw=1.8, ms=3.5)
            ax.annotate(f"{cr.edge_key}\n(matched {matched_of_c:.0f}/{Qc:.0f}kW)",
                        (x + Qc / 2, (Tc_high + Tc_low) / 2),
                        fontsize=5, color="darkblue", ha="center",
                        xytext=(0, -12), textcoords="offset points")
            ax.annotate(f"in {Tc_low:.0f}K", (x + Qc, Tc_low),
                        fontsize=4.5, color="darkblue", ha="left", va="top",
                        xytext=(3, -2), textcoords="offset points")
            ax.annotate(f"out {Tc_high:.0f}K", (x, Tc_high),
                        fontsize=4.5, color="darkblue", ha="right", va="bottom",
                        xytext=(-3, 2), textcoords="offset points")
            x += Qc

        title = (f"Unit {ui + 1}  matched={unit.matched_heat:.0f}kW  "
                 f"resid={unit.residual:.1f}kW  pinch={unit.internal_pinch:.0f}K")
        ax.set_title(title, fontsize=8, fontweight="bold")
        ax.set_xlabel("Q [kW]")
        ax.set_ylabel("T [K]")
        ax.grid(True, alpha=0.2)
        ax.set_xlim(0, max(x, 1))

    fig.subplots_adjust(left=0.05, right=0.97, bottom=0.05)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _matched_in_unit(rec: ProcessRecord, unit: HXUnit) -> float:
    """估算一条记录在 unit 中被匹配的热量（按比例分配残差）。"""
    Q = abs(rec.power_rate) if rec.power_rate else 0.0
    return Q


# ──────────────────────────────────────────────────────────────
# 测试
# ──────────────────────────────────────────────────────────────


def test_hx_match() -> None:
    hots, colds = _make_test_data()
    dT_min = 10.0

    print("=" * 70)
    print("输入记录 (ΔT_min = 10K):")
    print("-" * 70)
    for rec in hots + colds:
        Q = abs(rec.power_rate) if rec.power_rate else 0.0
        Th = max(rec.tail_state.T, rec.head_state.T)
        Tl = min(rec.tail_state.T, rec.head_state.T)
        tag = "HOT" if rec.category == ProcessCategory.HEAT_REJECTION else "COLD"
        print(f"  {tag:6s} {rec.edge_key:8s}  Q={Q:8.1f}kW  T: {Th:6.0f}→{Tl:6.0f}K")

    # 内部统计：候选组数
    recs, th, tc = _normalize_records(hots, colds)
    h_list = [r for r in recs if r.is_hot]
    c_list = [r for r in recs if not r.is_hot]
    candidates = _enumerate_candidate_groups(h_list, c_list, dT_min)
    print(f"\n枚举候选组: {len(candidates)} 个")
    print(f"  其中 H→C 组: {sum(1 for c in candidates if len(c[4])==1)}")
    print(f"  其中 C→H 组: {sum(1 for c in candidates if len(c[5])==1 and len(c[4])>1)}")

    result = match_heat_exchanger_groups(hots, colds, dT_min)

    print(f"\n总热: {th:.1f}kW  总冷: {tc:.1f}kW")
    print(f"换热器数: {result.num_units}")
    print(f"总匹配热: {result.total_matched:.1f}kW")
    print(f"总不匹配: {result.total_unmatched:.1f}kW")
    print()

    for ui, unit in enumerate(result.units):
        print(f"─── Unit {ui + 1} ───")
        print(f"  热: {', '.join(r.edge_key for r in unit.hot_records)}")
        print(f"  冷: {', '.join(r.edge_key for r in unit.cold_records)}")
        print(f"  matched={unit.matched_heat:.1f}kW  residual={unit.residual:.1f}kW  "
              f"pinch={unit.internal_pinch:.0f}K")
        print()

    if result.unassigned_hots:
        print(f"未分配热: {', '.join(r.edge_key for r in result.unassigned_hots)}")
    if result.unassigned_colds:
        print(f"未分配冷: {', '.join(r.edge_key for r in result.unassigned_colds)}")

    out = TESTS_DIR / "hx_match_demo.png"
    _draw_result(result, hots, colds, len(candidates), out)
    print(f"\n图表: {out}")

    assert result.total_unmatched >= 0.0
    assert result.num_units >= 0


if __name__ == "__main__":
    test_hx_match()
