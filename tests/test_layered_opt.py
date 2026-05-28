r"""分层优化测试：外层 LHS 采样分位空间 → 内层 CMA 优化子循环流量。

   外层: t_min, t_max, t_q[], p_q[] — LHS 200点
   内层: mf[N] 直接变量（无基函数）, interp He 快跑
   验证: top-5 用 CoolProp 重新评估

   Usage:  python tests/test_layered_opt.py
"""

from __future__ import annotations

import math
import random
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import pytest

from core import (
    ClosedCycleLayer,
    ClosedCycleTPInput,
    CycleConfig,
    ExternalSourceInput,
    ProcessCategory,
    PropertyRegistry,
    SystemInput,
)
from core.heat_exchanger import match_heat_exchanger_groups, HXMatchResult
from optimize import Optimizer

TESTS_DIR = Path(__file__).resolve().parent

# ──────────────────────────────────────────────────────────────────
# 系统输入
# ──────────────────────────────────────────────────────────────────


def _make_system_input(
    t_quantiles: tuple[float, ...] = (),
    p_quantiles: tuple[float, ...] = (0.33, 0.67),
    h2_mass_flow: float = 4.3,
) -> SystemInput:
    hot = ExternalSourceInput(
        fluid="Air", mass_flow=100.0,
        T_in=1250.0, P_in=200.0, T_out=500.0, P_out=180.0,
    )
    cold = ExternalSourceInput(
        fluid="Hydrogen", mass_flow=h2_mass_flow,
        T_in=20.0, P_in=5000.0, T_out=900.0, P_out=4500.0,
    )
    cycle_input = ClosedCycleTPInput(
        fluid="He", t_min=40.0, t_max=1000.0, p_min=2000.0, p_max=10000.0,
        t_quantiles=t_quantiles, p_quantiles=p_quantiles,
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


# ──────────────────────────────────────────────────────────────────
# LHS 采样
# ──────────────────────────────────────────────────────────────────


def _lhs(n: int, bounds: list[tuple[float, float]], seed: int = 42) -> list[list[float]]:
    """Latin Hypercube Sampling.

    :param n: 样本数
    :param bounds: 各维 (lo, hi) 边界
    :returns: n 个采样点
    """
    d = len(bounds)
    rng = random.Random(seed)
    samples: list[list[float]] = [[0.0] * d for _ in range(n)]
    for dim in range(d):
        lo, hi = bounds[dim]
        # 分成 n 个等距桶，每桶内随机采一点
        bucket_order = list(range(n))
        rng.shuffle(bucket_order)
        for i in range(n):
            b = bucket_order[i]
            lo_i = lo + (hi - lo) * b / n
            hi_i = lo + (hi - lo) * (b + 1) / n
            samples[i][dim] = rng.uniform(lo_i, hi_i)
    return samples


# ──────────────────────────────────────────────────────────────────
# 分层优化器
# ──────────────────────────────────────────────────────────────────


@dataclass
class LayerResult:
    n: int
    t_min: float
    t_max: float
    p_max: float
    t_q: tuple[float, ...]
    p_q: tuple[float, ...]
    flows: list[float]
    h2_mf: float
    obj: float
    n_subcycles: int
    n_evals: int
    runtime: float


def _eval_inner_direct(
    tp_input: ClosedCycleTPInput,
    mf_flows: list[float],
    sys_inp_template: SystemInput,
    props: PropertyRegistry,
    hx_dT: float,
    he_solver=None,
) -> float:
    """固定拓扑下，用给定子循环流量评估 hx_unmatched 目标。"""
    from optimize.objective import _HX_DT_MIN, _HX_MAX_GROUP_SIZE
    _HX_DT_MIN = hx_dT
    _HX_MAX_GROUP_SIZE = 3

    tp = tp_input
    try:
        layer = ClosedCycleLayer(tp, properties=he_solver)
    except Exception:
        return 1e9
    n_sc = len(layer.subcycles)
    if n_sc == 0:
        return 1e9
    trunc = mf_flows[:n_sc]
    h2_mf = mf_flows[n_sc] if len(mf_flows) > n_sc else 4.3

    # 用可变 H2 流量重建冷源
    cold_src = ExternalSourceInput(
        fluid="Hydrogen", mass_flow=h2_mf,
        T_in=20.0, P_in=5000.0, T_out=900.0, P_out=4500.0,
    )

    cfg = CycleConfig(
        input=tp,
        use_non_ideal=True,
        subcycle_mass_flows=trunc,
        delta_T_min=20.0,
        heat_method=None,
    )
    si = SystemInput(
        heat_sources=sys_inp_template.heat_sources,
        cold_sources=(cold_src,),
        cycles=(cfg,),
        delta_T_min=20.0,
        heat_method="system_pinch",
    )
    try:
        from core.system import SystemPipeline
        raw = SystemPipeline(si).run(props, cycle_properties=he_solver)
    except Exception:
        return 1e9

    # hx_unmatched with N+1 penalty
    all_hots = list(raw.heat_source_records)
    all_colds = list(raw.cold_source_records)
    for report in raw.cycle_reports:
        for _, rec in report.by_edge:
            if rec.kind == "heat" and rec.power_rate is not None:
                if rec.category == ProcessCategory.HEAT_REJECTION:
                    all_hots.append(rec)
                else:
                    all_colds.append(rec)
    total_q = sum(abs(float(r.power_rate)) for r in all_hots + all_colds if r.power_rate)
    if total_q < 1e-12:
        return 1.0
    hx = match_heat_exchanger_groups(all_hots, all_colds, dT_min=hx_dT, max_group_size=3)
    return hx.total_unmatched / total_q


def _inner_cma(
    tp_input: ClosedCycleTPInput,
    sys_inp: SystemInput,
    props: PropertyRegistry,
    hx_dT: float,
    *,
    maxiter: int = 20,
    restarts: int = 2,
    early_stop: int = 5,
    sigma0: float = 15.0,
    seed: int = 42,
) -> LayerResult:
    """内层 CMA：固定拓扑，优化子循环流量。"""
    import cma

    from optimize.solver import _round_and_dedup
    from core import InterpolatingHeliumSolver

    t0 = time.perf_counter()
    he_solver = InterpolatingHeliumSolver(
        T_min=tp_input.t_min, T_max=tp_input.t_max,
        P_min=tp_input.p_min, P_max=tp_input.p_max,
    )
    he_solver._build()
    probe = ClosedCycleLayer(tp_input, properties=he_solver)
    n_sc = len(probe.subcycles)
    dim = n_sc + 1  # 子循环流量 + H2 流量

    lo = [0.0] * n_sc + [3.0]   # 子循环 mf∈[0,50], H2 mf∈[3,6]
    hi = [50.0] * n_sc + [6.0]
    stds = [sigma0 * (hi[i] - lo[i]) for i in range(dim)]
    x0_init = [(lo[i] + hi[i]) / 2 for i in range(dim)]

    global_best_val = float("inf")
    global_best_x = x0_init[:]
    total_evals = 0

    for restart in range(restarts):
        if restart % 2 == 0:
            p = max(4, int(4 + 3 * math.log2(max(dim, 2))))
        else:
            p = max(10, min(dim * dim, 50))

        opts: dict = {
            "bounds": [lo, hi],
            "CMA_stds": stds,
            "popsize": p,
            "verbose": -9,
            "seed": seed + restart if seed is not None else None,
        }
        x0 = x0_init if restart == 0 else [random.Random(seed + restart).uniform(lo[j], hi[j]) for j in range(dim)]
        es = cma.CMAEvolutionStrategy(x0, sigma0, opts)
        stall = 0
        gen = 0

        while not es.stop():
            solutions = es.ask()
            scores = []
            for sol in solutions:
                scores.append(_eval_inner_direct(tp_input, list(sol), sys_inp, props, hx_dT, he_solver))
            es.tell(solutions, scores)
            total_evals += len(solutions)
            gen += 1

            improved = False
            for s, sc in zip(solutions, scores):
                if sc < global_best_val:
                    global_best_val = sc
                    global_best_x = list(s)
                    improved = True
                    stall = 0
            if not improved:
                stall += 1
            if gen >= maxiter or stall >= early_stop * len(solutions):
                break

        if callback_progress:
            rt = time.perf_counter() - t0
            sys.stdout.write(
                f"    r{restart}: gen={gen} obj={global_best_val:.4f} "
                f"evals={total_evals} n_sc={n_sc} dim={dim} "
                f"best_flow_sum={sum(global_best_x[:-1]):.0f} H2={global_best_x[-1]:.1f}\n"
            )
            sys.stdout.flush()

    flows = global_best_x[:n_sc]
    flows = [round(f / 0.5) * 0.5 for f in flows]
    h2_mf = round(global_best_x[-1] * 10) / 10  # 0.1 精度
    full_x = flows + [h2_mf]
    obj_final = _eval_inner_direct(tp_input, full_x, sys_inp, props, hx_dT, he_solver)
    runtime = time.perf_counter() - t0

    return LayerResult(
        n=n_sc,
        t_min=tp_input.t_min, t_max=tp_input.t_max, p_max=tp_input.p_max,
        t_q=tp_input.t_quantiles, p_q=tp_input.p_quantiles,
        flows=flows, h2_mf=h2_mf, obj=obj_final,
        n_subcycles=n_sc, n_evals=total_evals, runtime=runtime,
    )


# ──────────────────────────────────────────────────────────────────
# 图表
# ──────────────────────────────────────────────────────────────────


def _draw_layered_charts(
    result: LayerResult,
    sys_inp: SystemInput,
    props: PropertyRegistry,
    run_dir: Path,
    hx_dT: float,
) -> None:
    """输出 2 张图：HX 匹配 + 收敛散点。"""
    from optimize.solver import _round_and_dedup
    from core.system import SystemPipeline, convert_sources
    from core import InterpolatingHeliumSolver

    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    tq_dedup = _round_and_dedup(result.t_q, 0.01, 0.0)
    pq_dedup = _round_and_dedup(result.p_q, 0.01, 0.0)
    tp = ClosedCycleTPInput(
        fluid="He", t_min=result.t_min, t_max=result.t_max,
        p_min=2000.0, p_max=result.p_max,
        t_quantiles=tq_dedup, p_quantiles=pq_dedup,
    )
    he_solver = InterpolatingHeliumSolver(
        T_min=result.t_min, T_max=result.t_max,
        P_min=2000.0, P_max=result.p_max,
    )
    he_solver._build()
    layer = ClosedCycleLayer(tp, properties=he_solver)
    n_sc = len(layer.subcycles)
    # pad or truncate flows to match actual subcycle count
    if len(result.flows) > n_sc:
        flows = result.flows[:n_sc]
    else:
        flows = list(result.flows) + [0.0] * (n_sc - len(result.flows))
    layer.subcycle_mass_flows = flows
    layer.commit_subcycle_mass_flows_to_topology()
    ni = layer.ensure_non_ideal()
    ni.apply_offsets()

    # 直接从已构建 layer 收集换热记录（跳过 SystemPipeline 避免 CoolProp 重建）
    report = layer.performance_report()
    hots: list = []
    colds: list = []
    for _, rec in report.by_edge:
        if rec.kind == "heat" and rec.power_rate is not None:
            if rec.category == ProcessCategory.HEAT_REJECTION:
                hots.append(rec)
            else:
                colds.append(rec)
    # 冷热源转换 — 使用优化后的 H2 流量
    from core.system import convert_sources
    src_hots, _ = convert_sources(sys_inp.heat_sources, (), props)
    # 用优化的 H2 流量重建冷源
    opt_cold = ExternalSourceInput(
        fluid="Hydrogen", mass_flow=result.h2_mf,
        T_in=20.0, P_in=5000.0, T_out=900.0, P_out=4500.0,
    )
    _, src_colds = convert_sources((), (opt_cold,), props)
    hots[:0] = src_hots
    colds[:0] = src_colds
    hx = match_heat_exchanger_groups(hots, colds, dT_min=hx_dT, max_group_size=3)

    # ── 图 1: 理想骨架 ──
    fig1, (ax_ts, ax_ps) = plt.subplots(1, 2, figsize=(14, 6))
    fig1.suptitle(f"Layered Opt — Ideal Topology ({result.n_subcycles} subcycles)", fontsize=12, fontweight="bold")
    colors = plt.cm.tab10.colors
    for n in layer.nodes.values():
        ax_ts.scatter(n.S, n.T, s=10, color="0.4", zorder=1)
        ax_ps.scatter(n.S, n.P, s=10, color="0.4", zorder=1)
    for i, sc in enumerate(layer.subcycles):
        c = colors[i % len(colors)]
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
    fig1.tight_layout(rect=[0, 0, 1, 0.95])
    fig1.savefig(run_dir / "01_ideal_grid.png", dpi=150)
    plt.close(fig1)

    # ── 图 2: HX 匹配 ──
    n_units = max(len(hx.units), 1)
    n_rows = (n_units + 2) // 3
    fig2 = plt.figure(figsize=(max(14, 5 * min(n_units, 3)), 3.0 + 3.5 * n_rows))
    fig2.suptitle(
        f"Layered Best — matched={hx.total_matched:.0f}kW "
        f"unmatched={hx.total_unmatched:.0f}kW {hx.num_units} units | "
        f"obj={result.obj:.4f}",
        fontsize=12, fontweight="bold")
    gs = fig2.add_gridspec(n_rows + 1, 1, height_ratios=[1] + [2] * n_rows, hspace=0.45, top=0.92)
    ax_ov = fig2.add_subplot(gs[0])
    from core.heat_exchanger import _normalize_records
    recs, _, _ = _normalize_records(hots, colds)
    sorted_recs = sorted(recs, key=lambda r: (r.T_high, 0 if r.is_hot else 1), reverse=True)
    cum_q = 0.0
    for r in sorted_recs:
        color = "tab:red" if r.is_hot else "tab:blue"
        ax_ov.plot([cum_q, cum_q + r.Q], [r.T_high, r.T_low], color=color, lw=2.2, marker="o", ms=4)
        mid_q = cum_q + r.Q / 2
        ax_ov.annotate(f"{r.record.edge_key}", (mid_q, (r.T_high + r.T_low) / 2),
                       fontsize=5.5, ha="center", color=color)
        cum_q += r.Q
    ax_ov.set_xlabel("Cumulative Q [kW]"); ax_ov.set_ylabel("T [K]")
    ax_ov.grid(True, alpha=0.2)
    ax_ov.legend(handles=[mpatches.Patch(color="tab:red", label="hot"),
                          mpatches.Patch(color="tab:blue", label="cold")],
                 fontsize=8, loc="upper left")
    for ui, unit in enumerate(hx.units):
        row = 1 + ui // 3
        col = ui % 3
        n_sub = min(n_units - (ui // 3) * 3, 3)
        sub_gs = gs[row].subgridspec(1, n_sub, wspace=0.35)
        ax = fig2.add_subplot(sub_gs[0, col])
        h_sorted = sorted(unit.hot_records, key=lambda rr: max(rr.tail_state.T, rr.head_state.T), reverse=True)
        c_sorted = sorted(unit.cold_records, key=lambda rr: max(rr.tail_state.T, rr.head_state.T), reverse=True)
        x = 0.0
        for hr in h_sorted:
            Th = max(hr.tail_state.T, hr.head_state.T)
            Tl = min(hr.tail_state.T, hr.head_state.T)
            Qh = abs(hr.power_rate) if hr.power_rate else 0.0
            ax.plot([x, x + Qh], [Th, Tl], "o-", color="tab:red", lw=2.2, ms=3)
            x += Qh
        x = 0.0
        for cr in c_sorted:
            Tc_h = max(cr.tail_state.T, cr.head_state.T)
            Tc_l = min(cr.tail_state.T, cr.head_state.T)
            Qc = abs(cr.power_rate) if cr.power_rate else 0.0
            ax.plot([x, x + Qc], [Tc_h, Tc_l], "s--", color="tab:blue", lw=1.8, ms=3)
            x += Qc
        ax.set_title(f"U{ui + 1}  matched={unit.matched_heat:.0f}kW  resid={unit.residual:.1f}kW  pinch={unit.internal_pinch:.0f}K",
                     fontsize=7.5, fontweight="bold")
        ax.set_xlabel("Q [kW]"); ax.set_ylabel("T [K]")
        ax.grid(True, alpha=0.2)
    fig2.tight_layout()
    fig2.savefig(run_dir / "02_hx_match.png", dpi=150)
    plt.close(fig2)

    print(f"  Charts: {run_dir / '01_ideal_grid.png'}, {run_dir / '02_hx_match.png'}")


# ──────────────────────────────────────────────────────────────────
# 回调
# ──────────────────────────────────────────────────────────────────

callback_progress: bool = True


# ──────────────────────────────────────────────────────────────────
# 主流程
# ──────────────────────────────────────────────────────────────────


def test_layered_optimization() -> None:
    """分层优化测试：LHS 外层 + CMA 内层。"""
    # ── 配置 ──
    hp = {"n_t_q": 0, "n_p_q": 2, "n_lhs": 400, "hx_dT": 10.0,
          "maxiter_inner": 20, "restarts_inner": 2}
    print("=" * 60)
    print("分层优化 (LHS + CMA inner)")
    print(f"  T分位={hp['n_t_q']}, P分位={hp['n_p_q']}, LHS={hp['n_lhs']}点")
    print(f"  inner CMA: maxiter={hp['maxiter_inner']}, restarts={hp['restarts_inner']}")
    print(f"  H2 mf: 3-6 kg/s  |  p_max: 8000-15000 kPa")
    print("=" * 60)

    props = PropertyRegistry()
    base = _make_system_input()
    hx_dT = hp["hx_dT"]

    # ── 外层 LHS 采样 ──
    outer_bounds: list[tuple[float, float]] = []
    outer_bounds.append((40.0, 300.0))
    outer_bounds.append((800.0, 1100.0))
    for _ in range(hp["n_t_q"]):
        outer_bounds.append((0.01, 0.99))
    for _ in range(hp["n_p_q"]):
        outer_bounds.append((0.01, 0.99))
    outer_bounds.append((8000.0, 15000.0))  # p_max

    lhs_samples = _lhs(hp["n_lhs"], outer_bounds, seed=42)

    # ── 遍历外层 ──
    global_best = LayerResult(n=0, t_min=0, t_max=0, p_max=0, t_q=(), p_q=(), flows=[], h2_mf=0,
                              obj=float("inf"), n_subcycles=0, n_evals=0, runtime=0)
    all_results: list[LayerResult] = []

    t0_total = time.perf_counter()
    for i, row in enumerate(lhs_samples):
        t_min = round(row[0])
        t_max = round(row[1])
        idx = 2
        t_q_vals = tuple(round(row[idx + j], 4) for j in range(hp["n_t_q"]))
        idx += hp["n_t_q"]
        p_q_vals = tuple(round(row[idx + j], 4) for j in range(hp["n_p_q"]))
        idx += hp["n_p_q"]
        p_max = round(row[idx])

        tp_in = ClosedCycleTPInput(
            fluid="He", t_min=t_min, t_max=t_max,
            p_min=2000.0, p_max=p_max,
            t_quantiles=t_q_vals, p_quantiles=p_q_vals,
        )

        try:
            from core import InterpolatingHeliumSolver
            _he = InterpolatingHeliumSolver(T_min=t_min, T_max=t_max, P_min=2000.0, P_max=p_max)
            _he._build()
            probe = ClosedCycleLayer(tp_in, properties=_he)
        except Exception:
            continue
        if len(probe.subcycles) == 0:
            continue

        sys.stdout.write(f"\n[{i+1}/{hp['n_lhs']}] t_min={t_min} t_max={t_max} "
                         f"t_q={[f'{v:.3f}' for v in t_q_vals]} p_q={[f'{v:.3f}' for v in p_q_vals]} "
                         f"p_max={p_max:.0f} n_sc={len(probe.subcycles)} ")
        sys.stdout.flush()

        res = _inner_cma(tp_in, base, props, hx_dT,
                         maxiter=hp["maxiter_inner"],
                         restarts=hp["restarts_inner"],
                         seed=42 + i)
        all_results.append(res)

        if res.obj < global_best.obj:
            global_best = res
            sys.stdout.write(" ← NEW BEST\n")
        else:
            sys.stdout.write("\n")
        sys.stdout.flush()

    t_total = time.perf_counter() - t0_total

    # ── 汇总 ──
    print("\n" + "=" * 55)
    print(f"分层优化完成: {len(all_results)}/{hp['n_lhs']} 有效解, {t_total:.1f}s")
    best = global_best
    print(f"  t_max={best.t_max:.0f}K  t_min={best.t_min:.0f}K  p_max={best.p_max:.0f}kPa")
    print(f"  t_q={[f'{v:.3f}' for v in best.t_q]}  p_q={[f'{v:.3f}' for v in best.p_q]}")
    print(f"  obj={best.obj:.4f}  n_sc={best.n_subcycles}  flows_sum={sum(best.flows):.0f} kg/s  H2={best.h2_mf:.1f} kg/s")
    print(f"  内层总 eval: {sum(r.n_evals for r in all_results)}")

    # top-5
    sorted_res = sorted([r for r in all_results if r.obj < 1e8], key=lambda r: r.obj)[:5]
    print(f"\n  Top-5:")
    for rank, r in enumerate(sorted_res):
        print(f"  #{rank+1} obj={r.obj:.4f} t_min={r.t_min:.0f} t_max={r.t_max:.0f} "
              f"n_sc={r.n_subcycles} flows_sum={sum(r.flows):.0f} H2={r.h2_mf:.1f}")

    # ── 图表 ──
    run_dir = TESTS_DIR / "run_layered_opt"
    run_dir.mkdir(parents=True, exist_ok=True)
    _draw_layered_charts(best, base, props, run_dir, hx_dT)

    assert best.obj < 1e3
    assert best.n_subcycles > 0


if __name__ == "__main__":
    test_layered_optimization()
