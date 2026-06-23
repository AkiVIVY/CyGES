"""内层优化方法对比: CMA-ES vs L-BFGS-B (多起点LHS + 光滑软惩罚).

固定外层拓扑 (1P1S, p_q=0.5, s_q=0.5), 在同一拓扑上跑两种内层方法,
对比 obj, eta, n_evals, wall_time. 含可视化图表.

运行:
    pytest -s tests/test_inner_compare.py::test_compare_inner
"""

import time
from pathlib import Path

from core import ClosedCycleLayer, ClosedCycleTPInput, PropertyRegistry

from tests.test_layered_opt import (
    _eval_fast,
    _inner_cma_fast,
    _inner_lbfgsb_fast,
    _make_h2_source,
)

_TESTS_DIR = Path(__file__).resolve().parent
_OUT_DIR = _TESTS_DIR / "inner_compare"
_OUT_DIR.mkdir(parents=True, exist_ok=True)

_HP = {
    "n_t_q": 0, "n_p_q": 1, "n_s_q": 1,
    "t_min_lo": 50.0, "t_min_hi": 50.0,
    "t_max_lo": 1000.0, "t_max_hi": 1000.0,
    "p_min_lo": 2000.0, "p_min_hi": 2000.0,
    "p_max_lo": 10000.0, "p_max_hi": 10000.0,
    "t_q_vals": (), "p_q_vals": (0.5,), "s_q_vals": (0.5,),
    "cold_fluid": "Hydrogen",
    "h2_T_in": 20.0, "h2_P_in": 5000.0, "h2_P_out": 4500.0,
    "h2_T_out_lo": 100.0, "h2_T_out_hi": 1000.0,
    "h2_mf_lo": 3.5, "h2_mf_hi": 3.5,
    "air_mf": 100.0, "air_T_in": 1250.0, "air_P_in": 200.0,
    "air_T_out": 500.0, "air_P_out": 180.0,
    "obj_mode": "eff_pinch",
    "util_tol": 1.0, "penalty_w": 100.0,
    "hx_dT": 10.0, "penalty_k": 10,
    "use_non_ideal": False, "use_interp_he": False,
    "mf_lo": 0.0, "mf_hi": 50.0,
    "qstep": 0.001, "flow_step": 0.05, "h2_step": 0.01,
}

_CMA_BUDGETS = [100, 250, 500, 1000]
_LBFGSB_STARTS = [16]
_LBFGSB_WORKERS = 16


def _build_fixed_layer(hp: dict) -> ClosedCycleLayer:
    tp_in = ClosedCycleTPInput(
        fluid="He", t_min=hp["t_min_lo"], t_max=hp["t_max_hi"],
        p_min=hp["p_min_lo"], p_max=hp["p_max_hi"],
        t_quantiles=hp["t_q_vals"], p_quantiles=hp["p_q_vals"],
        s_quantiles=hp["s_q_vals"],
    )
    return ClosedCycleLayer(tp_in)


def test_compare_inner() -> None:
    """对比 CMA-ES 与 L-BFGS-B (多起点 LHS + softplus)。"""
    hp = dict(_HP)
    layer = _build_fixed_layer(hp)
    n_sc = len(layer.subcycles)
    dim = n_sc + (1 if hp["h2_T_out_lo"] != hp["h2_T_out_hi"] else 0)

    from core.system import SystemInput, CycleConfig
    from core import ExternalSourceInput
    hot = ExternalSourceInput(fluid="Air", mass_flow=hp["air_mf"],
                              T_in=hp["air_T_in"], P_in=hp["air_P_in"],
                              T_out=hp["air_T_out"], P_out=hp["air_P_out"])
    cold = ExternalSourceInput(fluid=hp["cold_fluid"], mass_flow=hp["h2_mf_lo"],
                               T_in=hp["h2_T_in"], P_in=hp["h2_P_in"],
                               T_out=hp["h2_T_out_lo"], P_out=hp["h2_P_out"])
    cycle_in = ClosedCycleTPInput(
        fluid="He", t_min=hp["t_min_lo"], t_max=hp["t_max_hi"],
        p_min=hp["p_min_lo"], p_max=hp["p_max_hi"],
        t_quantiles=hp["t_q_vals"], p_quantiles=hp["p_q_vals"],
        s_quantiles=hp["s_q_vals"], subcycle_mass_flow_initial=20.0,
    )
    sys_inp = SystemInput(
        heat_sources=(hot,), cold_sources=(cold,),
        cycles=(CycleConfig(input=cycle_in, use_non_ideal=hp["use_non_ideal"],
                             delta_T_min=20.0, heat_method=None),),
        delta_T_min=20.0, heat_method="system_pinch",
    )

    results: list[dict] = []

    # ── CMA-ES ──
    for budget in _CMA_BUDGETS:
        hp_cma = dict(hp)
        hp_cma["restarts_inner"] = 3
        hp_cma["early_stop"] = 10
        hp_cma["sigma0"] = 20.0
        hp_cma["inner_method"] = "cma"
        avg_pop = (max(4, int(4 + 3 * __import__("math").log2(max(dim, 2)))) +
                   max(10, min(dim * dim, 50))) / 2
        hp_cma["maxiter_inner"] = max(2, round(budget * dim / (avg_pop * 3)))
        label = f"CMA-b{budget}"
        t0 = time.perf_counter()
        result, _ = _inner_cma_fast(cycle_in, sys_inp, hp_cma, seed=42)
        elapsed = time.perf_counter() - t0
        results.append({
            "method": label, "kind": "CMA",
            "obj": result.obj, "eta": -result.obj,
            "n_evals": result.n_evals, "time_s": elapsed,
            "h2_T_out": result.h2_T_out, "flows": result.flows,
        })

    # ── L-BFGS-B ──
    best_lb_layer = None
    best_lb_result = None
    for n_starts in _LBFGSB_STARTS:
        hp_lb = dict(hp)
        hp_lb["lbfgsb_starts"] = n_starts
        hp_lb["lbfgsb_maxiter"] = 50
        hp_lb["lbfgsb_workers"] = _LBFGSB_WORKERS
        hp_lb["inner_method"] = "lbfgsb"
        label = f"LBFGSB-s{n_starts}"
        t0 = time.perf_counter()
        result, layer_out = _inner_lbfgsb_fast(cycle_in, sys_inp, hp_lb, seed=42)
        elapsed = time.perf_counter() - t0
        results.append({
            "method": label, "kind": "LBFGSB",
            "obj": result.obj, "eta": -result.obj,
            "n_evals": result.n_evals, "time_s": elapsed,
            "h2_T_out": result.h2_T_out, "flows": result.flows,
        })
        if n_starts == 15:
            best_lb_layer = layer_out
            best_lb_result = result

    # ── 控制台输出 ──
    print(f"\n{'='*70}")
    print(f"{'Method':<16} {'obj':>10} {'eta':>8} {'evals':>8} {'time_s':>8} {'h2_T_out':>10}")
    print("-" * 64)
    for r in results:
        print(f"{r['method']:<16} {r['obj']:>10.5f} {r['eta']:>8.4f} "
              f"{r['n_evals']:>8} {r['time_s']:>8.1f} {r['h2_T_out']:>8.0f}")
    print(f"{'='*70}\n")

    # ── 可视化 ──
    _draw_comparison(results, n_sc, dim, hp)

    # ── 最佳解 PS 图 + 夹点 T-Q ──
    if best_lb_layer is not None and best_lb_result is not None:
        _draw_best_solution(best_lb_layer, best_lb_result, sys_inp, hp)

    # ── 保存汇总 CSV ──
    _save_csv(results, n_sc, dim, hp)


def _draw_comparison(results: list[dict], n_sc: int, dim: int, hp: dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cmap_cma = "tab:blue"
    cmap_lb = "tab:orange"
    cma = [r for r in results if r["kind"] == "CMA"]
    lb = [r for r in results if r["kind"] == "LBFGSB"]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"Inner Optimizer Comparison | 1P1S H2=3.5kg/s "
                 f"n_sc={n_sc} dim={dim} | softplus k={hp['penalty_k']}",
                 fontsize=12, fontweight="bold")

    # ── 左上: eta vs evals ──
    ax1 = axes[0, 0]
    for r in cma:
        ax1.scatter(r["n_evals"], r["eta"], color=cmap_cma, s=80, zorder=3)
        ax1.annotate(r["method"].replace("CMA-", ""),
                     (r["n_evals"], r["eta"]),
                     textcoords="offset points", xytext=(6, 6), fontsize=8,
                     color=cmap_cma)
    for r in lb:
        ax1.scatter(r["n_evals"], r["eta"], color=cmap_lb, s=80, zorder=3)
        ax1.annotate(r["method"].replace("LBFGSB-", ""),
                     (r["n_evals"], r["eta"]),
                     textcoords="offset points", xytext=(6, 6), fontsize=8,
                     color=cmap_lb)
    ax1.set_xlabel("Total evaluations"); ax1.set_ylabel("eta")
    ax1.set_title("Efficiency vs evaluation budget"); ax1.grid(True, alpha=0.25)
    ax1.legend(["CMA-ES", "L-BFGS-B"], loc="lower right")

    # ── 右上: eta vs time ──
    ax2 = axes[0, 1]
    for r in cma:
        ax2.scatter(r["time_s"], r["eta"], color=cmap_cma, s=80, zorder=3)
        ax2.annotate(r["method"].replace("CMA-", ""),
                     (r["time_s"], r["eta"]),
                     textcoords="offset points", xytext=(6, 6), fontsize=8,
                     color=cmap_cma)
    for r in lb:
        ax2.scatter(r["time_s"], r["eta"], color=cmap_lb, s=80, zorder=3)
        ax2.annotate(r["method"].replace("LBFGSB-", ""),
                     (r["time_s"], r["eta"]),
                     textcoords="offset points", xytext=(6, 6), fontsize=8,
                     color=cmap_lb)
    ax2.set_xlabel("Wall time [s]"); ax2.set_ylabel("eta")
    ax2.set_title("Efficiency vs wall time"); ax2.grid(True, alpha=0.25)

    # ── 左下: bar chart - evals ──
    ax3 = axes[1, 0]
    labels = [r["method"] for r in results]
    etas = [r["eta"] for r in results]
    colors = [cmap_cma if r["kind"] == "CMA" else cmap_lb for r in results]
    bars = ax3.bar(range(len(results)), etas, color=colors)
    ax3.set_xticks(range(len(results)))
    ax3.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax3.set_ylabel("eta"); ax3.set_title("Efficiency per method")
    for bar, eta in zip(bars, etas):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                 f"{eta:.3f}", ha="center", fontsize=8)
    ax3.grid(True, axis="y", alpha=0.25)

    # ── 右下: h2_T_out vs eta ──
    ax4 = axes[1, 1]
    for r in cma:
        ax4.scatter(r["h2_T_out"], r["eta"], color=cmap_cma, s=80, zorder=3)
        ax4.annotate(r["method"].replace("CMA-", ""),
                     (r["h2_T_out"], r["eta"]),
                     textcoords="offset points", xytext=(6, 6), fontsize=8,
                     color=cmap_cma)
    for r in lb:
        ax4.scatter(r["h2_T_out"], r["eta"], color=cmap_lb, s=80, zorder=3)
        ax4.annotate(r["method"].replace("LBFGSB-", ""),
                     (r["h2_T_out"], r["eta"]),
                     textcoords="offset points", xytext=(6, 6), fontsize=8,
                     color=cmap_lb)
    ax4.set_xlabel("h2_T_out [K]"); ax4.set_ylabel("eta")
    ax4.set_title("Efficiency vs H2 outlet temperature"); ax4.grid(True, alpha=0.25)

    plt.tight_layout()
    out_path = _OUT_DIR / "inner_compare.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  图表: {out_path}")


def _save_csv(results: list[dict], n_sc: int, dim: int, hp: dict) -> None:
    import csv
    out_path = _OUT_DIR / "inner_compare.csv"
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["method", "kind", "obj", "eta", "n_evals", "time_s",
                     "h2_T_out", "n_sc", "dim", "penalty_k", "flows"])
        for r in results:
            w.writerow([r["method"], r["kind"], r["obj"], r["eta"],
                         r["n_evals"], r["time_s"], r["h2_T_out"],
                         n_sc, dim, hp["penalty_k"],
                         ",".join(f"{v:.1f}" for v in r["flows"])])
    print(f"  CSV: {out_path}")


def _draw_best_solution(layer, result, sys_inp, hp: dict) -> None:
    """对最佳 LBFGSB-s20 解绘制 TS/PS 子循环图 + 夹点 T-Q 图。"""
    from core import PropertyRegistry
    from core.cycle_performance import ProcessCategory
    from core.system import convert_sources
    from core.postprocess import analyze_pinch
    from tests.test_layered_opt import (
        _draw_cycle_ts_ps, _draw_pinch_tq, _make_h2_source,
    )

    props = PropertyRegistry()
    try:
        layer.subcycle_mass_flows = list(result.flows)
        layer.commit_subcycle_mass_flows_to_topology()
    except Exception:
        pass
    report = layer.performance_report()
    cold_src = _make_h2_source(hp, result.h2_mf, result.h2_T_out)

    hots: list = []
    colds: list = []
    for _, rec in report.by_edge:
        if rec.kind == "heat" and rec.power_rate is not None:
            if rec.category == ProcessCategory.HEAT_REJECTION:
                hots.append(rec)
            else:
                colds.append(rec)
    src_hots, _ = convert_sources(sys_inp.heat_sources, (), props)
    _, src_colds = convert_sources((), (cold_src,), props)
    hots[:0] = src_hots
    colds[:0] = src_colds

    pinch = analyze_pinch(colds, hots, hp["hx_dT"], props)

    _draw_cycle_ts_ps(report, layer, "LBFGSB-s20 best",
                       _OUT_DIR / "best_ts_ps.png")
    _draw_pinch_tq(hots, colds, pinch, hp["hx_dT"],
                    _OUT_DIR / "best_pinch_tq.png")
    print(f"  TS/PS: {_OUT_DIR / 'best_ts_ps.png'}")
    print(f"  Pinch: {_OUT_DIR / 'best_pinch_tq.png'}")


def test_1p0s_lbfgsb_multiseed() -> None:
    """1P0S L-BFGS-B 多seed确定性对比: S96, w=8, maxiter=80, ftol=1e-8, 9 seeds."""
    hp = dict(_HP)
    hp["n_p_q"] = 1
    hp["n_s_q"] = 0
    hp["p_q_vals"] = (0.5,)
    hp["s_q_vals"] = ()
    hp["lbfgsb_starts"] = 256
    hp["lbfgsb_maxiter"] = 80
    hp["lbfgsb_workers"] = 8
    hp["lbfgsb_sampler"] = "lhs"
    hp["inner_method"] = "lbfgsb"
    hp["penalty_w"] = 1000.0

    layer = _build_fixed_layer(hp)
    n_sc = len(layer.subcycles)
    dim = n_sc + (1 if hp["h2_T_out_lo"] != hp["h2_T_out_hi"] else 0)

    from core.system import SystemInput, CycleConfig
    from core import ExternalSourceInput
    hot = ExternalSourceInput(fluid="Air", mass_flow=hp["air_mf"],
                              T_in=hp["air_T_in"], P_in=hp["air_P_in"],
                              T_out=hp["air_T_out"], P_out=hp["air_P_out"])
    cold = ExternalSourceInput(fluid=hp["cold_fluid"], mass_flow=hp["h2_mf_lo"],
                               T_in=hp["h2_T_in"], P_in=hp["h2_P_in"],
                               T_out=hp["h2_T_out_lo"], P_out=hp["h2_P_out"])
    cycle_in = ClosedCycleTPInput(
        fluid="He", t_min=hp["t_min_lo"], t_max=hp["t_max_hi"],
        p_min=hp["p_min_lo"], p_max=hp["p_max_hi"],
        t_quantiles=hp["t_q_vals"], p_quantiles=hp["p_q_vals"],
        s_quantiles=hp["s_q_vals"], subcycle_mass_flow_initial=20.0,
    )
    sys_inp = SystemInput(
        heat_sources=(hot,), cold_sources=(cold,),
        cycles=(CycleConfig(input=cycle_in, use_non_ideal=hp["use_non_ideal"],
                             delta_T_min=20.0, heat_method=None),),
        delta_T_min=20.0, heat_method="system_pinch",
    )

    seeds = [0, 100, 200, 300, 400, 500, 600, 700, 800]
    all_results: list[dict] = []

    print(f"\n{'='*90}")
    print(f"1P0S L-BFGS-B S256 w=8 maxiter=80 | penalty_w={hp['penalty_w']} k={hp['penalty_k']} ftol=1e-8")
    print(f"n_sc={n_sc} dim={dim}  seeds={len(seeds)}")
    print(f"{'='*90}")
    print(f"{'seed':>5} {'obj':>10} {'eta':>8} {'evals':>8} {'time_s':>7} "
          f"{'h2_T_out':>9}  {'flows':>20}")
    print("-" * 88)

    for seed in seeds:
        hp_run = dict(hp)
        hp_run["inner_method"] = "lbfgsb"
        t0 = time.perf_counter()
        result, _ = _inner_lbfgsb_fast(cycle_in, sys_inp, hp_run, seed=seed)
        elapsed = time.perf_counter() - t0
        eta_approx = -result.obj
        flows_str = "[" + ",".join(f"{v:.1f}" for v in result.flows) + "]"
        all_results.append({
            "seed": seed, "obj": result.obj, "eta": eta_approx,
            "n_evals": result.n_evals, "time_s": elapsed,
            "h2_T_out": result.h2_T_out, "flows": result.flows,
        })
        print(f"{seed:>5} {result.obj:>10.5f} {eta_approx:>8.4f} "
              f"{result.n_evals:>8} {elapsed:>7.1f} {result.h2_T_out:>8.0f}K  "
              f"{flows_str}")

    etas = [r["eta"] for r in all_results]
    evals_list = [r["n_evals"] for r in all_results]
    times_list = [r["time_s"] for r in all_results]
    print("-" * 88)
    print(f"{'best':>5} {min(r['obj'] for r in all_results):>10.5f} "
          f"{max(etas):>8.4f} {sum(evals_list):>8} {sum(times_list):>7.1f}")
    print(f"{'mean':>5} {'':>10} {sum(etas)/len(etas):>8.4f} "
          f"{sum(evals_list)//len(evals_list):>8} {sum(times_list)/len(times_list):>7.1f}")
    print(f"{'spread':>5} {'':>10} {max(etas)-min(etas):>8.4f}")
    print(f"{'='*90}\n")

    import csv
    csv_path = _OUT_DIR / "1p0s_s256_multiseed.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["seed", "obj", "eta", "n_evals", "time_s", "h2_T_out", "flows"])
        for r in all_results:
            w.writerow([r["seed"], r["obj"], r["eta"], r["n_evals"],
                        r["time_s"], r["h2_T_out"],
                        ",".join(f"{v:.1f}" for v in r["flows"])])
    print(f"  CSV: {csv_path}")


def test_1p0s_consistency_h2_750() -> None:
    """1P0S 内层一致性诊断: h2_T_out=750K 固定, 扫描 LHS starts × seeds。

    目的: 退出 h2_T_out 窄峰后, 观察 flows 盆地结构是否仍是 spread 主因。
    若 750K 下内层纯 flows 的 spread 仍然较大 → flows 自身多盆地。
    若 750K 下 spread 很小 → h2_T_out 是唯一 spread 源 → 新框架方向正确。
    """
    from core.system import SystemInput, CycleConfig
    from core import ExternalSourceInput
    import csv

    hp = dict(_HP)
    hp["n_p_q"] = 1
    hp["n_s_q"] = 0
    hp["p_q_vals"] = (0.5,)
    hp["s_q_vals"] = ()
    hp["lbfgsb_maxiter"] = 80
    hp["lbfgsb_workers"] = 8
    hp["lbfgsb_sampler"] = "lhs"
    hp["inner_method"] = "lbfgsb"
    hp["penalty_w"] = 1000.0

    h2_fixed = 750.0
    hp["h2_T_out_lo"] = h2_fixed
    hp["h2_T_out_hi"] = h2_fixed   # → h2_T_fixed → dim=n_sc

    layer = _build_fixed_layer(hp)
    n_sc = len(layer.subcycles)
    dim = n_sc

    hot = ExternalSourceInput(fluid="Air", mass_flow=hp["air_mf"],
                              T_in=hp["air_T_in"], P_in=hp["air_P_in"],
                              T_out=hp["air_T_out"], P_out=hp["air_P_out"])
    cold = ExternalSourceInput(fluid=hp["cold_fluid"], mass_flow=hp["h2_mf_lo"],
                               T_in=hp["h2_T_in"], P_in=hp["h2_P_in"],
                               T_out=h2_fixed, P_out=hp["h2_P_out"])
    cycle_in = ClosedCycleTPInput(
        fluid="He", t_min=hp["t_min_lo"], t_max=hp["t_max_hi"],
        p_min=hp["p_min_lo"], p_max=hp["p_max_hi"],
        t_quantiles=hp["t_q_vals"], p_quantiles=hp["p_q_vals"],
        s_quantiles=hp["s_q_vals"], subcycle_mass_flow_initial=20.0,
    )
    sys_inp = SystemInput(
        heat_sources=(hot,), cold_sources=(cold,),
        cycles=(CycleConfig(input=cycle_in, use_non_ideal=hp["use_non_ideal"],
                             delta_T_min=20.0, heat_method=None),),
        delta_T_min=20.0, heat_method="system_pinch",
    )

    starts_levels = [16, 32, 48, 64, 96]
    seeds = [0, 100, 200, 300, 400, 500, 600, 700, 800]

    all_summaries: list[dict] = []
    all_seed_results: list[dict] = []

    for n_starts in starts_levels:
        hp_run = dict(hp)
        hp_run["lbfgsb_starts"] = n_starts

        print(f"\n{'='*85}")
        print(f"h2_T_out=750K | LHS S{n_starts} w=8 | dim={dim} | 9 seeds")
        print(f"{'='*85}")
        print(f"{'seed':>5} {'obj':>10} {'eta':>8} {'evals':>8} {'time_s':>7}  {'flows':>20}")
        print("-" * 77)

        seed_vars: list[dict] = []
        for seed in seeds:
            t0 = time.perf_counter()
            result, _ = _inner_lbfgsb_fast(cycle_in, sys_inp, dict(hp_run), seed=seed)
            elapsed = time.perf_counter() - t0
            eta_app = -result.obj
            flows_str = "[" + ",".join(f"{v:.1f}" for v in result.flows) + "]"
            seed_vars.append({
                "seed": seed, "obj": result.obj, "eta": eta_app,
                "n_evals": result.n_evals, "time_s": elapsed,
                "flows": result.flows,
            })
            all_seed_results.append({
                "starts": n_starts, "seed": seed, "obj": result.obj,
                "eta": eta_app, "n_evals": result.n_evals,
                "time_s": elapsed, "flows": result.flows,
            })
            print(f"{seed:>5} {result.obj:>10.5f} {eta_app:>8.4f} "
                  f"{result.n_evals:>8} {elapsed:>7.1f}  {flows_str}")

        etas = [s["eta"] for s in seed_vars]
        ev_list = [s["n_evals"] for s in seed_vars]
        best_eta = max(etas)
        mean_eta = sum(etas) / len(etas)
        spread = max(etas) - min(etas)
        all_summaries.append({
            "starts": n_starts, "best": best_eta, "mean": mean_eta,
            "spread": spread, "evals": sum(ev_list),
        })
        print("-" * 77)
        print(f"  best={best_eta:.4f}  mean={mean_eta:.4f}  spread={spread:.4f}  "
              f"evals={sum(ev_list)}")
        print(f"{'='*85}")

    # ── 汇总表 ──
    print(f"\n{'='*70}")
    print(f"  汇总: h2_T_out=750K 固定, dim={dim} (flows only)")
    print(f"{'='*70}")
    print(f"{'starts':>7} {'best_eta':>10} {'mean_eta':>10} {'spread':>10} {'evals':>10}")
    print("-" * 55)
    for s in all_summaries:
        print(f"{s['starts']:>7} {s['best']:>10.4f} {s['mean']:>10.4f} "
              f"{s['spread']:>10.4f} {s['evals']:>10}")
    print(f"{'='*70}")

    # ── 保存 CSV ──
    csv_path = _OUT_DIR / "consistency_h2_750.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["starts", "seed", "obj", "eta", "n_evals", "time_s", "flows"])
        for r in all_seed_results:
            w.writerow([r["starts"], r["seed"], r["obj"], r["eta"],
                        r["n_evals"], r["time_s"],
                        ",".join(f"{v:.1f}" for v in r["flows"])])
    print(f"\n  CSV: {csv_path}")

    # ── 可视化 ──
    _draw_consistency_sweep(all_summaries, all_seed_results, dim, h2_fixed)


def _draw_consistency_sweep(
    summaries: list[dict], all_records: list[dict],
    dim: int, h2_fixed: float,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    starts_list = [s["starts"] for s in summaries]
    best_list = [s["best"] for s in summaries]
    mean_list = [s["mean"] for s in summaries]
    spread_list = [s["spread"] for s in summaries]

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(
        f"Inner Consistency | 1P0S h2_T_out={h2_fixed:.0f}K dim={dim} "
        f"| LHS starts × seeds sweep",
        fontsize=13, fontweight="bold",
    )

    # ── 左上: best/mean ──
    ax0 = axes[0, 0]
    ax0.plot(starts_list, best_list, "b-o", label="best η")
    ax0.plot(starts_list, mean_list, color="orange", marker="s", label="mean η")
    ax0.set_xlabel("LHS starts"); ax0.set_ylabel("η")
    ax0.set_title("Best & Mean η vs starts")
    ax0.legend(); ax0.grid(True, alpha=0.25)

    # ── 中上: spread ──
    ax1 = axes[0, 1]
    ax1.plot(starts_list, spread_list, "r-D", linewidth=2, markersize=8)
    ax1.set_xlabel("LHS starts"); ax1.set_ylabel("spread (max η − min η)")
    ax1.set_title("Spread vs starts")
    ax1.grid(True, alpha=0.25)
    for x, y in zip(starts_list, spread_list):
        ax1.annotate(f"{y:.4f}", (x, y), textcoords="offset points",
                      xytext=(0, 8), fontsize=9, ha="center")

    # ── 右上: 箱线图 ──
    ax2 = axes[0, 2]
    box_data: list[list[float]] = []
    for n_starts in starts_list:
        box_data.append([r["eta"] for r in all_records if r["starts"] == n_starts])
    bp = ax2.boxplot(box_data, labels=[str(s) for s in starts_list], patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("lightblue")
    ax2.set_xlabel("LHS starts"); ax2.set_ylabel("η")
    ax2.set_title("η distribution per starts level")
    ax2.grid(True, axis="y", alpha=0.25)

    # ── 左下: evals ──
    ax3 = axes[1, 0]
    for n_starts in starts_list:
        recs = [r for r in all_records if r["starts"] == n_starts]
        ev = [r["n_evals"] for r in recs]
        ax3.scatter(np.full(len(ev), n_starts), ev, alpha=0.6, s=30)
    ax3.set_xlabel("LHS starts"); ax3.set_ylabel("n_evals")
    ax3.set_title("Evaluations per seed"); ax3.grid(True, alpha=0.25)

    # ── 中下: flow 分布 (取每个 starts 的 best seed) ──
    ax4 = axes[1, 1]
    flow_all: list[tuple[int, list[float]]] = []
    for n_starts in starts_list:
        recs = [r for r in all_records if r["starts"] == n_starts]
        best_r = max(recs, key=lambda r: r["eta"])
        flow_all.append((n_starts, best_r["flows"]))
    bar_width = 0.15
    x_pos = np.arange(len(flow_all))
    palette = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
    for fi in range(dim):
        vals = [f[1][fi] if fi < len(f[1]) else 0 for f in flow_all]
        ax4.bar(x_pos + fi * bar_width, vals, bar_width,
                label=f"flow{fi + 1}", color=palette[fi % 4])
    ax4.set_xticks(x_pos + bar_width * (dim - 1) / 2)
    ax4.set_xticklabels([f"S{s}" for s in starts_list])
    ax4.set_ylabel("flow [kg/s]"); ax4.set_title("Best-seed flows per starts")
    ax4.legend(fontsize=8); ax4.grid(True, axis="y", alpha=0.25)

    # ── 右下: 迭代散点 (obj vs starts, colored by seed) ──
    ax5 = axes[1, 2]
    from matplotlib.cm import get_cmap
    cmap = get_cmap("tab10")
    for si, seed in enumerate(sorted({r["seed"] for r in all_records})):
        recs = [(r["starts"], r["eta"]) for r in all_records if r["seed"] == seed]
        recs.sort()
        sx = [r[0] for r in recs]
        sy = [r[1] for r in recs]
        ax5.plot(sx, sy, "o-", color=cmap(si % 10), alpha=0.7, label=f"seed{seed}")
    ax5.set_xlabel("LHS starts"); ax5.set_ylabel("η")
    ax5.set_title("η per seed (line trace)"); ax5.grid(True, alpha=0.25)
    ax5.legend(fontsize=6, ncol=2, loc="lower right")

    plt.tight_layout()
    out_path = _OUT_DIR / "consistency_h2_750.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  可视化: {out_path}")


def test_1p0s_flow_precision() -> None:
    """验证内层 flows 的实际精度: h2_T_out=750K 固定, S48, 3 seeds, 打印6位小数。"""
    from core.system import SystemInput, CycleConfig
    from core import ExternalSourceInput

    hp = dict(_HP)
    hp["n_p_q"] = 1
    hp["n_s_q"] = 0
    hp["p_q_vals"] = (0.5,)
    hp["s_q_vals"] = ()
    hp["lbfgsb_starts"] = 48
    hp["lbfgsb_maxiter"] = 80
    hp["lbfgsb_workers"] = 8
    hp["lbfgsb_sampler"] = "lhs"
    hp["inner_method"] = "lbfgsb"
    hp["penalty_w"] = 1000.0
    hp["h2_T_out_lo"] = 750.0
    hp["h2_T_out_hi"] = 750.0

    layer = _build_fixed_layer(hp)
    n_sc = len(layer.subcycles)
    hot = ExternalSourceInput(fluid="Air", mass_flow=hp["air_mf"],
                              T_in=hp["air_T_in"], P_in=hp["air_P_in"],
                              T_out=hp["air_T_out"], P_out=hp["air_P_out"])
    cold = ExternalSourceInput(fluid=hp["cold_fluid"], mass_flow=hp["h2_mf_lo"],
                               T_in=hp["h2_T_in"], P_in=hp["h2_P_in"],
                               T_out=750.0, P_out=hp["h2_P_out"])
    cycle_in = ClosedCycleTPInput(
        fluid="He", t_min=hp["t_min_lo"], t_max=hp["t_max_hi"],
        p_min=hp["p_min_lo"], p_max=hp["p_max_hi"],
        t_quantiles=hp["t_q_vals"], p_quantiles=hp["p_q_vals"],
        s_quantiles=hp["s_q_vals"], subcycle_mass_flow_initial=20.0,
    )
    sys_inp = SystemInput(
        heat_sources=(hot,), cold_sources=(cold,),
        cycles=(CycleConfig(input=cycle_in, use_non_ideal=hp["use_non_ideal"],
                             delta_T_min=20.0, heat_method=None),),
        delta_T_min=20.0, heat_method="system_pinch",
    )

    print(f"\n  flows 精度验证: h2_T_out=750K, n_sc={n_sc}, S48 w=8")
    print(f"  {'seed':>5} {'obj':>12}  {'flows (6 decimal places)':>55}")
    print(f"  {'-'*75}")

    for seed in [0, 100, 500]:
        result, _ = _inner_lbfgsb_fast(cycle_in, sys_inp, dict(hp), seed=seed)
        flows_hi = ", ".join(f"{v:.6f}" for v in result.flows)
        print(f"  {seed:>5} {result.obj:>12.10f}  [{flows_hi}]")
        print(f"         sum={sum(result.flows):.6f}")


def test_1p0s_penalty_sweep() -> None:
    """1P0S 惩罚参数扫描: eff_pinch(w=1000,k=10) vs pinch vs eff_pinch(w=2000,k=5), h2_T_out=750K 固定。"""
    from core.system import SystemInput, CycleConfig
    from core import ExternalSourceInput

    configs = [
        ("eff_pinch w=1000 k=10", {"obj_mode": "eff_pinch", "penalty_w": 1000, "penalty_k": 10}),
        ("pinch direct",          {"obj_mode": "pinch"}),
        ("eff_pinch w=2000 k=5",  {"obj_mode": "eff_pinch", "penalty_w": 2000, "penalty_k": 5}),
    ]

    hp = dict(_HP)
    hp["n_p_q"] = 1
    hp["n_s_q"] = 0
    hp["p_q_vals"] = (0.5,)
    hp["s_q_vals"] = ()
    hp["lbfgsb_starts"] = 48
    hp["lbfgsb_maxiter"] = 80
    hp["lbfgsb_workers"] = 8
    hp["lbfgsb_sampler"] = "lhs"
    hp["inner_method"] = "lbfgsb"
    hp["h2_T_out_lo"] = 750.0
    hp["h2_T_out_hi"] = 750.0

    layer = _build_fixed_layer(hp)
    n_sc = len(layer.subcycles)

    hot = ExternalSourceInput(fluid="Air", mass_flow=hp["air_mf"],
                              T_in=hp["air_T_in"], P_in=hp["air_P_in"],
                              T_out=hp["air_T_out"], P_out=hp["air_P_out"])
    cycle_in = ClosedCycleTPInput(
        fluid="He", t_min=hp["t_min_lo"], t_max=hp["t_max_hi"],
        p_min=hp["p_min_lo"], p_max=hp["p_max_hi"],
        t_quantiles=hp["t_q_vals"], p_quantiles=hp["p_q_vals"],
        s_quantiles=hp["s_q_vals"], subcycle_mass_flow_initial=20.0,
    )
    seeds = [0, 100, 200, 300, 400, 500, 600, 700, 800]

    all_summaries: list[dict] = []

    for label, c_extra in configs:
        hp_c = dict(hp)
        hp_c.update(c_extra)
        cold = ExternalSourceInput(fluid=hp["cold_fluid"], mass_flow=hp["h2_mf_lo"],
                                   T_in=hp["h2_T_in"], P_in=hp["h2_P_in"],
                                   T_out=750.0, P_out=hp["h2_P_out"])
        sys_inp = SystemInput(
            heat_sources=(hot,), cold_sources=(cold,),
            cycles=(CycleConfig(input=cycle_in, use_non_ideal=hp["use_non_ideal"],
                                 delta_T_min=20.0, heat_method=None),),
            delta_T_min=20.0, heat_method="system_pinch",
        )

        print(f"\n{'='*80}")
        print(f"  {label}")
        print(f"{'='*80}")
        print(f"{'seed':>5} {'obj':>12} {'evals':>8}  {'flows':>25}")
        print("-" * 65)

        all_objs: list[float] = []
        for seed in seeds:
            result, _ = _inner_lbfgsb_fast(cycle_in, sys_inp, dict(hp_c), seed=seed)
            flows_str = ",".join(f"{v:.1f}" for v in result.flows)
            all_objs.append(result.obj)
            print(f"{seed:>5} {result.obj:>12.6f} {result.n_evals:>8}  [{flows_str}]")

        best = min(all_objs)
        mean = sum(all_objs) / len(all_objs)
        spread = max(all_objs) - min(all_objs)
        all_summaries.append({"label": label, "best": best, "mean": mean, "spread": spread})
        print("-" * 65)
        print(f"  best={best:.4f}  mean={mean:.4f}  spread={spread:.4f}")
        print(f"{'='*80}")

    print(f"\n{'='*80}")
    print(f"  对比汇总")
    print(f"{'='*80}")
    print(f"  {'config':<30} {'best':>10} {'mean':>10} {'spread':>10}")
    print(f"  {'-'*62}")
    for s in all_summaries:
        print(f"  {s['label']:<30} {s['best']:>10.4f} {s['mean']:>10.4f} {s['spread']:>10.4f}")
    print(f"{'='*80}\n")
    """1P0S 分拆惩罚对比: eff_pinch(合并) vs eff_pinch_split(分开冷热), h2_T_out=750K 固定。"""
    from core.system import SystemInput, CycleConfig
    from core import ExternalSourceInput, ProcessCategory, PropertyRegistry
    from core.postprocess import analyze_pinch
    from core.system import convert_sources

    hp = dict(_HP)
    hp["n_p_q"] = 1
    hp["n_s_q"] = 0
    hp["p_q_vals"] = (0.5,)
    hp["s_q_vals"] = ()
    hp["lbfgsb_starts"] = 48
    hp["lbfgsb_maxiter"] = 80
    hp["lbfgsb_workers"] = 8
    hp["lbfgsb_sampler"] = "lhs"
    hp["inner_method"] = "lbfgsb"
    hp["penalty_w"] = 1000.0
    hp["obj_mode"] = "eff_pinch"  # 合并（对照组）
    hp["h2_T_out_lo"] = 750.0
    hp["h2_T_out_hi"] = 750.0

    layer = _build_fixed_layer(hp)
    n_sc = len(layer.subcycles)
    hot = ExternalSourceInput(fluid="Air", mass_flow=hp["air_mf"],
                              T_in=hp["air_T_in"], P_in=hp["air_P_in"],
                              T_out=hp["air_T_out"], P_out=hp["air_P_out"])
    cycle_in = ClosedCycleTPInput(
        fluid="He", t_min=hp["t_min_lo"], t_max=hp["t_max_hi"],
        p_min=hp["p_min_lo"], p_max=hp["p_max_hi"],
        t_quantiles=hp["t_q_vals"], p_quantiles=hp["p_q_vals"],
        s_quantiles=hp["s_q_vals"], subcycle_mass_flow_initial=20.0,
    )

    seeds = [0, 100, 200, 300, 400, 500, 600, 700, 800]

    # ── 合并惩罚（对照组）──
    hp_merged = dict(hp)
    hp_merged["obj_mode"] = "eff_pinch"
    cold_m = ExternalSourceInput(fluid=hp["cold_fluid"], mass_flow=hp["h2_mf_lo"],
                                 T_in=hp["h2_T_in"], P_in=hp["h2_P_in"],
                                 T_out=750.0, P_out=hp["h2_P_out"])
    sys_m = SystemInput(
        heat_sources=(hot,), cold_sources=(cold_m,),
        cycles=(CycleConfig(input=cycle_in, use_non_ideal=hp["use_non_ideal"],
                             delta_T_min=20.0, heat_method=None),),
        delta_T_min=20.0, heat_method="system_pinch",
    )
    merged_results: list[dict] = []
    print(f"\n{'='*80}")
    print("  合并惩罚 eff_pinch (hot+cold)")
    print(f"{'='*80}")
    for seed in seeds:
        result, _ = _inner_lbfgsb_fast(cycle_in, sys_m, dict(hp_merged), seed=seed)
        merged_results.append({"seed": seed, "obj": result.obj, "eta": -result.obj,
                               "evals": result.n_evals, "flows": list(result.flows)})

    # ── 分拆惩罚 ──
    hp_split = dict(hp)
    hp_split["obj_mode"] = "eff_pinch_split"
    cold_s = ExternalSourceInput(fluid=hp["cold_fluid"], mass_flow=hp["h2_mf_lo"],
                                 T_in=hp["h2_T_in"], P_in=hp["h2_P_in"],
                                 T_out=750.0, P_out=hp["h2_P_out"])
    sys_s = SystemInput(
        heat_sources=(hot,), cold_sources=(cold_s,),
        cycles=(CycleConfig(input=cycle_in, use_non_ideal=hp["use_non_ideal"],
                             delta_T_min=20.0, heat_method=None),),
        delta_T_min=20.0, heat_method="system_pinch",
    )
    split_results: list[dict] = []
    print(f"\n{'='*80}")
    print("  分拆惩罚 eff_pinch_split (hot + cold 分开)")
    print(f"{'='*80}")
    for seed in seeds:
        result, _ = _inner_lbfgsb_fast(cycle_in, sys_s, dict(hp_split), seed=seed)
        split_results.append({"seed": seed, "obj": result.obj, "eta": -result.obj,
                              "evals": result.n_evals, "flows": list(result.flows)})

    # ── 对比 ──
    print(f"\n{'='*80}")
    print(f"  {'seed':>5} {'merged_obj':>12} {'split_obj':>12} {'Δobj':>10}  "
          f"{'merged_flows':>25}  {'split_flows':>25}")
    print(f"  {'-'*78}")
    for i, seed in enumerate(seeds):
        mo = merged_results[i]["obj"]
        so = split_results[i]["obj"]
        mf = ",".join(f"{v:.1f}" for v in merged_results[i]["flows"])
        sf = ",".join(f"{v:.1f}" for v in split_results[i]["flows"])
        print(f"  {seed:>5} {mo:>12.6f} {so:>12.6f} {mo-so:>10.6f}  "
              f"[{mf}]  [{sf}]")

    m_etas = [-r["obj"] for r in merged_results]
    s_etas = [-r["obj"] for r in split_results]
    print(f"  {'-'*78}")
    print(f"  {'':>5} merged: best={max(m_etas):.4f} mean={sum(m_etas)/len(m_etas):.4f} "
          f"spread={max(m_etas)-min(m_etas):.4f}")
    print(f"  {'':>5} split:  best={max(s_etas):.4f} mean={sum(s_etas)/len(s_etas):.4f} "
          f"spread={max(s_etas)-min(s_etas):.4f}")
    print(f"{'='*80}\n")
