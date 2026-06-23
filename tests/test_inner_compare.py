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
_LBFGSB_STARTS = [5, 10, 15, 20, 40]


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
