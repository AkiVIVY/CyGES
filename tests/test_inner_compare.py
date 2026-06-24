"""内层优化方法对比与诊断: L-BFGS-B 多seed, 一致性诊断, 2D 外层优化方法对比.

测试:
    test_1p0s_lbfgsb_multiseed    — 1P0S L-BFGS-B S16 w=16 dim=n_sc 9 seeds
    test_1p0s_consistency_h2_750  — h2_T_out=750K 固定, LHS starts × seeds sweep
    test_1p0s_2d_opt_compare      — 2D (h2_T_out, p_q) Grid/CMA/DE 方法对比

运行:
    pytest -s tests/test_inner_compare.py::<test_name>
"""

import sys
import time
from pathlib import Path

from core import ClosedCycleLayer, ClosedCycleTPInput, PropertyRegistry

from tests.test_layered_opt import (
    _eval_fast,
    _inner_lbfgsb_fast,
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


def _build_fixed_layer(hp: dict) -> ClosedCycleLayer:
    tp_in = ClosedCycleTPInput(
        fluid="He", t_min=hp["t_min_lo"], t_max=hp["t_max_hi"],
        p_min=hp["p_min_lo"], p_max=hp["p_max_hi"],
        t_quantiles=hp["t_q_vals"], p_quantiles=hp["p_q_vals"],
        s_quantiles=hp["s_q_vals"],
    )
    return ClosedCycleLayer(tp_in)


def test_1p0s_lbfgsb_multiseed() -> None:
    """1P0S L-BFGS-B 内层一致性: S32, w=16, maxiter=20, pinch, 9 seeds."""
    hp = dict(_HP)
    hp["n_p_q"] = 1; hp["n_s_q"] = 0
    hp["p_q_vals"] = (0.5,); hp["s_q_vals"] = ()
    hp["lbfgsb_starts"] = 32
    hp["lbfgsb_maxiter"] = 20; hp["lbfgsb_workers"] = 16
    hp["lbfgsb_sampler"] = "lhs"; hp["inner_method"] = "lbfgsb"
    hp["obj_mode"] = "pinch"
    hp["h2_T_out_lo"] = 800.0; hp["h2_T_out_hi"] = 800.0

    layer = _build_fixed_layer(hp)
    n_sc = len(layer.subcycles)
    dim = n_sc

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
    hp_run = dict(hp)

    print(f"\n{'='*85}")
    print(f"1P0S S32 w=16 maxiter=20 pinch ftol=1e-8 | dim={dim} | 9 seeds")
    print(f"{'='*85}")
    print(f"{'seed':>5} {'obj':>10} {'evals':>8} {'time_s':>7}  {'flows':>20}")
    print("-" * 70)

    for seed_val in seeds:
        t0 = time.perf_counter()
        result, _ = _inner_lbfgsb_fast(cycle_in, sys_inp, dict(hp_run), seed=seed_val)
        elapsed = time.perf_counter() - t0
        flows_str = "[" + ",".join(f"{v:.1f}" for v in result.flows) + "]"
        all_results.append({
            "seed": seed_val, "obj": result.obj, "n_evals": result.n_evals,
            "time_s": elapsed, "flows": list(result.flows),
        })
        print(f"{seed_val:>5} {result.obj:>10.5f} {result.n_evals:>8} {elapsed:>7.1f}  {flows_str}")

    objs = [r["obj"] for r in all_results]
    evals_list = [r["n_evals"] for r in all_results]
    best = min(objs); mean = sum(objs) / len(objs); spread = max(objs) - min(objs)
    print("-" * 70)
    print(f"  best={best:.5f}  mean={mean:.5f}  spread={spread:.5f}  evals={sum(evals_list)}")
    print(f"{'='*85}\n")

    import csv
    csv_path = _OUT_DIR / "1p0s_s32_multiseed.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["seed", "obj", "n_evals", "time_s", "flows"])
        for r in all_results:
            w.writerow([r["seed"], r["obj"], r["n_evals"], r["time_s"],
                        ",".join(f"{v:.1f}" for v in r["flows"])])
    print(f"  CSV: {csv_path}")
    """1P0S L-BFGS-B 多seed: S32, w=16, maxiter=20, pinch, ftol=1e-8, eps sweep."""
    hp = dict(_HP)
    hp["n_p_q"] = 1; hp["n_s_q"] = 0
    hp["p_q_vals"] = (0.5,); hp["s_q_vals"] = ()
    hp["lbfgsb_starts"] = 32
    hp["lbfgsb_maxiter"] = 20; hp["lbfgsb_workers"] = 16
    hp["lbfgsb_sampler"] = "lhs"; hp["inner_method"] = "lbfgsb"
    hp["penalty_w"] = 1000.0
    hp["h2_T_out_lo"] = 800.0; hp["h2_T_out_hi"] = 800.0
    hp["obj_mode"] = "pinch"

    layer = _build_fixed_layer(hp)
    n_sc = len(layer.subcycles)
    dim = n_sc

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
    all_summaries: list[dict] = []

    for eps_val in [None, 1e-3, 1e-5]:
        hp_run = dict(hp)
        hp_run["lbfgsb_starts"] = 32
        if eps_val is not None:
            hp_run["lbfgsb_eps"] = eps_val
        label = f"eps={eps_val}" if eps_val else "eps=default(~1e-8)"
        all_results: list[dict] = []

        print(f"\n{'='*80}")
        print(f"1P0S S32 w=16 maxiter=20 pinch | {label} | 9 seeds")
        print(f"{'='*80}")
        print(f"{'seed':>5} {'obj':>10} {'evals':>8} {'time_s':>7}  {'flows':>20}")
        print("-" * 70)

        for seed_val in seeds:
            t0 = time.perf_counter()
            result, _ = _inner_lbfgsb_fast(cycle_in, sys_inp, dict(hp_run), seed=seed_val)
            elapsed = time.perf_counter() - t0
            flows_str = "[" + ",".join(f"{v:.1f}" for v in result.flows) + "]"
            all_results.append({
                "seed": seed_val, "obj": result.obj, "n_evals": result.n_evals,
                "time_s": elapsed, "flows": list(result.flows),
            })
            print(f"{seed_val:>5} {result.obj:>10.5f} {result.n_evals:>8} {elapsed:>7.1f}  {flows_str}")

        objs = [r["obj"] for r in all_results]
        evals_list = [r["n_evals"] for r in all_results]
        best = min(objs); mean = sum(objs) / len(objs); spread = max(objs) - min(objs)
        all_summaries.append({"label": label, "best": best, "mean": mean,
                               "spread": spread, "evals": sum(evals_list)})
        print("-" * 70)
        print(f"  best={best:.5f}  mean={mean:.5f}  spread={spread:.5f}  evals={sum(evals_list)}")
        print(f"{'='*80}")

    print(f"\n{'='*60}")
    print(f"  eps 对比 (S32)")
    print(f"{'='*60}")
    print(f"  {'eps':>18} {'best':>10} {'mean':>10} {'spread':>10} {'evals':>10}")
    print(f"  {'-'*52}")
    for s in all_summaries:
        print(f"  {s['label']:>18} {s['best']:>10.5f} {s['mean']:>10.5f} "
              f"{s['spread']:>10.5f} {s['evals']:>10}")
    print(f"{'='*60}\n")


def test_1p0s_consistency_h2_750() -> None:
    """1P0S 内层一致性诊断: h2_T_out=750K 固定, 扫描 LHS starts × seeds.

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


def test_1p0s_2d_opt_compare() -> None:
    """2D 优化方法对比: Grid / CMA-ES / DE 搜索 (h2_T_out, s_q) — 0P1S。
    内层: pinch(收flows) / 外层 obj: eff_pinch(区分η)。"""
    import scipy.optimize as _sopt
    import cma

    _TOPO = "0P1S"
    hp = dict(_HP)
    hp["n_t_q"] = 0; hp["n_p_q"] = 0; hp["n_s_q"] = 1
    hp["t_min_lo"] = 50.0; hp["t_min_hi"] = 50.0
    hp["t_max_lo"] = 1000.0; hp["t_max_hi"] = 1000.0
    hp["p_min_lo"] = 2000.0; hp["p_min_hi"] = 2000.0
    hp["p_max_lo"] = 10000.0; hp["p_max_hi"] = 10000.0
    hp["t_q_vals"] = (); hp["p_q_vals"] = (); hp["s_q_vals"] = (0.5,)
    hp["lbfgsb_starts"] = 32; hp["lbfgsb_maxiter"] = 20; hp["lbfgsb_workers"] = 16
    hp["lbfgsb_sampler"] = "lhs"; hp["inner_method"] = "lbfgsb"
    hp["penalty_w"] = 1000.0; hp["penalty_k"] = 10; hp["util_tol"] = 1.0
    hp["hx_dT"] = 10.0; hp["use_non_ideal"] = False
    hp["mf_lo"] = 0.0; hp["mf_hi"] = 50.0

    from core import ExternalSourceInput, PropertyRegistry, ProcessCategory
    from core.system import SystemInput, CycleConfig, convert_sources

    out_dir = _OUT_DIR / "2d_opt_0p1s"
    out_dir.mkdir(parents=True, exist_ok=True)

    # ──── 评估函数 ────
    def _eval_outer(h2t: float, sq_val: float, seed: int = 42) -> dict:
        t0 = time.perf_counter()
        tp_in = ClosedCycleTPInput(
            fluid="He", t_min=50.0, t_max=1000.0,
            p_min=2000.0, p_max=10000.0,
            t_quantiles=(), p_quantiles=(), s_quantiles=(sq_val,),
        )
        layer = ClosedCycleLayer(tp_in)
        if len(layer.subcycles) == 0:
            return {"obj": 1e9, "eta": 0, "flows": [], "evals": 0,
                    "time_s": time.perf_counter() - t0}

        hp_inner = dict(hp)
        hp_inner["h2_T_out_lo"] = h2t; hp_inner["h2_T_out_hi"] = h2t
        hp_inner["s_q_vals"] = (sq_val,); hp_inner["obj_mode"] = "pinch"
        hp_inner["inner_method"] = "lbfgsb"

        hot = ExternalSourceInput(fluid="Air", mass_flow=hp["air_mf"],
                                  T_in=hp["air_T_in"], P_in=hp["air_P_in"],
                                  T_out=hp["air_T_out"], P_out=hp["air_P_out"])
        cold = ExternalSourceInput(fluid=hp.get("cold_fluid", "Hydrogen"),
                                   mass_flow=hp["h2_mf_lo"],
                                   T_in=hp["h2_T_in"], P_in=hp["h2_P_in"],
                                   T_out=h2t, P_out=hp["h2_P_out"])
        cycle_in = ClosedCycleTPInput(
            fluid="He", t_min=50.0, t_max=1000.0, p_min=2000.0, p_max=10000.0,
            t_quantiles=(), p_quantiles=(), s_quantiles=(sq_val,),
            subcycle_mass_flow_initial=20.0,
        )
        sys_inp = SystemInput(
            heat_sources=(hot,), cold_sources=(cold,),
            cycles=(CycleConfig(input=cycle_in, use_non_ideal=False,
                                 delta_T_min=20.0, heat_method=None),),
            delta_T_min=20.0, heat_method="system_pinch",
        )
        result, _ = _inner_lbfgsb_fast(tp_in, sys_inp, hp_inner, seed=seed)

        hp_outer = dict(hp)
        hp_outer["obj_mode"] = "eff_pinch"
        hp_outer["h2_T_out_lo"] = h2t; hp_outer["h2_T_out_hi"] = h2t
        hp_outer["s_q_vals"] = (sq_val,)
        outer_obj = 1e9
        try:
            layer2 = ClosedCycleLayer(tp_in)
            layer2.subcycle_mass_flows = list(result.flows)
            layer2.commit_subcycle_mass_flows_to_topology()
            outer_obj = _eval_fast(layer2, hp["h2_mf_lo"], h2t, sys_inp, hp_outer)
        except Exception:
            pass

        props = PropertyRegistry()
        _, lst_colds = convert_sources((), sys_inp.cold_sources, props)
        lst_hots, _ = convert_sources(sys_inp.heat_sources, (), props)
        q_source = sum(abs(float(r.power_rate)) for r in lst_hots if r.power_rate)
        q_cold = sum(abs(float(r.power_rate)) for r in lst_colds if r.power_rate)
        eta = (q_source - q_cold) / q_source if q_source > 1e-12 else 0.0

        return {"obj": outer_obj, "eta": eta, "flows": list(result.flows),
                "evals": result.n_evals, "time_s": time.perf_counter() - t0}

    # ═══════════════════════════════════════════════════════
    # 1. Grid
    # ═══════════════════════════════════════════════════════
    h2t_grid = list(range(600, 1025, 25))
    sq_grid = [round(0.05 + i * 0.1, 2) for i in range(10)]
    print(f"\n  Grid: h2_T_out x s_q = {len(h2t_grid)*len(sq_grid)}pts")
    grid_best = {"obj": float("inf"), "eta": 0.0}
    grid_data: list[dict] = []
    t_grid_start = time.perf_counter()
    for h2t in h2t_grid:
        for sq_i in sq_grid:
            r = _eval_outer(float(h2t), sq_i)
            grid_data.append({"h2t": h2t, "sq": sq_i, **r})
            if r["obj"] < grid_best["obj"] or (r["obj"] == grid_best["obj"] and r["eta"] > grid_best["eta"]):
                grid_best = {"obj": r["obj"], "eta": r["eta"], "h2t": h2t, "sq": sq_i}
            sys.stdout.write(f"\r    ({h2t}K,{sq_i:.2f}) obj={r['obj']:.5f} eta={r['eta']:.4f}")
    sys.stdout.write("\n")
    t_grid = time.perf_counter() - t_grid_start
    print(f"  Grid best: h2tout={grid_best['h2t']}K s_q={grid_best['sq']} obj={grid_best['obj']:.5f} eta={grid_best['eta']:.4f}  time={t_grid:.1f}s")

    # ═══════════════════════════════════════════════════════
    # 2. CMA-ES
    # ═══════════════════════════════════════════════════════
    print(f"\n  CMA-ES 0P1S (sigma=0.15, maxfevals=120)")
    cma_history: list[dict] = []
    cma_best = [float("inf"), 0.0, 800.0, 0.5]
    t_cma = time.perf_counter()

    def _cma_obj_cb(xx):
        r = _eval_outer(xx[0], xx[1])
        cma_history.append({"h2t": xx[0], "sq": xx[1], "obj": r["obj"], "eta": r["eta"],
                            "evals_add": r["evals"]})
        if r["obj"] < cma_best[0] or (r["obj"] == cma_best[0] and r["eta"] > cma_best[1]):
            cma_best[0] = r["obj"]; cma_best[1] = r["eta"]
            cma_best[2] = xx[0]; cma_best[3] = xx[1]
        return r["obj"]

    try:
        es = cma.CMAEvolutionStrategy([800.0, 0.5], 0.15,
                                       {"maxfevals": 120, "verbose": -9, "seed": 42})
        es.optimize(_cma_obj_cb)
    except Exception as e:
        print(f"    CMA error: {e}")
    t_cma = time.perf_counter() - t_cma
    cma_final = _eval_outer(cma_best[2], cma_best[3], seed=42)
    print(f"  CMA best: h2tout={cma_best[2]:.1f}K s_q={cma_best[3]:.4f} obj={cma_final['obj']:.5f} eta={cma_final['eta']:.4f}  calls={len(cma_history)} evals={sum(h['evals_add'] for h in cma_history) if cma_history else 0} time={t_cma:.1f}s")

    # ═══════════════════════════════════════════════════════
    # 3. DE
    # ═══════════════════════════════════════════════════════
    bounds = [(600.0, 1000.0), (0.05, 0.95)]
    import numpy as np
    import random as _rnd
    rng_de = _rnd.Random(42)
    de_init = np.array([[cma_best[2], cma_best[3]]])
    for _ in range(4):
        de_init = np.vstack([de_init, [rng_de.uniform(*b) for b in bounds]])
    popsize, maxiter_de = 12, 10
    de_best = [float("inf"), 0.0, 800.0, 0.5]
    de_history: list[dict] = []

    def _de_obj(xx):
        r = _eval_outer(xx[0], xx[1], seed=42)
        de_history.append({"h2t": xx[0], "sq": xx[1], "obj": r["obj"], "eta": r["eta"],
                           "evals_add": r["evals"]})
        if r["obj"] < de_best[0] or (r["obj"] == de_best[0] and r["eta"] > de_best[1]):
            de_best[0] = r["obj"]; de_best[1] = r["eta"]
            de_best[2] = xx[0]; de_best[3] = xx[1]
        return r["obj"]

    print(f"\n  DE pop={popsize} iter={maxiter_de} mut=(0.5,1.5) | seeded: [{cma_best[2]:.0f}K,{cma_best[3]:.3f}]")
    t_de = time.perf_counter()
    _sopt.differential_evolution(
        _de_obj, bounds, seed=42, popsize=popsize, maxiter=maxiter_de,
        strategy="best1bin", mutation=(0.5, 1.5), recombination=0.9,
        init=de_init, polish=False,
    )
    t_de = time.perf_counter() - t_de
    de_final = _eval_outer(de_best[2], de_best[3], seed=42)
    print(f"  DE best: h2tout={de_best[2]:.1f}K s_q={de_best[3]:.4f} obj={de_final['obj']:.5f} eta={de_final['eta']:.4f}  calls={len(de_history)} evals={sum(h['evals_add'] for h in de_history)} time={t_de:.1f}s")

    # ── 汇总 ──
    methods = [
        ("Grid", grid_best["h2t"], grid_best["sq"], grid_best["obj"], grid_best["eta"],
         len(grid_data), sum(d["evals"] for d in grid_data), t_grid),
        ("CMA-ES", cma_best[2], cma_best[3], cma_final["obj"], cma_final["eta"],
         len(cma_history), sum(h["evals_add"] for h in cma_history) if cma_history else 0, t_cma),
        ("DE", de_best[2], de_best[3], de_final["obj"], de_final["eta"],
         len(de_history), sum(h["evals_add"] for h in de_history), t_de),
    ]
    print(f"\n  {'='*80}")
    print(f"  0P1S 2D 对比 (h2_T_out, s_q)")
    print(f"  {'='*80}")
    print(f"  {'Method':<12} {'h2_T_out':>8} {'s_q':>7} {'obj':>10} {'eta':>7} {'calls':>6} {'evals':>9} {'time':>7}")
    for name, h2t, sq_v, obj, eta, fc, ev, t in methods:
        print(f"  {name:<12} {h2t:>6.1f}K {sq_v:>5.4f}    {obj:>9.5f} {eta:>7.4f} {fc:>6} {ev:>9} {t:>6.1f}s")
    print(f"  {'='*80}\n")

    cma_penalty = cma_final["obj"] + cma_final["eta"]
    print(f"  penalty check (CMA): obj+η = {cma_penalty:.5f}")

    import csv
    csv_path = out_dir / "grid_0p1s.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["h2_T_out", "s_q", "obj", "eta", "evals", "flows"])
        for d in grid_data:
            w.writerow([d["h2t"], d["sq"], d["obj"], d["eta"], d["evals"],
                        ",".join(f"{v:.1f}" for v in d.get("flows", []))])
    print(f"  CSV: {csv_path}")

    # 简易可视化
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt; import numpy as np
    h2t_v = sorted(set(d["h2t"] for d in grid_data))
    sq_v = sorted(set(d["sq"] for d in grid_data))
    Z = np.full((len(sq_v), len(h2t_v)), -1.0)
    for d in grid_data:
        Z[sq_v.index(d["sq"]), h2t_v.index(d["h2t"])] = d.get("eta", -1)
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.pcolormesh(h2t_v, sq_v, Z, cmap="RdYlGn", shading="auto", vmin=0, vmax=0.55)
    ax.scatter(grid_best["h2t"], grid_best["sq"], marker="*", color="black", s=200, zorder=5, label="Grid")
    ax.scatter(cma_best[2], cma_best[3], marker="D", color="tab:green", s=70, zorder=5, label="CMA-ES")
    ax.scatter(de_best[2], de_best[3], marker="D", color="tab:blue", s=70, zorder=5, label="DE")
    ax.set_xlabel("h2_T_out [K]"); ax.set_ylabel("s_q"); ax.set_title("0P1S 2D"); ax.legend()
    fig.colorbar(im, ax=ax, label="eta"); fig.tight_layout()
    fig.savefig(out_dir / "0p1s_2d.png", dpi=150); plt.close(fig)
    print(f"  可视化: {out_dir / '0p1s_2d.png'}")


def _draw_2d_opt_v2(
    grid_data, de_hist, lb_hist, cma_hist,
    grid_best, de_best, lb_best, cma_best, out_dir: Path,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    h2t_vals = sorted(set(d["h2t"] for d in grid_data))
    pq_vals = sorted(set(d["pq"] for d in grid_data))
    Z_eta = np.full((len(pq_vals), len(h2t_vals)), -1.0)
    Z_obj = np.full((len(pq_vals), len(h2t_vals)), 1e9)
    for d in grid_data:
        hi = h2t_vals.index(d["h2t"]); pi = pq_vals.index(d["pq"])
        Z_eta[pi, hi] = d.get("eta", -1)
        Z_obj[pi, hi] = d.get("obj", 1e9)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle("2D (h2_T_out, p_q) Optimization v2 | 1P0S S24 inner=pinch outer=eff_pinch",
                 fontsize=11, fontweight="bold")

    # 左: eta 热力图
    ax0 = axes[0]
    im0 = ax0.pcolormesh(h2t_vals, pq_vals, Z_eta, cmap="RdYlGn", shading="auto", vmin=0.0, vmax=0.55)
    ax0.scatter(grid_best["h2t"], grid_best["pq"], marker="*", color="black", s=200,
                edgecolors="white", linewidths=1.5, zorder=5, label="Grid")
    ax0.scatter(de_best[0], de_best[1], marker="D", color="tab:blue", s=70, zorder=5, label="DE")
    ax0.scatter(lb_best[0], lb_best[1], marker="D", color="tab:orange", s=70, zorder=5, label="L-BFGS-B")
    ax0.scatter(cma_best[0], cma_best[1], marker="D", color="tab:green", s=70, zorder=5, label="CMA-ES")
    ax0.set_xlabel("h2_T_out [K]"); ax0.set_ylabel("p_q"); ax0.set_title("eta heatmap")
    ax0.legend(fontsize=7); fig.colorbar(im0, ax=ax0, label="eta")

    # 右: 迭代路径
    ax1 = axes[1]
    if de_hist:
        dh = de_hist
        ax1.plot([d["h2t"] for d in dh], [d["pq"] for d in dh], "o-",
                 color="tab:blue", alpha=0.3, markersize=2, linewidth=0.5, label="DE")
    if lb_hist:
        ax1.plot([d["h2t"] for d in lb_hist], [d["pq"] for d in lb_hist], "s-",
                 color="tab:orange", alpha=0.3, markersize=2, linewidth=0.5, label="L-BFGS-B")
    if cma_hist:
        ax1.plot([d["h2t"] for d in cma_hist], [d["pq"] for d in cma_hist], "^-",
                 color="tab:green", alpha=0.3, markersize=2, linewidth=0.5, label="CMA-ES")
    ax1.scatter(grid_best["h2t"], grid_best["pq"], marker="*", color="black", s=200, zorder=5)
    ax1.set_xlabel("h2_T_out [K]"); ax1.set_ylabel("p_q"); ax1.set_title("Search paths")
    ax1.legend(fontsize=7); ax1.grid(True, alpha=0.2)

    plt.tight_layout()
    out_path = out_dir / "2d_opt_v2.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  可视化: {out_path}")


def test_2d_cma_penalty_sweep() -> None:
    """0P1S 2D CMA-ES penalty_w 扫描: 降低惩罚能否帮CMA找到更低h2_T_out? S16."""
    import cma

    hp_base = dict(_HP)
    hp_base["n_t_q"] = 0; hp_base["n_p_q"] = 0; hp_base["n_s_q"] = 1
    hp_base["t_min_lo"] = 50.0; hp_base["t_min_hi"] = 50.0
    hp_base["t_max_lo"] = 1000.0; hp_base["t_max_hi"] = 1000.0
    hp_base["p_min_lo"] = 2000.0; hp_base["p_min_hi"] = 2000.0
    hp_base["p_max_lo"] = 10000.0; hp_base["p_max_hi"] = 10000.0
    hp_base["lbfgsb_starts"] = 16; hp_base["lbfgsb_maxiter"] = 20; hp_base["lbfgsb_workers"] = 16
    hp_base["lbfgsb_sampler"] = "lhs"; hp_base["inner_method"] = "lbfgsb"
    hp_base["penalty_k"] = 10; hp_base["util_tol"] = 1.0
    hp_base["hx_dT"] = 10.0; hp_base["use_non_ideal"] = False
    hp_base["mf_lo"] = 0.0; hp_base["mf_hi"] = 50.0

    from core import ExternalSourceInput, PropertyRegistry
    from core.system import SystemInput, CycleConfig, convert_sources

    out_dir = _OUT_DIR / "penalty_sweep_0p1s"
    out_dir.mkdir(parents=True, exist_ok=True)

    def _eval_2d(h2t: float, sq: float, hp: dict) -> dict:
        tp_in = ClosedCycleTPInput(
            fluid="He", t_min=50.0, t_max=1000.0,
            p_min=2000.0, p_max=10000.0,
            t_quantiles=(), p_quantiles=(), s_quantiles=(sq,),
        )
        layer = ClosedCycleLayer(tp_in)
        if len(layer.subcycles) == 0:
            return {"obj": 1e9, "eta": 0, "n_sc": 0, "flows": [], "evals": 0}
        hp_inner = dict(hp)
        hp_inner["h2_T_out_lo"] = h2t; hp_inner["h2_T_out_hi"] = h2t
        hp_inner["s_q_vals"] = (sq,); hp_inner["obj_mode"] = "pinch"
        hot = ExternalSourceInput(fluid="Air", mass_flow=hp["air_mf"],
                                  T_in=hp["air_T_in"], P_in=hp["air_P_in"],
                                  T_out=hp["air_T_out"], P_out=hp["air_P_out"])
        cold = ExternalSourceInput(fluid=hp.get("cold_fluid", "Hydrogen"),
                                   mass_flow=hp["h2_mf_lo"],
                                   T_in=hp["h2_T_in"], P_in=hp["h2_P_in"],
                                   T_out=h2t, P_out=hp["h2_P_out"])
        cycle_in = ClosedCycleTPInput(
            fluid="He", t_min=50.0, t_max=1000.0, p_min=2000.0, p_max=10000.0,
            t_quantiles=(), p_quantiles=(), s_quantiles=(sq,),
            subcycle_mass_flow_initial=20.0,
        )
        sys_inp = SystemInput(
            heat_sources=(hot,), cold_sources=(cold,),
            cycles=(CycleConfig(input=cycle_in, use_non_ideal=False,
                                 delta_T_min=20.0, heat_method=None),),
            delta_T_min=20.0, heat_method="system_pinch",
        )
        result, _ = _inner_lbfgsb_fast(tp_in, sys_inp, hp_inner, seed=42)
        hp_outer = dict(hp)
        hp_outer["obj_mode"] = "eff_pinch"
        hp_outer["h2_T_out_lo"] = h2t; hp_outer["h2_T_out_hi"] = h2t
        hp_outer["s_q_vals"] = (sq,)
        outer_obj = 1e9
        try:
            layer2 = ClosedCycleLayer(tp_in)
            layer2.subcycle_mass_flows = list(result.flows)
            layer2.commit_subcycle_mass_flows_to_topology()
            outer_obj = _eval_fast(layer2, hp["h2_mf_lo"], h2t, sys_inp, hp_outer)
        except Exception:
            pass
        props = PropertyRegistry()
        _, lst_colds = convert_sources((), sys_inp.cold_sources, props)
        lst_hots, _ = convert_sources(sys_inp.heat_sources, (), props)
        q_source = sum(abs(float(r.power_rate)) for r in lst_hots if r.power_rate)
        q_cold = sum(abs(float(r.power_rate)) for r in lst_colds if r.power_rate)
        eta = (q_source - q_cold) / q_source if q_source > 1e-12 else 0.0
        return {"obj": outer_obj, "eta": eta, "n_sc": len(layer.subcycles),
                "flows": list(result.flows), "evals": result.n_evals}

    penalty_vals = [10, 50, 100, 500, 1000]
    all_summaries: list[dict] = []

    for pw in penalty_vals:
        hp = dict(hp_base)
        hp["penalty_w"] = float(pw)
        history: list[dict] = []
        best_w = {"obj": float("inf"), "eta": 0.0, "x": [800.0, 0.5]}

        def _cma_obj_pw(xx):
            r = _eval_2d(xx[0], xx[1], hp)
            history.append({"h2t": xx[0], "sq": xx[1], "obj": r["obj"], "eta": r["eta"],
                            "evals": r["evals"]})
            if r["obj"] < best_w["obj"] or (r["obj"] == best_w["obj"] and r["eta"] > best_w["eta"]):
                best_w["obj"] = r["obj"]; best_w["eta"] = r["eta"]
                best_w["x"] = [xx[0], xx[1]]
            return r["obj"]

        print(f"\n  penalty_w={pw}")
        es = cma.CMAEvolutionStrategy([800.0, 0.5], 0.15,
                                       {"maxfevals": 120, "verbose": -9, "seed": 42})
        es.optimize(_cma_obj_pw)
        final = _eval_2d(best_w["x"][0], best_w["x"][1], hp)
        penalty = best_w["obj"] + best_w["eta"]
        all_summaries.append({
            "penalty_w": pw, "h2tout": best_w["x"][0], "sq": best_w["x"][1],
            "obj": best_w["obj"], "eta": best_w["eta"], "penalty": penalty,
            "calls": len(history), "evals": sum(h["evals"] for h in history),
        })
        print(f"    h2tout={best_w['x'][0]:.1f}K s_q={best_w['x'][1]:.4f} "
              f"obj={best_w['obj']:.5f} eta={best_w['eta']:.4f} "
              f"penalty={penalty:.5f} calls={len(history)} evals={sum(h['evals'] for h in history)}")

    print(f"\n{'='*80}")
    print(f"  penalty_w 对比汇总 (0P1S 2D CMA-ES)")
    print(f"{'='*80}")
    print(f"  {'w':>8} {'h2_T':>8} {'s_q':>7} {'obj':>10} {'eta':>7} {'penalty':>10} {'calls':>6} {'evals':>9}")
    print(f"  {'-'*75}")
    for s in all_summaries:
        print(f"  {s['penalty_w']:>8.0f} {s['h2tout']:>7.1f}K {s['sq']:>5.4f} {s['obj']:>10.5f} {s['eta']:>7.4f} "
              f"{s['penalty']:>10.5f} {s['calls']:>6} {s['evals']:>9}")
    print(f"{'='*80}\n")


def test_1p0s_pinch_aligned() -> None:
    """0P1S 内层对比: pinch vs pinch_aligned, h2tout=800K固定, S32 w=16, 9 seeds."""
    hp = dict(_HP)
    hp["n_t_q"] = 0; hp["n_p_q"] = 0; hp["n_s_q"] = 1
    hp["s_q_vals"] = (0.5,); hp["p_q_vals"] = (); hp["t_q_vals"] = ()
    hp["t_min_lo"] = 50.0; hp["t_min_hi"] = 50.0
    hp["t_max_lo"] = 1000.0; hp["t_max_hi"] = 1000.0
    hp["p_min_lo"] = 2000.0; hp["p_min_hi"] = 2000.0
    hp["p_max_lo"] = 10000.0; hp["p_max_hi"] = 10000.0
    hp["lbfgsb_starts"] = 32; hp["lbfgsb_maxiter"] = 20; hp["lbfgsb_workers"] = 16
    hp["lbfgsb_sampler"] = "lhs"; hp["inner_method"] = "lbfgsb"
    hp["penalty_w"] = 1000.0; hp["penalty_k"] = 10; hp["util_tol"] = 1.0
    hp["pinch_reward"] = 1.0
    hp["hx_dT"] = 10.0; hp["use_non_ideal"] = False
    hp["mf_lo"] = 0.0; hp["mf_hi"] = 50.0
    hp["h2_T_out_lo"] = 800.0; hp["h2_T_out_hi"] = 800.0

    from core.system import SystemInput, CycleConfig
    from core import ExternalSourceInput
    layer = _build_fixed_layer(hp)
    n_sc = len(layer.subcycles)
    hot = ExternalSourceInput(fluid="Air", mass_flow=hp["air_mf"],
                              T_in=hp["air_T_in"], P_in=hp["air_P_in"],
                              T_out=hp["air_T_out"], P_out=hp["air_P_out"])
    cold = ExternalSourceInput(fluid=hp.get("cold_fluid", "Hydrogen"),
                               mass_flow=hp["h2_mf_lo"],
                               T_in=hp["h2_T_in"], P_in=hp["h2_P_in"],
                               T_out=800.0, P_out=hp["h2_P_out"])
    cycle_in = ClosedCycleTPInput(
        fluid="He", t_min=hp["t_min_lo"], t_max=hp["t_max_hi"],
        p_min=hp["p_min_lo"], p_max=hp["p_max_hi"],
        t_quantiles=hp["t_q_vals"], p_quantiles=hp["p_q_vals"],
        s_quantiles=hp["s_q_vals"], subcycle_mass_flow_initial=20.0,
    )
    sys_inp = SystemInput(
        heat_sources=(hot,), cold_sources=(cold,),
        cycles=(CycleConfig(input=cycle_in, use_non_ideal=False,
                             delta_T_min=20.0, heat_method=None),),
        delta_T_min=20.0, heat_method="system_pinch",
    )
    seeds = [0, 100, 200, 300, 400, 500, 600, 700, 800]

    for label, obj_mode in [("pinch", "pinch"), ("pinch_aligned", "pinch_aligned")]:
        hp_run = dict(hp)
        hp_run["obj_mode"] = obj_mode
        all_r: list[dict] = []
        print(f"\n{'='*85}")
        print(f"  0P1S | {label} | h2tout=800K | n_sc={n_sc} | S32 w=16 | 9 seeds")
        print(f"{'='*85}")
        print(f"{'seed':>5} {'obj':>12} {'evals':>8} {'time_s':>7}  {'flows':>20}")
        print("-" * 72)
        for seed in seeds:
            t0 = time.perf_counter()
            result, _ = _inner_lbfgsb_fast(cycle_in, sys_inp, dict(hp_run), seed=seed)
            elapsed = time.perf_counter() - t0
            flows_str = "[" + ",".join(f"{v:.1f}" for v in result.flows) + "]"
            all_r.append({"seed": seed, "obj": result.obj, "n_evals": result.n_evals,
                          "time_s": elapsed, "flows": list(result.flows)})
            print(f"{seed:>5} {result.obj:>12.5f} {result.n_evals:>8} {elapsed:>7.1f}  {flows_str}")
        objs = [r["obj"] for r in all_r]
        best = min(objs); mean = sum(objs) / len(objs); spread = max(objs) - min(objs)
        print("-" * 72)
        print(f"  best={best:.5f}  mean={mean:.5f}  spread={spread:.5f}  evals-tot={sum(r['n_evals'] for r in all_r)}")
        print(f"{'='*85}")

    print(f"\n  ← pinch_aligned 应该具有更低的 spread——min_dT 提供连续梯度")
