r"""CMA-ES 5D 外层优化: [t_min, t_max, p_q, p_max, h2_T_out] + 内层 pinch→eff_pinch。

运行: pytest -s tests/test_cma_5d.py
"""

import sys, time, csv
from pathlib import Path

from core import ClosedCycleLayer, ClosedCycleTPInput, ExternalSourceInput, PropertyRegistry
from core.system import SystemInput, CycleConfig, convert_sources
from tests.test_layered_opt import _eval_fast, _inner_lbfgsb_fast

_TESTS_DIR = Path(__file__).resolve().parent
_OUT_DIR = _TESTS_DIR / "inner_compare" / "cma_5d"
_OUT_DIR.mkdir(parents=True, exist_ok=True)

_HP = {
    "n_t_q": 0, "n_p_q": 1, "n_s_q": 0,
    "h2_mf_lo": 3.5, "h2_mf_hi": 3.5,
    "cold_fluid": "Hydrogen",
    "h2_T_in": 20.0, "h2_P_in": 5000.0, "h2_P_out": 4500.0,
    "air_mf": 100.0, "air_T_in": 1250.0, "air_P_in": 200.0,
    "air_T_out": 500.0, "air_P_out": 180.0,
    "obj_mode": "eff_pinch",
    "util_tol": 1.0, "penalty_w": 1000.0, "penalty_k": 10,
    "hx_dT": 10.0, "use_non_ideal": False, "use_interp_he": False,
    "lbfgsb_starts": 32,
    "lbfgsb_maxiter": 80,
    "lbfgsb_workers": 16,
    "lbfgsb_sampler": "lhs",
    "inner_method": "lbfgsb",
    "mf_lo": 0.0, "mf_hi": 50.0,
}


def test_cma_5d() -> None:
    import cma

    hp = dict(_HP)

    # 5D bounds: t_min, t_max, p_q, p_max, h2_T_out
    bounds_5d = [
        (50.0, 500.0),
        (800.0, 1100.0),
        (0.05, 0.95),
        (8000.0, 15000.0),
        (600.0, 1000.0),
    ]

    def _eval_5d(xx: list[float]) -> dict:
        t0 = time.perf_counter()
        t_min = round(xx[0]); t_max = round(xx[1]); pq = xx[2]
        p_max = round(xx[3]); h2t = xx[4]; p_min = 2000.0

        tp_in = ClosedCycleTPInput(
            fluid="He", t_min=t_min, t_max=t_max,
            p_min=p_min, p_max=p_max,
            t_quantiles=(), p_quantiles=(pq,), s_quantiles=(),
        )
        try:
            layer = ClosedCycleLayer(tp_in)
        except Exception:
            return {"obj": 1e9, "eta": 0, "n_sc": 0, "flows": [], "evals": 0,
                    "time_s": time.perf_counter() - t0}
        n_sc = len(layer.subcycles)
        if n_sc == 0:
            return {"obj": 1e9, "eta": 0, "n_sc": 0, "flows": [], "evals": 0,
                    "time_s": time.perf_counter() - t0}

        # 内层 pinch: 收 flows
        hp_inner = dict(hp)
        hp_inner["h2_T_out_lo"] = h2t; hp_inner["h2_T_out_hi"] = h2t
        hp_inner["p_q_vals"] = (pq,); hp_inner["obj_mode"] = "pinch"
        hp_inner["inner_method"] = "lbfgsb"

        hot = ExternalSourceInput(fluid="Air", mass_flow=hp["air_mf"],
                                  T_in=hp["air_T_in"], P_in=hp["air_P_in"],
                                  T_out=hp["air_T_out"], P_out=hp["air_P_out"])
        cold = ExternalSourceInput(fluid=hp.get("cold_fluid", "Hydrogen"),
                                   mass_flow=hp["h2_mf_lo"],
                                   T_in=hp["h2_T_in"], P_in=hp["h2_P_in"],
                                   T_out=h2t, P_out=hp["h2_P_out"])
        cycle_in = ClosedCycleTPInput(
            fluid="He", t_min=t_min, t_max=t_max, p_min=p_min, p_max=p_max,
            t_quantiles=(), p_quantiles=(pq,), s_quantiles=(),
            subcycle_mass_flow_initial=20.0,
        )
        sys_inp = SystemInput(
            heat_sources=(hot,), cold_sources=(cold,),
            cycles=(CycleConfig(input=cycle_in, use_non_ideal=False,
                                 delta_T_min=20.0, heat_method=None),),
            delta_T_min=20.0, heat_method="system_pinch",
        )
        result, _ = _inner_lbfgsb_fast(tp_in, sys_inp, hp_inner, seed=42)

        # 外层 eff_pinch: 重算 obj
        hp_outer = dict(hp)
        hp_outer["obj_mode"] = "eff_pinch"
        hp_outer["h2_T_out_lo"] = h2t; hp_outer["h2_T_out_hi"] = h2t
        hp_outer["p_q_vals"] = (pq,)
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

        return {"obj": outer_obj, "eta": eta, "n_sc": n_sc,
                "flows": list(result.flows), "evals": result.n_evals,
                "time_s": time.perf_counter() - t0}

    history: list[dict] = []
    best = {"obj": float("inf"), "eta": 0.0, "x": [50, 1000, 0.5, 10000, 800]}

    def _cma_obj(xx):
        r = _eval_5d(list(xx))
        history.append({"x": list(xx), "obj": r["obj"], "eta": r["eta"],
                        "n_sc": r["n_sc"], "evals": r["evals"], "time_s": r["time_s"]})
        if r["obj"] < best["obj"] or (best["obj"] == r["obj"] and r["eta"] > best["eta"]):
            best["obj"] = r["obj"]; best["eta"] = r["eta"]; best["x"] = list(xx)
            sys.stdout.write(
                f"\r    [{len(history):>4}] obj={r['obj']:.5f} eta={r['eta']:.4f} "
                f"T=[{round(xx[0])},{round(xx[1])}] pq={xx[2]:.3f} "
                f"pmax={round(xx[3])} h2={xx[4]:.0f}K nsc={r['n_sc']}  ")
            sys.stdout.flush()
        return r["obj"]

    print(f"\n{'='*90}")
    print(f"  CMA-ES 5D (maxfevals=300)")
    print(f"  t_min∈[50,500] t_max∈[800,1100] p_q∈[0.05,0.95] p_max∈[8000,15000] h2tout∈[600,1000]")
    print(f"  inner: S32 w=16 pinch → outer: eff_pinch")
    print(f"{'='*90}")
    t_total = time.perf_counter()

    es = None
    xopt = None
    try:
        xopt, es = cma.fmin(
            _cma_obj, [200.0, 950.0, 0.5, 11500.0, 800.0], 0.25,
            options={
                "maxfevals": 300, "verbose": -9, "seed": 42,
                "bounds": [[50, 800, 0.05, 8000, 600],
                           [500, 1100, 0.95, 15000, 1000]],
            },
        )
    except Exception as e:
        print(f"\n    CMA fmin exception: {e}")

    t_total = time.perf_counter() - t_total
    sys.stdout.write("\n")

    best_r = _eval_5d(best["x"])
    penalty = best["obj"] + best["eta"]

    print(f"\n{'='*80}")
    print(f"  CMA-ES 5D 结果")
    print(f"{'='*80}")
    print(f"  t_min     = {best['x'][0]:.0f} K")
    print(f"  t_max     = {best['x'][1]:.0f} K")
    print(f"  p_q       = {best['x'][2]:.4f}")
    print(f"  p_max     = {best['x'][3]:.0f} kPa")
    print(f"  h2_T_out  = {best['x'][4]:.0f} K")
    print(f"  n_sc      = {best_r['n_sc']}")
    print(f"  flows     = {[f'{v:.1f}' for v in best_r['flows']]}")
    print(f"  sum flow  = {sum(best_r['flows']):.1f} kg/s")
    print(f"  obj       = {best['obj']:.5f}")
    print(f"  eta       = {best['eta']:.4f}")
    print(f"  penalty   = {penalty:.6f}")
    print(f"  funcalls  = {len(history)}")
    print(f"  evals     = {sum(h['evals'] for h in history)}")
    print(f"  time      = {t_total:.1f}s")
    print(f"{'='*80}\n")

    csv_path = _OUT_DIR / "cma_5d_history.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["n", "obj", "eta", "t_min", "t_max", "p_q", "p_max", "h2tout",
                     "n_sc", "evals", "time_s"])
        for i, h in enumerate(history):
            x = h["x"]
            w.writerow([i + 1, h["obj"], h["eta"], round(x[0]), round(x[1]),
                        x[2], round(x[3]), round(x[4]), h["n_sc"], h["evals"], h["time_s"]])
    print(f"  CSV: {csv_path}")

    if history:
        _draw_cma_5d_convergence(history, best, _OUT_DIR)


def _draw_cma_5d_convergence(history: list, best: dict, out_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    objs = [h["obj"] for h in history]
    etas = [h["eta"] for h in history]
    xs = list(range(1, len(objs) + 1))
    feasible = [i for i, e in enumerate(etas) if e > 0]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(f"CMA-ES 5D | maxfevals=300 | best eta={best['eta']:.4f} obj={best['obj']:.5f}",
                 fontsize=11, fontweight="bold")

    ax0 = axes[0]
    ax0.plot(xs, objs, ".-", markersize=3, alpha=0.7, color="tab:blue")
    ax0.axhline(y=best["obj"], color="red", linestyle="--", alpha=0.5, label=f"best={best['obj']:.5f}")
    ax0.set_xlabel("Function evaluations"); ax0.set_ylabel("obj")
    ax0.set_title("obj convergence"); ax0.legend(); ax0.grid(True, alpha=0.2)

    ax1 = axes[1]
    ax1.plot(xs, etas, ".-", markersize=3, alpha=0.5, color="tab:green")
    if feasible:
        ax1.scatter([xs[i] for i in feasible], [etas[i] for i in feasible],
                    color="tab:green", s=10, label=f"feasible={len(feasible)}")
    ax1.axhline(y=best["eta"], color="red", linestyle="--", alpha=0.5, label=f"best={best['eta']:.4f}")
    ax1.set_xlabel("Function evaluations"); ax1.set_ylabel("eta")
    ax1.set_title("eta trace"); ax1.legend(); ax1.grid(True, alpha=0.2)

    plt.tight_layout()
    out_path = out_dir / "cma_5d_convergence.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  可视化: {out_path}")
