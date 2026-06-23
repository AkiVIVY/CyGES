"""对比: 内层 workers=1 vs workers=8, 外层从 (0.2,0.2) L-BFGS-B 优化 (p_q,s_q).

运行: pytest -s tests/test_outer_quantile_opt.py::test_compare_workers
"""

import time, sys
from pathlib import Path
from core import ClosedCycleTPInput
from tests.test_layered_opt import _inner_lbfgsb_fast

_TESTS_DIR = Path(__file__).resolve().parent
_OUT_DIR = _TESTS_DIR / "outer_quantile_opt"
_OUT_DIR.mkdir(parents=True, exist_ok=True)

_HP_BASE = {
    "n_t_q": 0, "n_p_q": 1, "n_s_q": 1,
    "t_min_lo": 50.0, "t_min_hi": 50.0,
    "t_max_lo": 1000.0, "t_max_hi": 1000.0,
    "p_min_lo": 2000.0, "p_min_hi": 2000.0,
    "p_max_lo": 10000.0, "p_max_hi": 10000.0,
    "cold_fluid": "Hydrogen", "h2_mf_lo": 3.5, "h2_mf_hi": 3.5,
    "h2_T_in": 20.0, "h2_P_in": 5000.0, "h2_P_out": 4500.0,
    "h2_T_out_lo": 100.0, "h2_T_out_hi": 1000.0,
    "air_mf": 100.0, "air_T_in": 1250.0, "air_P_in": 200.0,
    "air_T_out": 500.0, "air_P_out": 180.0,
    "obj_mode": "eff_pinch", "util_tol": 1.0, "penalty_w": 1000,
    "hx_dT": 10.0, "penalty_k": 10,
    "use_non_ideal": False, "use_interp_he": False,
    "mf_lo": 0.0, "mf_hi": 50.0, "qstep": 0.001,
    "flow_step": 0.05, "h2_step": 0.01,
    "lbfgsb_starts": 16, "lbfgsb_maxiter": 50,
    "inner_method": "lbfgsb",
}


def _outer_obj(x, sys_inp, hp):
    p_q = max(0.01, min(0.99, float(x[0])))
    s_q = max(0.01, min(0.99, float(x[1])))
    tp_in = ClosedCycleTPInput(
        fluid="He", t_min=hp["t_min_lo"], t_max=hp["t_max_hi"],
        p_min=hp["p_min_lo"], p_max=hp["p_max_hi"],
        t_quantiles=(), p_quantiles=(p_q,), s_quantiles=(s_q,),
    )
    result, _ = _inner_lbfgsb_fast(tp_in, sys_inp, hp, seed=42)
    return result.obj


def run_outer(workers, label):
    hp = dict(_HP_BASE)
    hp["lbfgsb_workers"] = workers
    from core.system import SystemInput, CycleConfig
    from core import ExternalSourceInput
    hot = ExternalSourceInput(fluid="Air", mass_flow=hp["air_mf"],
                              T_in=hp["air_T_in"], P_in=hp["air_P_in"],
                              T_out=hp["air_T_out"], P_out=hp["air_P_out"])
    cold = ExternalSourceInput(fluid=hp["cold_fluid"], mass_flow=hp["h2_mf_lo"],
                               T_in=hp["h2_T_in"], P_in=hp["h2_P_in"],
                               T_out=hp["h2_T_out_lo"], P_out=hp["h2_P_out"])
    sys_inp = SystemInput(
        heat_sources=(hot,), cold_sources=(cold,),
        cycles=(CycleConfig(input=ClosedCycleTPInput(
            fluid="He", t_min=hp["t_min_lo"], t_max=hp["t_max_hi"],
            p_min=hp["p_min_lo"], p_max=hp["p_max_hi"],
            t_quantiles=(), p_quantiles=(0.2,), s_quantiles=(0.2,)),
            use_non_ideal=False, delta_T_min=20.0, heat_method=None),),
        delta_T_min=20.0, heat_method="system_pinch",
    )
    import scipy.optimize as _sopt
    x0 = [0.2, 0.2]
    lo = [0.1, 0.1]; hi = [0.9, 0.9]
    history = []

    def _f(x):
        obj = _outer_obj(x, sys_inp, hp)
        history.append((float(x[0]), float(x[1]), obj))
        sys.stdout.write(f"  {label} eval{len(history)}: ({x[0]:.4f},{x[1]:.4f}) obj={obj:.5f}\n")
        sys.stdout.flush()
        return obj

    print(f"\n--- {label} (workers={workers}) ---")
    t0 = time.perf_counter()
    res = _sopt.minimize(_f, x0, method="L-BFGS-B",
                          bounds=[(lo[i], hi[i]) for i in range(2)],
                          options={"maxiter": 15, "ftol": 1e-6, "eps": 0.01})
    elapsed = time.perf_counter() - t0
    print(f"  Done: ({res.x[0]:.4f},{res.x[1]:.4f}) obj={res.fun:.5f} nfev={res.nfev} {elapsed:.1f}s")
    return history, res, elapsed


def test_compare_workers():
    print("=" * 60)
    print("对比: 内层 w=1 vs w=8, 外层从(0.2,0.2) L-BFGS-B")
    print("=" * 60)

    h1, r1, t1 = run_outer(1, "SER")
    h2, r2, t2 = run_outer(8, "PAR")

    print(f"\n{'='*60}")
    print(f"  SER: ({r1.x[0]:.4f},{r1.x[1]:.4f}) obj={r1.fun:.5f}  {len(h1)} evals  {t1:.1f}s")
    print(f"  PAR: ({r2.x[0]:.4f},{r2.x[1]:.4f}) obj={r2.fun:.5f}  {len(h2)} evals  {t2:.1f}s")
    print(f"{'='*60}")

    # 检查路径一致性
    same = True
    for i, (h1p, h2p) in enumerate(zip(h1, h2)):
        if abs(h1p[0]-h2p[0]) > 1e-6 or abs(h1p[1]-h2p[1]) > 1e-6 or abs(h1p[2]-h2p[2]) > 1e-6:
            print(f"  路径分歧 @ eval{i+1}: SER({h1p[0]:.4f},{h1p[1]:.4f})={h1p[2]:.5f}  PAR({h2p[0]:.4f},{h2p[1]:.4f})={h2p[2]:.5f}")
            same = False
    if same:
        print("  路径完全一致 — 并行未引入误差")
