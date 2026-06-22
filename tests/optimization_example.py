"""热效率优化示例 — 1P1S 拓扑，零公用工程约束下最大化循环效率。

外层 DE: t_min ∈ [50,400], t_max ∈ [600,1000], p_max ∈ [2000,12000], p_q ∈ (0,1), s_q ∈ (0,1)
内层 CMA: 子循环流量 + h2_T_out ∈ [100,1000]K
目标:   η = (Q_source − Q_cold) / Q_source，软惩罚 utility

输出: tests/optimization_example/result.csv + result.json + 图表

运行:
    pytest -s tests/optimization_example.py::test_h2_scan_efficiency
"""

import time
from pathlib import Path

from tests.test_layered_opt import _run_layered

_EXAMPLE_DIR = Path(__file__).resolve().parent / "optimization_example"

_HP_BASE = {
    # ── 拓扑 ──
    "n_t_q": 0, "n_p_q": 1, "n_s_q": 1,

    # ── 外层边界 ──
    "t_min_lo": 150.0, "t_min_hi": 450.0,
    "t_max_lo": 600.0, "t_max_hi": 950.0,
    "p_min_lo": 2000.0, "p_min_hi": 2000.0,       # p_min 固定
    "p_max_lo": 2000.0, "p_max_hi": 12000.0,

    # ── 内层冷源 (CH₄) ──
    "cold_fluid": "Methane",
    "h2_T_in": 120.0, "h2_P_in": 5000.0,
    "h2_P_out": 4500.0,
    "h2_T_out_lo": 300.0, "h2_T_out_hi": 900.0,    # 内层 CMA 自由优化

    # ── 热源 Air ──
    "air_mf": 100.0, "air_T_in": 1000.0, "air_P_in": 179.0,
    "air_T_out": 500.0, "air_P_out": 151.0,

    # ── 目标函数 ──
    "obj_mode": "eff_pinch",
    "util_tol": 1.0,
    "penalty_w": 100.0,
    "hx_dT": 10.0,

    # ── 优化器预算 ──
    "outer_method": "de",
    "n_lhs": 50,
    "de_popsize": 20, "de_maxiter": 35, "de_F": 0.8, "de_CR": 0.9,
    "n_workers": 8,
    "maxiter_inner": 45, "restarts_inner": 4, "sigma0": 20.0,

    # ── 其他 ──
    "out_dir": str(_EXAMPLE_DIR),
    "qstep": 0.001,
    "flow_step": 0.05, "h2_step": 0.01,
    "mf_lo": 0.0, "mf_hi": 50.0,
    "use_non_ideal": False,
    "use_interp_he": False,
    "hx_max_group_size": 3,
    "early_stop": 10,
}


def test_h2_scan_efficiency() -> None:
    r"""H2 流量扫描：固定 H2 流量，内层自由优化 h2_T_out，最大化热效率。

    扫描范围: H2 ∈ {2.5, 3.0, 3.5, 4.0} kg/s
    预期: 存在使 η 最大的最优 H2（权衡冷源用量与零公用约束）。
    """
    h2_levels = [11.0]  # CH4=11kg/s
    for h2_val in h2_levels:
        hp = {**_HP_BASE,
              "h2_mf_lo": float(h2_val), "h2_mf_hi": float(h2_val)}
        tag = f"eff1P1S_CH4-{h2_val}"
        print(f"\n{'#'*70}")
        print(f"# η-opt 1P1S H2 = {h2_val} kg/s")
        print(f"{'#'*70}")
        t0 = time.perf_counter()
        _run_layered(hp, tag=tag)
        print(f"  H2={h2_val}: {time.perf_counter()-t0:.1f}s")

    print("\n" + "=" * 80)
    print("η-opt H2 扫描完成 — 见上方各 run 的 η 值")
    print("=" * 80)


if __name__ == "__main__":
    test_h2_scan_efficiency()
