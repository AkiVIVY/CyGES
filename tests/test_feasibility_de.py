r"""h2_T_out 外层纳入 + 内层 flows-only 的优化框架。

核心思路: h2_T_out 是系统级变量 (决定冷源功率/η), 不应塞进内层 random-start L-BFGS-B。
  改为由外层直接搜索 h2_T_out, 内层仅优化子循环流量 (dim=n_sc)。

架构:
  外层 (DE/CMA/…) → 搜索 [t_min, t_max, p_max, p_min, t_q, p_q, s_q, h2_T_out]
  │  每候选: h2_T_out 固定 → 内层 L-BFGS-B(flows only, dim=n_sc)
  └─ obj = -η + penalty_w * softplus(util, tol, k) / q_source

运行:
  pytest -s tests/test_feasibility_de.py::test_1p0s_baseline
"""

from __future__ import annotations

import math
import random
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import pytest

from core import (
    ClosedCycleLayer,
    ClosedCycleTPInput,
    CycleConfig,
    ExternalSourceInput,
    PropertyRegistry,
    SystemInput,
)
from core.heat_exchanger import match_heat_exchanger_groups
from tests.test_layered_opt import (
    _eval_fast,
    _inner_lbfgsb_fast,
    _make_h2_source,
)

_TESTS_DIR = Path(__file__).resolve().parent
_OUT_DIR = _TESTS_DIR / "feasibility_de"
_OUT_DIR.mkdir(parents=True, exist_ok=True)

# ══════════════════════════════════════════════════════════════════
# 基准 HP 配置
# ══════════════════════════════════════════════════════════════════

_HP = {
    # ── 拓扑 ──
    "n_t_q": 0, "n_p_q": 1, "n_s_q": 0,
    # ── 拓扑搜索边界 ──
    "t_min_lo": 50.0, "t_min_hi": 500.0,
    "t_max_lo": 800.0, "t_max_hi": 1100.0,
    "p_min_lo": 2000.0, "p_min_hi": 2000.0,       # p_min 固定
    "p_max_lo": 8000.0, "p_max_hi": 15000.0,
    # ── h2_T_out 搜索边界 (外层变量) ──
    "h2_T_out_lo": 500.0, "h2_T_out_hi": 1000.0,
    "h2_mf_lo": 3.5, "h2_mf_hi": 3.5,              # H2 流量固定
    # ── 外部源 ──
    "cold_fluid": "Hydrogen",
    "h2_T_in": 20.0, "h2_P_in": 5000.0, "h2_P_out": 4500.0,
    "air_mf": 100.0, "air_T_in": 1250.0, "air_P_in": 200.0,
    "air_T_out": 500.0, "air_P_out": 180.0,
    # ── 目标函数 ──
    "obj_mode": "pinch",
    "util_tol": 1.0, "penalty_w": 1000.0, "penalty_k": 10,
    "hx_dT": 10.0,
    "use_non_ideal": False, "use_interp_he": False,
    # ── 内层 L-BFGS-B ──
    "inner_method": "lbfgsb",
    "lbfgsb_starts": 32,
    "lbfgsb_maxiter": 20,
    "lbfgsb_workers": 16,
    "lbfgsb_sampler": "lhs",
    "mf_lo": 0.0, "mf_hi": 50.0,
}


# ══════════════════════════════════════════════════════════════════
# 数据结构
# ══════════════════════════════════════════════════════════════════

@dataclass
class TrialResult:
    """外层单候选评估结果。"""
    t_min: float = 50.0
    t_max: float = 1000.0
    p_min: float = 2000.0
    p_max: float = 10000.0
    t_q: tuple[float, ...] = ()
    p_q: tuple[float, ...] = ()
    s_q: tuple[float, ...] = ()
    h2_T_out: float = 800.0
    flows: list[float] = field(default_factory=list)
    obj: float = 1e9
    n_sc: int = 0
    n_evals: int = 0
    time_s: float = 0.0


# ══════════════════════════════════════════════════════════════════
# 外层变量打包/解包
# ══════════════════════════════════════════════════════════════════

def _pack_bounds(hp: dict) -> list[tuple[float, float]]:
    """从 hp 构造外层搜索边界, 顺序: t_min, t_max, t_q…, p_q…, s_q…, p_max, p_min, h2_T_out."""
    b: list[tuple[float, float]] = [
        (hp["t_min_lo"], hp["t_min_hi"]),
        (hp["t_max_lo"], hp["t_max_hi"]),
    ]
    for _ in range(hp.get("n_t_q", 0)):
        b.append((0.0, 1.0))
    for _ in range(hp.get("n_p_q", 0)):
        b.append((0.0, 1.0))
    for _ in range(hp.get("n_s_q", 0)):
        b.append((0.0, 1.0))
    b.append((hp["p_max_lo"], hp["p_max_hi"]))
    b.append((hp["p_min_lo"], hp["p_min_hi"]))
    b.append((hp["h2_T_out_lo"], hp["h2_T_out_hi"]))
    return b


def _unpack_vector(x: list[float], hp: dict) -> tuple:
    """从 DE 向量解包拓扑参数 + h2_T_out。"""
    t_min = round(x[0])
    t_max = round(x[1])
    idx = 2
    n_tq = hp.get("n_t_q", 0)
    n_pq = hp.get("n_p_q", 0)
    n_sq = hp.get("n_s_q", 0)
    t_q = tuple(float(x[idx + j]) for j in range(n_tq))
    idx += n_tq
    p_q = tuple(float(x[idx + j]) for j in range(n_pq))
    idx += n_pq
    s_q = tuple(float(x[idx + j]) for j in range(n_sq))
    idx += n_sq
    p_max = round(x[idx])
    p_min = round(x[idx + 1])
    h2_T_out = x[idx + 2]
    return t_min, t_max, t_q, p_q, s_q, p_max, p_min, h2_T_out


def _vector_brief(x: list[float], hp: dict) -> str:
    t_min, t_max, t_q, p_q, s_q, p_max, p_min, h2_T_out = _unpack_vector(x, hp)
    parts = [f"T=[{t_min},{t_max}] P=[{p_min},{p_max}]"]
    if t_q: parts.append(f"tq={t_q}")
    if p_q: parts.append(f"pq={p_q}")
    if s_q: parts.append(f"sq={s_q}")
    parts.append(f"h2tout={h2_T_out:.0f}K")
    return " ".join(parts)


# ══════════════════════════════════════════════════════════════════
# 单候选评估
# ══════════════════════════════════════════════════════════════════

def _eval_candidate(x: list[float], hp: dict, seed: int = 42) -> TrialResult:
    """评估一个外层候选向量: 1) 构建拓扑 2) 内层 L-BFGS-B(flows only) 3) 返回 obj。"""
    t0 = time.perf_counter()
    t_min, t_max, t_q, p_q, s_q, p_max, p_min, h2_T_out = _unpack_vector(x, hp)

    # 1. 拓扑
    tp_in = ClosedCycleTPInput(
        fluid="He", t_min=t_min, t_max=t_max,
        p_min=p_min, p_max=p_max,
        t_quantiles=t_q, p_quantiles=p_q, s_quantiles=s_q,
    )
    try:
        probe = ClosedCycleLayer(tp_in)
    except Exception:
        return TrialResult(
            t_min=t_min, t_max=t_max, p_min=p_min, p_max=p_max,
            t_q=t_q, p_q=p_q, s_q=s_q, h2_T_out=h2_T_out,
            obj=1e9, n_sc=0, time_s=time.perf_counter() - t0,
        )
    n_sc = len(probe.subcycles)
    if n_sc == 0:
        return TrialResult(
            t_min=t_min, t_max=t_max, p_min=p_min, p_max=p_max,
            t_q=t_q, p_q=p_q, s_q=s_q, h2_T_out=h2_T_out,
            obj=1e9, n_sc=0, time_s=time.perf_counter() - t0,
        )

    # 2. 构造系统输入 (h2_T_out 已固定)
    sys_inp = _make_sys_inp_from_params(hp, h2_T_out)

    # 3. 内层 L-BFGS-B (h2_T_out 固定 → dim=n_sc)
    hp_inner = dict(hp)
    hp_inner["h2_T_out_lo"] = h2_T_out
    hp_inner["h2_T_out_hi"] = h2_T_out       # → h2_T_fixed=True
    hp_inner["inner_method"] = "lbfgsb"

    result, _ = _inner_lbfgsb_fast(tp_in, sys_inp, hp_inner, seed=seed)

    return TrialResult(
        t_min=t_min, t_max=t_max, p_min=p_min, p_max=p_max,
        t_q=t_q, p_q=p_q, s_q=s_q, h2_T_out=h2_T_out,
        flows=list(result.flows),
        obj=result.obj, n_sc=n_sc,
        n_evals=result.n_evals, time_s=time.perf_counter() - t0,
    )


def _make_sys_inp_from_params(hp: dict, h2_T_out: float) -> SystemInput:
    """用可变 h2_T_out 构造 SystemInput。"""
    hot = ExternalSourceInput(fluid="Air", mass_flow=hp["air_mf"],
                              T_in=hp["air_T_in"], P_in=hp["air_P_in"],
                              T_out=hp["air_T_out"], P_out=hp["air_P_out"])
    cold = ExternalSourceInput(fluid=hp.get("cold_fluid", "Hydrogen"),
                               mass_flow=hp["h2_mf_lo"],
                               T_in=hp["h2_T_in"], P_in=hp["h2_P_in"],
                               T_out=h2_T_out, P_out=hp["h2_P_out"])
    cycle_input = ClosedCycleTPInput(
        fluid="He", t_min=hp["t_min_lo"], t_max=hp["t_max_hi"],
        p_min=hp["p_min_lo"], p_max=hp["p_max_hi"],
        t_quantiles=(), p_quantiles=(),
        subcycle_mass_flow_initial=20.0,
    )
    return SystemInput(
        heat_sources=(hot,), cold_sources=(cold,),
        cycles=(CycleConfig(input=cycle_input, use_non_ideal=hp["use_non_ideal"],
                             delta_T_min=20.0, heat_method=None),),
        delta_T_min=20.0, heat_method="system_pinch",
    )


# ══════════════════════════════════════════════════════════════════
# DE 外层优化
# ══════════════════════════════════════════════════════════════════

def _run_de(
    hp: dict,
    popsize: int = 20,
    maxiter: int = 80,
    seed: int = 42,
    workers: int = 1,
) -> tuple[TrialResult, list[float], list[float]]:
    """运行 scipy DE 外层优化, 返回 best_result + 收敛曲线(每代最优+均值)。"""
    import scipy.optimize as _sopt

    bounds = _pack_bounds(hp)
    dim = len(bounds)
    convergence_best: list[float] = []
    convergence_mean: list[float] = []

    # DE 每次 eval 需要 seed, 每代递增
    class _EvalCounter:
        def __init__(self):
            self.count = 0
            self.generation_pop: list[float] = []

        def callback(self, xk, convergence):
            pass   # 每代结束后由 workers 参数控制

    best_obj = [float("inf")]
    best_result = [None]

    _eval_counter = _EvalCounter()

    # 用 workers=-1 让 scipy 多进程执行为 n_workers 的并行 task
    # 但 _eval_candidate 内部已经并行, 这里外层 worker=1 避免嵌套并行冲突
    # 也可以用 multiprocessing.pool 统一管理

    # 替代: 直接跑 DE 自定义

    def _objective_wrapped(xx: list[float]) -> float:
        _eval_counter.count += 1
        trial = _eval_candidate(xx, hp, seed=seed + _eval_counter.count)
        if trial.obj < best_obj[0]:
            best_obj[0] = trial.obj
            best_result[0] = trial
        return trial.obj

    # scipy DE 内置收敛追踪
    def _per_gen(intermediate_result: _sopt.OptimizeResult):
        gen_pop = [intermediate_result.fun]
        if hasattr(intermediate_result, "population_energies"):
            gen_pop = list(intermediate_result.population_energies)
        convergence_best.append(min(gen_pop))
        convergence_mean.append(sum(gen_pop) / len(gen_pop) if gen_pop else float("inf"))

    result = _sopt.differential_evolution(
        _objective_wrapped, bounds, seed=seed,
        popsize=popsize, maxiter=maxiter,
        strategy="best1bin", mutation=(0.5, 1.0), recombination=0.9,
        workers=workers, updating="deferred",
    )

    best_trial = best_result[0]
    if best_trial is None:
        best_trial = _eval_candidate(list(result.x), hp, seed=seed)

    return best_trial, convergence_best, convergence_mean


# ══════════════════════════════════════════════════════════════════
# 可视化
# ══════════════════════════════════════════════════════════════════

def _draw_de_convergence(
    conv_best: list[float], conv_mean: list[float],
    tag: str, out_path: Path,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5))
    gens = list(range(1, len(conv_best) + 1))
    ax.plot(gens, conv_best, "b-", label="best")
    if conv_mean:
        ax.plot(gens, conv_mean, "orange", alpha=0.5, label="mean")
    ax.set_xlabel("Generation"); ax.set_ylabel("obj")
    ax.set_title(f"DE Convergence | {tag}")
    ax.legend(); ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  收敛图: {out_path}")


# ══════════════════════════════════════════════════════════════════
# 总结输出
# ══════════════════════════════════════════════════════════════════

def _print_result_summary(result: TrialResult, hp: dict, tag: str) -> None:
    print(f"\n{'='*70}")
    print(f"  {tag}")
    print(f"{'='*70}")
    print(f"  T        = [{result.t_min:.0f}, {result.t_max:.0f}] K")
    print(f"  P        = [{result.p_min:.0f}, {result.p_max:.0f}] kPa")
    if result.t_q:
        print(f"  t_q      = {result.t_q}")
    if result.p_q:
        print(f"  p_q      = {result.p_q}")
    if result.s_q:
        print(f"  s_q      = {result.s_q}")
    print(f"  h2_T_out = {result.h2_T_out:.0f} K")
    print(f"  n_sc     = {result.n_sc}")
    print(f"  flows    = {result.flows}")
    print(f"  sum flow = {sum(result.flows):.1f} kg/s")
    print(f"  obj      = {result.obj:.5f}")
    print(f"  n_evals  = {result.n_evals}")
    print(f"  time     = {result.time_s:.1f}s")
    print(f"{'='*70}\n")


# ══════════════════════════════════════════════════════════════════
# 测试入口
# ══════════════════════════════════════════════════════════════════

def test_1p0s_baseline() -> None:
    """1P0S 基准测试: DE 外层 (含 h2_T_out) + 内层 L-BFGS-B (flows only)。"""
    hp = dict(_HP)
    hp["n_t_q"] = 0
    hp["n_p_q"] = 1
    hp["n_s_q"] = 0

    bounds = _pack_bounds(hp)
    dim = len(bounds)
    print(f"\n{'='*70}")
    print(f"  1P0S DE 外层 (dim={dim}) + 内层 L-BFGS-B S{hp['lbfgsb_starts']} w={hp['lbfgsb_workers']}")
    print(f"  obj_mode={hp['obj_mode']}  penalty_w={hp['penalty_w']}  k={hp['penalty_k']}")
    print(f"  h2_T_out bounds: [{hp['h2_T_out_lo']}, {hp['h2_T_out_hi']}]K")
    print(f"{'='*70}")

    result, conv_best, conv_mean = _run_de(
        hp, popsize=15, maxiter=60, seed=42, workers=1,
    )

    _print_result_summary(result, hp, "1P0S DE Outer + L-BFGS-B Inner")

    tag = f"1P0S_DE15x60_S{hp['lbfgsb_starts']}"
    out_csv = _OUT_DIR / f"{tag}.csv"
    import csv
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_min", "t_max", "p_min", "p_max", "t_q", "p_q", "s_q",
                     "h2_T_out", "n_sc", "obj", "n_evals", "time_s", "flows"])
        w.writerow([
            result.t_min, result.t_max, result.p_min, result.p_max,
            " ".join(f"{q:.3f}" for q in result.t_q),
            " ".join(f"{q:.3f}" for q in result.p_q),
            " ".join(f"{q:.3f}" for q in result.s_q),
            result.h2_T_out, result.n_sc, result.obj,
            result.n_evals, result.time_s,
            " ".join(f"{f:.1f}" for f in result.flows),
        ])
    print(f"  CSV: {out_csv}")

    if conv_best:
        _draw_de_convergence(conv_best, conv_mean, tag, _OUT_DIR / f"{tag}.png")


def test_0p1s_baseline() -> None:
    """0P1S 基准测试: DE 外层 (含 h2_T_out) + 内层 L-BFGS-B (flows only)。"""
    hp = dict(_HP)
    hp["n_t_q"] = 0
    hp["n_p_q"] = 0
    hp["n_s_q"] = 1

    bounds = _pack_bounds(hp)
    dim = len(bounds)
    print(f"\n{'='*70}")
    print(f"  0P1S DE 外层 (dim={dim}) + 内层 L-BFGS-B S{hp['lbfgsb_starts']} w={hp['lbfgsb_workers']}")
    print(f"  obj_mode={hp['obj_mode']}  penalty_w={hp['penalty_w']}  k={hp['penalty_k']}")
    print(f"  h2_T_out bounds: [{hp['h2_T_out_lo']}, {hp['h2_T_out_hi']}]K")
    print(f"{'='*70}")

    result, conv_best, conv_mean = _run_de(
        hp, popsize=15, maxiter=60, seed=42, workers=1,
    )

    _print_result_summary(result, hp, "0P1S DE Outer + L-BFGS-B Inner")

    tag = f"0P1S_DE15x60_S{hp['lbfgsb_starts']}"
    out_csv = _OUT_DIR / f"{tag}.csv"
    import csv
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_min", "t_max", "p_min", "p_max", "t_q", "p_q", "s_q",
                     "h2_T_out", "n_sc", "obj", "n_evals", "time_s", "flows"])
        w.writerow([
            result.t_min, result.t_max, result.p_min, result.p_max,
            " ".join(f"{q:.3f}" for q in result.t_q),
            " ".join(f"{q:.3f}" for q in result.p_q),
            " ".join(f"{q:.3f}" for q in result.s_q),
            result.h2_T_out, result.n_sc, result.obj,
            result.n_evals, result.time_s,
            " ".join(f"{f:.1f}" for f in result.flows),
        ])
    print(f"  CSV: {out_csv}")

    if conv_best:
        _draw_de_convergence(conv_best, conv_mean, tag, _OUT_DIR / f"{tag}.png")


def test_h2tout_scan_sweep() -> None:
    """h2_T_out 扫描: 给定拓扑下固定 h2_T_out 扫描, 验证 obj(h2_T_out) 地貌。"""
    hp = dict(_HP)
    hp["n_t_q"] = 0
    hp["n_p_q"] = 1
    hp["n_s_q"] = 0

    # 固定拓扑
    tp_in = ClosedCycleTPInput(
        fluid="He", t_min=50.0, t_max=1000.0,
        p_min=2000.0, p_max=10000.0,
        t_quantiles=(), p_quantiles=(0.5,), s_quantiles=(),
    )
    probe = ClosedCycleLayer(tp_in)
    n_sc = len(probe.subcycles)

    print(f"\n{'='*70}")
    print(f"  h2_T_out 单拓扑扫描 [500,1000]K step=25K")
    print(f"  n_sc={n_sc}  S{hp['lbfgsb_starts']} w={hp['lbfgsb_workers']}")
    print(f"{'='*70}")
    print(f"{'h2_T_out':>9} {'obj':>10} {'eta':>8} {'evals':>8}  {'flows':>20}")
    print("-" * 65)

    all_pts: list[dict] = []
    for h2t in range(500, 1025, 25):
        sys_inp = _make_sys_inp_from_params(hp, float(h2t))
        hp_i = dict(hp)
        hp_i["h2_T_out_lo"] = float(h2t)
        hp_i["h2_T_out_hi"] = float(h2t)
        hp_i["inner_method"] = "lbfgsb"

        result, _ = _inner_lbfgsb_fast(tp_in, sys_inp, hp_i, seed=42)
        eta = -result.obj
        flows_str = "[" + ",".join(f"{v:.1f}" for v in result.flows) + "]"
        all_pts.append({"h2t": h2t, "obj": result.obj, "eta": eta,
                        "evals": result.n_evals, "flows": result.flows})
        print(f"  {h2t:>5}K   {result.obj:>10.5f} {eta:>8.4f} "
              f"{result.n_evals:>8}  {flows_str}")

    best_pt = min(all_pts, key=lambda p: p["obj"])
    print("-" * 65)
    print(f"  best_h2tout={best_pt['h2t']}K  obj={best_pt['obj']:.5f}  eta={best_pt['eta']:.4f}")
    print(f"{'='*70}\n")

    import csv
    csv_path = _OUT_DIR / "h2tout_scan_sweep.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["h2_T_out", "obj", "eta", "evals", "flows"])
        for pt in all_pts:
            w.writerow([pt["h2t"], pt["obj"], pt["eta"], pt["evals"],
                        " ".join(f"{v:.1f}" for v in pt["flows"])])
    print(f"  CSV: {csv_path}")
