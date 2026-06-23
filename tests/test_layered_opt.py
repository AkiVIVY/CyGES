r"""分层循环拓扑+流量联合优化。

  外层: LHS/DE 搜索拓扑参数 (t_min/t_max/p_min/p_max/t_q/p_q/s_q),
       多进程并行评估各候选拓扑。
  内层: 三种方法可选 — CMA-ES / L-BFGS-B(LHS多起点+best-track+并行) / hybrid两阶段,
       在固定拓扑下优化子循环质量流量 + H₂冷源流量 + H₂出口温度。
  目标: 4 种 obj_mode — group(星型HX) / series(串联夹点HX) /
        pinch(直接最小化公用工程) / eff_pinch(最大化热效率η + 软惩罚)。

  输出 TS/PS 拓扑图、HX T-Q 复合曲线、DE 收敛曲线、最优参数 CSV/JSON。
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
    ProcessCategory,
    PropertyRegistry,
    SystemInput,
)
from core.heat_exchanger import (
    _RecInfo,
    _normalize_records,
    HXMatchResult,
    HXUnit,
    match_heat_exchanger_groups,
)

TESTS_DIR = Path(__file__).resolve().parent

_DEFAULT_HP = {
    "sigma0": 15.0, "early_stop": 5,
    "p_min": 2000.0, "hx_max_group_size": 3,
    "flow_step": 0.05, "h2_step": 0.01,
    "h2_mf_lo": 3.0, "h2_mf_hi": 6.0,
    "mf_lo": 0.0, "mf_hi": 50.0,
    "n_workers": 6,
    "use_interp_he": False,
    "n_s_q": 0,
    # 外层搜索边界
    "t_min_lo": 40.0, "t_min_hi": 500.0,
    "t_max_lo": 800.0, "t_max_hi": 1100.0,
    "p_min_lo": 2000.0, "p_min_hi": 5000.0,
    "p_max_lo": 8000.0, "p_max_hi": 15000.0,
    # 外部源参数
    "air_mf": 100.0, "air_T_in": 1250.0, "air_P_in": 200.0,
    "air_T_out": 500.0, "air_P_out": 180.0,
    "h2_T_in": 20.0, "h2_P_in": 5000.0,
    "h2_T_out": 900.0, "h2_P_out": 4500.0,
    "h2_T_out_lo": 400.0, "h2_T_out_hi": 900.0,
    "use_non_ideal": True,
    # L-BFGS-B 内层参数
    "inner_method": "cma",
    "lbfgsb_starts": 20,
    "lbfgsb_maxiter": 50,
}


def _make_h2_source(hp: dict, h2_mf: float, h2_T_out: float | None = None) -> ExternalSourceInput:
    """用 hp 中的冷源参数 + 可变流量/出口温构造冷源。"""
    return ExternalSourceInput(
        fluid=hp.get("cold_fluid", "Hydrogen"), mass_flow=h2_mf,
        T_in=hp["h2_T_in"], P_in=hp["h2_P_in"],
        T_out=h2_T_out if h2_T_out is not None else hp["h2_T_out"],
        P_out=hp["h2_P_out"],
    )

# ──────────────────────────────────────────────────────────────────
# 系统输入
# ──────────────────────────────────────────────────────────────────


def _make_system_input(
    hp: dict,
    t_quantiles: tuple[float, ...] = (),
    p_quantiles: tuple[float, ...] = (0.33, 0.67),
    h2_mass_flow: float = 4.3,
) -> SystemInput:
    hot = ExternalSourceInput(fluid="Air",
                              mass_flow=hp["air_mf"],
                              T_in=hp["air_T_in"], P_in=hp["air_P_in"],
                              T_out=hp["air_T_out"], P_out=hp["air_P_out"])
    cold = ExternalSourceInput(fluid=hp.get("cold_fluid", "Hydrogen"), mass_flow=h2_mass_flow,
                               T_in=hp["h2_T_in"], P_in=hp["h2_P_in"],
                               T_out=hp["h2_T_out"], P_out=hp["h2_P_out"])
    cycle_input = ClosedCycleTPInput(
        fluid="He", t_min=hp["t_min_lo"], t_max=hp["t_max_hi"],
        p_min=hp["p_min"], p_max=hp["p_max_hi"],
        t_quantiles=t_quantiles, p_quantiles=p_quantiles,
        subcycle_mass_flow_initial=20.0,
    )
    return SystemInput(
        heat_sources=(hot,), cold_sources=(cold,),
        cycles=(CycleConfig(input=cycle_input, use_non_ideal=hp["use_non_ideal"],
                            delta_T_min=20.0, heat_method=None),),
        delta_T_min=20.0, heat_method="system_pinch",
    )


# ──────────────────────────────────────────────────────────────────
# LHS 采样
# ──────────────────────────────────────────────────────────────────


def _lhs(n: int, bounds: list[tuple[float, float]], seed: int = 42) -> list[list[float]]:
    d = len(bounds)
    rng = random.Random(seed)
    samples: list[list[float]] = [[0.0] * d for _ in range(n)]
    for dim in range(d):
        lo, hi = bounds[dim]
        bucket_order = list(range(n))
        rng.shuffle(bucket_order)
        for i in range(n):
            b = bucket_order[i]
            lo_i = lo + (hi - lo) * b / n
            hi_i = lo + (hi - lo) * (b + 1) / n
            samples[i][dim] = rng.uniform(lo_i, hi_i)
    return samples


# ──────────────────────────────────────────────────────────────────
# 数据结构
# ──────────────────────────────────────────────────────────────────


@dataclass
class LayerResult:
    n: int
    t_min: float
    t_max: float
    p_min: float
    p_max: float
    t_q: tuple[float, ...]
    p_q: tuple[float, ...]
    s_q: tuple[float, ...] = ()
    flows: list[float] = field(default_factory=list)
    h2_mf: float = 0.0
    h2_T_out: float = 900.0
    obj: float = 1e9
    n_subcycles: int = 0
    n_evals: int = 0
    runtime: float = 0.0


# ──────────────────────────────────────────────────────────────────
# A. 快速内层 eval: 复用已构建 layer, 只更新流量
# ──────────────────────────────────────────────────────────────────


def _eval_fast(layer: ClosedCycleLayer, h2_mf: float, h2_T_out: float,
               sys_inp_template: SystemInput, hp: dict) -> float:
    """固定拓扑下：更新流量 → 非理想 → 性能 → HX。不重建 ClosedCycleLayer。"""
    from core.system import SystemPipeline
    from core import InterpolatingHeliumSolver

    props = PropertyRegistry()
    if hp["use_non_ideal"]:
        try:
            ni = layer.ensure_non_ideal()
            ni.apply_offsets()
        except Exception:
            return 1e9

    cold_src = _make_h2_source(hp, h2_mf, h2_T_out)

    report = layer.performance_report()
    hots: list = []
    colds: list = []
    for _, rec in report.by_edge:
        if rec.kind == "heat" and rec.power_rate is not None:
            if rec.category == ProcessCategory.HEAT_REJECTION:
                hots.append(rec)
            else:
                colds.append(rec)

    from core.system import convert_sources
    src_hots, _ = convert_sources(sys_inp_template.heat_sources, (), props)
    _, src_colds = convert_sources((), (cold_src,), props)
    hots[:0] = src_hots
    colds[:0] = src_colds

    q_source = sum(abs(float(r.power_rate)) for r in src_hots if r.power_rate)
    if q_source < 1e-12:
        return 1.0

    obj_mode = hp.get("obj_mode", "series" if hp.get("hx_series") else "group")
    if obj_mode == "pinch":
        from core.postprocess import analyze_pinch
        pinch = analyze_pinch(colds, hots, hp["hx_dT"], props)
        obj = (pinch.hot_utility_demand + pinch.cold_utility_demand) / q_source
        return obj
    elif obj_mode == "eff_pinch":
        from core.postprocess import analyze_pinch
        import math
        pinch = analyze_pinch(colds, hots, hp["hx_dT"], props)
        util = pinch.hot_utility_demand + pinch.cold_utility_demand
        q_cold = sum(abs(float(r.power_rate)) for r in src_colds if r.power_rate)
        eta = (q_source - q_cold) / q_source if q_source > 1e-12 else 0.0
        penalty_w = hp.get("penalty_w", 100.0)
        util_tol = hp.get("util_tol", 1.0)
        k = hp.get("penalty_k", None)
        if k is not None:
            excess = k * (util - util_tol)
            if excess > 50:
                penalty = util - util_tol
            elif excess < -10:
                penalty = 0.0
            else:
                penalty = math.log(1 + math.exp(excess)) / k
        else:
            penalty = max(0, util - util_tol)
        obj = -eta + penalty_w * penalty / q_source
        return obj
    elif obj_mode == "eff_pinch_split":
        from core.postprocess import analyze_pinch
        import math
        pinch = analyze_pinch(colds, hots, hp["hx_dT"], props)
        q_cold = sum(abs(float(r.power_rate)) for r in src_colds if r.power_rate)
        eta = (q_source - q_cold) / q_source if q_source > 1e-12 else 0.0
        penalty_w = hp.get("penalty_w", 100.0)
        util_tol = hp.get("util_tol", 1.0)
        k = hp.get("penalty_k", 10)
        def _sp(excess: float) -> float:
            if excess > 50: return excess / k
            if excess < -10: return 0.0
            return math.log(1 + math.exp(excess)) / k
        p_hot = _sp(k * (pinch.hot_utility_demand - util_tol))
        p_cold = _sp(k * (pinch.cold_utility_demand - util_tol))
        obj = -eta + penalty_w * (p_hot + p_cold) / q_source
        return obj
    elif obj_mode == "series":
        from core.heat_exchanger import match_series_pinch
        hx = match_series_pinch(hots, colds, dT_min=hp["hx_dT"])
        unmatched_ratio = hx.total_unmatched / q_source
        obj = unmatched_ratio + 1e-2 * hx.pinch_violation
        return obj
    else:
        hx = match_heat_exchanger_groups(hots, colds, dT_min=hp["hx_dT"],
                                          max_group_size=hp["hx_max_group_size"])
        unmatched_ratio = hx.total_unmatched / q_source
        num_unmatched = len(hx.unassigned_hots) + len(hx.unassigned_colds)
        obj = unmatched_ratio + 1e-3 * num_unmatched
        return obj


# ──────────────────────────────────────────────────────────────────
# B. 内层 CMA: 复用拓扑
# ──────────────────────────────────────────────────────────────────


def _inner_cma_fast(
    tp_in: ClosedCycleTPInput,
    sys_inp: SystemInput,
    hp: dict,
    *,
    seed: int = 42,
    layer: ClosedCycleLayer | None = None,
) -> tuple[LayerResult, ClosedCycleLayer]:
    """内层 CMA 优化子循环流量 + H2 流量；可选复用预建 layer。"""
    import cma
    from core import InterpolatingHeliumSolver

    t0 = time.perf_counter()
    use_he = hp.get("use_interp_he", True)
    if use_he:
        from core import InterpolatingHeliumSolver
        he_solver = InterpolatingHeliumSolver(
            T_min=tp_in.t_min, T_max=tp_in.t_max,
            P_min=tp_in.p_min, P_max=tp_in.p_max,
        )
        he_solver._build()

    if layer is not None:
        base_layer = layer
    elif use_he:
        try:
            base_layer = ClosedCycleLayer(tp_in, properties=he_solver)
            if len(base_layer.subcycles) == 0:
                base_layer = ClosedCycleLayer(tp_in)
            else:
                ref = ClosedCycleLayer(tp_in)
                if len(ref.subcycles) != len(base_layer.subcycles):
                    base_layer = ref
        except Exception:
            try:
                base_layer = ClosedCycleLayer(tp_in)
            except Exception:
                runtime = time.perf_counter() - t0
                return LayerResult(
                    n=0, t_min=tp_in.t_min, t_max=tp_in.t_max,
                    p_min=tp_in.p_min, p_max=tp_in.p_max,
                    t_q=tp_in.t_quantiles, p_q=tp_in.p_quantiles,
                    flows=[], h2_mf=hp["h2_mf_lo"], h2_T_out=hp["h2_T_out_lo"],
                    obj=1.0, n_subcycles=0, n_evals=0, runtime=runtime,
                ), None
    else:
        try:
            base_layer = ClosedCycleLayer(tp_in)
        except Exception:
            runtime = time.perf_counter() - t0
            return LayerResult(
                n=0, t_min=tp_in.t_min, t_max=tp_in.t_max,
                p_min=tp_in.p_min, p_max=tp_in.p_max,
                t_q=tp_in.t_quantiles, p_q=tp_in.p_quantiles,
                flows=[], h2_mf=hp["h2_mf_lo"], h2_T_out=hp["h2_T_out_lo"],
                obj=1.0, n_subcycles=0, n_evals=0, runtime=runtime,
            ), None
    n_sc = len(base_layer.subcycles)
    if n_sc == 0:
        runtime = time.perf_counter() - t0
        return LayerResult(
            n=0, t_min=tp_in.t_min, t_max=tp_in.t_max,
            p_min=tp_in.p_min, p_max=tp_in.p_max,
            t_q=tp_in.t_quantiles, p_q=tp_in.p_quantiles,
            flows=[], h2_mf=hp["h2_mf_lo"], h2_T_out=hp["h2_T_out_lo"],
            obj=1.0, n_subcycles=0, n_evals=0, runtime=runtime,
        ), base_layer
    h2_fixed = hp["h2_mf_lo"] == hp["h2_mf_hi"]
    h2_T_fixed = hp["h2_T_out_lo"] == hp["h2_T_out_hi"]

    def _extract_extra(sol: list[float], idx: int) -> tuple[float, float]:
        ei = idx
        h2_mf = hp["h2_mf_lo"] if h2_fixed else sol[ei]; ei += 0 if h2_fixed else 1
        h2_T_out = hp["h2_T_out_lo"] if h2_T_fixed else sol[ei]
        return h2_mf, h2_T_out

    dim = n_sc
    lo = [hp["mf_lo"]] * n_sc
    hi = [hp["mf_hi"]] * n_sc
    if not h2_fixed:
        dim += 1
        lo.append(hp["h2_mf_lo"])
        hi.append(hp["h2_mf_hi"])
    if not h2_T_fixed:
        dim += 1
        lo.append(hp["h2_T_out_lo"])
        hi.append(hp["h2_T_out_hi"])
    import numpy as np
    stds = [hp["sigma0"] * (hi[i] - lo[i]) for i in range(dim)]
    stds = np.array(stds, dtype=float) if dim == 1 else stds
    x0_init = [(lo[i] + hi[i]) / 2 for i in range(dim)]

    global_best_val = float("inf")
    global_best_x = x0_init[:]
    total_evals = 0

    for restart in range(hp["restarts_inner"]):
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
        x0 = x0_init if restart == 0 else [
            random.Random(seed + restart).uniform(lo[j], hi[j]) for j in range(dim)]
        es = cma.CMAEvolutionStrategy(x0, hp["sigma0"], opts)
        stall = 0
        gen = 0

        while not es.stop():
            solutions = es.ask()
            scores = []
            for sol in solutions:
                flows = list(sol[:n_sc])
                h2_mf, h2_T_out = _extract_extra(sol, n_sc)
                try:
                    base_layer.subcycle_mass_flows = flows
                    base_layer.commit_subcycle_mass_flows_to_topology()
                except Exception:
                    scores.append(1e9)
                    continue
                scores.append(_eval_fast(base_layer, h2_mf, h2_T_out, sys_inp, hp))
            es.tell(solutions, scores)
            total_evals += len(solutions)
            gen += 1

            if es.sigma < 0.5:
                break

            improved = False
            for s, sc in zip(solutions, scores):
                if sc < global_best_val:
                    global_best_val = sc
                    global_best_x = list(s)
                    improved = True
                    stall = 0
            if not improved:
                stall += 1
            if gen >= hp["maxiter_inner"] or stall >= 3 * len(solutions):
                break

        h2_val, h2tout_val = _extract_extra(global_best_x, n_sc)
        sys.stdout.write(
            f"    r{restart}: gen={gen} obj={global_best_val:.5f} "
            f"evals={total_evals} dim={dim} n_sc={n_sc} "
            f"best_flow_sum={sum(global_best_x[:n_sc]):.0f} "
            f"H2={h2_val:.1f} h2tout={h2tout_val:.0f}\n"
        )
        sys.stdout.flush()

    flows = list(global_best_x[:n_sc])
    h2_mf, h2_T_out = _extract_extra(global_best_x, n_sc)

    try:
        base_layer.subcycle_mass_flows = flows
        base_layer.commit_subcycle_mass_flows_to_topology()
    except Exception:
        pass
    obj_final = _eval_fast(base_layer, h2_mf, h2_T_out, sys_inp, hp)
    runtime = time.perf_counter() - t0

    return LayerResult(
        n=n_sc, t_min=tp_in.t_min, t_max=tp_in.t_max,
        p_min=tp_in.p_min, p_max=tp_in.p_max,
        t_q=tp_in.t_quantiles, p_q=tp_in.p_quantiles,
        s_q=tp_in.s_quantiles,
        flows=flows, h2_mf=h2_mf, h2_T_out=h2_T_out, obj=obj_final,
        n_subcycles=n_sc, n_evals=total_evals, runtime=runtime,
    ), base_layer


# ──────────────────────────────────────────────────────────────────
# B2. 内层 L-BFGS-B + 多起点
# ──────────────────────────────────────────────────────────────────


def _build_base_layer(tp_in, hp) -> ClosedCycleLayer | None:
    """构建拓扑基座（供 CMA 和 L-BFGS-B 共用）。"""
    t0 = time.perf_counter()
    use_he = hp.get("use_interp_he", True)
    if use_he:
        from core import InterpolatingHeliumSolver
        he_solver = InterpolatingHeliumSolver(
            T_min=tp_in.t_min, T_max=tp_in.t_max,
            P_min=tp_in.p_min, P_max=tp_in.p_max,
        )
        he_solver._build()
        solver = he_solver
    else:
        solver = None

    try:
        base_layer = ClosedCycleLayer(tp_in, properties=solver)
        if len(base_layer.subcycles) == 0:
            base_layer = ClosedCycleLayer(tp_in)
        else:
            ref = ClosedCycleLayer(tp_in)
            if len(ref.subcycles) != len(base_layer.subcycles):
                base_layer = ref
    except Exception:
        try:
            base_layer = ClosedCycleLayer(tp_in)
        except Exception:
            return None
    return base_layer


def _lbfgsb_single_worker(args: tuple) -> tuple[float, list[float], int]:
    """单起点 L-BFGS-B worker（模块级，供多进程调用）。"""
    (tp_in, sys_inp, hp, x0, lo, hi, n_sc, dim, h2_fixed, h2_T_fixed, maxiter_per) = args
    import scipy.optimize as _sopt

    layer = _build_base_layer(tp_in, hp)
    if layer is None or len(layer.subcycles) == 0:
        return (1e9, [], 0)

    _best_obj = [float("inf")]
    _best_x: list[float] = []

    def _obj(x: list[float]) -> float:
        flows = list(x[:n_sc])
        ei = n_sc
        _h2_mf = hp["h2_mf_lo"] if h2_fixed else x[ei]; ei += 0 if h2_fixed else 1
        _h2_T_out = hp["h2_T_out_lo"] if h2_T_fixed else x[ei]
        try:
            layer.subcycle_mass_flows = flows
            layer.commit_subcycle_mass_flows_to_topology()
        except Exception:
            return 1e9
        val = _eval_fast(layer, _h2_mf, _h2_T_out, sys_inp, hp)
        if val < _best_obj[0]:
            _best_obj[0] = val
            _best_x[:] = list(x)
        return val

    opts: dict = {"maxiter": maxiter_per, "ftol": 1e-8}
    _eps = hp.get("lbfgsb_eps", None)
    if _eps is not None:
        opts["eps"] = float(_eps)
    try:
        res = _sopt.minimize(
            _obj, x0, method="L-BFGS-B",
            bounds=[(lo[i], hi[i]) for i in range(dim)],
            options=opts,
        )
        n_evals = getattr(res, "nfev", maxiter_per * dim * 2) if hasattr(res, "nfev") else maxiter_per * dim * 2
        best_val = min(res.fun, _best_obj[0])
        best_x_out = list(_best_x) if _best_x and _best_obj[0] < res.fun else list(res.x)
        return (best_val, best_x_out, n_evals)
    except Exception:
        return (1e9, [], 0)


def _inner_lbfgsb_fast(
    tp_in: ClosedCycleTPInput,
    sys_inp: SystemInput,
    hp: dict,
    *,
    seed: int = 42,
    layer: ClosedCycleLayer | None = None,
) -> tuple[LayerResult, ClosedCycleLayer]:
    """内层 L-BFGS-B + 多起点优化子循环流量 + H2 出口温度。"""
    import random as _rnd

    t0 = time.perf_counter()

    if layer is not None:
        base_layer = layer
    else:
        base_layer = _build_base_layer(tp_in, hp)
        if base_layer is None:
            runtime = time.perf_counter() - t0
            return LayerResult(
                n=0, t_min=tp_in.t_min, t_max=tp_in.t_max,
                p_min=tp_in.p_min, p_max=tp_in.p_max,
                t_q=tp_in.t_quantiles, p_q=tp_in.p_quantiles,
                flows=[], h2_mf=hp["h2_mf_lo"], h2_T_out=hp["h2_T_out_lo"],
                obj=1.0, n_subcycles=0, n_evals=0, runtime=runtime,
            ), None

    n_sc = len(base_layer.subcycles)
    if n_sc == 0:
        runtime = time.perf_counter() - t0
        return LayerResult(
            n=0, t_min=tp_in.t_min, t_max=tp_in.t_max,
            p_min=tp_in.p_min, p_max=tp_in.p_max,
            t_q=tp_in.t_quantiles, p_q=tp_in.p_quantiles,
            flows=[], h2_mf=hp["h2_mf_lo"], h2_T_out=hp["h2_T_out_lo"],
            obj=1.0, n_subcycles=0, n_evals=0, runtime=runtime,
        ), base_layer

    h2_fixed = hp["h2_mf_lo"] == hp["h2_mf_hi"]
    h2_T_fixed = hp["h2_T_out_lo"] == hp["h2_T_out_hi"]

    dim = n_sc
    lo = [hp["mf_lo"]] * n_sc
    hi = [hp["mf_hi"]] * n_sc
    if not h2_fixed:
        dim += 1
        lo.append(hp["h2_mf_lo"])
        hi.append(hp["h2_mf_hi"])
    if not h2_T_fixed:
        dim += 1
        lo.append(hp["h2_T_out_lo"])
        hi.append(hp["h2_T_out_hi"])

    n_starts = hp.get("lbfgsb_starts", 20)
    maxiter_per = hp.get("lbfgsb_maxiter", 80)
    n_workers = hp.get("lbfgsb_workers", 1)
    rng = _rnd.Random(seed)

    sampler = hp.get("lbfgsb_sampler", "lhs")
    sobol_samples: list[list[float]] = []
    lhs_samples: list[list[float]] = []
    if n_starts > 1 and sampler == "sobol":
        try:
            from scipy.stats.qmc import Sobol
            sbl = Sobol(d=dim, scramble=True, seed=seed + 1)
            raw = sbl.random(n=n_starts)
            sobol_samples = []
            for i in range(n_starts):
                sobol_samples.append([lo[j] + (hi[j] - lo[j]) * float(raw[i, j])
                                      for j in range(dim)])
        except Exception:
            pass
    if not sobol_samples and n_starts > 1:
        lhs_samples = _lhs(n_starts, [(lo[i], hi[i]) for i in range(dim)], seed=seed + 1)

    start_points: list[list[float]] = []
    for i in range(n_starts):
        if i == 0:
            x0 = [(lo[j] + hi[j]) / 2 for j in range(dim)]
        elif sobol_samples:
            x0 = sobol_samples[i - 1]
        elif lhs_samples:
            x0 = lhs_samples[i - 1]
        else:
            x0 = [rng.uniform(lo[j], hi[j]) for j in range(dim)]
        start_points.append(x0)

    worker_args = [
        (tp_in, sys_inp, hp, x0, lo, hi, n_sc, dim, h2_fixed, h2_T_fixed, maxiter_per)
        for x0 in start_points
    ]

    best_val = float("inf")
    best_x: list[float] = []
    total_evals = 0

    if n_workers > 1:
        import concurrent.futures, multiprocessing
        ctx = multiprocessing.get_context("spawn")
        with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as ex:
            futures = [ex.submit(_lbfgsb_single_worker, a) for a in worker_args]
            for fut in futures:
                try:
                    val, x, n_ev = fut.result(timeout=600)
                    total_evals += n_ev
                    if val < best_val:
                        best_val = val
                        best_x = x
                except Exception:
                    pass
    else:
        for args_i in worker_args:
            val, x, n_ev = _lbfgsb_single_worker(args_i)
            total_evals += n_ev
            if val < best_val:
                best_val = val
                best_x = x

    if not best_x:
        runtime = time.perf_counter() - t0
        return LayerResult(
            n=n_sc, t_min=tp_in.t_min, t_max=tp_in.t_max,
            p_min=tp_in.p_min, p_max=tp_in.p_max,
            t_q=tp_in.t_quantiles, p_q=tp_in.p_quantiles,
            flows=[], h2_mf=hp["h2_mf_lo"], h2_T_out=hp["h2_T_out_lo"],
            obj=1.0, n_subcycles=n_sc, n_evals=0, runtime=runtime,
        ), base_layer

    flows = list(best_x[:n_sc])
    h2_mf = hp["h2_mf_lo"]
    h2_T_out = hp["h2_T_out_lo"]
    if not h2_fixed and not h2_T_fixed:
        h2_mf = best_x[n_sc]
        h2_T_out = best_x[n_sc + 1]
    elif not h2_fixed:
        h2_mf = best_x[n_sc]
    elif not h2_T_fixed:
        h2_T_out = best_x[n_sc]

    try:
        base_layer.subcycle_mass_flows = flows
        base_layer.commit_subcycle_mass_flows_to_topology()
    except Exception:
        pass
    obj_final = _eval_fast(base_layer, h2_mf, h2_T_out, sys_inp, hp)
    runtime = time.perf_counter() - t0

    sys.stdout.write(
        f"    lbfgsb: starts={n_starts} obj={obj_final:.5f} "
        f"evals={total_evals} dim={dim} n_sc={n_sc} "
        f"best_flow_sum={sum(flows):.0f} "
        f"H2={h2_mf:.1f} h2tout={h2_T_out:.0f}\n"
    )
    sys.stdout.flush()

    return LayerResult(
        n=n_sc, t_min=tp_in.t_min, t_max=tp_in.t_max,
        p_min=tp_in.p_min, p_max=tp_in.p_max,
        t_q=tp_in.t_quantiles, p_q=tp_in.p_quantiles,
        s_q=tp_in.s_quantiles,
        flows=flows, h2_mf=h2_mf, h2_T_out=h2_T_out, obj=obj_final,
        n_subcycles=n_sc, n_evals=total_evals, runtime=runtime,
    ), base_layer


# ──────────────────────────────────────────────────────────────────
# B3. 内层 CMA(粗) + L-BFGS-B(精) 两阶段
# ──────────────────────────────────────────────────────────────────


def _inner_hybrid_fast(
    tp_in: ClosedCycleTPInput,
    sys_inp: SystemInput,
    hp: dict,
    *,
    seed: int = 42,
    layer: ClosedCycleLayer | None = None,
) -> tuple[LayerResult, ClosedCycleLayer]:
    """两阶段内层优化：CMA 粗搜 + L-BFGS-B 精搜索。"""
    import cma
    import math
    import random as _rnd
    import scipy.optimize as _sopt

    t0 = time.perf_counter()

    if layer is not None:
        base_layer = layer
    else:
        base_layer = _build_base_layer(tp_in, hp)
        if base_layer is None:
            runtime = time.perf_counter() - t0
            return LayerResult(
                n=0, t_min=tp_in.t_min, t_max=tp_in.t_max,
                p_min=tp_in.p_min, p_max=tp_in.p_max,
                t_q=tp_in.t_quantiles, p_q=tp_in.p_quantiles,
                flows=[], h2_mf=hp["h2_mf_lo"], h2_T_out=hp["h2_T_out_lo"],
                obj=1.0, n_subcycles=0, n_evals=0, runtime=runtime,
            ), None

    n_sc = len(base_layer.subcycles)
    if n_sc == 0:
        runtime = time.perf_counter() - t0
        return LayerResult(
            n=0, t_min=tp_in.t_min, t_max=tp_in.t_max,
            p_min=tp_in.p_min, p_max=tp_in.p_max,
            t_q=tp_in.t_quantiles, p_q=tp_in.p_quantiles,
            flows=[], h2_mf=hp["h2_mf_lo"], h2_T_out=hp["h2_T_out_lo"],
            obj=1.0, n_subcycles=0, n_evals=0, runtime=runtime,
        ), base_layer

    h2_fixed = hp["h2_mf_lo"] == hp["h2_mf_hi"]
    h2_T_fixed = hp["h2_T_out_lo"] == hp["h2_T_out_hi"]

    dim = n_sc
    lo = [hp["mf_lo"]] * n_sc
    hi = [hp["mf_hi"]] * n_sc
    if not h2_fixed:
        dim += 1
        lo.append(hp["h2_mf_lo"])
        hi.append(hp["h2_mf_hi"])
    if not h2_T_fixed:
        dim += 1
        lo.append(hp["h2_T_out_lo"])
        hi.append(hp["h2_T_out_hi"])

    def _objective(x: list[float]) -> float:
        flows = list(x[:n_sc])
        ei = n_sc
        h2_mf = hp["h2_mf_lo"] if h2_fixed else x[ei]; ei += 0 if h2_fixed else 1
        h2_T_out = hp["h2_T_out_lo"] if h2_T_fixed else x[ei]
        try:
            base_layer.subcycle_mass_flows = flows
            base_layer.commit_subcycle_mass_flows_to_topology()
        except Exception:
            return 1e9
        return _eval_fast(base_layer, h2_mf, h2_T_out, sys_inp, hp)

    total_evals = 0
    best_x: list[float] = []
    best_val = float("inf")

    # ── 阶段 1: CMA 粗搜 ──
    cma_maxiter = hp.get("hybrid_cma_iter", 10)
    cma_restarts = hp.get("hybrid_cma_rst", 1)
    n_keep = hp.get("hybrid_n_keep", 3)
    sigma0_val = hp.get("hybrid_sigma0", 15.0)

    stds = [sigma0_val * (hi[i] - lo[i]) for i in range(dim)]
    x0_init = [(lo[i] + hi[i]) / 2 for i in range(dim)]
    top_k: list[tuple[float, list[float]]] = []  # (obj, x)

    for restart in range(cma_restarts):
        if restart % 2 == 0:
            p = max(4, int(4 + 3 * math.log2(max(dim, 2))))
        else:
            p = max(10, min(dim * dim, 50))

        opts = {
            "bounds": [lo, hi], "CMA_stds": stds, "popsize": p,
            "verbose": -9, "seed": seed + restart if seed is not None else None,
        }
        x0 = x0_init if restart == 0 else [
            _rnd.Random(seed + restart).uniform(lo[i], hi[i]) for i in range(dim)]
        es = cma.CMAEvolutionStrategy(x0, sigma0_val, opts)
        gen = 0

        while gen < cma_maxiter and not es.stop():
            solutions = es.ask()
            scores = [_objective(list(s)) for s in solutions]
            es.tell(solutions, scores)
            total_evals += len(solutions)
            gen += 1
            for s, sc in zip(solutions, scores):
                top_k.append((sc, list(s)))

        top_k.sort(key=lambda t: t[0])
        top_k = top_k[:n_keep]

    # ── 阶段 2: L-BFGS-B 精搜 ──
    lbfgsb_maxiter = hp.get("hybrid_lbfgsb_iter", 50)

    for si in range(min(n_keep, len(top_k))):
        obj_k, xk = top_k[si]
        try:
            res = _sopt.minimize(
                _objective, xk, method="L-BFGS-B",
                bounds=[(lo[i], hi[i]) for i in range(dim)],
                options={"maxiter": lbfgsb_maxiter, "ftol": 1e-12},
            )
        except Exception:
            continue
        n_evals = getattr(res, "nfev", 0)
        total_evals += n_evals
        if res.fun < best_val:
            best_val = res.fun
            best_x = list(res.x)

        sys.stdout.write(
            f"    hybrid_fine[{si}]: start_obj={obj_k:.5f} "
            f"→ fine_obj={res.fun:.5f} n_fev={n_evals}\n"
        )
        sys.stdout.flush()

    if not best_x:
        # fallback: use best CMA result
        if top_k:
            best_val, best_x = top_k[0]
            best_x = list(best_x)
        else:
            runtime = time.perf_counter() - t0
            return LayerResult(
                n=n_sc, t_min=tp_in.t_min, t_max=tp_in.t_max,
                p_min=tp_in.p_min, p_max=tp_in.p_max,
                t_q=tp_in.t_quantiles, p_q=tp_in.p_quantiles,
                flows=[], h2_mf=hp["h2_mf_lo"], h2_T_out=hp["h2_T_out_lo"],
                obj=1.0, n_subcycles=n_sc, n_evals=total_evals, runtime=runtime,
            ), base_layer

    flows = list(best_x[:n_sc])
    h2_mf = hp["h2_mf_lo"]
    h2_T_out = hp["h2_T_out_lo"]
    if not h2_fixed and not h2_T_fixed:
        h2_mf = best_x[n_sc]
        h2_T_out = best_x[n_sc + 1]
    elif not h2_fixed:
        h2_mf = best_x[n_sc]
    elif not h2_T_fixed:
        h2_T_out = best_x[n_sc]

    try:
        base_layer.subcycle_mass_flows = flows
        base_layer.commit_subcycle_mass_flows_to_topology()
    except Exception:
        pass
    obj_final = _eval_fast(base_layer, h2_mf, h2_T_out, sys_inp, hp)
    runtime = time.perf_counter() - t0

    sys.stdout.write(
        f"    hybrid: cma_iter={cma_maxiter} cma_rst={cma_restarts} "
        f"n_keep={n_keep} lbfgsb_iter={lbfgsb_maxiter} "
        f"obj={obj_final:.5f} evals={total_evals} dim={dim} n_sc={n_sc}\n"
    )
    sys.stdout.flush()

    return LayerResult(
        n=n_sc, t_min=tp_in.t_min, t_max=tp_in.t_max,
        p_min=tp_in.p_min, p_max=tp_in.p_max,
        t_q=tp_in.t_quantiles, p_q=tp_in.p_quantiles,
        s_q=tp_in.s_quantiles,
        flows=flows, h2_mf=h2_mf, h2_T_out=h2_T_out, obj=obj_final,
        n_subcycles=n_sc, n_evals=total_evals, runtime=runtime,
    ), base_layer


# ──────────────────────────────────────────────────────────────────
# C. 外层样本评估 (供 ProcessPoolExecutor 并行调用)
# ──────────────────────────────────────────────────────────────────


def _sample_worker(args: tuple) -> LayerResult:
    """模块级 worker: 评估一个外层 LHS 样本。"""
    (t_min, t_max, t_q_vals, p_q_vals, s_q_vals, p_max, p_min, sys_inp, hp, seed) = args

    tp_in = ClosedCycleTPInput(
        fluid="He", t_min=t_min, t_max=t_max,
        p_min=p_min, p_max=p_max,
        t_quantiles=t_q_vals, p_quantiles=p_q_vals,
        s_quantiles=s_q_vals,
    )

    from core import InterpolatingHeliumSolver
    if hp.get("use_interp_he", True):
        try:
            _he = InterpolatingHeliumSolver(
                T_min=t_min, T_max=t_max, P_min=p_min, P_max=p_max)
            _he._build()
            probe = ClosedCycleLayer(tp_in, properties=_he)
            if len(probe.subcycles) == 0:
                probe = ClosedCycleLayer(tp_in)
            else:
                ref = ClosedCycleLayer(tp_in)
                if len(ref.subcycles) != len(probe.subcycles):
                    probe = ref
        except Exception:
            try:
                probe = ClosedCycleLayer(tp_in)
            except Exception:
                return LayerResult(n=0, t_min=t_min, t_max=t_max,
                                   p_min=p_min, p_max=p_max,
                                   t_q=t_q_vals, p_q=p_q_vals, flows=[], h2_mf=0,
                                   obj=1e9, n_subcycles=0, n_evals=0, runtime=0)
    else:
        try:
            probe = ClosedCycleLayer(tp_in)
        except Exception:
            return LayerResult(n=0, t_min=t_min, t_max=t_max,
                               p_min=p_min, p_max=p_max,
                               t_q=t_q_vals, p_q=p_q_vals, flows=[], h2_mf=0,
                               obj=1e9, n_subcycles=0, n_evals=0, runtime=0)
    if len(probe.subcycles) == 0:
        return LayerResult(n=0, t_min=t_min, t_max=t_max,
                           p_min=p_min, p_max=p_max,
                           t_q=t_q_vals, p_q=p_q_vals, flows=[], h2_mf=0,
                           obj=1e9, n_subcycles=0, n_evals=0, runtime=0)

    method = hp.get("inner_method", "cma")
    if method == "lbfgsb":
        result, _ = _inner_lbfgsb_fast(tp_in, sys_inp, hp, seed=seed)
    elif method == "hybrid":
        result, _ = _inner_hybrid_fast(tp_in, sys_inp, hp, seed=seed)
    else:
        result, _ = _inner_cma_fast(tp_in, sys_inp, hp, seed=seed)
    return result


def _de_trial_worker(args: tuple) -> float:
    """模块级 worker: 外层参数 → 内层 CMA → obj（仅回传标量）。"""
    (x_outer, sys_inp, hp, seed) = args

    t_min = round(x_outer[0])
    t_max = round(x_outer[1])
    idx = 2
    tq = tuple(float(x_outer[idx + j]) for j in range(hp["n_t_q"]))
    idx += hp["n_t_q"]
    pq = tuple(float(x_outer[idx + j]) for j in range(hp["n_p_q"]))
    idx += hp["n_p_q"]
    sq = tuple(float(x_outer[idx + j]) for j in range(hp["n_s_q"]))
    idx += hp["n_s_q"]
    p_max = round(x_outer[idx])
    p_min = round(x_outer[idx + 1])

    tp_in = ClosedCycleTPInput(
        fluid="He", t_min=t_min, t_max=t_max,
        p_min=p_min, p_max=p_max,
        t_quantiles=tq, p_quantiles=pq,
        s_quantiles=sq,
    )
    from core import InterpolatingHeliumSolver
    try:
        if hp.get("use_interp_he", True):
            try:
                _he = InterpolatingHeliumSolver(
                    T_min=t_min, T_max=t_max, P_min=p_min, P_max=p_max)
                _he._build()
                probe = ClosedCycleLayer(tp_in, properties=_he)
                if len(probe.subcycles) == 0:
                    probe = ClosedCycleLayer(tp_in)
                else:
                    ref = ClosedCycleLayer(tp_in)
                    if len(ref.subcycles) != len(probe.subcycles):
                        probe = ref
            except Exception:
                probe = ClosedCycleLayer(tp_in)
        else:
            probe = ClosedCycleLayer(tp_in)
    except Exception:
        return 1e9
    if len(probe.subcycles) == 0:
        return 1e9

    method = hp.get("inner_method", "cma")
    if method == "lbfgsb":
        result, _ = _inner_lbfgsb_fast(tp_in, sys_inp, hp, seed=seed)
    elif method == "hybrid":
        result, _ = _inner_hybrid_fast(tp_in, sys_inp, hp, seed=seed)
    else:
        result, _ = _inner_cma_fast(tp_in, sys_inp, hp, seed=seed)
    return result.obj


# ──────────────────────────────────────────────────────────────────
# E. 最优解重建 + CoolProp 验证
# ──────────────────────────────────────────────────────────────────


def _rebuild_result(x_outer: list[float], sys_inp: SystemInput,
                    hp: dict, *, use_coolprop: bool = False) -> tuple:
    """用最优参数重建层、获取 report + HX；可选 CoolProp 模式。"""
    props = PropertyRegistry()
    t_min = round(x_outer[0])
    t_max = round(x_outer[1])
    idx = 2
    tq = tuple(round(round(x_outer[idx + j] / hp["qstep"]) * hp["qstep"], 3)
               for j in range(hp["n_t_q"]))
    idx += hp["n_t_q"]
    pq = tuple(round(round(x_outer[idx + j] / hp["qstep"]) * hp["qstep"], 3)
               for j in range(hp["n_p_q"]))
    idx += hp["n_p_q"]
    sq = tuple(round(round(x_outer[idx + j] / hp["qstep"]) * hp["qstep"], 3)
               for j in range(hp["n_s_q"]))
    idx += hp["n_s_q"]
    p_max = round(x_outer[idx])
    p_min = round(x_outer[idx + 1])

    tp_in = ClosedCycleTPInput(
        fluid="He", t_min=t_min, t_max=t_max,
        p_min=p_min, p_max=p_max,
        t_quantiles=tq, p_quantiles=pq,
        s_quantiles=sq,
    )
    from core import InterpolatingHeliumSolver
    try:
        if hp.get("use_interp_he", True):
            he_solver = InterpolatingHeliumSolver(
                T_min=t_min, T_max=t_max, P_min=p_min, P_max=p_max)
            he_solver._build()
            layer = ClosedCycleLayer(tp_in, properties=he_solver)
            if len(layer.subcycles) == 0:
                layer = ClosedCycleLayer(tp_in)
            else:
                ref = ClosedCycleLayer(tp_in)
                if len(ref.subcycles) != len(layer.subcycles):
                    layer = ref
        else:
            layer = ClosedCycleLayer(tp_in)
    except Exception:
        layer = ClosedCycleLayer(tp_in)

    if hp.get("inner_method", "cma") == "lbfgsb":
        result, layer = _inner_lbfgsb_fast(tp_in, sys_inp, hp, seed=0, layer=layer)
    elif hp.get("inner_method", "cma") == "hybrid":
        result, layer = _inner_hybrid_fast(tp_in, sys_inp, hp, seed=0, layer=layer)
    else:
        result, layer = _inner_cma_fast(tp_in, sys_inp, hp, seed=0, layer=layer)

    if use_coolprop:
        try:
            layer_cp = ClosedCycleLayer(tp_in)
        except Exception:
            return result, None, None, None, layer, [], []
        try:
            if len(layer_cp.subcycles) == len(result.flows):
                layer_cp.subcycle_mass_flows = list(result.flows)
                layer_cp.commit_subcycle_mass_flows_to_topology()
        except Exception:
            pass
        ideal_rep_cp = layer_cp.performance_report()
        if hp["use_non_ideal"]:
            ni_cp = layer_cp.ensure_non_ideal()
            ni_cp.apply_offsets()
            ni_rep_cp = layer_cp.performance_report()
        else:
            ni_rep_cp = ideal_rep_cp
        return result, None, ideal_rep_cp, ni_rep_cp, layer_cp, [], []

    ideal_report = layer.performance_report()
    if hp["use_non_ideal"]:
        ni = layer.ensure_non_ideal()
        ni.apply_offsets()
        ni_report = layer.performance_report()
    else:
        ni_report = ideal_report
    hots: list = []
    colds: list = []
    for _, rec in ni_report.by_edge:
        if rec.kind == "heat" and rec.power_rate is not None:
            if rec.category == ProcessCategory.HEAT_REJECTION:
                hots.append(rec)
            else:
                colds.append(rec)
    from core.system import convert_sources
    cold_src = _make_h2_source(hp, result.h2_mf, result.h2_T_out)
    src_hots, _ = convert_sources(sys_inp.heat_sources, (), props)
    _, src_colds = convert_sources((), (cold_src,), props)
    hots[:0] = src_hots
    colds[:0] = src_colds
    obj_mode = hp.get("obj_mode", "series" if hp.get("hx_series") else "group")
    if obj_mode in ("pinch", "eff_pinch"):
        from core.postprocess import analyze_pinch
        pinch_result = analyze_pinch(colds, hots, hp["hx_dT"], props)
        hx = pinch_result
    elif hp.get("hx_series", False):
        from core.heat_exchanger import match_series_pinch
        hx = match_series_pinch(hots, colds, dT_min=hp["hx_dT"])
    else:
        hx = match_heat_exchanger_groups(hots, colds, dT_min=hp["hx_dT"],
                                          max_group_size=hp["hx_max_group_size"])
    return result, hx, ideal_report, ni_report, layer, hots, colds


# ──────────────────────────────────────────────────────────────────
# 绘图函数
# ──────────────────────────────────────────────────────────────────


def _draw_cycle_ts_ps(report, layer, label: str, out_path: Path,
                       ni_nodes: dict | None = None) -> None:
    """绘制单个循环 TS/PS 图（带子循环多边形覆盖 + 流量标注）。
    
    非理想图传入 ``ni_nodes`` 以正确对齐子循环多边形位置。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    overlay_nodes = ni_nodes if ni_nodes is not None else layer.nodes

    fig, (ax_ts, ax_ps) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(f"Cycle T\u2013S & P\u2013S [{label}]", fontsize=11)

    for _, snap in report.nodes:
        ax_ts.scatter(snap.S, snap.T, s=10, color="0.3", zorder=2)
        ax_ps.scatter(snap.S, snap.P, s=10, color="0.3", zorder=2)

    for _, rec in report.by_edge:
        c = "tab:orange" if rec.kind == "mechanical" else "tab:green"
        ax_ts.plot(
            [rec.tail_state.S, rec.head_state.S],
            [rec.tail_state.T, rec.head_state.T],
            color=c, linewidth=1.0, alpha=0.7,
        )
        ax_ps.plot(
            [rec.tail_state.S, rec.head_state.S],
            [rec.tail_state.P, rec.head_state.P],
            color=c, linewidth=1.0, alpha=0.7,
        )
        ts_fr, ts_to = 0.38, 0.62
        for ax, y1, y2 in [
            (ax_ts, rec.tail_state.T, rec.head_state.T),
            (ax_ps, rec.tail_state.P, rec.head_state.P),
        ]:
            ax.annotate(
                "", xytext=(
                    rec.tail_state.S + ts_fr * (rec.head_state.S - rec.tail_state.S),
                    y1 + ts_fr * (y2 - y1),
                ), xy=(
                    rec.tail_state.S + ts_to * (rec.head_state.S - rec.tail_state.S),
                    y1 + ts_to * (y2 - y1),
                ), arrowprops=dict(arrowstyle="->", color="0.3", lw=1.2, alpha=0.7))
        mx_s = (rec.tail_state.S + rec.head_state.S) / 2
        ax_ts.annotate(rec.edge_key, (
            mx_s, (rec.tail_state.T + rec.head_state.T) / 2),
            textcoords="offset points", xytext=(3, 3), fontsize=5, color="0.4")

    colors_sc = plt.cm.tab10.colors
    for i, sc in enumerate(layer.subcycles):
        c = colors_sc[i % len(colors_sc)]
        q = sc.mass_flow if sc.mass_flow is not None else 0.0
        n0, n1, n2, n3 = [overlay_nodes[idx] for idx in sc.nodes]
        for ax, y_get in [(ax_ts, lambda n: n.T), (ax_ps, lambda n: n.P)]:
            xs = [n0.S, n1.S, n2.S, n3.S, n0.S]
            ys = [y_get(n0), y_get(n1), y_get(n2), y_get(n3), y_get(n0)]
            ax.plot(xs, ys, "-", color=c, linewidth=1.2, alpha=0.5)
            cx = (n0.S + n2.S) / 2
            cy = (y_get(n0) + y_get(n2)) / 2
            ax.annotate(f"SC{i}\n{q:.1f}", (cx, cy), ha="center", va="center",
                        fontsize=6, color=c, fontweight="bold")

    ax_ts.set_xlabel("S [kJ/(kg·K)]"); ax_ts.set_ylabel("T [K]")
    ax_ts.set_title(f"{label} T\u2013S"); ax_ts.grid(True, alpha=0.25)
    ax_ps.set_xlabel("S [kJ/(kg·K)]"); ax_ps.set_ylabel("P [kPa]")
    ax_ps.set_title(f"{label} P\u2013S"); ax_ps.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _draw_hx_tq(hots: list, colds: list, hx_result: HXMatchResult,
                out_path: Path) -> None:
    """绘制 HX 匹配 T-Q 图：概览 + 每单元详情。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    all_recs, _, _ = _normalize_records(list(hots), list(colds))
    if not all_recs:
        return
    sorted_recs = sorted(all_recs, key=lambda r: (r.T_high, 0 if r.is_hot else 1),
                         reverse=True)

    unit_rec_keys: list[set[str]] = []
    for u in hx_result.units:
        keys: set[str] = set()
        for r in u.hot_records:
            keys.add(r.edge_key)
        for r in u.cold_records:
            keys.add(r.edge_key)
        unit_rec_keys.append(keys)
    ua_keys: set[str] = set()
    for r in hx_result.unassigned_hots:
        ua_keys.add(r.edge_key)
    for r in hx_result.unassigned_colds:
        ua_keys.add(r.edge_key)

    unit_colors = ["#FFD0D0", "#D0D0FF", "#D0FFD0", "#FFFFD0",
                   "#FFD0FF", "#D0FFFF", "#FFE8D0", "#E8D0FF"]

    n_units = max(len(hx_result.units), 1)
    n_rows = (n_units + 2) // 3
    fig = plt.figure(figsize=(max(14, 5 * min(n_units, 3)), 3.5 + 3.5 * n_rows))
    fig.suptitle(
        f"HX Optimal Matching | {len(hots) + len(colds)} records, "
        f"{hx_result.num_units} units | "
        f"matched={hx_result.total_matched:.0f}kW, "
        f"unmatched={hx_result.total_unmatched:.0f}kW",
        fontsize=12, fontweight="bold")

    gs = fig.add_gridspec(n_rows + 1, 1, height_ratios=[1] + [2] * n_rows,
                          hspace=0.45, top=0.90)

    # ── Row 0: 概览 ──
    ax_ov = fig.add_subplot(gs[0])
    cum_q = 0.0
    for r in sorted_recs:
        color = "tab:red" if r.is_hot else "tab:blue"
        ax_ov.plot([cum_q, cum_q + r.Q], [r.T_high, r.T_low],
                   color=color, linewidth=2.2, marker="o", markersize=4)
        label = r.record.edge_key
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
        else:
            ax_ov.annotate(f"{label}\n{r.Q:.0f}kW",
                           (cum_q + r.Q / 2, (r.T_high + r.T_low) / 2),
                           fontsize=6, ha="center", color="gray",
                           bbox=dict(boxstyle="round,pad=0.1",
                                     facecolor="#EEEEEE", alpha=0.7))
        cum_q += r.Q
    ax_ov.set_xlabel("Cumulative Q [kW]"); ax_ov.set_ylabel("T [K]")
    ax_ov.grid(True, alpha=0.2)
    ax_ov.legend(handles=[
        mpatches.Patch(color="tab:red", label="hot"),
        mpatches.Patch(color="tab:blue", label="cold"),
        mpatches.Patch(color="gray", label="unassigned"),
    ], fontsize=8, loc="upper left", ncol=3)
    ax_ov.set_title("Overview — all records (sorted by T_high)", fontsize=9, loc="left")

    # ── 每单元详情 ──
    for ui, unit in enumerate(hx_result.units):
        row = 1 + ui // 3
        col = ui % 3
        sub_gs = gs[row].subgridspec(1, min(n_units - (ui // 3) * 3, 3), wspace=0.35)
        ax = fig.add_subplot(sub_gs[0, col])

        h_sorted = sorted(unit.hot_records,
                          key=lambda rr: max(rr.tail_state.T, rr.head_state.T), reverse=True)
        c_sorted = sorted(unit.cold_records,
                          key=lambda rr: max(rr.tail_state.T, rr.head_state.T), reverse=True)

        x = 0.0
        for hr in h_sorted:
            Th_h, Th_l = max(hr.tail_state.T, hr.head_state.T), min(hr.tail_state.T, hr.head_state.T)
            Qh = abs(hr.power_rate) if hr.power_rate else 0.0
            ax.plot([x, x + Qh], [Th_h, Th_l], "o-", color="tab:red", lw=2.2, ms=3.5)
            ax.annotate(f"{hr.edge_key}\n({Qh:.0f}kW)",
                        (x + Qh / 2, (Th_h + Th_l) / 2),
                        fontsize=5, color="darkred", ha="center",
                        xytext=(0, 8), textcoords="offset points")
            x += Qh
        q_max = x
        x = 0.0
        for cr in c_sorted:
            Tc_l, Tc_h = min(cr.tail_state.T, cr.head_state.T), max(cr.tail_state.T, cr.head_state.T)
            Qc = abs(cr.power_rate) if cr.power_rate else 0.0
            ax.plot([x, x + Qc], [Tc_h, Tc_l], "s--", color="tab:blue", lw=1.8, ms=3.5)
            ax.annotate(f"{cr.edge_key}\n({Qc:.0f}kW)",
                        (x + Qc / 2, (Tc_h + Tc_l) / 2),
                        fontsize=5, color="darkblue", ha="center",
                        xytext=(0, -12), textcoords="offset points")
            x += Qc

        ax.set_title(f"Unit {ui + 1}  matched={unit.matched_heat:.0f}kW  "
                     f"pinch={unit.internal_pinch:.0f}K", fontsize=8, fontweight="bold")
        ax.set_xlabel("Q [kW]"); ax.set_ylabel("T [K]")
        ax.grid(True, alpha=0.2)
        ax.set_xlim(0, max(max(q_max, x), 1))

    fig.subplots_adjust(left=0.05, right=0.97, bottom=0.05)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _draw_pinch_tq(hots: list, colds: list, pinch_result,
                   dT_min: float, out_path: Path) -> None:
    """绘制夹点分析 T-Q 图：逐个过程线段 + 复合曲线 + 夹点标记。

    热源(Air)和冷源(H2)用粗线高亮，循环换热过程用细线；
    叠加 ``PinchResult`` 中的复合曲线和夹点标记。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax_proc, ax_comp) = plt.subplots(2, 1, figsize=(14, 12),
                                            gridspec_kw={"height_ratios": [1, 1], "hspace": 0.35})
    fig.suptitle(f"Pinch Analysis | dT_min={dT_min:.0f}K | "
                 f"Hot util={pinch_result.hot_utility_demand:.0f}kW  "
                 f"Cold util={pinch_result.cold_utility_demand:.0f}kW  "
                 f"ΔQ={pinch_result.delta_Q:.0f}kW",
                 fontsize=11, fontweight="bold")

    # ── 上子图: 逐个过程 T-Q 线段 ──
    from core.heat_exchanger import _normalize_records
    all_recs, _, _ = _normalize_records(list(hots), list(colds))
    sorted_recs = sorted(all_recs, key=lambda r: r.T_high, reverse=True)

    ax_proc.set_title("Individual heat processes (sorted by T_high)", fontsize=9, loc="left")
    cum_q = 0.0
    for r in sorted_recs:
        is_hot = r.is_hot
        edge_key = r.record.edge_key
        # 高亮外部源
        is_source = edge_key in ("Air_source", "H2_source")  # convert_sources 生成的 key
        lw = 3.0 if is_source else 1.2
        alpha_val = 1.0 if is_source else 0.55
        color = "darkred" if (is_hot and is_source) else \
                "darkblue" if (not is_hot and is_source) else \
                "tab:red" if is_hot else "tab:blue"

        Th = r.T_high; Tl = r.T_low
        ax_proc.plot([cum_q, cum_q + r.Q], [Th, Tl],
                     color=color, linewidth=lw, alpha=alpha_val, marker="o", markersize=3)

        # 温度标注
        ax_proc.annotate(f"{Th:.0f}K", (cum_q, Th),
                         textcoords="offset points", xytext=(-2, 4),
                         fontsize=5, color=color, ha="right")
        ax_proc.annotate(f"{Tl:.0f}K", (cum_q + r.Q, Tl),
                         textcoords="offset points", xytext=(2, -8),
                         fontsize=5, color=color, ha="left")

        # 过程标签
        if is_source:
            ax_proc.annotate(f"{edge_key}\n{r.Q:.0f} kW",
                             (cum_q + r.Q / 2, (Th + Tl) / 2),
                             fontsize=7, color=color, ha="center", fontweight="bold",
                             bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.85))
        else:
            ax_proc.annotate(f"{edge_key}",
                             (cum_q + r.Q / 2, (Th + Tl) / 2),
                             fontsize=5, color=color, ha="center", alpha=0.7)
        cum_q += r.Q

    ax_proc.set_xlabel("Cumulative Q [kW]")
    ax_proc.set_ylabel("T [K]")
    ax_proc.grid(True, alpha=0.2)
    from matplotlib.patches import Patch
    ax_proc.legend(handles=[
        Patch(color="darkred", label="Air source"),
        Patch(color="darkblue", label="H₂ source"),
        Patch(color="tab:red", label="Cycle rejection"),
        Patch(color="tab:blue", label="Cycle absorption"),
    ], fontsize=8, loc="upper left", ncol=2)

    # ── 下子图: 复合曲线 + 夹点 ──
    ax = ax_comp
    ax.set_title("Composite curves (pinch analysis)", fontsize=9, loc="left")

    delta_Q = pinch_result.delta_Q

    # 热复合曲线（放热侧）
    if pinch_result.matched_rejection.q_points:
        q_h = list(pinch_result.matched_rejection.q_points)
        t_h = list(pinch_result.matched_rejection.t_points)
        ax.plot(q_h, t_h, "-", color="darkred", lw=2.5, label="Hot composite")
        if pinch_result.extra_rejection is not None and pinch_result.extra_rejection.q_points:
            q_hr = list(pinch_result.extra_rejection.q_points)
            t_hr = list(pinch_result.extra_rejection.t_points)
            ax.plot(q_hr, t_hr, "-", color="darkred", lw=2.5)

    # 冷复合曲线（吸热侧，平移后）
    if pinch_result.matched_absorption.q_points:
        q_c = [q + delta_Q for q in pinch_result.matched_absorption.q_points]
        t_c = list(pinch_result.matched_absorption.t_points)
        ax.plot(q_c, t_c, "-", color="darkblue", lw=2.5, label="Cold composite (shifted)")
        if pinch_result.extra_absorption is not None and pinch_result.extra_absorption.q_points:
            q_cr = [q + delta_Q for q in pinch_result.extra_absorption.q_points]
            t_cr = list(pinch_result.extra_absorption.t_points)
            ax.plot(q_cr, t_cr, "-", color="darkblue", lw=2.5)

    # 夹点标记
    p_th = pinch_result.pinch_T_hot
    p_tc = pinch_result.pinch_T_cold
    p_qh = pinch_result.pinch_Q_hot
    if p_qh > 0 and p_th > 0:
        ax.plot([p_qh], [p_th], "o", color="black", ms=12, zorder=5)
        ax.annotate(f"Pinch\nT_h={p_th:.0f}K\nT_c={p_tc:.0f}K\nΔT_min={dT_min:.0f}K",
                    (p_qh, p_th), textcoords="offset points", xytext=(12, 12),
                    fontsize=9, color="black",
                    bbox=dict(boxstyle="round", facecolor="yellow", alpha=0.8))

    # 公用工程阴影标注
    all_q_h = list(pinch_result.matched_rejection.q_points)
    all_q_c_shifted = [q + delta_Q for q in pinch_result.matched_absorption.q_points]
    q_max_comp = max(max(all_q_h, default=0), max(all_q_c_shifted, default=0))

    hot_util = pinch_result.hot_utility_demand
    cold_util = pinch_result.cold_utility_demand
    q_plot_max = q_max_comp + delta_Q + max(hot_util, 1000)

    # ── 源温度段高亮：从记录中提取外部源温度范围 ──
    src_air_T_hi = max(float(hots[0].tail_state.T), float(hots[0].head_state.T)) if hots else 1250.0
    src_air_T_lo = min(float(hots[0].tail_state.T), float(hots[0].head_state.T)) if hots else 500.0
    src_h2_T_hi  = max(float(colds[0].tail_state.T), float(colds[0].head_state.T)) if colds else 800.0
    src_h2_T_lo  = min(float(colds[0].tail_state.T), float(colds[0].head_state.T)) if colds else 20.0

    ax.hlines([src_air_T_lo, src_air_T_hi], 0, q_plot_max, colors="darkred", lw=1.5, ls="--", alpha=0.5)
    ax.fill_between([0, q_plot_max], src_air_T_lo, src_air_T_hi,
                     color="darkred", alpha=0.06)
    ax.annotate(f"Air source\n{src_air_T_hi:.0f}->{src_air_T_lo:.0f}K",
                (q_plot_max * 0.97, (src_air_T_lo + src_air_T_hi) / 2),
                ha="right", va="center", fontsize=8, color="darkred",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

    ax.hlines([src_h2_T_lo, src_h2_T_hi], 0, q_plot_max, colors="darkblue", lw=1.5, ls="--", alpha=0.5)
    ax.fill_between([0, q_plot_max], src_h2_T_lo, src_h2_T_hi,
                     color="darkblue", alpha=0.06)
    ax.annotate(f"H2 source\n{src_h2_T_lo:.0f}->{src_h2_T_hi:.0f}K",
                (q_plot_max * 0.97, (src_h2_T_lo + src_h2_T_hi) / 2),
                ha="right", va="center", fontsize=8, color="darkblue",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

    if cold_util > 1e-6:
        ax.axvspan(0, delta_Q, alpha=0.10, color="darkred")
        ax.annotate(f"Cold utility\n{cold_util:.0f} kW",
                    (delta_Q / 2, 250), ha="center", fontsize=9, color="darkred",
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
    if hot_util > 1e-6:
        ax.axvspan(q_max_comp + delta_Q, q_max_comp + delta_Q + hot_util,
                    alpha=0.10, color="darkblue")
        ax.annotate(f"Hot utility\n{hot_util:.0f} kW",
                    (q_max_comp + delta_Q + hot_util / 2, 250),
                    ha="center", fontsize=9, color="darkblue",
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

    ax.set_xlabel("Q [kW]")
    ax.set_ylabel("T [K]")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.2)

    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _draw_perf_bars(report, label: str, out_path: Path) -> None:
    """绘制循环性能柱状图（机械功率 + 类别汇总）。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax_mech, ax_cat) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"Cycle performance [{label}]", fontsize=11)

    mech_data = []
    cat_sums = {pc: 0.0 for pc in (
        ProcessCategory.COMPRESSION, ProcessCategory.EXPANSION,
        ProcessCategory.HEAT_ABSORPTION, ProcessCategory.HEAT_REJECTION)}
    for _, rec in report.by_edge:
        if rec.power_rate is None:
            continue
        cat_sums[rec.category] += rec.power_rate
        if rec.kind == "mechanical":
            signed = abs(rec.power_rate) if rec.category == ProcessCategory.COMPRESSION else -abs(rec.power_rate)
            mech_data.append((rec.edge_key, signed, abs(rec.power_rate)))
    mech_data.sort(key=lambda x: x[2])
    if mech_data:
        labels, vals_m = zip(*[(r[0], r[1]) for r in mech_data])
        colors = ["tab:blue" if v >= 0 else "tab:orange" for v in vals_m]
        ax_mech.bar(range(len(labels)), vals_m, color=colors, alpha=0.9)
        ax_mech.set_xticks(range(len(labels)))
        ax_mech.set_xticklabels(labels, rotation=60, ha="right", fontsize=6)
    ax_mech.axhline(0, color="0.2"); ax_mech.set_title("mechanical [kW]")
    ax_mech.grid(True, axis="y", alpha=0.25)

    cn = ["compression", "expansion", "absorption", "rejection"]
    cv = [cat_sums[pc] for pc in (
        ProcessCategory.COMPRESSION, ProcessCategory.EXPANSION,
        ProcessCategory.HEAT_ABSORPTION, ProcessCategory.HEAT_REJECTION)]
    ax_cat.bar(cn, cv, color=["tab:orange", "tab:blue", "tab:red", "tab:green"], alpha=0.9)
    ax_cat.axhline(0, color="0.2"); ax_cat.set_title("category sums [kW]")
    ax_cat.grid(True, axis="y", alpha=0.25)
    ax_cat.text(3.5, max(abs(v) for v in cv) * 0.9,
                f"sum mech={cv[0] + cv[1]:.0f}\nsum heat={cv[2] + cv[3]:.0f}",
                fontsize=8, ha="center",
                bbox=dict(boxstyle="round", alpha=0.2))

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────
# 主测试: LHS 种子 + DE 外层
# ──────────────────────────────────────────────────────────────────


def test_layered_optimization() -> None:
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 2, "n_lhs": 50, "hx_dT": 10.0,
          "maxiter_inner": 10, "restarts_inner": 2,
          "de_popsize": 15, "de_maxiter": 20, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001}
    _run_layered(hp, tag="2P0T")


def test_layered_1t2p() -> None:
    """外层加 1 个 T 分位 → 更多拓扑自由度。"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 1, "n_p_q": 2, "n_lhs": 80, "hx_dT": 10.0,
          "maxiter_inner": 10, "restarts_inner": 2,
          "de_popsize": 20, "de_maxiter": 25, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001}
    _run_layered(hp, tag="1T2P")


def test_layered_2p0t_long() -> None:
    """2P+0T 长跑: 更多 LHS 种子 + 更多 DE 代数。"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 2, "n_lhs": 100, "hx_dT": 10.0,
          "maxiter_inner": 10, "restarts_inner": 2,
          "de_popsize": 20, "de_maxiter": 40, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001}
    _run_layered(hp, tag="2P0T_long")


def test_layered_0p0t() -> None:
    """0P+0T 基线: 无分位，仅调 t_min/t_max/p_max。"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 0, "n_lhs": 30, "hx_dT": 10.0,
          "maxiter_inner": 10, "restarts_inner": 2,
          "de_popsize": 10, "de_maxiter": 15, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001}
    _run_layered(hp, tag="0P0T")


def test_layered_0p0t_wide() -> None:
    """0P+0T 宽边界 + 无非理想：t_min∈[50,500] t_max∈[700,1100] p_min∈[1000,4000] p_max∈[6000,15000] H2∈[2,6]"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 0, "n_lhs": 50, "hx_dT": 10.0,
          "maxiter_inner": 10, "restarts_inner": 2,
          "de_popsize": 15, "de_maxiter": 40, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 500.0,
          "t_max_lo": 700.0, "t_max_hi": 1100.0,
          "p_min_lo": 1000.0, "p_min_hi": 4000.0,
          "p_max_lo": 6000.0, "p_max_hi": 15000.0,
          "h2_mf_lo": 2.0, "h2_mf_hi": 6.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="0P0T_wide")


def test_layered_1p0t_wide() -> None:
    """1P+0T 宽边界 + 无非理想 + CoolProp + 变量 H2。"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 1, "n_lhs": 80, "hx_dT": 10.0,
          "maxiter_inner": 10, "restarts_inner": 2,
          "de_popsize": 20, "de_maxiter": 40, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 500.0,
          "t_max_lo": 700.0, "t_max_hi": 1100.0,
          "p_min_lo": 1000.0, "p_min_hi": 4000.0,
          "p_max_lo": 6000.0, "p_max_hi": 15000.0,
          "h2_mf_lo": 2.0, "h2_mf_hi": 6.0,
          "use_non_ideal": False,
          "use_interp_he": False,
          "seed": 123,
          }
    _run_layered(hp, tag="1P0T_cp_s123")


def test_layered_2p0t_ideal() -> None:
    """2P+0T 理想 t_max∈[800,1100] t_min∈[50,400] p_min∈[2000,4000] p_max∈[8000,12000] H2∈[3,6] H2_T_out∈[400,900]"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 2, "n_lhs": 50, "hx_dT": 10.0,
          "hx_max_group_size": 2,
          "maxiter_inner": 10, "restarts_inner": 2,
          "de_popsize": 15, "de_maxiter": 40, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 4000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 3.0, "h2_mf_hi": 6.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 900.0,
          "use_non_ideal": False,
           }
    _run_layered(hp, tag="2P0T_HX2")


def test_layered_2p1t_ideal() -> None:
    """2P+1T 理想 t_max∈[800,1100] t_min∈[50,400] p_min∈[2000,4000] p_max∈[8000,12000] H2∈[3,6] HXmax=2"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 1, "n_p_q": 2, "n_lhs": 60, "hx_dT": 10.0,
          "hx_max_group_size": 2,
          "maxiter_inner": 10, "restarts_inner": 2,
          "de_popsize": 20, "de_maxiter": 40, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 4000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 3.0, "h2_mf_hi": 6.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 900.0,
          "use_non_ideal": False,
           }
    _run_layered(hp, tag="2P1T_HX2")


def test_layered_1p1t_ideal() -> None:
    """1P+1T 理想 t_max∈[800,1100] t_min∈[50,400] p_min∈[2000,4000] p_max∈[8000,12000] H2∈[3,6] HXmax=2"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 1, "n_p_q": 1, "n_lhs": 60, "hx_dT": 10.0,
          "hx_max_group_size": 2,
          "maxiter_inner": 10, "restarts_inner": 2,
          "de_popsize": 20, "de_maxiter": 40, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 4000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 3.0, "h2_mf_hi": 6.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 900.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="1P1T_HX2")


def test_layered_0p0t_hx1() -> None:
    """0P0T0S 理想 HXmax=1 H2固定4.3 h2_T_out∈[400,900] p_min=2000kPa t_max∈[800,1100] t_min∈[50,400] p_max∈[8000,12000]"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 0, "n_s_q": 0, "n_lhs": 30, "hx_dT": 10.0,
          "hx_max_group_size": 1,
          "maxiter_inner": 10, "restarts_inner": 2,
          "de_popsize": 10, "de_maxiter": 25, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 2000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 4.3, "h2_mf_hi": 4.3,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 900.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="0P0T0S_HX1")


def test_layered_0p1t_hx1() -> None:
    """1P+0T+1S 理想 串联夹点 H2∈[2,5.5] h2_T_out∈[400,900] p_min=2000kPa t_max∈[800,1100] t_min∈[50,400] p_max∈[8000,12000]"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 1, "n_s_q": 1, "n_lhs": 40, "hx_dT": 10.0,
          "hx_max_group_size": 2,
          "maxiter_inner": 15, "restarts_inner": 3,
          "de_popsize": 12, "de_maxiter": 60, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001,
          "hx_series": True,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 2000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 2.0, "h2_mf_hi": 5.5,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 900.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="1P0T1S_series")


def test_layered_0p0t0s_h2_4() -> None:
    """0P0T0S 理想 串联夹点 H2=4固定 h2_T_out∈[400,900] p_min=2000kPa t_max∈[800,1100] t_min∈[50,400] p_max∈[8000,12000]"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 0, "n_s_q": 0, "n_lhs": 30, "hx_dT": 10.0,
          "hx_series": True,
          "maxiter_inner": 15, "restarts_inner": 3,
          "de_popsize": 10, "de_maxiter": 60, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 2000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 4.0, "h2_mf_hi": 4.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 900.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="0P0T0S_H2-4")


def test_layered_0p1t0s_h2_4() -> None:
    """0T+1P+0S 理想 串联夹点 H2=4固定 h2_T_out∈[400,900] p_min=2000kPa t_max∈[800,1100] t_min∈[50,400] p_max∈[8000,12000]"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 1, "n_s_q": 0, "n_lhs": 40, "hx_dT": 10.0,
          "hx_series": True,
          "maxiter_inner": 20, "restarts_inner": 4,
          "de_popsize": 15, "de_maxiter": 100, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 2000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 4.0, "h2_mf_hi": 4.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 900.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="1P0T0S_H2-4")


def test_layered_2p1s_h2_4() -> None:
    """0T+2P+1S 理想 串联夹点 H2=4固定 h2_T_out∈[400,900] p_min=2000kPa t_max∈[800,1100] t_min∈[50,400] p_max∈[8000,12000]"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 2, "n_s_q": 1, "n_lhs": 40, "hx_dT": 10.0,
          "hx_series": True,
          "maxiter_inner": 20, "restarts_inner": 4,
          "de_popsize": 15, "de_maxiter": 100, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 2000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 4.0, "h2_mf_hi": 4.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 900.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="2P1S_H2-4")


def test_layered_0p1s_h2_4() -> None:
    """0T+0P+1S 理想 串联夹点 H2=4固定 h2_T_out∈[400,1000] p_min=2000kPa t_max∈[800,1100] t_min∈[50,400] p_max∈[8000,12000]"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 0, "n_s_q": 1, "n_lhs": 40, "hx_dT": 10.0,
          "hx_series": True,
          "maxiter_inner": 20, "restarts_inner": 4,
          "de_popsize": 15, "de_maxiter": 100, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 2000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 4.0, "h2_mf_hi": 4.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 1000.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="0P1S_H2-4")


def test_layered_budget_sweep():
    """两层计算量逐步调整: 0P1S TPE外层+CMA内层, 4档对比。"""
    print("=" * 70)
    print("两层计算量扫描 — 0P1S TPE+CMA")
    print("=" * 70)

    combos = [
        ( 10, 1,  50, 10.0, 15),
        ( 20, 2, 100, 15.0, 20),
        ( 40, 3, 150, 20.0, 25),
        ( 60, 6, 200, 20.0, 30),
    ]
    labels = ["轻量 10×1 T50", "标准 20×2 T100", "加强 40×3 T150", "重装 60×6 T200"]

    for label, (mi, ri, tt, sig, nlhs) in zip(labels, combos):
        hp = {**_DEFAULT_HP,
              "n_t_q": 0, "n_p_q": 0, "n_s_q": 1, "n_lhs": nlhs, "hx_dT": 10.0,
              "hx_series": True,
              "outer_method": "tpe",
              "tpe_trials": tt,
              "maxiter_inner": mi, "restarts_inner": ri, "sigma0": sig,
              "qstep": 0.001,
              "t_min_lo": 50.0, "t_min_hi": 400.0,
              "t_max_lo": 800.0, "t_max_hi": 1100.0,
              "p_min_lo": 2000.0, "p_min_hi": 2000.0,
              "p_max_lo": 8000.0, "p_max_hi": 12000.0,
              "h2_mf_lo": 4.0, "h2_mf_hi": 4.0,
              "h2_T_out_lo": 400.0, "h2_T_out_hi": 1240.0,
              "use_non_ideal": False,
              }
        t0 = time.perf_counter()
        _run_layered(hp, tag=f"bd_{label.split()[0]}")
        t1 = time.perf_counter()
        print(f"\n  >>> {label}: {t1 - t0:.1f}s <<<\n")

    print("=" * 70)
    print("  扫描完成，各结果见上面的汇总输出")
    print("=" * 70)


def test_layered_0p1s_h2_4_lbfgsb() -> None:
    """0T+0P+1S L-BFGS-B 内层 串联夹点 H2=4固定 h2_T_out∈[400,1000]"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 0, "n_s_q": 1, "n_lhs": 40, "hx_dT": 10.0,
          "hx_series": True,
          "inner_method": "lbfgsb",
          "lbfgsb_starts": 15, "lbfgsb_maxiter": 40,
          "de_popsize": 15, "de_maxiter": 100, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 2000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 4.0, "h2_mf_hi": 4.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 1000.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="0P1S_LBFGSB")


def test_layered_0p1s_h2_4_lbfgsb() -> None:
    """0T+0P+1S L-BFGS-B 内层 串联夹点 H2=4固定 h2_T_out∈[400,1000]"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 0, "n_s_q": 1, "n_lhs": 40, "hx_dT": 10.0,
          "hx_series": True,
          "inner_method": "lbfgsb",
          "lbfgsb_starts": 15, "lbfgsb_maxiter": 40,
          "de_popsize": 15, "de_maxiter": 100, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 2000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 4.0, "h2_mf_hi": 4.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 1000.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="0P1S_LBFGSB")


def test_layered_0p0t0s_h2_var():
    """0T+0P+0S H2流量可优化∈[2,6] — 模板配置。"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 0, "n_s_q": 0,
          "n_lhs": 30, "hx_dT": 10.0,
          "hx_series": True,
          "outer_method": "de",
          "de_popsize": 20, "de_maxiter": 40, "de_F": 0.8, "de_CR": 0.9,
          "n_workers": 16,
          "maxiter_inner": 40, "restarts_inner": 3, "sigma0": 20.0,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 2000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 2.0, "h2_mf_hi": 6.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 1240.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="0P0T0S_H2var")


def test_layered_0p0t1s_h2_var():
    """0T+0P+1S H2流量可优化∈[2,6] — 模板配置+1个S分位。"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 0, "n_s_q": 1,
          "n_lhs": 30, "hx_dT": 10.0,
          "hx_series": True,
          "outer_method": "de",
          "de_popsize": 20, "de_maxiter": 40, "de_F": 0.8, "de_CR": 0.9,
          "n_workers": 16,
          "maxiter_inner": 40, "restarts_inner": 3, "sigma0": 20.0,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 2000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 2.0, "h2_mf_hi": 6.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 1240.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="0P0T1S_H2var")


def test_layered_0p1s_h2_4_de_large():
    """0T+0P+1S DE大种群外层+CMA内层 — popsize=20, 16进程并行。"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 0, "n_s_q": 1, "n_lhs": 30, "hx_dT": 10.0,
          "hx_series": True,
          "outer_method": "de",
          "de_popsize": 20, "de_maxiter": 40, "de_F": 0.8, "de_CR": 0.9,
          "n_workers": 16,
          "maxiter_inner": 40, "restarts_inner": 3, "sigma0": 20.0,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 2000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 4.0, "h2_mf_hi": 4.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 1240.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="0P1S_DE20x40_16w")


def test_layered_0p1s_h2_4_tpe() -> None:
    """0T+0P+1S TPE外层+CMA内层 H2=4固定 — 16进程并行。"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 0, "n_s_q": 1, "n_lhs": 30, "hx_dT": 10.0,
          "hx_series": True,
          "outer_method": "tpe",
          "tpe_trials": 200,
          "maxiter_inner": 40, "restarts_inner": 3, "sigma0": 20.0,
          "n_workers": 16,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 2000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 4.0, "h2_mf_hi": 4.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 1240.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="0P1S_TPE_16w")


def test_layered_1p1s_h2_4_lbfgsb() -> None:
    """0T+1P+1S L-BFGS-B 内层 串联夹点 H2=4固定 h2_T_out∈[400,1000]"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 1, "n_s_q": 1, "n_lhs": 40, "hx_dT": 10.0,
          "hx_series": True,
          "inner_method": "lbfgsb",
          "lbfgsb_starts": 25, "lbfgsb_maxiter": 40,
          "de_popsize": 15, "de_maxiter": 100, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 2000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 4.0, "h2_mf_hi": 4.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 1000.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="1P1S_LBFGSB")


def test_layered_1p1s_h2_4() -> None:
    """0T+1P+1S 理想 串联夹点 H2=4固定 h2_T_out∈[400,900] p_min=2000kPa t_max∈[800,1100] t_min∈[50,400] p_max∈[8000,12000]"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 1, "n_s_q": 1, "n_lhs": 40, "hx_dT": 10.0,
          "hx_series": True,
          "maxiter_inner": 20, "restarts_inner": 4,
          "de_popsize": 15, "de_maxiter": 100, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 2000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 4.0, "h2_mf_hi": 4.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 900.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="1P1S_H2-4")


def test_inner_global_search_0p1s():
    """0P1S 内层网格穷举: flow∈[0,50]步长1.0, h2_T_out∈[400,1000]步长10K → L-BFGS-B 精搜。"""
    import scipy.optimize as _sopt

    print("=" * 70)
    print("0P1S 内层网格穷举 (flow 步长 1.0, h2_T_out 步长 10K)")
    print("=" * 70)

    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 0, "n_s_q": 1, "hx_dT": 10.0,
          "hx_series": True, "h2_mf_lo": 4.0, "h2_mf_hi": 4.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 1240.0,
          "use_non_ideal": False,
          "maxiter_inner": 5, "restarts_inner": 1, "sigma0": 10.0,  # 极轻量快速
          }
    base = _make_system_input(hp)

    # 用 CMA 最优拓扑: t_max=1074, p_max=11588, s_q=0.610
    tp_in = ClosedCycleTPInput(
        fluid="He", t_min=50.0, t_max=1074.0,
        p_min=2000.0, p_max=11588.0,
        t_quantiles=(), p_quantiles=(),
        s_quantiles=(0.610,),
    )
    layer = ClosedCycleLayer(tp_in)
    n_sc = len(layer.subcycles)
    dim = n_sc + 1
    print(f"  n_sc={n_sc} dim={dim}  t_max=1074 p_max=11588 s_q=0.610")
    h2_mf = hp["h2_mf_lo"]

    t0 = time.perf_counter()

    # ── 网格穷举 ──
    flows_grid = [f for f in range(0, 51)]       # 0,1,...,50
    ht_grid = [t for t in range(400, 1001, 10)]   # 400,410,...,1000
    total_grid = len(flows_grid) ** n_sc * len(ht_grid)
    print(f"  网格: {len(flows_grid)}^{n_sc} × {len(ht_grid)} = {total_grid} 点")
    t_start = time.perf_counter()

    all_r: list[tuple[float, float, float, float]] = []  # (obj, f1, f2, ht)

    for f1 in flows_grid:
        for f2 in flows_grid:
            flows = [float(f1), float(f2)]
            try:
                layer.subcycle_mass_flows = flows
                layer.commit_subcycle_mass_flows_to_topology()
            except Exception:
                continue
            for ht in ht_grid:
                obj = _eval_fast(layer, h2_mf, float(ht), base, hp)
                all_r.append((obj, float(f1), float(f2), float(ht)))

    all_r.sort(key=lambda t: t[0])
    t_grid = time.perf_counter()
    print(f"  网格完成: {len(all_r)} 有效, {t_grid - t_start:.1f}s")
    print(f"  Top-10 网格:")
    for rank in range(min(10, len(all_r))):
        obj, f1, f2, ht = all_r[rank]
        print(f"    {rank + 1:2d}: obj={obj:.5f}  flows=({f1:.0f},{f2:.0f})  "
              f"Σ={f1 + f2:.0f}  h2tout={ht:.0f}K")

    # ── L-BFGS-B 精搜 top-500 ──
    n_refine = min(500, len(all_r))
    print(f"\n  L-BFGS-B 精搜 top-{n_refine} (maxiter=20)...")
    lo = [0.0] * n_sc + [400.0]
    hi = [50.0] * n_sc + [1000.0]
    best_grid_obj, best_f1, best_f2, best_ht = all_r[0]
    best_refined_obj = best_grid_obj
    best_refined_x = [best_f1, best_f2, best_ht]

    improved = 0
    for rank in range(n_refine):
        _, f1, f2, ht = all_r[rank]
        x0 = [f1, f2, ht]
        try:
            res = _sopt.minimize(
                lambda xx: _eval_fast_inner(xx, layer, n_sc, hp, base),
                x0, method="L-BFGS-B",
                bounds=[(lo[i], hi[i]) for i in range(dim)],
                options={"maxiter": 20, "ftol": 1e-12},
            )
        except Exception:
            continue
        if res.fun < best_refined_obj:
            improved += 1
            old = best_refined_obj
            best_refined_obj = res.fun
            best_refined_x = list(res.x)
            print(f"    refine[{rank:3d}]: grid_obj={all_r[rank][0]:.5f} → "
                  f"refined={res.fun:.5f}  flows=({res.x[0]:.1f},{res.x[1]:.1f})  "
                  f"Σ={res.x[0] + res.x[1]:.1f}  ht={res.x[2]:.0f}K  "
                  f"Δ={old - res.fun:.5f}  nfev={res.nfev}")

    t_end = time.perf_counter()

    print("\n" + "=" * 55)
    print("  0P1S 内层网格穷举结果:")
    print(f"    最佳网格 obj     = {best_grid_obj:.5f}")
    print(f"    精搜后 obj       = {best_refined_obj:.5f}")
    print(f"    flows            = ({best_refined_x[0]:.1f}, {best_refined_x[1]:.1f})")
    print(f"    Σflow            = {best_refined_x[0] + best_refined_x[1]:.1f}")
    print(f"    h2_T_out         = {best_refined_x[2]:.0f} K")
    print(f"    精搜改进次数      = {improved}/{n_refine}")
    print(f"    t_grid={t_grid - t_start:.1f}s  refine={t_end - t_grid:.1f}s  总={t_end - t0:.1f}s")
    print("=" * 55)


def _eval_fast_inner(x, layer, n_sc, hp, base):
    flows = list(x[:n_sc])
    ht = x[n_sc]
    try:
        layer.subcycle_mass_flows = flows
        layer.commit_subcycle_mass_flows_to_topology()
    except Exception:
        return 1e9
    return _eval_fast(layer, hp["h2_mf_lo"], ht, base, hp)


def test_inner_h2tout_1240():
    """0P1S 内层 CMA: h2_T_out 上限放至 1240K，验证是否继续改善。"""
    print("=" * 60)
    print("0P1S 内层 CMA — h2_T_out 上限 1240K")
    print("=" * 60)

    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 0, "n_s_q": 1, "hx_dT": 10.0,
          "hx_series": True, "h2_mf_lo": 4.0, "h2_mf_hi": 4.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 1240.0,
          "use_non_ideal": False,
          "maxiter_inner": 60, "restarts_inner": 6, "sigma0": 20.0,
          }
    base = _make_system_input(hp)

    tp_in = ClosedCycleTPInput(
        fluid="He", t_min=50.0, t_max=1074.0,
        p_min=2000.0, p_max=11588.0,
        t_quantiles=(), p_quantiles=(),
        s_quantiles=(0.610,),
    )

    best_overall = float("inf")
    for seed in [42, 123, 456, 789, 1024, 2048]:
        cfg = ClosedCycleTPInput(
            fluid="He", t_min=50.0, t_max=1074.0,
            p_min=2000.0, p_max=11588.0,
            t_quantiles=(), p_quantiles=(),
            s_quantiles=(0.610,),
        )
        result, _ = _inner_cma_fast(cfg, base, hp, seed=seed)
        if result.obj < best_overall:
            best_overall = result.obj
            best_flows = result.flows
            best_h2tout = result.h2_T_out
            best_evals = result.n_evals

    print(f"\n  best obj={best_overall:.5f}  h2tout={best_h2tout:.0f}K  "
          f"flows={[f'{f:.1f}' for f in best_flows]}  Σ={sum(best_flows):.0f}  evals={best_evals}")
    print(f"  对比上限1000K: obj=0.03743 (穷举) / 0.03784 (CMA)")


def test_outer_method_compare():
    """外层算法对比：TPE vs GP-BayesOpt vs CMA-ES 多重启，0P1S 工况。"""
    print("=" * 70)
    print("外层优化算法对比: TPE vs GP vs CMA-ES — 0P1S")
    print("=" * 70)

    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 0, "n_s_q": 1, "hx_dT": 10.0,
          "hx_series": True, "h2_mf_lo": 4.0, "h2_mf_hi": 4.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 1240.0,
          "use_non_ideal": False,
          "maxiter_inner": 60, "restarts_inner": 6, "sigma0": 20.0,
          }
    base = _make_system_input(hp)

    bounds = [
        (50.0, 400.0),     # t_min
        (800.0, 1100.0),   # t_max
        (0.01, 0.99),      # s_q
        (8000.0, 12000.0), # p_max
    ]

    def _outer_obj(x: list[float]) -> float:
        t_min, t_max, s_q, p_max = x
        try:
            tp_in = ClosedCycleTPInput(
                fluid="He", t_min=round(t_min), t_max=round(t_max),
                p_min=2000.0, p_max=round(p_max),
                t_quantiles=(), p_quantiles=(),
                s_quantiles=(round(s_q / 0.001) * 0.001,),
            )
        except Exception:
            return 1.0
        try:
            probe_layer = ClosedCycleLayer(tp_in)
        except Exception:
            return 1.0
        if len(probe_layer.subcycles) == 0:
            return 1.0
        result, _ = _inner_cma_fast(tp_in, base, hp, seed=42)
        return result.obj

    budget = 50

    # ── 1. TPE (optuna) ──
    print("\n  [1/3] TPE (optuna)...")
    t0 = time.perf_counter()
    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)

    def _obj_tpe(trial):
        x = [trial.suggest_float(f"x{i}", lo, hi) for i, (lo, hi) in enumerate(bounds)]
        return _outer_obj(x)

    study = optuna.create_study(sampler=optuna.samplers.TPESampler(seed=42))
    study.optimize(_obj_tpe, n_trials=budget)
    tpe_x = [study.best_params[f"x{i}"] for i in range(4)]
    tpe_val = study.best_value
    tpe_t = time.perf_counter() - t0

    # ── 2. GP-BayesOpt (scikit-optimize) ──
    print("  [2/3] GP-BayesOpt (skopt)...")
    t0 = time.perf_counter()
    import skopt

    def _obj_gp(x):
        return _outer_obj(list(x))

    gp_res = skopt.gp_minimize(
        _obj_gp, bounds, n_calls=budget,
        n_initial_points=30, random_state=42,
        noise=1e-6, n_jobs=1,
    )
    gp_x = gp_res.x
    gp_val = gp_res.fun
    gp_t = time.perf_counter() - t0

    # ── 3. CMA-ES 多重启 ──
    print("  [3/3] CMA-ES 多重启...")
    t0 = time.perf_counter()
    import cma as cma_mod
    import random as _rnd_cma

    dim = 4
    lo = [b[0] for b in bounds]
    hi = [b[1] for b in bounds]
    x0_cma = [(lo[i] + hi[i]) / 2 for i in range(dim)]
    sigma0_cma = 0.3
    total_evals_cma = 0
    best_x_cma = x0_cma[:]
    best_val_cma = float("inf")

    for restart in range(6):
        p = max(8, int(4 + 3 * __import__('math').log2(dim)))
        opts = {
            "bounds": [lo, hi],
            "CMA_stds": [sigma0_cma * (hi[i] - lo[i]) for i in range(dim)],
            "popsize": p,
            "verbose": -9,
            "seed": 42 + restart,
        }
        x0 = x0_cma if restart == 0 else [
            _rnd_cma.Random(42 + restart).uniform(lo[i], hi[i]) for i in range(dim)]
        es = cma_mod.CMAEvolutionStrategy(x0, sigma0_cma, opts)
        gen = 0
        while total_evals_cma < budget and gen < 100 and not es.stop():
            solutions = es.ask()
            scores = [_outer_obj(list(s)) for s in solutions]
            es.tell(solutions, scores)
            total_evals_cma += len(solutions)
            gen += 1
            for s, sc in zip(solutions, scores):
                if sc < best_val_cma:
                    best_val_cma = sc
                    best_x_cma = list(s)
        if total_evals_cma >= budget:
            break

    cma_val = best_val_cma
    cma_x = best_x_cma
    cma_t = time.perf_counter() - t0

    # ── 汇总 ──
    print("\n" + "=" * 70)
    print(f"  外层算法对比结果 (budget={budget} eval, 0P1S):")
    print(f"  {'方法':<20s} {'obj':>8s}  {'t':>6s}  "
          f"{'t_min':>6s} {'t_max':>6s} {'s_q':>6s} {'p_max':>7s}")
    print("  " + "-" * 65)
    for label, val, xt, t_elapsed in [
        ("TPE (optuna)", tpe_val, tpe_x, tpe_t),
        ("GP-BayesOpt (skopt)", gp_val, gp_x, gp_t),
        ("CMA-ES 多重启(6)", cma_val, cma_x, cma_t),
    ]:
        print(f"  {label:<20s} {val:>8.5f}  {t_elapsed:>5.0f}s  "
              f"{xt[0]:>6.0f} {xt[1]:>6.0f} {xt[2]:>6.3f} {xt[3]:>7.0f}")
    print("=" * 70)
    """0P1S 内层 CMA: h2_T_out 上限放至 1240K，验证是否继续改善。"""
    print("=" * 60)
    print("0P1S 内层 CMA — h2_T_out 上限 1240K")
    print("=" * 60)

    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 0, "n_s_q": 1, "hx_dT": 10.0,
          "hx_series": True, "h2_mf_lo": 4.0, "h2_mf_hi": 4.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 1240.0,
          "use_non_ideal": False,
          "maxiter_inner": 60, "restarts_inner": 6, "sigma0": 20.0,
          }
    base = _make_system_input(hp)

    # 用 CMA 最优拓扑
    tp_in = ClosedCycleTPInput(
        fluid="He", t_min=50.0, t_max=1074.0,
        p_min=2000.0, p_max=11588.0,
        t_quantiles=(), p_quantiles=(),
        s_quantiles=(0.610,),
    )

    best_overall = float("inf")
    for seed in [42, 123, 456, 789, 1024, 2048]:
        cfg = ClosedCycleTPInput(
            fluid="He", t_min=50.0, t_max=1074.0,
            p_min=2000.0, p_max=11588.0,
            t_quantiles=(), p_quantiles=(),
            s_quantiles=(0.610,),
        )
        result, _ = _inner_cma_fast(cfg, base, hp, seed=seed)
        if result.obj < best_overall:
            best_overall = result.obj
            best_flows = result.flows
            best_h2tout = result.h2_T_out
            best_evals = result.n_evals

    print(f"\n  best obj={best_overall:.5f}  h2tout={best_h2tout:.0f}K  "
          f"flows={[f'{f:.1f}' for f in best_flows]}  Σ={sum(best_flows):.0f}  evals={best_evals}")
    print(f"  对比上限1000K: obj=0.03743 (穷举) / 0.03784 (CMA)")


def test_inner_method_compare() -> None:
    """CMA vs L-BFGS-B vs Hybrid 对比：固定 1P1S + 0P1S 拓扑。"""
    print("=" * 70)
    print("内层优化方法对比: CMA-ES vs L-BFGS-B vs Hybrid")
    print("=" * 70)

    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 1, "n_s_q": 1, "hx_dT": 10.0,
          "hx_series": True, "h2_mf_lo": 4.0, "h2_mf_hi": 4.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 1000.0,
          "use_non_ideal": False,
          "maxiter_inner": 20, "restarts_inner": 3, "sigma0": 15.0,
          "lbfgsb_starts": 30, "lbfgsb_maxiter": 50,
          "hybrid_cma_iter": 10, "hybrid_cma_rst": 2,
          "hybrid_n_keep": 5, "hybrid_lbfgsb_iter": 50, "hybrid_sigma0": 15.0,
          }

    base = _make_system_input(hp)

    # 1P1S 拓扑
    print("\n  1P1S [n_sc~6, dim=7]:")
    for method, label in [
        ("cma", "CMA 20x3"),
        ("lbfgsb", "L-BFGS-B 30x50"),
        ("hybrid", "Hybrid CMA10x2+L3x50"),
    ]:
        hp_m = {**hp, "inner_method": method}
        best_obj = float("inf")
        best_flows = []
        best_h2tout = 0.0
        best_evals = 0
        best_t = 0.0
        for seed in [42, 123, 456]:
            cfg = ClosedCycleTPInput(
                fluid="He", t_min=51.0, t_max=1041.0,
                p_min=2000.0, p_max=10922.0,
                t_quantiles=(), p_quantiles=(0.990,),
                s_quantiles=(0.585,),
            )
            t_s = time.perf_counter()
            if method == "lbfgsb":
                result, _ = _inner_lbfgsb_fast(cfg, base, hp_m, seed=seed)
            elif method == "hybrid":
                result, _ = _inner_hybrid_fast(cfg, base, hp_m, seed=seed)
            else:
                result, _ = _inner_cma_fast(cfg, base, hp_m, seed=seed)
            t_e = time.perf_counter()
            if result.obj < best_obj:
                best_obj = result.obj
                best_flows = result.flows
                best_h2tout = result.h2_T_out
                best_evals = result.n_evals
                best_t = t_e - t_s
        print(f"    {label:20s}  best={best_obj:.5f}  h2tout={best_h2tout:.0f}K  "
              f"Σflow={sum(best_flows):.0f}  evals={best_evals}  t={best_t:.1f}s")

    # 0P1S 拓扑
    print("\n  0P1S [n_sc~2, dim=3]:")
    for method, label in [
        ("cma", "CMA 20x3"),
        ("lbfgsb", "L-BFGS-B 30x50"),
        ("hybrid", "Hybrid CMA10x2+L3x50"),
    ]:
        hp_m = {**hp, "inner_method": method}
        best_obj = float("inf")
        best_flows = []
        best_h2tout = 0.0
        best_evals = 0
        best_t = 0.0
        for seed in [42, 123, 456]:
            cfg = ClosedCycleTPInput(
                fluid="He", t_min=50.0, t_max=1100.0,
                p_min=2000.0, p_max=8000.0,
                t_quantiles=(), p_quantiles=(),
                s_quantiles=(0.674,),
            )
            t_s = time.perf_counter()
            if method == "lbfgsb":
                result, _ = _inner_lbfgsb_fast(cfg, base, hp_m, seed=seed)
            elif method == "hybrid":
                result, _ = _inner_hybrid_fast(cfg, base, hp_m, seed=seed)
            else:
                result, _ = _inner_cma_fast(cfg, base, hp_m, seed=seed)
            t_e = time.perf_counter()
            if result.obj < best_obj:
                best_obj = result.obj
                best_flows = result.flows
                best_h2tout = result.h2_T_out
                best_evals = result.n_evals
                best_t = t_e - t_s
        print(f"    {label:20s}  best={best_obj:.5f}  h2tout={best_h2tout:.0f}K  "
              f"Σflow={sum(best_flows):.0f}  evals={best_evals}  t={best_t:.1f}s")

    print("\n  0P1S [n_sc~2, dim=3]:")
    for method, label in [("cma", "CMA-ES 20×3"), ("lbfgsb", "L-BFGS-B 30×50")]:
        hp_m = {**hp, "inner_method": method}
        best_obj = float("inf")
        best_flows = []
        best_h2tout = 0.0
        best_evals = 0
        best_t = 0.0
        for seed in [42, 123, 456]:
            cfg = ClosedCycleTPInput(
                fluid="He", t_min=50.0, t_max=1100.0,
                p_min=2000.0, p_max=8000.0,
                t_quantiles=(), p_quantiles=(),
                s_quantiles=(0.674,),
            )
            t_s = time.perf_counter()
            if method == "lbfgsb":
                result, _ = _inner_lbfgsb_fast(cfg, base, hp_m, seed=seed)
            else:
                result, _ = _inner_cma_fast(cfg, base, hp_m, seed=seed)
            t_e = time.perf_counter()
            if result.obj < best_obj:
                best_obj = result.obj
                best_flows = result.flows
                best_h2tout = result.h2_T_out
                best_evals = result.n_evals
                best_t = t_e - t_s
        print(f"    {label:16s}  best={best_obj:.5f}  h2tout={best_h2tout:.0f}K  "
              f"Σflow={sum(best_flows):.0f}  evals={best_evals}  t={best_t:.1f}s")
    """内层 CMA 参数扫描：固定 1P1S 拓扑，测试不同 maxiter/restarts/sigma0 组合。"""
    print("=" * 70)
    print("内层 CMA 参数扫描 — 1P1S 固定拓扑")
    print("=" * 70)

    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 1, "n_s_q": 1, "hx_dT": 10.0,
          "hx_series": True, "h2_mf_lo": 4.0, "h2_mf_hi": 4.0,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 1000.0,
          "use_non_ideal": False,
          }

    tp_in = ClosedCycleTPInput(
        fluid="He", t_min=51.0, t_max=1041.0,
        p_min=2000.0, p_max=10922.0,
        t_quantiles=(), p_quantiles=(0.990,),
        s_quantiles=(0.585,),
    )

    # 预建拓扑（仅供流量信息，不参与 CMA 内数据传递）
    base_layer = ClosedCycleLayer(tp_in)
    n_sc = len(base_layer.subcycles)
    print(f"  n_sc={n_sc}")

    base = _make_system_input(hp)

    configs = [
        # (label, maxiter, restarts, sigma0)
        ("A: 20×3 σ15", 20, 3, 15.0),
        ("B: 20×5 σ15", 20, 5, 15.0),
        ("C: 40×3 σ15", 40, 3, 15.0),
        ("D: 40×5 σ15", 40, 5, 15.0),
        ("E: 60×5 σ15", 60, 5, 15.0),
        ("F: 40×5 σ5",  40, 5, 5.0),
        ("G: 40×5 σ30", 40, 5, 30.0),
        ("H: 100×5 σ15", 100, 5, 15.0),
    ]

    results: list[tuple] = []

    for label, maxiter, restarts, sigma0_val in configs:
        hp_cfg = {**hp, "maxiter_inner": maxiter,
                  "restarts_inner": restarts, "sigma0": sigma0_val}
        best_obj = float("inf")
        best_flows = []
        best_h2tout = 0.0
        total_t = 0.0

        for seed in [42, 123, 456]:
            cfg_tp = ClosedCycleTPInput(
                fluid="He", t_min=51.0, t_max=1041.0,
                p_min=2000.0, p_max=10922.0,
                t_quantiles=(), p_quantiles=(0.990,),
                s_quantiles=(0.585,),
            )
            result, _ = _inner_cma_fast(cfg_tp, base, hp_cfg, seed=seed)
            if result.obj < best_obj:
                best_obj = result.obj
                best_flows = result.flows
                best_h2tout = result.h2_T_out
            total_t += result.runtime

        avg_t = total_t / 3.0
        results.append((label, maxiter, restarts, sigma0_val,
                        best_obj, best_h2tout, sum(best_flows), avg_t))
        sys.stdout.write(
            f"  {label:14s}  best={best_obj:.5f}  h2tout={best_h2tout:.0f}K  "
            f"Σflow={sum(best_flows):.0f}  avg_t={avg_t:.1f}s\n"
        )
        sys.stdout.flush()

    print("\n" + "-" * 70)
    print(f"{'配置':<14s} {'iter':>5s} {'rst':>4s} {'σ':>5s} "
          f"{'obj':>8s} {'h2tout':>7s} {'Σflow':>6s} {'avg_t':>7s}")
    print("-" * 70)
    for row in results:
        label, mi, rst, sig, obj, ht, sf, at = row
        print(f"{label:<14s} {mi:>5d} {rst:>4d} {sig:>5.0f} "
              f"{obj:>8.5f} {ht:>7.0f} {sf:>6.0f} {at:>7.1f}")
    print("-" * 70)


_de_best_x: list[float] = []
_tpe_best_x: list[float] = []
_cma_multi_best_x: list[float] = []
_tpe_state: dict = {}  # 供 _optuna_objective 使用的全局状态


def _optuna_objective(trial) -> float:
    """模块级目标函数（供 optuna 多进程并行）。通过 _tpe_state 获取参数。"""
    st = _tpe_state
    dim = st["dim"]
    x = [trial.suggest_float(f"x{i}", st["lo"][i], st["hi"][i]) for i in range(dim)]
    args = (x, st["base"], st["hp"], 0)
    return _de_trial_worker(args)


def _run_layered_de_phase(hp, bounds, lhs_objs, base, ex, conv_history, seed):
    """外层 DE 阶段（原有逻辑，提取为独立函数）。"""
    import random as _rnd_
    rng = _rnd_.Random(seed)
    lo = [b[0] for b in bounds]
    hi = [b[1] for b in bounds]
    dim = len(bounds)
    popsize = hp["de_popsize"]

    pop: list[list[float]] = []
    for r in lhs_objs[:popsize]:
        x = [r.t_min, r.t_max]
        for q in r.t_q: x.append(q)
        for q in r.p_q: x.append(q)
        for q in getattr(r, 's_q', ()): x.append(q)
        x.append(r.p_max); x.append(r.p_min)
        pop.append(x)
    while len(pop) < popsize:
        pop.append([rng.uniform(lo[d], hi[d]) for d in range(dim)])

    g0_args = [(p, base, hp, 0) for p in pop]
    g0_futures = [ex.submit(_de_trial_worker, a) for a in g0_args]
    scores = [f.result(timeout=600) for f in g0_futures]

    best_idx = min(range(len(scores)), key=lambda i: scores[i])
    global _de_best_x
    _de_best_x = pop[best_idx][:]
    best_val = scores[best_idx]
    stall = 0
    n_evals = len(scores)
    conv_history.append(best_val)
    print(f"  DE gen 0: best={best_val:.5f}")

    for gen in range(hp["de_maxiter"]):
        trials: list[tuple[int, list[float]]] = []
        for i in range(popsize):
            pool = [j for j in range(popsize) if j != i]
            a, b, c = rng.sample(pool, 3)
            v = [pop[a][d] + hp["de_F"] * (pop[b][d] - pop[c][d]) for d in range(dim)]
            for d in range(dim): v[d] = max(lo[d], min(hi[d], v[d]))
            jr = rng.randrange(dim)
            u = [v[d] if rng.random() < hp["de_CR"] or d == jr else pop[i][d] for d in range(dim)]
            trials.append((i, u))

        trial_args = [(t, base, hp, 0) for _, t in trials]
        t_futures = [ex.submit(_de_trial_worker, a) for a in trial_args]
        trial_scores = [f.result(timeout=600) for f in t_futures]
        n_evals += len(trial_scores)

        improved = False
        for (i, u), us in zip(trials, trial_scores):
            if us < scores[i]:
                pop[i] = u; scores[i] = us
                if us < best_val:
                    best_val = us; _de_best_x = u[:]; improved = True; stall = 0
        if not improved:
            stall += popsize
        conv_history.append(best_val)

        sys.stdout.write(f"  DE gen {gen+1:2d}: best={best_val:.5f} "
                         f"stall={stall}/{popsize*hp['de_maxiter']} "
                         f"t_max={_de_best_x[1]:.0f} t_min={_de_best_x[0]:.0f} "
                         f"pmax={_de_best_x[-2]:.0f} pmin={_de_best_x[-1]:.0f}\n")
        sys.stdout.flush()
        if stall >= popsize * hp["de_maxiter"] * 2:
            break


def _run_layered_tpe_phase(hp, bounds, lhs_objs, base, ex, conv_history, seed):
    """外层 TPE 阶段（optuna，多进程并行）。"""
    import optuna
    import tempfile, os
    optuna.logging.set_verbosity(optuna.logging.WARNING)

    dim = len(bounds)
    lo = [b[0] for b in bounds]
    hi = [b[1] for b in bounds]
    n_trials = hp.get("tpe_trials", 200)
    n_workers = hp["n_workers"]

    # 设置全局状态供 worker 进程使用
    global _tpe_state
    _tpe_state = {"dim": dim, "lo": lo, "hi": hi, "base": base, "hp": hp}

    # 文件存储以支持多进程
    tmpdir = tempfile.gettempdir()
    storage_path = os.path.join(tmpdir, f"optuna_tpe_{os.getpid()}.log")
    storage = optuna.storages.JournalStorage(
        optuna.storages.JournalFileStorage(storage_path)
    )

    study = optuna.create_study(
        storage=storage,
        sampler=optuna.samplers.TPESampler(seed=seed, n_startup_trials=min(20, n_trials)),
        direction="minimize",
    )

    # 喂入 LHS 结果
    for r in lhs_objs[:min(30, len(lhs_objs))]:
        x = [r.t_min, r.t_max]
        for q in r.t_q: x.append(q)
        for q in r.p_q: x.append(q)
        for q in getattr(r, 's_q', ()): x.append(q)
        x.append(r.p_max); x.append(r.p_min)
        trial = optuna.trial.create_trial(
            params={f"x{i}": x[i] for i in range(dim)},
            distributions={f"x{i}": optuna.distributions.FloatDistribution(lo[i], hi[i])
                           for i in range(dim)},
            values=[r.obj],
        )
        study.add_trial(trial)

    print(f"\n  [阶段2] TPE ({n_trials} trials + {len(lhs_objs[:30])} LHS seeds, "
          f"{n_workers} 进程并行)...")

    if lhs_objs:
        conv_history.append(lhs_objs[0].obj)

    def _callback(study, trial):
        if study.best_value < conv_history[-1]:
            conv_history.append(study.best_value)
            bp = study.best_params
            global _tpe_best_x
            _tpe_best_x = [bp[f"x{i}"] for i in range(dim)]
            sys.stdout.write(
                f"  TPE trial {trial.number:3d}: best={study.best_value:.5f} "
                f"t_max={_tpe_best_x[1]:.0f} t_min={_tpe_best_x[0]:.0f} "
                f"s_q={_tpe_best_x[2]:.3f} pmax={_tpe_best_x[3]:.0f}\n")
            sys.stdout.flush()

    study.optimize(_optuna_objective, n_trials=n_trials, n_jobs=n_workers,
                   callbacks=[_callback], show_progress_bar=False)

    # 确保 _tpe_best_x 已设置
    global _tpe_best_x
    if not _tpe_best_x:
        bp = study.best_params
        _tpe_best_x = [bp[f"x{i}"] for i in range(dim)]

    # 清理临时文件
    try:
        os.unlink(storage_path)
    except Exception:
        pass


def _cma_multi_worker(args):
    """模块级 worker: 一个完整 CMA 搜索。bounds 含 lo==hi 的固定维度需填充。"""
    seed_i, full_lo, full_hi, base, hp = args
    import cma as cma_mod
    import random as _rnd
    import math
    # 只对 lo < hi 的维度做 CMA
    active_idx = [i for i in range(len(full_lo)) if full_hi[i] - full_lo[i] > 1e-6]
    alo = [full_lo[i] for i in active_idx]
    ahi = [full_hi[i] for i in active_idx]
    adim = len(alo)
    rng = _rnd.Random(seed_i)
    sigma0 = 0.3
    popsize = max(10, int(4 + 3 * math.log2(max(adim, 2))))
    x0 = [(alo[i] + ahi[i]) / 2 for i in range(adim)]
    es = cma_mod.CMAEvolutionStrategy(x0, sigma0, {
        "bounds": [alo, ahi], "popsize": popsize, "verbose": -9, "seed": seed_i,
    })
    best_full_x = x0[:]
    best_val = float("inf")
    stall = 0
    gen = 0
    while gen < 50 and not es.stop() and stall < popsize * 3:
        solutions = es.ask()
        scores = []
        for s in solutions:
            # 填充回完整维度
            full_x = list(full_lo)
            for ai, si in enumerate(active_idx):
                full_x[si] = s[ai]
            scores.append(_de_trial_worker((full_x, base, hp, 0)))
        es.tell(solutions, scores)
        gen += 1
        improved = False
        for s, sc in zip(solutions, scores):
            if sc < best_val:
                best_val = sc
                # 重建完整 x
                bx = list(full_lo)
                for ai, si in enumerate(active_idx):
                    bx[si] = s[ai]
                best_full_x = bx
                improved = True; stall = 0
        if not improved:
            stall += len(solutions)
    return (best_val, best_full_x, gen, seed_i)


def _run_layered_cma_multi_phase(hp, bounds, lhs_objs, base, ex, conv_history, seed):
    # 过滤掉 lo==hi 的固定维度（如 p_min=2000 固定）
    active_bounds = [(lo, hi) for lo, hi in bounds if hi - lo > 1e-6]
    dim = len(active_bounds)
    lo = [b[0] for b in active_bounds]
    hi = [b[1] for b in active_bounds]
    n_seeds = hp.get("cma_multi_seeds", 16)

    sys.stdout.write(f"\n  [阶段2] CMA多重启 ({n_seeds} 独立种子)...\n")
    sys.stdout.flush()

    full_lo = [b[0] for b in bounds]
    full_hi = [b[1] for b in bounds]
    worker_args = [(seed + i, list(full_lo), list(full_hi), base, hp) for i in range(n_seeds)]
    futures = [ex.submit(_cma_multi_worker, a) for a in worker_args]

    best_val = float("inf")
    best_x = None
    gen_sum = 0
    for f in futures:
        val, x, g, si = f.result(timeout=600)
        gen_sum += g
        sys.stdout.write(f"    seed {si:2d}: gen={g} obj={val:.5f}\n")
        sys.stdout.flush()
        if val < best_val:
            best_val = val
            best_x = x
            conv_history.append(val)

    global _cma_multi_best_x
    global _cma_multi_best_x
    _cma_multi_best_x = best_x
    sys.stdout.write(f"    avg gen={gen_sum / n_seeds:.0f}, best={best_val:.5f}\n")
    sys.stdout.flush()


def _run_layered(hp: dict, tag: str = "") -> None:
    hp = {**_DEFAULT_HP, **hp}
    obj_mode_label = hp.get("obj_mode", "series" if hp.get("hx_series") else "group")
    tq_tag = f"{hp['n_t_q']}T{hp['n_p_q']}P"
    if hp.get("n_s_q", 0) > 0:
        tq_tag += f"{hp['n_s_q']}S"
    print("=" * 60)
    olabel = ("LHS+TPE" if hp.get("outer_method") == "tpe"
              else ("LHS+CMA×N" if hp.get("outer_method") == "cma_multi"
                    else "LHS+DE并行"))
    print(f"分层优化 ({olabel}) [{tag}] {tq_tag}")
    outer_label = (f"TPE={hp.get('tpe_trials', '?')}trials"
                   if hp.get("outer_method") == "tpe"
                   else (f"CMA×{hp.get('cma_multi_seeds',16)}"
                         if hp.get("outer_method") == "cma_multi"
                         else f"DE={hp['de_popsize']}x{hp['de_maxiter']}"))
    print(f"  LHS={hp['n_lhs']}点 {outer_label} "
          f"t_min∈[{hp['t_min_lo']},{hp['t_min_hi']}] "
          f"t_max∈[{hp['t_max_lo']},{hp['t_max_hi']}] "
          f"p_min∈[{hp['p_min_lo']},{hp['p_min_hi']}] "
          f"p_max∈[{hp['p_max_lo']},{hp['p_max_hi']}] "
          f"H2∈[{hp['h2_mf_lo']},{hp['h2_mf_hi']}] "
           f"interp_He={hp.get('use_interp_he', False)} non_ideal={hp['use_non_ideal']} "
           f"inner={hp.get('inner_method','cma')} hx_dT={hp['hx_dT']} "
           f"obj={obj_mode_label}")
    print("=" * 60)

    base = _make_system_input(hp)

    outer_bounds = [
        (hp["t_min_lo"], hp["t_min_hi"]),
        (hp["t_max_lo"], hp["t_max_hi"]),
    ]
    for _ in range(hp["n_t_q"]):
        outer_bounds.append((0.01, 0.99))
    for _ in range(hp["n_p_q"]):
        outer_bounds.append((0.01, 0.99))
    for _ in range(hp["n_s_q"]):
        outer_bounds.append((0.01, 0.99))
    outer_bounds.append((hp["p_max_lo"], hp["p_max_hi"]))
    outer_bounds.append((hp["p_min_lo"], hp["p_min_hi"]))
    dim = len(outer_bounds)

    import concurrent.futures, multiprocessing, random as _rnd_
    ctx = multiprocessing.get_context("spawn")
    n_workers = hp["n_workers"]

    t0_all = time.perf_counter()
    best_inner: LayerResult | None = None
    conv_history: list[float] = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as ex:
        # ── 阶段1: LHS 并行采样 ──
        print(f"\n  [阶段1] LHS {hp['n_lhs']} 点并行采样...")
        seed = hp.get("seed", 42)
        # ── 熵预检：跳过 S(t_max,p_max) < S(t_min,p_min) 的无效样本 ──
        from core import PropertyRegistry
        _pre_props = PropertyRegistry()
        _pre_fluid = "He"
        tasks = []
        _lhs_seed = seed
        _target_tasks = max(hp['n_lhs'], hp.get('de_popsize', 15))
        while len(tasks) < _target_tasks:
            lhs_samples = _lhs(hp["n_lhs"], outer_bounds, seed=_lhs_seed)
            _lhs_seed += 1
            for i, row in enumerate(lhs_samples):
                if len(tasks) >= _target_tasks:
                    break
                t_min = round(row[0]); t_max = round(row[1])
                idx = 2
                t_q_vals = tuple(float(row[idx + j]) for j in range(hp["n_t_q"]))
                idx += hp["n_t_q"]
                p_q_vals = tuple(float(row[idx + j]) for j in range(hp["n_p_q"]))
                idx += hp["n_p_q"]
                s_q_vals = tuple(float(row[idx + j]) for j in range(hp["n_s_q"]))
                idx += hp["n_s_q"]
                p_max = round(row[idx]); p_min = round(row[idx + 1])
                try:
                    s_top = _pre_props(_pre_fluid, "TP", float(t_max), float(p_max))["S"]
                    s_bot = _pre_props(_pre_fluid, "TP", float(t_min), float(p_min))["S"]
                    if s_top < s_bot:
                        continue
                except Exception:
                    continue
                tasks.append((t_min, t_max, t_q_vals, p_q_vals, s_q_vals, p_max, p_min,
                              base, hp, seed + len(tasks)))

            lhs_results: list[LayerResult] = []
            futures = [ex.submit(_sample_worker, t) for t in tasks]
            for fut in futures:
                r = fut.result(timeout=600)
                if r.n_subcycles > 0:
                    lhs_results.append(r)

        lhs_objs = sorted(lhs_results, key=lambda r: r.obj)
        if not lhs_objs:
            print("  LHS: 0 有效 — 拓扑无子循环，终止")
            assert False, "所有 LHS 样本均无有效拓扑"
        print(f"  LHS: {len(lhs_results)} 有效, best={lhs_objs[0].obj:.5f}")

        # ── 阶段2: 外层搜索 (DE / TPE / CMA多重启) ──
        outer_method = hp.get("outer_method", "de")
        if outer_method == "tpe":
            _run_layered_tpe_phase(hp, outer_bounds, lhs_objs, base, ex,
                                   conv_history, seed)
        elif outer_method == "cma_multi":
            _run_layered_cma_multi_phase(hp, outer_bounds, lhs_objs, base, ex,
                                          conv_history, seed)
        else:
            _run_layered_de_phase(hp, outer_bounds, lhs_objs, base, ex,
                                  conv_history, seed)

    t_total = time.perf_counter() - t0_all

    best_idx = min(range(len(conv_history)), key=lambda i: conv_history[i])
    best_val = conv_history[best_idx]
    best_x = None
    if outer_method == "tpe":
        best_x = _tpe_best_x
    elif outer_method == "cma_multi":
        best_x = _cma_multi_best_x
    else:
        best_x = _de_best_x

    # ── 汇总 ──
    print("\n" + "=" * 55)
    outer_label = (f"LHS{len(lhs_results)}+TPE{hp.get('tpe_trials', 200)}"
                   if outer_method == "tpe"
                   else (f"LHS{len(lhs_results)}+CMA×{hp.get('cma_multi_seeds',16)}"
                         if outer_method == "cma_multi"
                         else f"LHS{len(lhs_results)}+DE{hp.get('de_popsize',15)}x{hp.get('de_maxiter',100)}gen"))
    print(f"优化完成: {outer_label} "
          f"= {len(conv_history)} 次改进, {t_total:.1f}s")
    idx = 2
    tq_out = tuple(round(best_x[idx + j], 4) for j in range(hp["n_t_q"]))
    idx += hp["n_t_q"]
    pq_out = tuple(round(best_x[idx + j], 4) for j in range(hp["n_p_q"]))
    idx += hp["n_p_q"]
    sq_out = tuple(round(best_x[idx + j], 4) for j in range(hp["n_s_q"]))
    print(f"  t_max={best_x[1]:.0f}K  t_min={best_x[0]:.0f}K  "
          f"p_max={best_x[-2]:.0f}kPa  p_min={best_x[-1]:.0f}kPa")
    print(f"  t_q={[f'{v:.3f}' for v in tq_out]}  p_q={[f'{v:.3f}' for v in pq_out]}")
    if hp["n_s_q"] > 0:
        print(f"  s_q={[f'{v:.3f}' for v in sq_out]}")
    print(f"  obj={best_val:.5f}")

    # ── 最终重建（唯一一次，绘图复用）──
    print("\n  [最终重建]")
    try:
        opt_result, full_hx, full_ideal_rep, full_ni_rep, full_layer, full_hots, full_colds = \
            _rebuild_result(best_x, base, hp, use_coolprop=False)
        best_inner = opt_result
    except Exception:
        best_inner = LayerResult(n=0, t_min=0, t_max=0, p_min=0, p_max=0,
                                  t_q=(), p_q=(), flows=[], h2_mf=0,
                                  obj=best_val, n_subcycles=0, n_evals=0, runtime=0)
        full_hx = full_ideal_rep = full_ni_rep = full_layer = full_hots = full_colds = None

    print(f"    flows={[f'{v:.1f}' for v in best_inner.flows]} "
          f"H2={best_inner.h2_mf:.2f}kg/s  n_sc={best_inner.n_subcycles}")

    # ── CoolProp 验证（仅 use_interp_he=True 时需要）──
    if hp.get("use_interp_he"):
        print("\n  [CoolProp 验证]")
        try:
            cp_result, _, cp_ir, cp_nr, cp_layer, _, _ = \
                _rebuild_result(best_x, base, hp, use_coolprop=True)
            if cp_nr.nodes:
                cp_p_min = min(snap.P for _, snap in cp_nr.nodes)
                cp_p_max = max(snap.P for _, snap in cp_nr.nodes)
                print(f"    interp He: obj={best_val:.5f} n_sc={best_inner.n_subcycles} "
                      f"H2={best_inner.h2_mf:.2f}")
                print(f"    CoolProp:  n_sc={len(cp_layer.subcycles)} "
                      f"P∈[{cp_p_min:.0f},{cp_p_max:.0f}]kPa "
                      f"T∈[{min(snap.T for _, snap in cp_nr.nodes):.0f},"
                      f"{max(snap.T for _, snap in cp_nr.nodes):.0f}]K")
            else:
                print("    (CoolProp 与 InterpHe 拓扑不一致，跳过详细对比)")
        except Exception as e:
            print(f"    CoolProp 验证失败: {e}")

    # ── 绘图 ──
    run_dir = Path(hp.get("out_dir", TESTS_DIR / "run_layered_fast"))
    run_dir.mkdir(parents=True, exist_ok=True)

    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig1, ax1 = plt.subplots(figsize=(8, 4))
    ax1.plot(range(len(conv_history)), conv_history, "o-", ms=4, color="tab:blue")
    ax1.axhline(y=best_val, color="tab:red", lw=1, ls="--", label=f"best={best_val:.5f}")
    if obj_mode_label == "eff_pinch":
        ylabel = "obj = -η"
    elif obj_mode_label == "pinch":
        ylabel = "obj (utility_demand/q_source)"
    else:
        ylabel = "obj (unmatched/total_Q)"
    ax1.set_xlabel("DE generation"); ax1.set_ylabel(ylabel)
    ax1.set_title(f"{'TPE' if outer_method == 'tpe' else 'DE'} convergence [{tag}] {tq_tag}")
    ax1.legend(); ax1.grid(True, alpha=0.2)
    fig1.tight_layout(); fig1.savefig(run_dir / f"conv_{tag}.png", dpi=150); plt.close(fig1)
    print(f"  收敛图: {run_dir / f'conv_{tag}.png'}")

    try:
        # 复用最终重建的结果（不再重建）
        best_hx, ideal_rep, ni_rep, best_layer, _hots, _colds = \
            full_hx, full_ideal_rep, full_ni_rep, full_layer, full_hots, full_colds

        # 理想 TS/PS
        _draw_cycle_ts_ps(ideal_rep, best_layer,
                          f"{tag} ideal", run_dir / f"ideal_ts_ps_{tag}.png")
        print(f"  理想 TS/PS: {run_dir / f'ideal_ts_ps_{tag}.png'}")

        # 非理想 TS/PS
        if hp["use_non_ideal"] and best_layer.non_ideal is not None:
            _draw_cycle_ts_ps(ni_rep, best_layer,
                              f"{tag} non-ideal", run_dir / f"nonideal_ts_ps_{tag}.png",
                              ni_nodes=best_layer.non_ideal.nodes)
            print(f"  非理想 TS/PS: {run_dir / f'nonideal_ts_ps_{tag}.png'}")
        else:
            print(f"  非理想 TS/PS: 跳过 (use_non_ideal=False)")

        obj_mode = hp.get("obj_mode", "series" if hp.get("hx_series") else "group")
        if obj_mode == "pinch":
            _draw_pinch_tq(_hots, _colds, best_hx, hp["hx_dT"],
                           run_dir / f"pinch_tq_{tag}.png")
            print(f"  夹点 T-Q: {run_dir / f'pinch_tq_{tag}.png'}")
            print(f"  热公用={best_hx.hot_utility_demand:.0f}kW  "
                  f"冷公用={best_hx.cold_utility_demand:.0f}kW")
        elif obj_mode == "eff_pinch":
            _draw_pinch_tq(_hots, _colds, best_hx, hp["hx_dT"],
                           run_dir / f"pinch_tq_{tag}.png")
            print(f"  夹点 T-Q: {run_dir / f'pinch_tq_{tag}.png'}")
            print(f"  热公用={best_hx.hot_utility_demand:.0f}kW  "
                  f"冷公用={best_hx.cold_utility_demand:.0f}kW  "
                  f"η={-best_val:.4f}")
        else:
            _draw_hx_tq(_hots, _colds, best_hx, run_dir / f"hx_tq_{tag}.png")
            print(f"  HX T-Q: {run_dir / f'hx_tq_{tag}.png'}")
            print(f"  HX单元={len(best_hx.units)} "
                  f"未匹配={best_hx.total_unmatched:.0f}kW "
                  f"总匹配={best_hx.total_matched:.0f}kW")
    except Exception as e:
        import traceback
        print(f"  绘图异常: {e}")
        traceback.print_exc()

    # ── 持久化结果 ──
    if hp.get("save_csv", True):
        _save_optimization_result(
            run_dir=run_dir,
            tag=tag,
            obj_mode=obj_mode_label,
            best_val=best_val,
            best_x=best_x,
            tq_out=tq_out,
            pq_out=pq_out,
            sq_out=sq_out,
            best_inner=best_inner,
            best_hx=best_hx,
            t_total=t_total,
            ni_rep=ni_rep,
            hots=_hots,
            colds=_colds,
        )

    assert best_val < 1e3


def _save_optimization_result(
    run_dir: Path,
    tag: str,
    obj_mode: str,
    best_val: float,
    best_x: list[float],
    tq_out: tuple,
    pq_out: tuple,
    sq_out: tuple,
    best_inner: LayerResult,
    best_hx,
    t_total: float,
    ni_rep,
    hots: list,
    colds: list,
) -> None:
    """将优化结果追加写入 CSV 和 JSON，同时保存全部换热/机械过程。"""
    import csv, json

    hot_util = getattr(best_hx, "hot_utility_demand", None)
    cold_util = getattr(best_hx, "cold_utility_demand", None)
    if hot_util is not None:
        hot_util = round(hot_util, 1)
    if cold_util is not None:
        cold_util = round(cold_util, 1)

    row = {
        "tag": tag,
        "obj_mode": obj_mode,
        "H2_kg_s": round(best_inner.h2_mf, 2),
        "h2_T_out_K": round(best_inner.h2_T_out, 1),
        "t_min_K": round(best_x[0], 1),
        "t_max_K": round(best_x[1], 1),
        "p_max_kPa": round(best_x[-2], 1),
        "p_min_kPa": round(best_x[-1], 1),
        "p_q": ",".join(f"{v:.3f}" for v in pq_out) if pq_out else "",
        "s_q": ",".join(f"{v:.3f}" for v in sq_out) if sq_out else "",
        "t_q": ",".join(f"{v:.3f}" for v in tq_out) if tq_out else "",
        "n_sc": best_inner.n_subcycles,
        "sum_flow_kg_s": round(sum(best_inner.flows), 1),
        "flows": ",".join(f"{v:.1f}" for v in best_inner.flows),
        "hot_util_kW": hot_util,
        "cold_util_kW": cold_util,
        "obj": round(best_val, 6),
        "time_s": round(t_total, 1),
        "n_evals": best_inner.n_evals,
    }

    # ── 收集所有过程（换热 + 机械）──
    processes: list[dict] = []
    process_rows: list[dict] = []

    def _rec_from(rec, is_source: bool = False):
        tail = rec.tail_state; head = rec.head_state
        t_lo = min(tail.T, head.T); t_hi = max(tail.T, head.T)
        mf = abs(float(rec.mass_flow)) if rec.mass_flow else None
        pr = float(rec.power_rate) if rec.power_rate else None
        pr_abs = abs(pr) if pr else None
        return {
            "edge_key": rec.edge_key,
            "kind": rec.kind,
            "category": rec.category.name if rec.category else "",
            "fluid": getattr(rec, "fluid", ""),
            "T_lo_K": round(t_lo, 1), "T_hi_K": round(t_hi, 1),
            "P_tail_kPa": round(tail.P, 1),
            "P_head_kPa": round(head.P, 1),
            "S_tail_kJkgK": round(tail.S, 4),
            "S_head_kJkgK": round(head.S, 4),
            "H_tail_kJkg": round(tail.H, 4),
            "H_head_kJkg": round(head.H, 4),
            "power_kW": round(pr_abs, 1) if pr_abs else None,
            "power_signed_kW": round(pr, 1) if pr else None,
            "mass_flow_kg_s": round(mf, 3) if mf else None,
            "is_source": is_source,
        }

    # 循环内部过程
    for _, rec in ni_rep.by_edge:
        d = _rec_from(rec)
        processes.append(d)
        process_rows.append(d)

    # 外部源
    if hots and hots[0] != ni_rep.by_edge[0][1] if ni_rep.by_edge else True:
        d = _rec_from(hots[0], is_source=True)
        processes.append(d)
        process_rows.append(d)
    if colds and hots and len(colds) > 0:
        cold_src = colds[0]
        if not any(cold_src.edge_key == p.get("edge_key") for p in processes):
            d = _rec_from(cold_src, is_source=True)
            processes.append(d)
            process_rows.append(d)

    # ── 写 processes.csv ──
    proc_csv_path = run_dir / "processes.csv"
    write_proc_header = not proc_csv_path.exists()
    with open(proc_csv_path, "a", newline="") as f:
        if process_rows:
            fieldnames = ["tag"] + list(process_rows[0].keys())
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            if write_proc_header:
                writer.writeheader()
            for prw in process_rows:
                writer.writerow({"tag": tag, **prw})

    # ── 写 results.csv ──
    csv_path = run_dir / "results.csv"
    write_header = not csv_path.exists()
    with open(csv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()))
        if write_header:
            writer.writeheader()
        writer.writerow(row)

    # ── 写 results.json ──
    json_path = run_dir / "results.json"
    existing = []
    if json_path.exists():
        with open(json_path, "r") as f:
            existing = json.load(f)
    row["processes"] = processes
    existing.append(row)
    with open(json_path, "w") as f:
        json.dump(existing, f, indent=2, ensure_ascii=False, default=str)

    print(f"  结果已保存: {csv_path}, {proc_csv_path}")


def test_h2_scan() -> None:
    """H2 硬上限扫描: 固定 H2∈{2,3,4,5} 分别优化 0P0T0S 串联夹点。
    
    对比不同 H2 用量的未匹配功率 vs 匹配效率，找性价比最优解。
    """
    h2_levels = [2, 3, 4, 5]

    for h2_val in h2_levels:
        hp = {**_DEFAULT_HP,
              "n_t_q": 0, "n_p_q": 0, "n_s_q": 0,
              "n_lhs": 30, "hx_dT": 10.0,
              "hx_series": True,
              "outer_method": "de",
              "de_popsize": 20, "de_maxiter": 40, "de_F": 0.8, "de_CR": 0.9,
              "n_workers": 16,
              "maxiter_inner": 40, "restarts_inner": 3, "sigma0": 20.0,
              "qstep": 0.001,
              "t_min_lo": 50.0, "t_min_hi": 400.0,
              "t_max_lo": 800.0, "t_max_hi": 1100.0,
              "p_min_lo": 2000.0, "p_min_hi": 2000.0,
              "p_max_lo": 8000.0, "p_max_hi": 12000.0,
              "h2_mf_lo": float(h2_val), "h2_mf_hi": float(h2_val),
              "h2_T_out_lo": 400.0, "h2_T_out_hi": 1240.0,
              "use_non_ideal": False,
              }
        tag = f"H2scan_{h2_val}"
        print(f"\n{'#' * 70}")
        print(f"# H2 = {h2_val} kg/s 固定扫描")
        print(f"{'#' * 70}")
        t0 = time.perf_counter()
        _run_layered(hp, tag=tag)
        elapsed = time.perf_counter() - t0
        print(f"\n  H2={h2_val}: 耗时 {elapsed:.1f}s")

    print("\n" + "=" * 80)
    print("H2 硬上限扫描汇总")
    print("=" * 80)
    print(f"{'H2 kg/s':>8s}  {'obj':>8s}  {'未匹配kW':>10s}  "
          f"{'匹配kW':>9s}  {'h2tout K':>8s}  "
          f"{'t_maxK':>7s} {'t_minK':>7s} {'p_maxkPa':>9s} "
          f"{'n_sc':>4s} {'Σflow':>7s}")
    print("-" * 80)
    print("  (结果见上方各 run 的输出)")
    print("=" * 80)


def test_h2_scan_1p1s() -> None:
    """H2 硬上限扫描 1P1S: 固定 H2∈{3,4,5} 分别优化 串联夹点。
    
    比 0P0T0S 多了 1 个 P 分位 + 1 个 S 分位，内外层迭代加重。
    """
    h2_levels = [3, 4, 5]

    for h2_val in h2_levels:
        hp = {**_DEFAULT_HP,
              "n_t_q": 0, "n_p_q": 1, "n_s_q": 1,
              "n_lhs": 30, "hx_dT": 10.0,
              "hx_series": True,
              "outer_method": "de",
              "de_popsize": 20, "de_maxiter": 30, "de_F": 0.8, "de_CR": 0.9,
              "n_workers": 8,
              "maxiter_inner": 25, "restarts_inner": 2, "sigma0": 20.0,
              "qstep": 0.001,
              "t_min_lo": 50.0, "t_min_hi": 400.0,
              "t_max_lo": 800.0, "t_max_hi": 1100.0,
              "p_min_lo": 2000.0, "p_min_hi": 2000.0,
              "p_max_lo": 8000.0, "p_max_hi": 12000.0,
              "h2_mf_lo": float(h2_val), "h2_mf_hi": float(h2_val),
              "h2_T_out_lo": 400.0, "h2_T_out_hi": 1240.0,
              "use_non_ideal": False,
              }
        tag = f"H2scan_1P1S_{h2_val}"
        print(f"\n{'#' * 70}")
        print(f"# 1P1S H2 = {h2_val} kg/s 固定扫描")
        print(f"{'#' * 70}")
        t0 = time.perf_counter()
        _run_layered(hp, tag=tag)
        elapsed = time.perf_counter() - t0
        print(f"\n  1P1S H2={h2_val}: 耗时 {elapsed:.1f}s")

    print("\n" + "=" * 80)
    print("1P1S H2 硬上限扫描汇总")
    print("=" * 80)
    print(f"{'H2 kg/s':>8s}  {'obj':>8s}  {'未匹配kW':>10s}  "
          f"{'匹配kW':>9s}  {'h2tout K':>8s}  "
          f"{'t_maxK':>7s} {'t_minK':>7s} {'p_maxkPa':>9s} "
          f"{'n_sc':>4s} {'Σflow':>7s}  {'s_q':>6s} {'p_q':>6s}")
    print("-" * 80)
    print("  (结果见上方各 run 的输出)")
    print("=" * 80)


def test_h2_1p1s_h2fix() -> None:
    """1P1S H2=5固定 h2_T_out=800K固定 — 只在循环拓扑+流量上优化。

    验证当 H2 完全受限时，P/S 分位是否提供有价值的额外自由度。
    """
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 1, "n_s_q": 1,
          "n_lhs": 40, "hx_dT": 10.0,
          "hx_series": True,
          "outer_method": "de",
          "de_popsize": 20, "de_maxiter": 50, "de_F": 0.8, "de_CR": 0.9,
          "n_workers": 8,
          "maxiter_inner": 50, "restarts_inner": 4, "sigma0": 20.0,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 2000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 5.0, "h2_mf_hi": 5.0,
          "h2_T_out_lo": 800.0, "h2_T_out_hi": 800.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="1P1S_H2-5_T800")


def test_h2_1p1s_h2fix_pinch() -> None:
    """1P1S H2=5固定 h2_T_out=800K + 全局夹点目标函数 — 轻量版。

    验证夹点法下 P/S 分位是否仍有价值。LHS 已能找到最优解，DE 仅做保险。
    """
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 1, "n_s_q": 1,
          "n_lhs": 20, "hx_dT": 10.0,
          "obj_mode": "pinch",
          "outer_method": "de",
          "de_popsize": 15, "de_maxiter": 10, "de_F": 0.8, "de_CR": 0.9,
          "n_workers": 8,
          "maxiter_inner": 25, "restarts_inner": 2, "sigma0": 20.0,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 2000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 5.0, "h2_mf_hi": 5.0,
          "h2_T_out_lo": 800.0, "h2_T_out_hi": 800.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="1P1S_H2-5_T800_pinch")


def test_h2_1p1s_h2_4_pinch() -> None:
    """1P1S H2=4固定 h2_T_out=800K + pinch 目标 — 轻量版。

    H2 从 5→4 kg/s，验证更严约束下 P/S 分位表现。
    """
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 1, "n_s_q": 1,
          "n_lhs": 20, "hx_dT": 10.0,
          "obj_mode": "pinch",
          "outer_method": "de",
          "de_popsize": 15, "de_maxiter": 10, "de_F": 0.8, "de_CR": 0.9,
          "n_workers": 8,
          "maxiter_inner": 25, "restarts_inner": 2, "sigma0": 20.0,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 2000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 4.0, "h2_mf_hi": 4.0,
          "h2_T_out_lo": 800.0, "h2_T_out_hi": 800.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="1P1S_H2-4_T800_pinch")


def test_h2_1p1s_h2_scan_pinch() -> None:
    """1P1S pinch 目标，逐级降低 H2 流量，找最低可行 H2。

    H2 = 3.5, 3.0, 2.5, 2.0, 1.5 kg/s，h2_T_out=800K 固定。
    记录各 H2 下的 obj 和公用工程需求。
    """
    h2_levels = [3.5, 3.0, 2.5, 2.0, 1.5]

    for h2_val in h2_levels:
        hp = {**_DEFAULT_HP,
              "n_t_q": 0, "n_p_q": 1, "n_s_q": 1,
              "n_lhs": 20, "hx_dT": 10.0,
              "obj_mode": "pinch",
              "outer_method": "de",
              "de_popsize": 15, "de_maxiter": 10, "de_F": 0.8, "de_CR": 0.9,
              "n_workers": 8,
              "maxiter_inner": 25, "restarts_inner": 2, "sigma0": 20.0,
              "qstep": 0.001,
              "t_min_lo": 50.0, "t_min_hi": 400.0,
              "t_max_lo": 800.0, "t_max_hi": 1100.0,
              "p_min_lo": 2000.0, "p_min_hi": 2000.0,
              "p_max_lo": 8000.0, "p_max_hi": 12000.0,
              "h2_mf_lo": float(h2_val), "h2_mf_hi": float(h2_val),
              "h2_T_out_lo": 800.0, "h2_T_out_hi": 800.0,
              "use_non_ideal": False,
              }
        print(f"\n{'#'*70}")
        print(f"# Pinch 1P1S H2 = {h2_val} kg/s")
        print(f"{'#'*70}")
        t0 = time.perf_counter()
        _run_layered(hp, tag=f"1P1S_H2-{h2_val}_T800_pinch")
        print(f"  H2={h2_val}: {time.perf_counter()-t0:.1f}s")

    print("\n" + "=" * 80)
    print("Pinch 1P1S H2 扫描完成 — 见上方各 run 输出")
    print("=" * 80)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "1t2p":
        test_layered_1t2p()
    elif len(sys.argv) > 1 and sys.argv[1] == "long":
        test_layered_2p0t_long()
    elif len(sys.argv) > 1 and sys.argv[1] == "0p0t":
        test_layered_0p0t()
    elif len(sys.argv) > 1 and sys.argv[1] == "wide":
        test_layered_0p0t_wide()
    elif len(sys.argv) > 1 and sys.argv[1] == "h2scan":
        test_h2_scan()
    else:
        test_layered_optimization()
