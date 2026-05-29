r"""分层优化 (加速版): 外层 LHS → 内层 CMA 复用拓扑 (maxiter=10, restarts=1), 外层多进程并行。

   A. 内层: 固定拓扑后只更新流量+非理想+性能+HX, ~50μs/eval
   B. 外层: ProcessPoolExecutor 并行 6 进程
   C. 减代数: maxiter=10, restarts=1
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
}


def _make_h2_source(hp: dict, h2_mf: float, h2_T_out: float | None = None) -> ExternalSourceInput:
    """用 hp 中的 H₂ 参数 + 可变流量/出口温构造冷源。"""
    return ExternalSourceInput(
        fluid="Hydrogen", mass_flow=h2_mf,
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
    cold = ExternalSourceInput(fluid="Hydrogen", mass_flow=h2_mass_flow,
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

    total_q = sum(abs(float(r.power_rate)) for r in hots + colds if r.power_rate)
    if total_q < 1e-12:
        return 1.0
    hx = match_heat_exchanger_groups(hots, colds, dT_min=hp["hx_dT"],
                                      max_group_size=hp["hx_max_group_size"])
    unmatched_ratio = hx.total_unmatched / total_q
    num_unmatched = len(hx.unassigned_hots) + len(hx.unassigned_colds)
    # 分阶段 merit: 最小化未匹配比(主) + 最少未匹配过程数(权重 1e-3)
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
    stds = [hp["sigma0"] * (hi[i] - lo[i]) for i in range(dim)]
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

            improved = False
            for s, sc in zip(solutions, scores):
                if sc < global_best_val:
                    global_best_val = sc
                    global_best_x = list(s)
                    improved = True
                    stall = 0
            if not improved:
                stall += 1
            if gen >= hp["maxiter_inner"] or stall >= hp["early_stop"] * len(solutions):
                break

        h2_val, h2tout_val = _extract_extra(global_best_x, n_sc)
        sys.stdout.write(
            f"    r{restart}: gen={gen} obj={global_best_val:.5f} "
            f"evals={total_evals} dim={dim} n_sc={n_sc} "
            f"best_flow_sum={sum(global_best_x[:n_sc]):.0f} "
            f"H2={h2_val:.1f} h2tout={h2tout_val:.0f}\n"
        )
        sys.stdout.flush()

    flows = global_best_x[:n_sc]
    flow_step = hp["flow_step"]
    flows = [round(f / flow_step) * flow_step for f in flows]
    h2_mf, h2_T_out_unrounded = _extract_extra(global_best_x, n_sc)
    h2_mf = round(h2_mf / hp["h2_step"]) * hp["h2_step"] if not h2_fixed else h2_mf
    h2_T_out = round(h2_T_out_unrounded) if not h2_T_fixed else h2_T_out_unrounded

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

    result, _ = _inner_cma_fast(tp_in, sys_inp, hp, seed=seed)
    return result


# ──────────────────────────────────────────────────────────────────
# D. 外层 DE 并行 worker
# ──────────────────────────────────────────────────────────────────


def _de_trial_worker(args: tuple) -> float:
    """模块级 worker: 外层参数 → 内层 CMA → obj（仅回传标量）。"""
    (x_outer, sys_inp, hp, seed) = args

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
    """0T+1P+1S 理想 HXmax=1 H2固定4.3 h2_T_out∈[400,900] p_min=2000kPa t_max∈[800,1100] t_min∈[50,400] p_max∈[8000,12000]"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 1, "n_s_q": 1, "n_lhs": 40, "hx_dT": 10.0,
          "hx_max_group_size": 2,
          "maxiter_inner": 15, "restarts_inner": 3,
          "de_popsize": 12, "de_maxiter": 60, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 2000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 2.0, "h2_mf_hi": 5.5,
          "h2_T_out_lo": 400.0, "h2_T_out_hi": 900.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="0P1T_S1_HX2")


def test_layered_0p0t_sq1() -> None:
    """0P0T + 1S HXmax=2 h2_T_out=600K p_min=2000kPa t_max∈[800,1100] t_min∈[50,400] p_max∈[8000,12000] H2∈[3,6]"""
    hp = {**_DEFAULT_HP,
          "n_t_q": 0, "n_p_q": 0, "n_s_q": 1, "n_lhs": 40, "hx_dT": 10.0,
          "hx_max_group_size": 2,
          "maxiter_inner": 10, "restarts_inner": 2,
          "de_popsize": 12, "de_maxiter": 40, "de_F": 0.8, "de_CR": 0.9,
          "qstep": 0.001,
          "t_min_lo": 50.0, "t_min_hi": 400.0,
          "t_max_lo": 800.0, "t_max_hi": 1100.0,
          "p_min_lo": 2000.0, "p_min_hi": 2000.0,
          "p_max_lo": 8000.0, "p_max_hi": 12000.0,
          "h2_mf_lo": 3.0, "h2_mf_hi": 6.0,
          "h2_T_out_lo": 600.0, "h2_T_out_hi": 600.0,
          "use_non_ideal": False,
          }
    _run_layered(hp, tag="0P0T_S1")


def _run_layered(hp: dict, tag: str = "") -> None:
    hp = {**_DEFAULT_HP, **hp}
    tq_tag = f"{hp['n_t_q']}T{hp['n_p_q']}P"
    if hp.get("n_s_q", 0) > 0:
        tq_tag += f"{hp['n_s_q']}S"
    print("=" * 60)
    print(f"分层优化 (LHS+DE并行) [{tag}] {tq_tag}")
    print(f"  LHS={hp['n_lhs']}点 DE={hp['de_popsize']}x{hp['de_maxiter']} "
          f"t_min∈[{hp['t_min_lo']},{hp['t_min_hi']}] "
          f"t_max∈[{hp['t_max_lo']},{hp['t_max_hi']}] "
          f"p_min∈[{hp['p_min_lo']},{hp['p_min_hi']}] "
          f"p_max∈[{hp['p_max_lo']},{hp['p_max_hi']}] "
          f"H2∈[{hp['h2_mf_lo']},{hp['h2_mf_hi']}] "
          f"interp_He={hp.get('use_interp_he', False)} non_ideal={hp['use_non_ideal']} hx_dT={hp['hx_dT']}")
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

    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as ex:
        # ── 阶段1: LHS 并行采样 ──
        print(f"\n  [阶段1] LHS {hp['n_lhs']} 点并行采样...")
        seed = hp.get("seed", 42)
        lhs_samples = _lhs(hp["n_lhs"], outer_bounds, seed=seed)
        tasks = []
        for i, row in enumerate(lhs_samples):
            t_min = round(row[0]); t_max = round(row[1])
            idx = 2
            t_q_vals = tuple(round(round(row[idx + j] / hp['qstep']) * hp['qstep'], 3)
                             for j in range(hp["n_t_q"]))
            idx += hp["n_t_q"]
            p_q_vals = tuple(round(round(row[idx + j] / hp['qstep']) * hp['qstep'], 3)
                             for j in range(hp["n_p_q"]))
            idx += hp["n_p_q"]
            s_q_vals = tuple(round(round(row[idx + j] / hp['qstep']) * hp['qstep'], 3)
                             for j in range(hp["n_s_q"]))
            idx += hp["n_s_q"]
            p_max = round(row[idx]); p_min = round(row[idx + 1])
            tasks.append((t_min, t_max, t_q_vals, p_q_vals, s_q_vals, p_max, p_min,
                          base, hp, seed + i))

        lhs_results: list[LayerResult] = []
        futures = [ex.submit(_sample_worker, t) for t in tasks]
        for fut in futures:
            r = fut.result(timeout=600)
            if r.obj < 1e8 and r.n_subcycles > 0:
                lhs_results.append(r)

        lhs_objs = sorted(lhs_results, key=lambda r: r.obj)
        if not lhs_objs:
            print("  LHS: 0 有效 — 拓扑无子循环，终止")
            assert False, "所有 LHS 样本均无有效拓扑"
        print(f"  LHS: {len(lhs_results)} 有效, best={lhs_objs[0].obj:.5f}")

        # ── 阶段2: DE 并行评估 ──
        print(f"\n  [阶段2] DE 并行 (popsize={hp['de_popsize']}, maxiter={hp['de_maxiter']})...")
        rng = _rnd_.Random(seed)
        lo = [b[0] for b in outer_bounds]
        hi = [b[1] for b in outer_bounds]

        pop: list[list[float]] = []
        for r in lhs_objs[:hp['de_popsize']]:
            x = [r.t_min, r.t_max]
            for q in r.t_q:
                x.append(q)
            for q in r.p_q:
                x.append(q)
            for q in getattr(r, 's_q', ()):
                x.append(q)
            x.append(r.p_max)
            x.append(r.p_min)
            pop.append(x)
        while len(pop) < hp['de_popsize']:
            pop.append([rng.uniform(lo[d], hi[d]) for d in range(dim)])

        g0_args = [(p, base, hp, 0) for p in pop]
        g0_futures = [ex.submit(_de_trial_worker, a) for a in g0_args]
        scores = [f.result(timeout=600) for f in g0_futures]

        best_idx = min(range(len(scores)), key=lambda i: scores[i])
        best_x = pop[best_idx][:]
        best_val = scores[best_idx]
        stall = 0
        n_evals = len(scores)
        conv_history: list[float] = [best_val]
        print(f"  DE gen 0: best={best_val:.5f}")

        for gen in range(hp["de_maxiter"]):
            trials: list[tuple[int, list[float]]] = []
            for i in range(hp["de_popsize"]):
                pool = [j for j in range(hp["de_popsize"]) if j != i]
                a, b, c = rng.sample(pool, 3)
                v = [pop[a][d] + hp["de_F"] * (pop[b][d] - pop[c][d]) for d in range(dim)]
                for d in range(dim):
                    v[d] = max(lo[d], min(hi[d], v[d]))
                jr = rng.randrange(dim)
                u = [v[d] if rng.random() < hp["de_CR"] or d == jr else pop[i][d]
                     for d in range(dim)]
                trials.append((i, u))

            trial_args = [(t, base, hp, 0) for _, t in trials]
            t_futures = [ex.submit(_de_trial_worker, a) for a in trial_args]
            trial_scores = [f.result(timeout=600) for f in t_futures]
            n_evals += len(trial_scores)

            improved = False
            for (i, u), us in zip(trials, trial_scores):
                if us < scores[i]:
                    pop[i] = u
                    scores[i] = us
                    if us < best_val:
                        best_val = us
                        best_x = u[:]
                        improved = True
                        stall = 0
            if not improved:
                stall += hp["de_popsize"]
            conv_history.append(best_val)

            sys.stdout.write(f"  DE gen {gen+1:2d}: best={best_val:.5f} "
                             f"stall={stall}/{hp['de_popsize']*hp['de_maxiter']} "
                             f"t_max={best_x[1]:.0f} t_min={best_x[0]:.0f} "
                             f"pmax={best_x[-2]:.0f} pmin={best_x[-1]:.0f}\n")
            sys.stdout.flush()

            if stall >= hp["de_popsize"] * hp["de_maxiter"] * 2:
                break

    t_total = time.perf_counter() - t0_all

    # ── 汇总 ──
    print("\n" + "=" * 55)
    print(f"优化完成: LHS{len(lhs_results)} + DE{hp['de_popsize']}x{gen+1}gen "
          f"= {n_evals} eval, {t_total:.1f}s")
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

    # ── 最终重建 + CoolProp 验证 ──
    print("\n  [最终重建]")
    try:
        # 用最优参数在主进程重建一次，获取完整结果（无需 pickle 传输）
        opt_result, _, _, _, _, _, _ = \
            _rebuild_result(best_x, base, hp, use_coolprop=False)
        best_inner = opt_result
    except Exception:
        best_inner = LayerResult(n=0, t_min=0, t_max=0, p_min=0, p_max=0,
                                 t_q=(), p_q=(), flows=[], h2_mf=0,
                                 obj=best_val, n_subcycles=0, n_evals=0, runtime=0)

    print(f"    flows={[f'{v:.1f}' for v in best_inner.flows]} "
          f"H2={best_inner.h2_mf:.2f}kg/s  n_sc={best_inner.n_subcycles}")

    # ── CoolProp 验证 ──
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
    run_dir = TESTS_DIR / "run_layered_fast"
    run_dir.mkdir(parents=True, exist_ok=True)

    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig1, ax1 = plt.subplots(figsize=(8, 4))
    ax1.plot(range(len(conv_history)), conv_history, "o-", ms=4, color="tab:blue")
    ax1.axhline(y=best_val, color="tab:red", lw=1, ls="--", label=f"best={best_val:.5f}")
    ax1.set_xlabel("DE generation"); ax1.set_ylabel("obj (unmatched/total_Q)")
    ax1.set_title(f"DE convergence [{tag}] {tq_tag}")
    ax1.legend(); ax1.grid(True, alpha=0.2)
    fig1.tight_layout(); fig1.savefig(run_dir / f"conv_{tag}.png", dpi=150); plt.close(fig1)
    print(f"  收敛图: {run_dir / f'conv_{tag}.png'}")

    try:
        _, best_hx, ideal_rep, ni_rep, best_layer, _hots, _colds = \
            _rebuild_result(best_x, base, hp)

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

        # HX T-Q
        _draw_hx_tq(_hots, _colds, best_hx, run_dir / f"hx_tq_{tag}.png")
        print(f"  HX T-Q: {run_dir / f'hx_tq_{tag}.png'}")
        print(f"  HX单元={len(best_hx.units)} "
              f"未匹配={best_hx.total_unmatched:.0f}kW "
              f"总匹配={best_hx.total_matched:.0f}kW")
    except Exception as e:
        import traceback
        print(f"  绘图异常: {e}")
        traceback.print_exc()

    assert best_val < 1e3


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
    else:
        test_layered_optimization()
