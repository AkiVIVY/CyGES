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
# B. 拓扑基座构建（供 CMA 和 L-BFGS-B 共用）
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


# ──────────────────────────────────────────────────────────────────
# B2. 内层 L-BFGS-B + 多起点 (Sobol 支持)
# ──────────────────────────────────────────────────────────────────


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
