"""差分进化/CMA-ES 求解器，支持多进程并行评估与高斯基编码。"""

from __future__ import annotations

import math
import random
import warnings
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

from core.closed_cycle_layer import ClosedCycleLayer, ClosedCycleTPInput
from core.fluid_property_solver import PropertyRegistry
from core.system import (
    CycleConfig,
    SystemInput,
    SystemPipeline,
    SystemResult,
    analyze_system_heat,
)
from optimize.objective import OBJECTIVES
from optimize.types import OptimizationResult


def _differential_evolution(
    fn: Callable[[tuple[float, ...]], float],
    bounds: list[tuple[float, float]],
    *,
    popsize: int = 15,
    maxiter: int = 200,
    F: float = 0.8,
    CR: float = 0.9,
    seed: int | None = 42,
    early_stop: int = 30,
    eval_state: _EvalState | None = None,
    n_workers: int = 1,
    callback: Callable[[int, int, list[float], float, int], None] | None = None,
) -> tuple[list[float], float, int]:
    """DE/rand/1/bin 同步型差分进化（每代批量评估，支持多进程并行）。

    :param eval_state: 多进程评估状态（``None`` = 用主进程 fn 串行）。
    :param callback: 每代末调用 ``callback(gen, restart, best_x, best_val, n_evals)``；``None`` 不回调。
    :returns: ``(best_x, best_val, n_evaluations)``。
    """
    rng = random.Random(seed)
    dim = len(bounds)
    lo = [b[0] for b in bounds]
    hi = [b[1] for b in bounds]

    # 持久 executor
    _executor = None
    if eval_state is not None and n_workers > 1:
        import concurrent.futures, multiprocessing
        ctx = multiprocessing.get_context("spawn")
        _executor = concurrent.futures.ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx)
    executor = _executor

    try:
        # 随机初始化种群（并行评估初代）
        pop = [
            [rng.uniform(lo[d], hi[d]) for d in range(dim)]
            for _ in range(popsize)
        ]
        scores = _eval_batch(fn, [list(p) for p in pop], eval_state, n_workers, executor)
        n_evals = len(scores)

        best_idx = min(range(popsize), key=lambda i: scores[i])
        best_x = pop[best_idx][:]
        best_val = scores[best_idx]
        stall = 0

        for gen in range(maxiter):
            # ── 同步型 DE：先生成所有尝试向量，再批量评估，最后统一替换 ──
            trials: list[tuple[int, list[float]]] = []
            for i in range(popsize):
                pool = [j for j in range(popsize) if j != i]
                a_idx, b_idx, c_idx = rng.sample(pool, 3)
                v = [
                    pop[a_idx][d] + F * (pop[b_idx][d] - pop[c_idx][d])
                    for d in range(dim)
                ]
                for d in range(dim):
                    v[d] = max(lo[d], min(hi[d], v[d]))
                j_rand = rng.randrange(dim)
                u = [
                    v[d] if rng.random() < CR or d == j_rand else pop[i][d]
                    for d in range(dim)
                ]
                trials.append((i, u))

            trial_scores = _eval_batch(fn, [t for _, t in trials], eval_state, n_workers, executor)
            n_evals += len(trial_scores)

            improved = False
            for (i, u), u_score in zip(trials, trial_scores):
                if u_score < scores[i]:
                    pop[i] = u
                    scores[i] = u_score
                    if u_score < best_val:
                        best_val = u_score
                        best_x = u[:]
                        improved = True
                        stall = 0
            if not improved:
                stall += popsize

            if callback is not None:
                callback(gen, 0, best_x, best_val, n_evals)

            if stall >= early_stop * popsize:
                break
    finally:
        if _executor is not None:
            _executor.shutdown(wait=False)

    return best_x, best_val, n_evals


def _cma_evolution(
    fn: Callable[[tuple[float, ...]], float],
    bounds: list[tuple[float, float]],
    *,
    sigma0: float = 0.3,
    maxiter: int = 200,
    seed: int | None = 42,
    early_stop: int = 30,
    restarts: int = 5,
    eval_state: _EvalState | None = None,
    n_workers: int = 1,
    callback: Callable[[int, int, list[float], float, int], None] | None = None,
) -> tuple[list[float], float, int]:
    """CMA-ES + 多起点重启（BIPOP 风格），自动交替大小种群避免局部最优。

    :param eval_state: 多进程评估状态（``None`` = 用主进程 fn 串行）。
    :param n_workers: 并行 worker 进程数。
    :param callback: 每代末调用 ``callback(gen, restart, best_x, best_val, n_evals)``。
    :returns: ``(best_x, best_val, n_evaluations)``。
    """
    import cma

    dim = len(bounds)
    lo = [b[0] for b in bounds]
    hi = [b[1] for b in bounds]
    stds = [sigma0 * (hi[i] - lo[i]) for i in range(dim)]
    x0_init = [(lo[i] + hi[i]) / 2 for i in range(dim)]

    executor = None
    if eval_state is not None and n_workers > 1:
        import concurrent.futures
        import multiprocessing
        ctx = multiprocessing.get_context("spawn")
        executor = concurrent.futures.ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx)

    try:
        global_best_val = float("inf")
        global_best_x = x0_init[:]
        total_evals = 0
        gen_offset = 0

        for restart in range(restarts):
            # BIPOP: 交替小种群（4+3*log(dim)）和大种群（dim²）
            if restart % 2 == 0:
                p = max(4, int(4 + 3 * math.log2(dim)))
            else:
                p = max(10, min(dim * dim, 50))

            options: dict[str, Any] = {
                "bounds": [lo, hi],
                "CMA_stds": stds,
                "popsize": p,
                "verbose": -9,
                "seed": seed + restart if seed is not None else None,
            }

            if restart == 0:
                x0 = x0_init[:]
            else:
                rng_x0 = random.Random(seed + restart if seed is not None else None)
                x0 = [rng_x0.uniform(lo[i], hi[i]) for i in range(dim)]

            es = cma.CMAEvolutionStrategy(x0, sigma0, options)
            stall = 0
            gen = 0

            while not es.stop():
                solutions = es.ask()
                scores = _eval_batch(fn, [list(s) for s in solutions], eval_state, n_workers, executor)
                es.tell(solutions, scores)
                total_evals += len(solutions)
                gen += 1

                improved = False
                for s, sc in zip(solutions, scores):
                    if sc < global_best_val:
                        global_best_val = sc
                        global_best_x = s
                        improved = True
                        stall = 0
                if not improved:
                    stall += 1

                if callback is not None:
                    callback(gen, restart, global_best_x, global_best_val, total_evals)

                if gen >= maxiter or stall >= early_stop * len(solutions):
                    break

            gen_offset += gen
    finally:
        if executor is not None:
            executor.shutdown(wait=False)

    return global_best_x, global_best_val, total_evals


def _round_and_dedup(values: tuple[float, ...], step: float, merge_ratio: float = 1.0) -> tuple[float, ...]:
    """舍入到 step 整数倍，去重，夹在 [step, 1-step] 内。

    ``step <= 0`` 时不做离散化，原值返回。
    ``merge_ratio > 0`` 时，舍入后相邻间距 ≤ ratio×step 的分位点合并为一个。
    """
    if step <= 0.0 or not values:
        return values
    lo, hi = 0.0, 1.0
    rounded = sorted(set(max(lo, min(hi, round(v / step) * step)) for v in values))
    if merge_ratio <= 0.0:
        return tuple(rounded)
    # 合并相邻过近点
    merged: list[float] = [rounded[0]]
    for v in rounded[1:]:
        if v - merged[-1] <= merge_ratio * step:
            continue
        merged.append(v)
    return tuple(merged)


# ────────────────────────────────────────────
# 多进程评估状态与 worker（可 pickle 序列化）
# ────────────────────────────────────────────


@dataclass
class _EvalState:
    """可序列化的评估配置，供 worker 进程独立重建管线。"""
    base_input: SystemInput
    objective_name: str
    mf_step: float
    qstep: float
    tstep: float
    qmerge: float
    mf_lo: float
    mf_hi: float
    n_t_q: int
    n_p_q: int
    basis: bool
    basis_centers: list[tuple[float, float]] | None = None
    basis_sigma: tuple[float, float] | None = None
    basis_mf_dim: int = 0


def _eval_worker(x: tuple[float, ...], state: _EvalState) -> float:
    """Worker 进程入口：重建 props 后执行完整评估管线。"""
    props = PropertyRegistry()
    return _eval_impl(tuple(x), state, props)


def _eval_chunk(solutions_and_state: tuple[list[list[float]], _EvalState]) -> list[float]:
    """模块级 worker：评估一个 chunk 内的所有解（供 ProcessPoolExecutor spawn 调用）。"""
    import pickle
    # reconstruct from tuple (ProcessPoolExecutor pickles args)
    solutions, state = solutions_and_state
    return [_eval_worker(tuple(s), state) for s in solutions]


def _eval_impl(x: tuple[float, ...], state: _EvalState, props: PropertyRegistry) -> float:
    """核心评估逻辑（复用主进程和 worker 进程）。"""
    import math

    idx = 0
    t_max = float(x[idx]); idx += 1
    t_min = float(x[idx]); idx += 1
    if state.tstep > 0:
        t_max = round(t_max / state.tstep) * state.tstep
        t_min = round(t_min / state.tstep) * state.tstep
    t_q = _round_and_dedup(tuple(x[idx + i] for i in range(state.n_t_q)), state.qstep, state.qmerge)
    idx += state.n_t_q
    p_q = _round_and_dedup(tuple(x[idx + i] for i in range(state.n_p_q)), state.qstep, state.qmerge)
    idx += state.n_p_q
    mf_all = x[idx:]

    if state.mf_step > 0:
        mf_all = [round(m / state.mf_step) * state.mf_step for m in mf_all]

    tp = ClosedCycleTPInput(
        fluid=state.base_input.cycles[0].input.fluid,
        t_min=t_min, t_max=t_max,
        p_min=state.base_input.cycles[0].input.p_min,
        p_max=state.base_input.cycles[0].input.p_max,
        t_quantiles=t_q, p_quantiles=p_q,
        subcycle_mass_flow_initial=state.base_input.cycles[0].input.subcycle_mass_flow_initial,
    )

    try:
        probe = ClosedCycleLayer(tp)
    except Exception:
        return 1.0
    n_sc = len(probe.subcycles)
    if n_sc == 0:
        return 1.0

    if state.basis and state.basis_centers is not None:
        weights = list(mf_all[:state.basis_mf_dim])
        mf_truncated = _decode_basis_flows_impl(probe, weights, state.basis_centers, state.basis_sigma)
        if state.mf_step > 0:
            mf_truncated = [round(f / state.mf_step) * state.mf_step for f in mf_truncated]
    else:
        mf_truncated = [float(mf_all[i]) for i in range(min(n_sc, len(mf_all)))]

    cfg = CycleConfig(
        input=tp,
        use_non_ideal=state.base_input.cycles[0].use_non_ideal,
        subcycle_mass_flows=mf_truncated,
        delta_T_min=state.base_input.cycles[0].delta_T_min,
        heat_method=state.base_input.cycles[0].heat_method,
    )
    sys_inp = SystemInput(
        heat_sources=state.base_input.heat_sources,
        cold_sources=state.base_input.cold_sources,
        cycles=(cfg,),
        delta_T_min=state.base_input.delta_T_min,
        heat_method=state.base_input.heat_method,
    )
    try:
        raw = SystemPipeline(sys_inp).run(props)
        result = analyze_system_heat(raw, sys_inp, props)
    except Exception:
        return 1.0
    return OBJECTIVES[state.objective_name](result)


def _decode_basis_flows_impl(layer, weights, centers, sigma) -> list[float]:
    """模块级基函数解码（供 worker 调用）。"""
    import math
    sig_s, sig_p = sigma
    denom_s = 2.0 * sig_s * sig_s
    denom_p = 2.0 * sig_p * sig_p
    flows: list[float] = []
    for sc in layer.subcycles:
        nodes = [layer.nodes[i] for i in sc.nodes]
        sc_s = sum(n.S for n in nodes) / len(nodes)
        sc_p = sum(n.P for n in nodes) / len(nodes)
        mf = 0.0
        for w, (cs, cp) in zip(weights, centers):
            mf += w * math.exp(
                -((sc_s - cs) ** 2) / denom_s
                - ((sc_p - cp) ** 2) / denom_p
            )
        flows.append(mf)
    return flows


def _eval_batch(
    fn_serial: Callable[[tuple[float, ...]], float],
    solutions: list[list[float]],
    state: _EvalState | None,
    n_workers: int,
    executor=None,
) -> list[float]:
    """评估一组解：有 executor 时用并行 chunk，否则串行。"""
    if executor is not None and n_workers > 1 and len(solutions) >= n_workers * 2:
        chunk_size = max(1, len(solutions) // n_workers)
        chunks = [solutions[i : i + chunk_size] for i in range(0, len(solutions), chunk_size)]
        futures = [executor.submit(_eval_chunk, (ch, state)) for ch in chunks]
        results: list[float] = []
        for f in futures:
            results.extend(f.result(timeout=600))
        return results
    return [fn_serial(tuple(s)) for s in solutions]


class Optimizer:
    """单目标优化器，支持 DE 和 CMA-ES。

    用法::

        opt = Optimizer(base_input=sys_inp, props=props, objective="min_max_utility")
        result = opt.run(popsize=15, maxiter=200)
        print(result.objective, result.x_opt)
    """

    def __init__(
        self,
        base_input: SystemInput,
        props: PropertyRegistry,
        objective: str | Callable[[SystemResult], float] = "min_max_utility",
        max_subcycles: int | None = None,
        mf_step_fraction: float = 0.001,
        quantile_step: float = 0.02,
        t_step: float = 10.0,
        quantile_merge_ratio: float = 1.0,
        mf_bounds: tuple[float, float] = (-10.0, 50.0),
        basis_encoding: bool = False,
        basis_s: int = 3,
        basis_p: int = 3,
        n_workers: int = 1,
    ):
        if not base_input.cycles:
            raise ValueError("base_input 须包含至少一个 CycleConfig")

        self._base = base_input
        self._props = props
        self._base_cfg = base_input.cycles[0]
        self._base_tp = self._base_cfg.input
        self._mf_step = mf_step_fraction
        self._qstep = quantile_step
        self._tstep = t_step
        self._qmerge = quantile_merge_ratio
        self._mf_lo, self._mf_hi = mf_bounds
        self._basis = basis_encoding
        self._n_workers = n_workers

        if isinstance(objective, str):
            if objective not in OBJECTIVES:
                raise ValueError(f"未知目标 '{objective}'，已知: {list(OBJECTIVES)}")
            self._objective_fn = OBJECTIVES[objective]
            self._objective_name = objective
        else:
            self._objective_fn = objective
            self._objective_name = getattr(objective, "__name__", "custom")

        if basis_encoding:
            self._n_basis_s = basis_s
            self._n_basis_p = basis_p
            self._mf_dim = basis_s * basis_p
            self._build_basis_grid()
        else:
            if max_subcycles is None:
                max_subcycles = self._probe_max_sc()
            self._mf_dim = max_subcycles

        self._max_sc = max_subcycles if not basis_encoding else 0
        self._n_t_q = len(self._base_tp.t_quantiles)
        self._n_p_q = len(self._base_tp.p_quantiles)

    def _probe_max_sc(self) -> int:
        """用初始分位探测子循环数，×1.5 安全系数作为 mf 维度上限。"""
        probe = ClosedCycleLayer(self._base_tp)
        return max(1, int(len(probe.subcycles) * 1.5))

    def _to_eval_state(self) -> _EvalState:
        """导出可 pickle 的评估状态（供多进程 worker 使用）。"""
        return _EvalState(
            base_input=self._base,
            objective_name=self._objective_name,
            mf_step=self._mf_step * (self._mf_hi - self._mf_lo) if self._mf_step > 0 else 0.0,
            qstep=self._qstep,
            tstep=self._tstep,
            qmerge=self._qmerge,
            mf_lo=self._mf_lo,
            mf_hi=self._mf_hi,
            n_t_q=self._n_t_q,
            n_p_q=self._n_p_q,
            basis=self._basis,
            basis_centers=list(self._basis_centers) if self._basis else None,
            basis_sigma=self._basis_sigma if self._basis else None,
            basis_mf_dim=self._mf_dim if self._basis else 0,
        )

    def _build_basis_grid(self) -> None:
        """构建高斯基网格：在 PS 空间的 [S_min, S_max] × [P_min, P_max] 上均匀放置核中心。"""
        import math
        probe = ClosedCycleLayer(self._base_tp)
        s_vals = [n.S for n in probe.nodes.values()]
        s_min, s_max = min(s_vals), max(s_vals)
        p_min, p_max = self._base_tp.p_min, self._base_tp.p_max

        # 网格间距
        ds = (s_max - s_min) / (self._n_basis_s + 1)
        dp = (p_max - p_min) / (self._n_basis_p + 1)
        self._basis_sigma = (ds / 2, dp / 2)
        # 核中心坐标 (S, P)
        self._basis_centers: list[tuple[float, float]] = []
        for si in range(1, self._n_basis_s + 1):
            for pi in range(1, self._n_basis_p + 1):
                self._basis_centers.append((s_min + si * ds, p_min + pi * dp))

    def _decode_basis_flows(
        self, layer,
        weights: list[float],
    ) -> list[float]:
        """将高斯基权重解码为每个子循环的流量值（含步长离散化）。"""
        import math
        centers = self._basis_centers
        sig_s, sig_p = self._basis_sigma
        denom_s = 2.0 * sig_s * sig_s
        denom_p = 2.0 * sig_p * sig_p

        flows: list[float] = []
        for sc in layer.subcycles:
            nodes = [layer.nodes[i] for i in sc.nodes]
            sc_s = sum(n.S for n in nodes) / len(nodes)
            sc_p = sum(n.P for n in nodes) / len(nodes)
            mf = 0.0
            for w, (cs, cp) in zip(weights, centers):
                mf += w * math.exp(
                    -((sc_s - cs) ** 2) / denom_s
                    - ((sc_p - cp) ** 2) / denom_p
                )
            flows.append(mf)

        # 步长离散化
        mf_step = self._mf_step * (self._mf_hi - self._mf_lo)
        if mf_step > 0:
            flows = [round(f / mf_step) * mf_step for f in flows]
        return flows

    @property
    def bounds(self) -> list[tuple[float, float]]:
        """参数边界：t_max [500,1000] + t_min [40,100] + 分位 [0, 1] + 流量。"""
        b: list[tuple[float, float]] = []
        b.append((800.0, 1100.0))
        b.append((40.0, 200.0))
        for _ in range(self._n_t_q):
            b.append((0.0, 1.0))
        for _ in range(self._n_p_q):
            b.append((0.0, 1.0))

        mf_lo = self._mf_lo
        mf_hi = self._mf_hi
        for _ in range(self._mf_dim):
            b.append((mf_lo, mf_hi))

        return b

    def _evaluate(self, x: tuple[float, ...]) -> float:
        """解包 x → 构造系统 → 运行管线 → 返回目标值。"""
        idx = 0
        t_max = float(x[idx]); idx += 1
        t_min = float(x[idx]); idx += 1
        # T 边界离散化
        if self._tstep > 0:
            t_max = round(t_max / self._tstep) * self._tstep
            t_min = round(t_min / self._tstep) * self._tstep
        t_q = _round_and_dedup(tuple(x[idx + i] for i in range(self._n_t_q)), self._qstep, self._qmerge)
        idx += self._n_t_q
        p_q = _round_and_dedup(tuple(x[idx + i] for i in range(self._n_p_q)), self._qstep, self._qmerge)
        idx += self._n_p_q
        mf_all = x[idx:]

        # 流量离散化
        mf_step = self._mf_step * (self._mf_hi - self._mf_lo)
        if mf_step > 0:
            mf_all = [round(m / mf_step) * mf_step for m in mf_all]

        tp = ClosedCycleTPInput(
            fluid=self._base_tp.fluid,
            t_min=t_min,
            t_max=t_max,
            p_min=self._base_tp.p_min,
            p_max=self._base_tp.p_max,
            t_quantiles=t_q,
            p_quantiles=p_q,
            subcycle_mass_flow_initial=self._base_tp.subcycle_mass_flow_initial,
        )

        # 试跑拓扑获取实际子循环数
        try:
            probe = ClosedCycleLayer(tp)
        except Exception:
            return 1.0
        n_sc = len(probe.subcycles)
        if n_sc == 0:
            return 1.0

        if self._basis:
            mf_truncated = self._decode_basis_flows(probe, list(mf_all[:self._mf_dim]))
        else:
            mf_truncated = [float(mf_all[i]) for i in range(min(n_sc, len(mf_all)))]

        cfg = CycleConfig(
            input=tp,
            use_non_ideal=self._base_cfg.use_non_ideal,
            subcycle_mass_flows=mf_truncated,
            delta_T_min=self._base_cfg.delta_T_min,
            heat_method=self._base_cfg.heat_method,
        )

        sys_inp = SystemInput(
            heat_sources=self._base.heat_sources,
            cold_sources=self._base.cold_sources,
            cycles=(cfg,),
            delta_T_min=self._base.delta_T_min,
            heat_method=self._base.heat_method,
        )

        try:
            result = SystemPipeline(sys_inp).run(self._props)
        except Exception:
            return 1.0

        return self._objective_fn(result)

    def run(
        self,
        *,
        method: str = "de",
        popsize: int = 15,
        maxiter: int = 200,
        F: float = 0.8,
        CR: float = 0.9,
        sigma0: float = 0.3,
        seed: int | None = 42,
        early_stop: int = 30,
        callback: Callable[[int, int, list[float], float, int], None] | None = None,
    ) -> OptimizationResult:
        """执行优化。

        :param method: ``"de"`` = 差分进化，``"cma"`` = CMA-ES。
        :param callback: 每代末回调 ``(gen, best_x, best_val, n_evals)``。
        :returns: ``OptimizationResult``（含最优参数与完整系统结果）。
        """
        bounds = self.bounds
        if method == "cma":
            best_x, best_val, n_evals = _cma_evolution(
                self._evaluate,
                bounds,
                sigma0=sigma0,
                maxiter=maxiter,
                seed=seed,
                early_stop=early_stop,
                restarts=5,
                eval_state=self._to_eval_state() if self._n_workers > 1 else None,
                n_workers=self._n_workers,
                callback=callback,
            )
        else:
            best_x, best_val, n_evals = _differential_evolution(
                self._evaluate,
                bounds,
                popsize=popsize,
                maxiter=maxiter,
                F=F,
                CR=CR,
                seed=seed,
                early_stop=early_stop,
                eval_state=self._to_eval_state() if self._n_workers > 1 else None,
                n_workers=self._n_workers,
                callback=callback,
            )

        # 用最优参数再运行一次，获取完整 SystemResult
        best_val_final = self._evaluate(tuple(best_x))
        # 重建一次获取 result（_evaluate 内部不返回 result，需单独跑）
        sys_result = self._run_with_params(tuple(best_x))

        return OptimizationResult(
            x_opt=tuple(best_x),
            objective=best_val_final,
            system_result=sys_result,
            n_evaluations=n_evals,
            bounds=tuple(bounds),
            objective_name=self._objective_name,
            success=True,
        )

    def _run_with_params(self, x: tuple[float, ...]) -> SystemResult:
        """用指定参数运行完整系统，返回 SystemResult。"""
        idx = 0
        t_max = float(x[idx]); idx += 1
        t_min = float(x[idx]); idx += 1
        if self._tstep > 0:
            t_max = round(t_max / self._tstep) * self._tstep
            t_min = round(t_min / self._tstep) * self._tstep
        t_q = _round_and_dedup(tuple(x[idx + i] for i in range(self._n_t_q)), self._qstep, self._qmerge)
        idx += self._n_t_q
        p_q = _round_and_dedup(tuple(x[idx + i] for i in range(self._n_p_q)), self._qstep, self._qmerge)
        idx += self._n_p_q
        mf_all = x[idx:]

        mf_step = self._mf_step * (self._mf_hi - self._mf_lo)
        if mf_step > 0:
            mf_all = [round(m / mf_step) * mf_step for m in mf_all]

        tp = ClosedCycleTPInput(
            fluid=self._base_tp.fluid,
            t_min=t_min,
            t_max=t_max,
            p_min=self._base_tp.p_min,
            p_max=self._base_tp.p_max,
            t_quantiles=t_q,
            p_quantiles=p_q,
            subcycle_mass_flow_initial=self._base_tp.subcycle_mass_flow_initial,
        )

        probe = ClosedCycleLayer(tp)
        n_sc = len(probe.subcycles)

        if self._basis:
            mf_truncated = self._decode_basis_flows(probe, list(mf_all[:self._mf_dim]))
        else:
            mf_truncated = [float(mf_all[i]) for i in range(min(n_sc, len(mf_all)))]

        cfg = CycleConfig(
            input=tp,
            use_non_ideal=self._base_cfg.use_non_ideal,
            subcycle_mass_flows=mf_truncated,
            delta_T_min=self._base_cfg.delta_T_min,
            heat_method=self._base_cfg.heat_method,
        )

        sys_inp = SystemInput(
            heat_sources=self._base.heat_sources,
            cold_sources=self._base.cold_sources,
            cycles=(cfg,),
            delta_T_min=self._base.delta_T_min,
            heat_method=self._base.heat_method,
        )

        raw = SystemPipeline(sys_inp).run(self._props)
        return analyze_system_heat(raw, sys_inp, self._props)
