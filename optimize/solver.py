"""差分进化求解器：单目标 DE/rand/1/bin，纯 Python 实现。"""

from __future__ import annotations

import math
import random
import warnings
from collections.abc import Callable
from typing import Any

from core.closed_cycle_layer import ClosedCycleLayer, ClosedCycleTPInput
from core.fluid_property_solver import PropertyRegistry
from core.system import CycleConfig, SystemInput, SystemPipeline, SystemResult
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
    callback: Callable[[int, int, list[float], float, int], None] | None = None,
) -> tuple[list[float], float, int]:
    """DE/rand/1/bin 单目标差分进化。

    :param callback: 每代末调用 ``callback(gen, restart, best_x, best_val, n_evals)``；``None`` 不回调。
    :returns: ``(best_x, best_val, n_evaluations)``。
    """
    rng = random.Random(seed)
    dim = len(bounds)
    lo = [b[0] for b in bounds]
    hi = [b[1] for b in bounds]

    # 随机初始化种群
    pop = [
        [rng.uniform(lo[d], hi[d]) for d in range(dim)]
        for _ in range(popsize)
    ]
    scores = [fn(tuple(p)) for p in pop]
    n_evals = len(scores)

    best_idx = min(range(popsize), key=lambda i: scores[i])
    best_x = pop[best_idx][:]
    best_val = scores[best_idx]
    stall = 0

    for gen in range(maxiter):
        for i in range(popsize):
            # 选三个不同的个体 a, b, c（不含 i）
            pool = [j for j in range(popsize) if j != i]
            a_idx, b_idx, c_idx = rng.sample(pool, 3)

            # 变异: v = a + F * (b - c)
            v = [
                pop[a_idx][d] + F * (pop[b_idx][d] - pop[c_idx][d])
                for d in range(dim)
            ]
            # clip to bounds
            for d in range(dim):
                v[d] = max(lo[d], min(hi[d], v[d]))

            # 交叉: binomial
            j_rand = rng.randrange(dim)
            u = [
                v[d] if rng.random() < CR or d == j_rand else pop[i][d]
                for d in range(dim)
            ]

            # 选择
            u_score = fn(tuple(u))
            n_evals += 1
            if u_score < scores[i]:
                pop[i] = u
                scores[i] = u_score
                if u_score < best_val:
                    best_val = u_score
                    best_x = u[:]
                    stall = 0
                else:
                    stall += 1
            else:
                stall += 1

        if callback is not None:
            callback(gen, 0, best_x, best_val, n_evals)

        if stall >= early_stop * popsize:
            break

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
    callback: Callable[[int, int, list[float], float, int], None] | None = None,
) -> tuple[list[float], float, int]:
    """CMA-ES + 多起点重启（BIPOP 风格），自动交替大小种群避免局部最优。

    :param callback: 每代末调用 ``callback(gen, restart, best_x, best_val, n_evals)``。
    :returns: ``(best_x, best_val, n_evaluations)``。
    """
    import cma

    dim = len(bounds)
    lo = [b[0] for b in bounds]
    hi = [b[1] for b in bounds]
    stds = [sigma0 * (hi[i] - lo[i]) for i in range(dim)]
    x0_init = [(lo[i] + hi[i]) / 2 for i in range(dim)]

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

        # 随机起点（restart=0 用中点，后续均匀散布）
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
            scores = [fn(tuple(s)) for s in solutions]
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

    return global_best_x, global_best_val, total_evals


def _round_and_dedup(values: tuple[float, ...], step: float) -> tuple[float, ...]:
    """舍入到 step 整数倍，去重，夹在 [step, 1-step] 内。

    ``step <= 0`` 时不做离散化，原值返回。
    """
    if step <= 0.0:
        return values
    lo, hi = step, 1.0 - step
    rounded = sorted(set(max(lo, min(hi, round(v / step) * step)) for v in values))
    return tuple(rounded)


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
        mf_bounds: tuple[float, float] = (-10.0, 50.0),
    ):
        if not base_input.cycles:
            raise ValueError("base_input 须包含至少一个 CycleConfig")

        self._base = base_input
        self._props = props
        self._base_cfg = base_input.cycles[0]
        self._base_tp = self._base_cfg.input
        self._mf_step = mf_step_fraction
        self._qstep = quantile_step
        self._mf_lo, self._mf_hi = mf_bounds

        if isinstance(objective, str):
            if objective not in OBJECTIVES:
                raise ValueError(f"未知目标 '{objective}'，已知: {list(OBJECTIVES)}")
            self._objective_fn = OBJECTIVES[objective]
            self._objective_name = objective
        else:
            self._objective_fn = objective
            self._objective_name = getattr(objective, "__name__", "custom")

        if max_subcycles is None:
            max_subcycles = self._probe_max_sc()

        self._max_sc = max_subcycles
        self._n_t_q = len(self._base_tp.t_quantiles)
        self._n_p_q = len(self._base_tp.p_quantiles)

    def _probe_max_sc(self) -> int:
        """用初始分位探测子循环数，×1.5 安全系数作为 mf 维度上限。"""
        probe = ClosedCycleLayer(self._base_tp)
        return max(1, int(len(probe.subcycles) * 1.5))

    @property
    def bounds(self) -> list[tuple[float, float]]:
        """参数边界：t_max [500,1000] + t_min [40,100] + 分位 [0.01, 0.99] + 流量。"""
        b: list[tuple[float, float]] = []
        b.append((500.0, 1000.0))
        b.append((40.0, 100.0))
        for _ in range(self._n_t_q):
            b.append((0.01, 0.99))
        for _ in range(self._n_p_q):
            b.append((0.01, 0.99))

        mf_lo = self._mf_lo
        mf_hi = self._mf_hi
        for _ in range(self._max_sc):
            b.append((mf_lo, mf_hi))

        return b

    def _evaluate(self, x: tuple[float, ...]) -> float:
        """解包 x → 构造系统 → 运行管线 → 返回目标值。"""
        idx = 0
        t_max = float(x[idx]); idx += 1
        t_min = float(x[idx]); idx += 1
        t_q = _round_and_dedup(tuple(x[idx + i] for i in range(self._n_t_q)), self._qstep)
        idx += self._n_t_q
        p_q = _round_and_dedup(tuple(x[idx + i] for i in range(self._n_p_q)), self._qstep)
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
        t_q = _round_and_dedup(tuple(x[idx + i] for i in range(self._n_t_q)), self._qstep)
        idx += self._n_t_q
        p_q = _round_and_dedup(tuple(x[idx + i] for i in range(self._n_p_q)), self._qstep)
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

        return SystemPipeline(sys_inp).run(self._props)
