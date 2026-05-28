"""目标函数注册表。

通过 ``@register("name")`` 装饰器注册新目标函数，在 ``Optimizer`` 中按名称调用。
"""

from __future__ import annotations

from collections.abc import Callable

from core.cycle_performance import ProcessCategory
from core.system import SystemResult

OBJECTIVES: dict[str, Callable[[SystemResult], float]] = {}

# 参数化目标的模块级配置（Optimizer 创建 / 每次评估前写入)
_HX_DT_MIN: float = 10.0
_HX_MAX_GROUP_SIZE: int | None = None


def register(name: str):
    """装饰器：将目标函数注册到 ``OBJECTIVES`` 表。

    ::

        @register("min_max_utility")
        def min_max_utility_ratio(result: SystemResult) -> float:
            ...
    """

    def decorator(fn: Callable[[SystemResult], float]):
        OBJECTIVES[name] = fn
        return fn

    return decorator


# ────────────────────────────────────────────
# 内置目标函数
# ────────────────────────────────────────────


@register("min_max_utility")
def min_max_utility_ratio(result: SystemResult) -> float:
    """最小化最大公用工程比例——闭环目标。

    取 ``system_pinch.hot_utility_ratio`` 与 ``cold_utility_ratio`` 中的较大值；
    值为 0 表示完全闭环，无需外部冷热源。
    """
    sp = result.system_pinch
    if sp is None:
        return 1.0
    return max(sp.hot_utility_ratio, sp.cold_utility_ratio)


@register("heat_balance")
def heat_balance_ratio(result: SystemResult) -> float:
    """最小化系统总放热与总吸热的不平衡比——纯热量守恒目标。

    .. math::
        obj = \\frac{|\\Sigma Q_\\text{rej} - \\Sigma Q_\\text{abs}|}
                    {\\max(\\Sigma Q_\\text{rej}, \\Sigma Q_\\text{abs})}

    值为 0 表示全系统放热与吸热完美匹配，无需外部冷热源；
    不依赖夹点分析。
    """
    all_recs: list = list(result.heat_source_records) + list(result.cold_source_records)
    for report in result.cycle_reports:
        for _, rec in report.by_edge:
            if rec.kind == "heat" and rec.power_rate is not None:
                all_recs.append(rec)

    Q_rej = sum(abs(float(r.power_rate)) for r in all_recs
                if r.category == ProcessCategory.HEAT_REJECTION)
    Q_abs = sum(abs(float(r.power_rate)) for r in all_recs
                if r.category == ProcessCategory.HEAT_ABSORPTION)
    max_q = max(Q_rej, Q_abs)
    if max_q < 1e-12:
        return 1.0
    return abs(Q_rej - Q_abs) / max_q


@register("hx_unmatched")
def hx_unmatched_ratio(result: SystemResult) -> float:
    """最小化换热器匹配的未匹配功率比——全局平衡 + 逐过程配对目标。

    内部调用 ``match_heat_exchanger_groups`` 进行星形逆流匹配，
    温差约束由模块变量 ``_HX_DT_MIN`` 控制（默认 10 K）。

    .. math::
        obj = \\frac{\\text{total\\_unmatched} \\times (N+1)}{\\Sigma |Q|}

    其中 *N* 为未匹配换热边总数（热+冷）。``N+1`` 倍的惩罚项
    使优化器优先消除碎片化未匹配边，而非仅最小化未匹配功率总和。

    值为 0 表示全部换热过程均可配对成热平衡换热器。
    """
    from core.heat_exchanger import match_heat_exchanger_groups

    all_hots: list = list(result.heat_source_records)
    all_colds: list = list(result.cold_source_records)
    for report in result.cycle_reports:
        for _, rec in report.by_edge:
            if rec.kind == "heat" and rec.power_rate is not None:
                if rec.category == ProcessCategory.HEAT_REJECTION:
                    all_hots.append(rec)
                else:
                    all_colds.append(rec)

    if not all_hots and not all_colds:
        return 1.0

    hx = match_heat_exchanger_groups(all_hots, all_colds, dT_min=_HX_DT_MIN,
                                      max_group_size=_HX_MAX_GROUP_SIZE)
    total_q = sum(abs(float(r.power_rate)) for r in all_hots + all_colds if r.power_rate)
    if total_q < 1e-12:
        return 1.0
    return hx.total_unmatched / total_q
