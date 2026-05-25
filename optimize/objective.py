"""目标函数注册表。

通过 ``@register("name")`` 装饰器注册新目标函数，在 ``Optimizer`` 中按名称调用。
"""

from __future__ import annotations

from collections.abc import Callable

from core.system import SystemResult

OBJECTIVES: dict[str, Callable[[SystemResult], float]] = {}


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
