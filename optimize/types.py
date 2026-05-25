"""优化器数据类型。"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

from core.system import SystemResult


@dataclass(frozen=True)
class OptimizationResult:
    """一次优化运行的冻结结果。"""

    x_opt: tuple[float, ...]
    """最优参数向量（分位点 + 子循环流量）。"""
    objective: float
    """最优目标函数值。"""
    system_result: SystemResult
    """最优参数下的完整系统结果。"""
    n_evaluations: int
    """总评估次数。"""
    bounds: tuple[tuple[float, float], ...]
    """参数边界。"""
    objective_name: str
    """使用的目标函数名称。"""
    success: bool = True
    """是否正常收敛。"""
