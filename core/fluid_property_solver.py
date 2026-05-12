"""
闭式循环层物性求解器：按工质实例化，对外统一单位
T[K]、P[kPa]、H[kJ/kg]、S[kJ/(kg·K)]；内部 CoolProp 使用 SI。

通过单一 ``state`` 接口指定已知量组合（"HP" / "TP" / "HS" / "PS"），
返回包含 T、P、H、S 的字典。

物性计算使用 ``AbstractState("HEOS", substance)`` 并长期复用同一实例，
相对反复 ``PropsSI`` 有利于降低每次状态更新的开销。
"""

from __future__ import annotations

from typing import Any, Literal, Protocol, TypedDict, runtime_checkable

import CoolProp.CoolProp as CP


# 已知强度量组合：首字母为 CoolProp 常用输入对（如 TP = 温度+压力）
PropertyPair = Literal["HP", "TP", "HS", "PS"]


class ThermoStateTPHS(TypedDict):
    """
    平衡态热力学状态，四个键对应 T、P、H、S（单位与全项目约定一致）。

    - T：温度 [K]
    - P：压力 [kPa]
    - H：比焓 [kJ/kg]
    - S：比熵 [kJ/(kg·K)]
    """

    T: float
    P: float
    H: float
    S: float


def _tphs(t_k: float, p_kpa: float, h_kjkg: float, s_kj_per_kgk: float) -> ThermoStateTPHS:
    """将四个标量打包为对外统一返回的字典结构。"""
    return ThermoStateTPHS(T=t_k, P=p_kpa, H=h_kjkg, S=s_kj_per_kgk)


@runtime_checkable
class FluidPropertySolver(Protocol):
    """物性求解抽象：由具体后端（如 CoolProp）实现 ``state``。"""

    @property
    def fluid(self) -> str:
        """当前求解器绑定的工质标识（与 CoolProp 流体名一致）。"""
        ...

    def state(self, pair: PropertyPair, x: float, y: float) -> ThermoStateTPHS:
        """
        按 ``pair`` 解释 ``(x, y)`` 为一对已知量，补全 T、P、H、S。

        - ``"HP"``：x = H [kJ/kg]，y = P [kPa]
        - ``"TP"``：x = T [K]，y = P [kPa]
        - ``"HS"``：x = H [kJ/kg]，y = S [kJ/(kg·K)]
        - ``"PS"``：x = P [kPa]，y = S [kJ/(kg·K)]
        """
        ...


class CoolPropFluidPropertySolver:
    """
    基于 CoolProp 的工质物性实现。
    每个闭式循环层按工质创建一个实例，避免在调用处重复传 fluid。
    内部持有一个 ``AbstractState("HEOS", fluid)``，多次 ``state`` 调用时复用。
    """

    __slots__ = ("_fluid", "_as")

    def __init__(self, fluid: str) -> None:
        if not fluid or not str(fluid).strip():
            raise ValueError("工质名称 fluid 不能为空")
        self._fluid = str(fluid)
        # 首次调用 state 时再创建 AbstractState，避免仅构造求解器就初始化物性后端
        self._as: Any = None

    @property
    def fluid(self) -> str:
        return self._fluid

    def _abstract_state(self) -> Any:
        """懒加载 Helmholtz 状态对象，全实例共享、反复 update。"""
        if self._as is None:
            self._as = CP.AbstractState("HEOS", self._fluid)
        return self._as

    def state(self, pair: PropertyPair, x: float, y: float) -> ThermoStateTPHS:
        """
        根据已知对更新 ``AbstractState`` 并读取 T、P、H、S（SI 读数再换为约定单位）。
        多相区或包线外可能由 CoolProp 抛错，由调用方捕获处理。
        """
        AS = self._abstract_state()

        # 已知 T、P：直接求 H、S（update 顺序为 P, T，单位 Pa、K）
        if pair == "TP":
            t_k, p_kpa = x, y
            p_pa = p_kpa * 1e3  # kPa → Pa
            AS.update(CP.PT_INPUTS, p_pa, t_k)
            return _tphs(float(AS.T()), p_kpa, float(AS.hmass()) / 1e3, float(AS.smass()) / 1e3)

        # 已知 P、S（如等熵过程指定终点压力）：求 T、H
        if pair == "PS":
            p_kpa, s_kj_per_kgk = x, y
            p_pa = p_kpa * 1e3
            s_si = s_kj_per_kgk * 1e3  # kJ/(kg·K) → J/(kg·K)
            AS.update(CP.PSmass_INPUTS, p_pa, s_si)
            return _tphs(float(AS.T()), p_kpa, float(AS.hmass()) / 1e3, s_kj_per_kgk)

        # 已知 H、P：求 T、S（update 顺序为 h [J/kg]、p [Pa]）
        if pair == "HP":
            h_kjkg, p_kpa = x, y
            h_si = h_kjkg * 1e3
            p_pa = p_kpa * 1e3
            AS.update(CP.HmassP_INPUTS, h_si, p_pa)
            return _tphs(float(AS.T()), p_kpa, h_kjkg, float(AS.smass()) / 1e3)

        # 已知 H、S：求 T、P（update 顺序为 h、s，单位 J/kg、J/(kg·K)）
        if pair == "HS":
            h_kjkg, s_kj_per_kgk = x, y
            h_si = h_kjkg * 1e3
            s_si = s_kj_per_kgk * 1e3
            AS.update(CP.HmassSmass_INPUTS, h_si, s_si)
            p_kpa = float(AS.p()) / 1e3
            return _tphs(float(AS.T()), p_kpa, h_kjkg, s_kj_per_kgk)

        # 理论上 Literal 已约束；保留分支便于扩展或运行时校验
        raise ValueError(f"不支持的 pair：{pair!r}，应为 HP、TP、HS、PS 之一")
