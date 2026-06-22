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
        """懒加载 Helmholtz 状态对象，本实例内复用，反复 update。"""
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


class InterpolatingHeliumSolver:
    """He 插值物性求解器：预计算 (T,P) 网格 → 双线性插值 / 二分查找。

    替代 CoolProp 迭代求根，精度与 CoolProp 一致（插值误差 < ΔT 步长对应偏差）。
    支持 TP、PS、HP 三种输入对；HS 未实现（代码库中从未调用）。

    用法::

        solver = InterpolatingHeliumSolver(dT=10.0, dP=200.0)
        st = solver.state("TP", 500.0, 3000.0)  # → ThermoStateTPHS
    """

    __slots__ = ("_dT", "_dP", "_T_min", "_T_max", "_P_min", "_P_max",
                 "_T_grid", "_P_grid", "_H_grid", "_S_grid", "_built")

    _T_DEFAULT = (40.0, 1200.0, 10.0)
    _P_DEFAULT = (200.0, 10000.0, 200.0)

    @property
    def fluid(self) -> str:
        return "He"

    def __init__(
        self,
        dT: float = 10.0,
        dP: float = 200.0,
        *,
        T_min: float | None = None,
        T_max: float | None = None,
        P_min: float | None = None,
        P_max: float | None = None,
    ) -> None:
        self._dT = dT
        self._dP = dP
        self._T_min = T_min if T_min is not None else self._T_DEFAULT[0]
        self._T_max = T_max if T_max is not None else self._T_DEFAULT[1]
        self._P_min = P_min if P_min is not None else self._P_DEFAULT[0]
        self._P_max = P_max if P_max is not None else self._P_DEFAULT[1]
        self._T_grid: list[float] = []
        self._P_grid: list[float] = []
        self._H_grid: list[list[float]] = []
        self._S_grid: list[list[float]] = []
        self._built = False

    def _build(self) -> None:
        if self._built:
            return
        import math
        as_cp = CP.AbstractState("HEOS", "He")
        nT = int((self._T_max - self._T_min) / self._dT) + 1
        nP = int((self._P_max - self._P_min) / self._dP) + 1
        T_vals = [self._T_min + i * self._dT for i in range(nT)]
        P_vals = [self._P_min + j * self._dP for j in range(nP)]
        self._H_grid = [[0.0] * nP for _ in range(nT)]
        self._S_grid = [[0.0] * nP for _ in range(nT)]
        for i, T in enumerate(T_vals):
            for j, P_kpa in enumerate(P_vals):
                try:
                    as_cp.update(CP.PT_INPUTS, P_kpa * 1e3, T)
                    self._H_grid[i][j] = float(as_cp.hmass()) / 1e3
                    self._S_grid[i][j] = float(as_cp.smass()) / 1e3
                except Exception:
                    self._H_grid[i][j] = float("nan")
                    self._S_grid[i][j] = float("nan")
        self._T_grid = T_vals
        self._P_grid = P_vals
        self._built = True

    # ── 双线性插值 ───────────────────────────────────────
    @staticmethod
    def _lerp(x: float, x0: float, x1: float, v0: float, v1: float) -> float:
        if x1 == x0:
            return v0
        return v0 + (v1 - v0) * (x - x0) / (x1 - x0)

    def _bilerp(self, T: float, P: float) -> tuple[float, float]:

        nT, nP = len(self._T_grid), len(self._P_grid)
        i = int((T - self._T_min) / self._dT)
        j = int((P - self._P_min) / self._dP)
        if i < 0:
            i = 0
        if j < 0:
            j = 0
        if i >= nT - 1:
            i = nT - 2
        if j >= nP - 1:
            j = nP - 2
        T0, T1 = self._T_grid[i], self._T_grid[i + 1]
        P0, P1 = self._P_grid[j], self._P_grid[j + 1]
        t_frac = (T - T0) / (T1 - T0) if T1 != T0 else 0.0
        p_frac = (P - P0) / (P1 - P0) if P1 != P0 else 0.0
        H00, H01 = self._H_grid[i][j], self._H_grid[i][j + 1]
        H10, H11 = self._H_grid[i + 1][j], self._H_grid[i + 1][j + 1]
        S00, S01 = self._S_grid[i][j], self._S_grid[i][j + 1]
        S10, S11 = self._S_grid[i + 1][j], self._S_grid[i + 1][j + 1]
        H = (1 - t_frac) * (1 - p_frac) * H00 + (1 - t_frac) * p_frac * H01 + \
            t_frac * (1 - p_frac) * H10 + t_frac * p_frac * H11
        S = (1 - t_frac) * (1 - p_frac) * S00 + (1 - t_frac) * p_frac * S01 + \
            t_frac * (1 - p_frac) * S10 + t_frac * p_frac * S11
        return H, S

    # ── 等压线二分查找（S/H → T） ─────────────────────────
    def _find_T_by_S(self, P: float, S_target: float) -> float:
        nT, nP = len(self._T_grid), len(self._P_grid)
        j = int((P - self._P_min) / self._dP)
        if j < 0:
            j = 0
        if j >= nP - 1:
            j = nP - 1
        P0, P1 = self._P_grid[j], self._P_grid[j + 1]

        def _T_on_isobar(jj: int, s: float) -> float:
            lo, hi = 0, nT - 1
            while lo < hi:
                mid = (lo + hi) // 2
                if self._S_grid[mid][jj] < s:
                    lo = mid + 1
                else:
                    hi = mid
            if lo == 0:
                lo = 1
            if lo >= nT:
                lo = nT - 1
            t0, t1 = self._T_grid[lo - 1], self._T_grid[lo]
            s0, s1 = self._S_grid[lo - 1][jj], self._S_grid[lo][jj]
            if abs(s1 - s0) < 1e-15:
                return (t0 + t1) / 2.0
            return t0 + (t1 - t0) * (s - s0) / (s1 - s0)

        T_j0 = _T_on_isobar(j, S_target)
        T_j1 = _T_on_isobar(j + 1, S_target)
        p_frac = (P - P0) / (P1 - P0) if P1 != P0 else 0.0
        return T_j0 + (T_j1 - T_j0) * p_frac

    def _find_T_by_H(self, P: float, H_target: float) -> float:
        nT, nP = len(self._T_grid), len(self._P_grid)
        j = int((P - self._P_min) / self._dP)
        if j < 0:
            j = 0
        if j >= nP - 1:
            j = nP - 1
        P0, P1 = self._P_grid[j], self._P_grid[j + 1]

        def _T_on_isobar(jj: int, h: float) -> float:
            lo, hi = 0, nT - 1
            while lo < hi:
                mid = (lo + hi) // 2
                if self._H_grid[mid][jj] < h:
                    lo = mid + 1
                else:
                    hi = mid
            if lo == 0:
                lo = 1
            if lo >= nT:
                lo = nT - 1
            t0, t1 = self._T_grid[lo - 1], self._T_grid[lo]
            h0, h1 = self._H_grid[lo - 1][jj], self._H_grid[lo][jj]
            if abs(h1 - h0) < 1e-15:
                return (t0 + t1) / 2.0
            return t0 + (t1 - t0) * (h - h0) / (h1 - h0)

        T_j0 = _T_on_isobar(j, H_target)
        T_j1 = _T_on_isobar(j + 1, H_target)
        p_frac = (P - P0) / (P1 - P0) if P1 != P0 else 0.0
        return T_j0 + (T_j1 - T_j0) * p_frac

    # ── 主接口 ────────────────────────────────────────────
    def state(self, pair: PropertyPair, x: float, y: float) -> ThermoStateTPHS:
        self._build()

        if pair == "TP":
            T, P = x, y
            H, S = self._bilerp(T, P)
            return _tphs(T, P, H, S)

        if pair == "PS":
            P, S_target = x, y
            T = self._find_T_by_S(P, S_target)
            H, _ = self._bilerp(T, P)
            return _tphs(T, P, H, S_target)

        if pair == "HP":
            H_target, P = x, y
            T = self._find_T_by_H(P, H_target)
            _, S = self._bilerp(T, P)
            return _tphs(T, P, H_target, S)

        if pair == "HS":
            raise NotImplementedError("InterpolatingHeliumSolver 不支持 HS 输入对")

        raise ValueError(f"不支持的 pair：{pair!r}，应为 HP、TP、HS、PS 之一")


class PropertyRegistry:
    """多工质物性注册表，按 fluid 缓存 ``CoolPropFluidPropertySolver``。

    统一入口，替代 ``ThermoLookup`` 与 ``EnthalpyLookup`` 类型别名：
    ``__call__`` 用于完整状态查询（convert_sources），
    ``enthalpy`` 用于焓值查询（T-Q 构建）。

    用法::

        props = PropertyRegistry()
        state = props("He", "TP", 300, 101.325)       # → ThermoStateTPHS
        h = props.enthalpy("He", 300, 101.325)         # → float
    """

    __slots__ = ("_solvers", "_cache")

    def __init__(self) -> None:
        self._solvers: dict[str, CoolPropFluidPropertySolver] = {}
        self._cache: dict[tuple[Any, ...], ThermoStateTPHS] = {}

    def __call__(self, fluid: str, pair: PropertyPair, x: float, y: float) -> ThermoStateTPHS:
        """``(fluid, pair, x, y) → ThermoStateTPHS``，匹配原 ``ThermoLookup`` 签名。"""
        key = (fluid, pair, round(x, 4), round(y, 4))
        cached = self._cache.get(key)
        if cached is not None:
            return cached
        if fluid not in self._solvers:
            self._solvers[fluid] = CoolPropFluidPropertySolver(fluid)
        val = self._solvers[fluid].state(pair, x, y)  # type: ignore[arg-type]
        self._cache[key] = val
        return val

    def enthalpy(self, fluid: str, T: float, P: float) -> float:
        """``(fluid, T, P) → H``，匹配原 ``EnthalpyLookup`` 签名。"""
        return self(fluid, "TP", T, P)["H"]
