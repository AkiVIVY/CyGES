"""
闭式循环层：由 ``ClosedCycleTPInput``（含工质）创建物性求解器，并提供 TP 拓扑分析。

TP 拓扑：由温压上下限与分位生成轴，一级为 TP 网格上的物性节点；二级在同一压力轴上对
一级状态做等熵延伸，仅筛温度是否仍在边界内。调用 ``ClosedCycleLayer.analyze_topology()`` 后，
全部节点写入 ``nodes``（``dict[int, Node]``，键为 ``index``）：``parent is None`` 为一级，
``parent`` 为一级节点 ``index`` 时为二级。不要求节点顺序。

``nodes`` 与 ``mechanical_edges`` / ``heat_edges`` 为 **拓扑快照**；边上 ``mass_flow`` 初始为
``None``，方向由 ``tail → head``（有向边尾端 / 头端节点序号）表示。后续可在此基础上再生成用于实际计算的边/流。

单位约定（字段名不缀单位）：T[K]、P[kPa]、H[kJ/kg]、S[kJ/(kg·K)]。
分位序列由调用方保证互不重复且落在合理区间，本模块不对分位做额外校验。
一级节点 ``index`` 从 0 起连续编号；二级节点 ``index`` 紧接一级最大编号之后连续编号。
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Literal, Sequence

from core.fluid_property_solver import CoolPropFluidPropertySolver, FluidPropertySolver


@dataclass(frozen=True)
class ClosedCycleTPInput:
    """闭式循环层拓扑输入：含工质，可单独作为创建该层的完整信息。"""

    fluid: str
    """工质标识，与物性后端（如 CoolProp）流体名一致。"""
    t_min: float
    t_max: float
    p_min: float
    p_max: float
    t_quantiles: Sequence[float] = field(default_factory=tuple)
    """[0,1] 分位，在 [t_min, t_max] 内线性插值；须互不重复（与端点亦不重复）。"""
    p_quantiles: Sequence[float] = field(default_factory=tuple)
    """[0,1] 分位，在 [p_min, p_max] 内线性插值；须互不重复（与端点亦不重复）。"""
    max_mass_flow: float | None = None
    """可选：全层质量流量上界等，拓扑阶段仅占位，供后续设备或边约束使用。"""


@dataclass(frozen=True)
class Node:
    """
    拓扑状态节点；一级（TP 网格）与二级（等熵延伸）共用本结构，用 ``parent`` 区分。

    - ``index``：全局唯一序号，亦作 ``nodes`` 字典键。
    - ``T, P, H, S``：状态量（项目约定单位）。
    - ``parent``：一级为 ``None``；二级为所属一级节点的 ``index``。
    """

    index: int
    T: float
    P: float
    H: float
    S: float
    parent: int | None = None


@dataclass
class Edge:
    """
    拓扑快照中的有向边（非物性计算边）。

    - ``kind``：``"mechanical"`` 为等熵链上压力离散相邻；``"heat"`` 为等压线上温度离散相邻。
    - ``tail`` / ``head``：端点节点序号，方向为尾 → 头。
    - ``mass_flow``：质量流 [kg/s] 等，初始 ``None``，由后续流程填写。
    """

    kind: Literal["mechanical", "heat"]
    tail: int
    head: int
    mass_flow: float | None = None


def build_axis(min_v: float, max_v: float, quantiles: Sequence[float]) -> list[float]:
    """由单轴上下限与分位生成一维采样序列（两端与各分位内插点排序，不做去重）。"""
    span = max_v - min_v
    pts = [min_v, max_v, *[min_v + q * span for q in quantiles]]
    return sorted(pts)


def build_grid_nodes(solver: FluidPropertySolver, inp: ClosedCycleTPInput) -> list[Node]:
    """
    在 T、P 轴笛卡尔积上逐点做 ``TP`` 物性查询，成功则生成一级 ``Node``。
    包线外或发散点跳过，不占位，故 ``index`` 在成功点上连续递增。
    """
    t_list = build_axis(inp.t_min, inp.t_max, inp.t_quantiles)
    p_list = build_axis(inp.p_min, inp.p_max, inp.p_quantiles)
    nodes: list[Node] = []
    idx = 0
    # 外层 T、内层 P：与历史行为一致；换热边在后续按 P 分桶统一生成，不依赖此处遍历顺序
    for Tk in t_list:
        for Pk in p_list:
            try:
                st = solver.state("TP", Tk, Pk)
            except Exception:
                continue
            nodes.append(Node(index=idx, T=Tk, P=Pk, H=st["H"], S=st["S"], parent=None))
            idx += 1
    return nodes


def expand_isentropic_nodes(
    solver: FluidPropertySolver,
    inp: ClosedCycleTPInput,
    grid: Sequence[Node],
) -> tuple[list[Node], list[Edge]]:
    """
    对每个一级节点，沿压力轴 ``p_axis`` 做 ``PS`` 等熵延伸得到二级 ``Node``；
    再将「该一级 + 本支全部二级」按 ``P`` 排序，在排序相邻点之间生成机械边（含一级↔二级、二级↔二级）。
    """
    p_axis = build_axis(inp.p_min, inp.p_max, inp.p_quantiles)
    out: list[Node] = []
    mechanical: list[Edge] = []
    # 已有 (T,P) 状态，避免二级与一级或其它支重复
    seen: set[tuple[float, float]] = {(n.T, n.P) for n in grid}
    idx = len(grid)

    for prim in grid:
        batch: list[Node] = []
        for p2 in p_axis:
            # 与一级同压则跳过（等熵线上该点即自身）
            if abs(p2 - prim.P) < 1e-9 * max(1.0, abs(p2)):
                continue
            try:
                st = solver.state("PS", p2, prim.S)
                t2 = st["T"]
            except Exception:
                continue
            # 二级仅校核温度窗口；压力已在轴上
            if not (inp.t_min <= t2 <= inp.t_max):
                continue
            key = (t2, p2)
            if key in seen:
                continue
            seen.add(key)
            node = Node(
                index=idx,
                T=t2,
                P=p2,
                H=st["H"],
                S=prim.S,
                parent=prim.index,
            )
            out.append(node)
            batch.append(node)
            idx += 1

        # 同熵线离散链：按 P（及 index 平局）排序后，相邻即「机械」相邻
        chain = sorted([prim, *batch], key=lambda n: (n.P, n.index))
        for k in range(len(chain) - 1):
            a, b = chain[k], chain[k + 1]
            mechanical.append(Edge(kind="mechanical", tail=a.index, head=b.index, mass_flow=None))

    return out, mechanical


def build_heat_edges(nodes: dict[int, Node]) -> list[Edge]:
    """
    对全部节点（一级 + 二级）按等压分桶，桶内按温度排序；
    排序相邻节点对之间生成换热边，方向沿温度增大（``tail`` 低温侧，``head`` 高温侧）。
    """
    by_p: dict[float, list[Node]] = defaultdict(list)
    for n in nodes.values():
        by_p[n.P].append(n)
    heat: list[Edge] = []
    for lst in by_p.values():
        lst.sort(key=lambda n: (n.T, n.index))
        for k in range(len(lst) - 1):
            a, b = lst[k], lst[k + 1]
            heat.append(Edge(kind="heat", tail=a.index, head=b.index, mass_flow=None))
    return heat


class ClosedCycleLayer:
    """
    闭式循环层。仅依赖 ``ClosedCycleTPInput`` 即可持有工质并创建默认物性后端；
    也可注入 ``properties`` 以便测试或替换实现。

    ``analyze_topology()`` 填充 ``nodes`` 与快照边 ``mechanical_edges``、``heat_edges``。
    """

    def __init__(self, inp: ClosedCycleTPInput, *, properties: FluidPropertySolver | None = None) -> None:
        self.input = inp
        self.fluid = inp.fluid
        self.properties = properties if properties is not None else CoolPropFluidPropertySolver(inp.fluid)
        # 分析前为空；``analyze_topology`` 一次性写入基准拓扑
        self.nodes: dict[int, Node] = {}
        self.mechanical_edges: list[Edge] = []
        self.heat_edges: list[Edge] = []

    def analyze_topology(self) -> None:
        """先一级网格、再二级与机械边，最后汇总 ``nodes`` 并生成换热边；结果写入本层属性。"""
        grid = build_grid_nodes(self.properties, self.input)
        secondary, mechanical = expand_isentropic_nodes(self.properties, self.input, grid)
        nodes: dict[int, Node] = {}
        for n in grid:
            nodes[n.index] = n
        for n in secondary:
            nodes[n.index] = n
        self.nodes = nodes
        self.mechanical_edges = mechanical
        self.heat_edges = build_heat_edges(nodes)
