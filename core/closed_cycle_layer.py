"""
闭式循环层：由 ``ClosedCycleTPInput``（含工质）创建物性求解器，并提供拓扑分析。

单位约定（字段名不缀单位）：T[K]、P[kPa]、H[kJ/kg]、S[kJ/(kg·K)]。
分位序列由调用方保证互不重复且落在合理区间，本模块不对分位做额外校验。

节点和边拓扑由 ``build_node_edge_topology`` 一次完成：
1.温压轴与分位生成一级节点；
2.在同一压力轴上对一级节点做等熵延伸，生成二级节点，并进一步生成机械边；
3.汇总节点后生成换热边；
4.机械边和换热边按生成顺序编号，写入 ``edges`` 字典，键名 ``M*`` / ``H*``；
5.将边键写入节点，具体规则参考PS平面约定。

build_node_edge_topology生成的节点和边为拓扑快照，后续具体计算时可在此基础上再生成用于实际计算的边/流。

最小 **子循环**（``SubCycle``）为沿 ``edge_up`` → ``edge_right`` → ``edge_down`` → ``edge_left`` 闭合的 **4 节点、4 边** 单元；
节点顺序为顺时针 **左下 → 左上 → 右上 → 右下**，边顺序为 **左、上、右、下**。
``SubCycle`` 含子循环 **单一** 质量流 ``mass_flow``（单位同 ``Edge.mass_flow``，可为负，初值 ``None``），供后续求解填写；该类非 ``frozen`` 以便就地更新。

由 ``build_subcycles`` 在拓扑快照上枚举（四节点互异；``frozenset`` 去重），**不再**对 ``P``/``S`` 几何一致性做额外校验。

一级节点 ``index`` 从 0 起连续编号；二级节点 ``index`` 紧接一级最大编号之后连续编号。

边上 ``mass_flow`` 初始为``None``，方向由 ``tail → head``（有向边尾端 / 头端节点序号）表示。

PS 平面约定：``P`` 由小到大为自下而上，``S`` 由小到大为自左而右；
有向边在生成时即取向为仅 **向上**（``P`` 增大或不变）且 **向右**（``S`` 增大或不变），
即 ``nodes[tail].P <= nodes[head].P`` 且 ``nodes[tail].S <= nodes[head].S``。

全边生成后，按 ``kind`` 将边键写入端点：机械边 ``tail`` 记 ``edge_up``、``head`` 记 ``edge_down``；
换热边 ``tail`` 记 ``edge_right``、``head`` 记 ``edge_left``（值为 ``edges`` 字典键，无邻边为 ``None``）。
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field, replace
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
    - ``edge_up`` / ``edge_down`` / ``edge_left`` / ``edge_right``：PS 平面上邻接边在 ``edges`` 中的键；拓扑生成末段由 ``tail``/``head`` 与边 ``kind`` 填入，初值为 ``None``。
    """

    index: int
    T: float
    P: float
    H: float
    S: float
    parent: int | None = None
    edge_up: str | None = None
    edge_down: str | None = None
    edge_left: str | None = None
    edge_right: str | None = None


@dataclass
class Edge:
    """
    拓扑快照中的有向边（非物性计算边）。

    - ``kind``：``"mechanical"`` 为等熵链上压力离散相邻；``"heat"`` 为等压线上熵 ``S`` 离散相邻。
    - ``tail`` / ``head``：端点节点序号，方向为尾 → 头。
    - ``mass_flow``：质量流 [kg/s] 等，初始 ``None``，由后续流程填写。
    """

    kind: Literal["mechanical", "heat"]
    tail: int
    head: int
    mass_flow: float | None = None


@dataclass
class SubCycle:
    """
    最小子循环（4 节点、4 边）。非 ``frozen``，便于后续写入流量等。

    - ``nodes``：``(左下, 左上, 右上, 右下)`` 顺时针，各元素为 ``Node.index``。
    - ``edges``：``(左, 上, 右, 下)``，各为 ``edges`` 字典键——左/右须为机械边，上/下须为换热边（由 ``build_subcycles`` 走边保证）。
    - ``mass_flow``：子循环质量流 [kg/s] 等标量，初值 ``None``；可为负（与约定正方向相反）。
    """

    nodes: tuple[int, int, int, int]
    edges: tuple[str, str, str, str]
    mass_flow: float | None = None


def _attach_edges_to_nodes_ps(nodes: dict[int, Node], edges: dict[str, Edge]) -> dict[int, Node]:
    """
    机械边视为沿 **P** 向上：``tail`` 的 ``edge_up``、``head`` 的 ``edge_down``。
    换热边视为沿 **S** 向右：``tail`` 的 ``edge_right``、``head`` 的 ``edge_left``。
    """
    slot: dict[int, dict[str, str | None]] = {
        i: {"up": None, "down": None, "left": None, "right": None} for i in nodes
    }

    def put(ni: int, side: str, ek: str) -> None:
        cur = slot[ni][side]
        if cur is not None:
            raise RuntimeError(f"node {ni} duplicate {side} edge: {cur!r} and {ek!r}")
        slot[ni][side] = ek

    for ek, e in edges.items():
        if e.kind == "mechanical":
            put(e.tail, "up", ek)
            put(e.head, "down", ek)
        else:
            put(e.tail, "right", ek)
            put(e.head, "left", ek)

    return {
        i: replace(
            n,
            edge_up=slot[i]["up"],
            edge_down=slot[i]["down"],
            edge_left=slot[i]["left"],
            edge_right=slot[i]["right"],
        )
        for i, n in nodes.items()
    }


def build_axis(min_v: float, max_v: float, quantiles: Sequence[float]) -> list[float]:
    """由单轴上下限与分位生成一维采样序列（两端与各分位内插点排序，不做去重）。"""
    span = max_v - min_v
    pts = [min_v, max_v, *[min_v + q * span for q in quantiles]]
    return sorted(pts)


def build_node_edge_topology(
    solver: FluidPropertySolver,
    inp: ClosedCycleTPInput,
) -> tuple[dict[int, Node], dict[str, Edge]]:
    """
    节点与边拓扑的一次性生成：一级 TP 网格、等熵二级与机械边（键 ``M*``）、
    全节点等压分桶换热边（键 ``H*``），边方向满足 PS 向上/向右约定；最后将边键写回各 ``Node`` 的 ``edge_*`` 字段。
    """

    def oriented_edge(
        kind: Literal["mechanical", "heat"],
        na: Node,
        nb: Node,
        *,
        rtol_p: float = 1e-9,
        rtol_s: float = 1e-12,
    ) -> Edge:
        eps_p = rtol_p * max(1.0, abs(na.P), abs(nb.P))
        eps_s = rtol_s * max(1.0, abs(na.S), abs(nb.S))
        a_sw = na.P <= nb.P + eps_p and na.S <= nb.S + eps_s
        b_sw = nb.P <= na.P + eps_p and nb.S <= na.S + eps_s
        if a_sw and not b_sw:
            return Edge(kind=kind, tail=na.index, head=nb.index, mass_flow=None)
        if b_sw and not a_sw:
            return Edge(kind=kind, tail=nb.index, head=na.index, mass_flow=None)
        if na.index <= nb.index:
            return Edge(kind=kind, tail=na.index, head=nb.index, mass_flow=None)
        return Edge(kind=kind, tail=nb.index, head=na.index, mass_flow=None)

    # —— 一级 TP 网格 ——
    t_list = build_axis(inp.t_min, inp.t_max, inp.t_quantiles)
    p_list = build_axis(inp.p_min, inp.p_max, inp.p_quantiles)
    grid: list[Node] = []
    idx = 0
    for Tk in t_list:
        for Pk in p_list:
            try:
                st = solver.state("TP", Tk, Pk)
            except Exception:
                continue
            grid.append(Node(index=idx, T=Tk, P=Pk, H=st["H"], S=st["S"], parent=None))
            idx += 1

    # —— 二级等熵 + 机械边 ——
    p_axis = build_axis(inp.p_min, inp.p_max, inp.p_quantiles)
    secondary: list[Node] = []
    mechanical: dict[str, Edge] = {}
    m = 0
    seen: set[tuple[float, float]] = {(n.T, n.P) for n in grid}
    idx = len(grid)

    for prim in grid:
        batch: list[Node] = []
        for p2 in p_axis:
            if abs(p2 - prim.P) < 1e-9 * max(1.0, abs(p2)):
                continue
            try:
                st = solver.state("PS", p2, prim.S)
                t2 = st["T"]
            except Exception:
                continue
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
            secondary.append(node)
            batch.append(node)
            idx += 1

        chain = sorted([prim, *batch], key=lambda n: (n.P, n.index))
        for k in range(len(chain) - 1):
            a, b = chain[k], chain[k + 1]
            m += 1
            mechanical[f"M{m}"] = oriented_edge("mechanical", a, b)

    nodes: dict[int, Node] = {}
    for n in grid:
        nodes[n.index] = n
    for n in secondary:
        nodes[n.index] = n

    # —— 换热边 ——
    by_p: dict[float, list[Node]] = defaultdict(list)
    for n in nodes.values():
        by_p[n.P].append(n)
    heat: dict[str, Edge] = {}
    h = 0
    for lst in by_p.values():
        lst.sort(key=lambda n: (n.S, n.index))
        for k in range(len(lst) - 1):
            a, b = lst[k], lst[k + 1]
            h += 1
            heat[f"H{h}"] = oriented_edge("heat", a, b)

    all_edges = {**mechanical, **heat}
    nodes = _attach_edges_to_nodes_ps(nodes, all_edges)
    return nodes, all_edges


def build_subcycles(nodes: dict[int, Node], edges: dict[str, Edge]) -> list[SubCycle]:
    """
    枚举最小子循环：对每个起始节点尝试 ``edge_up`` → 对端 ``edge_right`` → 对端 ``edge_down``（取该边 ``tail``）
    → 对端 ``edge_left``（取该边 ``tail``）须回到起始节点；四节点互异。

    使用 ``frozenset`` 对四角去重，每个物理单元至多输出一次。不对 ``P``/``S`` 几何做一致性校验。
    """
    seen_cells: set[frozenset[int]] = set()
    out: list[SubCycle] = []

    for n0 in nodes.values():
        ku = n0.edge_up
        if ku is None:
            continue
        e_left = edges[ku]
        if e_left.kind != "mechanical":
            continue
        n1 = nodes[e_left.head]
        kr = n1.edge_right
        if kr is None:
            continue
        e_top = edges[kr]
        if e_top.kind != "heat":
            continue
        n2 = nodes[e_top.head]
        kd = n2.edge_down
        if kd is None:
            continue
        e_right = edges[kd]
        if e_right.kind != "mechanical":
            continue
        n3 = nodes[e_right.tail]
        kl = n3.edge_left
        if kl is None:
            continue
        e_bottom = edges[kl]
        if e_bottom.kind != "heat":
            continue
        if e_bottom.tail != n0.index:
            continue

        idxs = (n0.index, n1.index, n2.index, n3.index)
        if len(set(idxs)) != 4:
            continue

        cell = frozenset(idxs)
        if cell in seen_cells:
            continue
        seen_cells.add(cell)

        out.append(
            SubCycle(
                nodes=(n0.index, n1.index, n2.index, n3.index),
                edges=(ku, kr, kd, kl),
            )
        )

    return out


class ClosedCycleLayer:
    """
    闭式循环层。仅依赖 ``ClosedCycleTPInput`` 即可持有工质并创建默认物性后端；
    也可注入 ``properties`` 以便测试或替换实现。

    ``analyze_topology()`` 填充 ``nodes``、``edges``（``dict[str, Edge]``，键 ``M*`` / ``H*``）与 ``subcycles``。
    """

    def __init__(self, inp: ClosedCycleTPInput, *, properties: FluidPropertySolver | None = None) -> None:
        self.input = inp
        self.fluid = inp.fluid
        self.properties = properties if properties is not None else CoolPropFluidPropertySolver(inp.fluid)
        # 分析前为空；``analyze_topology`` 一次性写入基准拓扑
        self.nodes: dict[int, Node] = {}
        self.edges: dict[str, Edge] = {}
        self.subcycles: list[SubCycle] = []

    def analyze_topology(self) -> None:
        """调用 ``build_node_edge_topology`` 写入 ``nodes`` 与 ``edges``，再 ``build_subcycles`` 写入 ``subcycles``。"""
        self.nodes, self.edges = build_node_edge_topology(self.properties, self.input)
        self.subcycles = build_subcycles(self.nodes, self.edges)
