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

build_node_edge_topology 生成的节点和边为拓扑快照；``ClosedCycleLayer.assign_edge_mass_flows_from_subcycles`` 将子循环有符号流量按顺时针约定汇聚到各 ``Edge.mass_flow``（``tail → head`` 方向）。

最小 **子循环**（``SubCycle``）为沿 ``edge_up`` → ``edge_right`` → ``edge_down`` → ``edge_left`` 闭合的 **4 节点、4 边** 单元；
节点顺序为顺时针 **左下 → 左上 → 右上 → 右下**，边顺序为 **左、上、右、下**。
``SubCycle`` 含子循环 **单一** 质量流 ``mass_flow``（单位同 ``Edge.mass_flow``，可为负）；
拓扑分析后由 ``ClosedCycleLayer.subcycle_mass_flows`` 经 ``sync_subcycle_mass_flows_to_subcycles`` 写入，外部优化器优先改列表再量化、同步。

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
from typing import TYPE_CHECKING, Literal, Sequence

if TYPE_CHECKING:
    from core.non_ideal_closed_cycle_layer import NonIdealClosedCycleLayer

import config as cyges_config

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
    """全层质量流量上界 [kg/s] 等；用于子循环默认流量（系数见根目录 ``config.SUBCYCLE_INITIAL_MASS_FLOW_FRACTION_OF_MAX``）与量化步长基准。"""
    subcycle_mass_flow_step_fraction: float = field(
        default_factory=lambda: float(cyges_config.SUBCYCLE_MASS_FLOW_STEP_FRACTION_DEFAULT)
    )
    """子循环质量流量化步长占 ``max_mass_flow`` 的比例；未传入时取根目录 ``config.SUBCYCLE_MASS_FLOW_STEP_FRACTION_DEFAULT``。步长 ``step = subcycle_mass_flow_step_fraction * max_mass_flow``。"""


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
    - ``mass_flow``：质量流 [kg/s] 等，初始 ``None``；经 ``ClosedCycleLayer.assign_edge_mass_flows_from_subcycles`` 后为沿 ``tail → head`` 的代数和（未被子循环覆盖的边保持 ``None``）。
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
    - ``mass_flow``：子循环质量流 [kg/s] 等标量，可为负；由 ``ClosedCycleLayer`` 的 ``subcycle_mass_flows`` 同步写入（见 ``sync_subcycle_mass_flows_to_subcycles``）。
    """

    nodes: tuple[int, int, int, int]
    edges: tuple[str, str, str, str]
    mass_flow: float | None = None


@dataclass(frozen=True)
class SkippedPoint:
    """
    拓扑构建中被物性求解器异常静默跳过的候选点；仅用于诊断，不参与拓扑/流量计算。

    - ``stage``：``"primary_TP"`` 为一级网格点构建阶段（按 ``(T, P)`` 求 ``H, S``）；
      ``"secondary_PS"`` 为沿一级节点等熵延伸阶段（按 ``(P, S)`` 求 ``T, H``）。
    - ``T``：``stage == "primary_TP"`` 时为请求的 T；``stage == "secondary_PS"`` 时为 ``None``（T 由 PS 解出，失败时不存在）。
    - ``P``：请求时的压力。
    - ``S``：``stage == "secondary_PS"`` 时为父节点 S；``stage == "primary_TP"`` 时为 ``None``。
    - ``reason``：触发跳过的异常字符串表示，便于复盘是否落在包线/相区外等。
    """

    stage: Literal["primary_TP", "secondary_PS"]
    T: float | None
    P: float
    S: float | None
    reason: str


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
            raise ValueError(f"node {ni} duplicate {side} edge: {cur!r} and {ek!r}")
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
) -> tuple[dict[int, Node], dict[str, Edge], list[SkippedPoint]]:
    """
    节点与边拓扑的一次性生成：一级 TP 网格、等熵二级与机械边（键 ``M*``）、
    全节点等压分桶换热边（键 ``H*``），边方向满足 PS 向上/向右约定；最后将边键写回各 ``Node`` 的 ``edge_*`` 字段。

    返回 ``(nodes, edges, skipped_points)``；``skipped_points`` 收集物性求解器异常静默跳过的候选点，仅供诊断。
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
        raise RuntimeError(
            f"PS 方向不可判定：节点 {na.index!r}/{nb.index!r} 互不严格小于；"
            "上游算法（机械边按 P 升序成链、换热边按 S 升序分桶）不应抵达此分支"
        )

    skipped_points: list[SkippedPoint] = []

    # —— 一级 TP 网格 ——
    t_list = build_axis(inp.t_min, inp.t_max, inp.t_quantiles)
    p_list = build_axis(inp.p_min, inp.p_max, inp.p_quantiles)
    grid: list[Node] = []
    idx = 0
    for Tk in t_list:
        for Pk in p_list:
            try:
                st = solver.state("TP", Tk, Pk)
            except Exception as exc:
                skipped_points.append(
                    SkippedPoint(stage="primary_TP", T=Tk, P=Pk, S=None, reason=repr(exc))
                )
                continue
            grid.append(Node(index=idx, T=Tk, P=Pk, H=st["H"], S=st["S"], parent=None))
            idx += 1

    # —— 二级等熵 + 机械边 ——
    p_axis = build_axis(inp.p_min, inp.p_max, inp.p_quantiles)
    secondary: list[Node] = []
    mechanical: dict[str, Edge] = {}
    m = 0
    idx = len(grid)

    for prim in grid:
        batch: list[Node] = []
        for p2 in p_axis:
            if abs(p2 - prim.P) < 1e-9 * max(1.0, abs(p2)):
                continue
            try:
                st = solver.state("PS", p2, prim.S)
                t2 = st["T"]
            except Exception as exc:
                skipped_points.append(
                    SkippedPoint(stage="secondary_PS", T=None, P=p2, S=prim.S, reason=repr(exc))
                )
                continue
            if not (inp.t_min <= t2 <= inp.t_max):
                continue
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
    return nodes, all_edges, skipped_points


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

    默认 ``auto_analyze=True`` 时，构造末尾会调用一次 ``analyze_topology()``；否则保持空拓扑直至显式调用。
    外部修改 ``subcycle_mass_flows`` 后若需离散化并写回子循环与边，调用 ``commit_subcycle_mass_flows_to_topology()``。
    非理想分析通过 ``ensure_non_ideal()`` 得到 ``non_ideal``；每次 ``analyze_topology()`` 或 ``commit_subcycle_mass_flows_to_topology()`` 会清空 ``non_ideal``（闭式循环层稳定后再重建非理想层）。

    ``analyze_topology()`` 填充 ``nodes``、``edges``、``subcycles``、``skipped_points``，以及 ``subcycle_mass_flows``（默认每项为 ``config.SUBCYCLE_INITIAL_MASS_FLOW_FRACTION_OF_MAX * max_mass_flow``）
    并同步到各 ``SubCycle.mass_flow``，最后调用 ``assign_edge_mass_flows_from_subcycles`` 写入 ``edges[*].mass_flow``。
    """

    def __init__(
        self,
        inp: ClosedCycleTPInput,
        *,
        properties: FluidPropertySolver | None = None,
        auto_analyze: bool = True,
    ) -> None:
        self.input = inp
        self.fluid = inp.fluid
        self.properties = properties if properties is not None else CoolPropFluidPropertySolver(inp.fluid)
        # 分析前为空；``analyze_topology`` 一次性写入基准拓扑
        self.nodes: dict[int, Node] = {}
        self.edges: dict[str, Edge] = {}
        self.subcycles: list[SubCycle] = []
        self.subcycle_mass_flows: list[float] = []
        self.skipped_points: list[SkippedPoint] = []
        self.non_ideal: NonIdealClosedCycleLayer | None = None
        if auto_analyze:
            self.analyze_topology()

    def _invalidate_non_ideal(self) -> None:
        """理想层基准发生变化时清空非理想层挂载；外部不直接调用。"""
        self.non_ideal = None

    def quantize_subcycle_mass_flows(self) -> None:
        """
        将 ``subcycle_mass_flows`` 就地舍入到 ``step`` 的整数倍：
        ``step = subcycle_mass_flow_step_fraction * max_mass_flow``（``subcycle_mass_flow_step_fraction`` 未在输入中指定时由 ``ClosedCycleTPInput`` 使用 ``config.SUBCYCLE_MASS_FLOW_STEP_FRACTION_DEFAULT``）。
        列表为空则不做任何事；调用方应保证 ``max_mass_flow`` 与 ``subcycle_mass_flow_step_fraction`` 在此前为有效正数，本层不再做防御性校验。
        """
        mf = self.input.max_mass_flow
        frac = self.input.subcycle_mass_flow_step_fraction
        step = frac * mf
        self.subcycle_mass_flows = [round(x / step) * step for x in self.subcycle_mass_flows]

    def sync_subcycle_mass_flows_to_subcycles(self) -> None:
        """将 ``subcycle_mass_flows[i]`` 写入 ``subcycles[i].mass_flow``；长度须一致。"""
        for i, sc in enumerate(self.subcycles):
            sc.mass_flow = self.subcycle_mass_flows[i]

    def assign_edge_mass_flows_from_subcycles(self) -> None:
        """
        将各子循环 ``mass_flow`` 汇聚到 ``self.edges[*].mass_flow``。

        约定：``SubCycle.nodes`` 为顺时针 ``(左下, 左上, 右上, 右下)``；``mass_flow > 0`` 表示沿该顺时针回路，
        ``< 0`` 表示逆时针。每条拓扑边 ``tail → head`` 上，``mass_flow`` 为沿该方向的代数和；``None`` 子循环按 ``0``。

        若某子循环的边键与顺时针段 ``(u→v)`` 和 ``(tail, head)`` 均不一致，抛出 ``ValueError``。
        先清空所有边的 ``mass_flow``，再仅对至少出现一次的边写入累加结果（含 ``0.0``）。
        """
        for e in self.edges.values():
            e.mass_flow = None
        totals: dict[str, float] = defaultdict(float)

        for sc in self.subcycles:
            q = 0.0 if sc.mass_flow is None else float(sc.mass_flow)
            n0, n1, n2, n3 = sc.nodes
            segment_edge: tuple[tuple[tuple[int, int], str], ...] = (
                ((n0, n1), sc.edges[0]),
                ((n1, n2), sc.edges[1]),
                ((n2, n3), sc.edges[2]),
                ((n3, n0), sc.edges[3]),
            )
            for (u, v), ek in segment_edge:
                edge = self.edges[ek]
                if edge.tail == u and edge.head == v:
                    totals[ek] += q
                elif edge.tail == v and edge.head == u:
                    totals[ek] -= q
                else:
                    raise ValueError(
                        f"子循环边 {ek!r} 的 tail/head {(edge.tail, edge.head)!r} 与顺时针段 {(u, v)!r} 不一致"
                    )

        for ek, tot in totals.items():
            self.edges[ek].mass_flow = tot

    def commit_subcycle_mass_flows_to_topology(self) -> None:
        """
        对 ``subcycle_mass_flows`` 就地量化到步长整数倍，再 ``sync`` 到各 ``SubCycle``，并 ``assign`` 到各边。

        供外部在写入任意浮点初值后、后续分析前保证离散化与拓扑一致。``subcycle_mass_flows`` 与 ``subcycles`` 长度须一致。
        无子循环时仅将各边 ``mass_flow`` 置为 ``None``（与 ``assign_edge_mass_flows_from_subcycles`` 行为一致），不调用量化。
        调用开始时清空 ``non_ideal``，与 ``analyze_topology`` 一致，保证仅在理想层流量稳定后再 ``ensure_non_ideal``。
        """
        self._invalidate_non_ideal()
        n_sc = len(self.subcycles)
        n_q = len(self.subcycle_mass_flows)
        if n_sc != n_q:
            raise ValueError(f"subcycle_mass_flows 长度 {n_q} 与 subcycles 长度 {n_sc} 不一致")
        if n_sc == 0:
            self.assign_edge_mass_flows_from_subcycles()
            return
        self.quantize_subcycle_mass_flows()
        self.sync_subcycle_mass_flows_to_subcycles()
        self.assign_edge_mass_flows_from_subcycles()

    def ensure_non_ideal(self) -> NonIdealClosedCycleLayer:
        """若尚无 ``non_ideal`` 则基于当前层构造并挂载；不修改理想层基准字段。须在 ``analyze_topology`` / ``commit`` 之后、理想层稳定时调用。"""
        from core.non_ideal_closed_cycle_layer import NonIdealClosedCycleLayer

        if self.non_ideal is None:
            self.non_ideal = NonIdealClosedCycleLayer.from_closed_cycle_layer(self)
        return self.non_ideal

    def analyze_topology(self) -> None:
        """构建拓扑与子循环；初始化 ``subcycle_mass_flows``（系数见根目录 ``config``）并同步，再赋边流量。"""
        self._invalidate_non_ideal()
        self.nodes, self.edges, self.skipped_points = build_node_edge_topology(self.properties, self.input)
        self.subcycles = build_subcycles(self.nodes, self.edges)
        n = len(self.subcycles)
        if n == 0:
            self.subcycle_mass_flows = []
            return
        mf = self.input.max_mass_flow
        q0 = cyges_config.SUBCYCLE_INITIAL_MASS_FLOW_FRACTION_OF_MAX * (0.0 if mf is None else float(mf))
        self.subcycle_mass_flows = [q0] * n
        self.sync_subcycle_mass_flows_to_subcycles()
        self.assign_edge_mass_flows_from_subcycles()
