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

import math
import warnings
from collections import defaultdict
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING, Literal, Sequence

if TYPE_CHECKING:
    from core.non_ideal_closed_cycle_layer import NonIdealClosedCycleLayer

import config as cyges_config

from core.fluid_property_solver import CoolPropFluidPropertySolver, FluidPropertySolver


MERGED_ISOLATED_NODE_EDGE_KEY = "__CyGES_ISOLATED__"
"""精简拓扑中四邻边槽全空、且未上链的节点占位键：在 ``SimplifiedTopology.merged_into`` 中作为「不参与活跃子图」的标记，无对应 ``SimplifiedEdge``。"""


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


@dataclass(frozen=True)
class SimplifiedEdge:
    """
    精简边：一段连续同类型原始边链合并而成的有向边。

    - ``kind``：``"mechanical"`` 或 ``"heat"``。
    - ``tail`` / ``head``：精简边两个保留端点的原始 ``Node.index``；沿 ``mass_flow`` 方向 ``tail → head``。
    - ``constituent_edges``：覆盖的原始边键，按 ``tail → head`` 顺序。
    - ``merged_nodes``：被合并掉的中间节点 ``index``，按 ``tail → head`` 顺序；与 ``constituent_edges`` 长度相差 1。
    - ``mass_flow``：链上统一流量沿 ``tail → head`` 方向（恒非负）；若链上全部为 ``None`` 则保持 ``None``。
    """

    kind: Literal["mechanical", "heat"]
    tail: int
    head: int
    constituent_edges: tuple[str, ...]
    merged_nodes: tuple[int, ...]
    mass_flow: float | None


@dataclass(frozen=True)
class SimplifiedTopology:
    """
    精简后的 PS 拓扑快照；全部字段为不可变值。

    - ``kept_nodes``：保留下来的原始节点 ``index`` 集合。
    - ``simplified_edges``：键 ``"SM*"``（机械）/ ``"SH*"``（换热） → ``SimplifiedEdge``；以 tuple-of-items 形式持有，便于 ``frozen=True``。
    - ``merged_into``：被合并或视为移出活跃子图的原始节点 ``index`` → 精简边键或占位键
      ``MERGED_ISOLATED_NODE_EDGE_KEY``（过滤后四邻边槽全空且未出现在任何 ``SimplifiedEdge.merged_nodes`` 中的节点）；
      以 tuple-of-items 形式持有。
    """

    kept_nodes: frozenset[int]
    simplified_edges: tuple[tuple[str, SimplifiedEdge], ...]
    merged_into: tuple[tuple[int, str], ...]

    def edges_dict(self) -> dict[str, SimplifiedEdge]:
        """便捷访问：把精简边以 ``dict`` 形式返回（构造时为 ``tuple`` 以维持冻结性）。"""
        return dict(self.simplified_edges)

    def merged_dict(self) -> dict[int, str]:
        """便捷访问：节点 ``index`` → 所在精简边键，或 ``MERGED_ISOLATED_NODE_EDGE_KEY``（过滤后四邻边全空且未上链）。"""
        return dict(self.merged_into)


def _mass_flow_close(a: float | None, b: float | None) -> bool:
    """链上 ``mass_flow`` 一致性比较：``None`` 视为缺测，与任意值兼容；两个非 ``None`` 用 ``math.isclose``。"""
    if a is None or b is None:
        return True
    return math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12)


def _subcycle_edge_key_set(subcycles: Sequence[SubCycle]) -> frozenset[str]:
    """所有子循环 ``SubCycle.edges`` 中出现过的原始边键并集。"""
    out: set[str] = set()
    for sc in subcycles:
        for ek in sc.edges:
            out.add(ek)
    return frozenset(out)


def _edge_has_nonzero_mass_flow(e: Edge) -> bool:
    """边是否视为「有流量」：``mass_flow`` 非 ``None`` 且绝对值在容差意义下非零。"""
    if e.mass_flow is None:
        return False
    return not math.isclose(float(e.mass_flow), 0.0, rel_tol=0.0, abs_tol=1e-12)


def filter_topology_for_non_ideal(
    nodes: dict[int, Node],
    edges: dict[str, Edge],
    subcycles: Sequence[SubCycle],
) -> tuple[dict[int, Node], dict[str, Edge]]:
    """
    在合并精简拓扑**之前**使用的拓扑拷贝：不修改传入的 ``nodes`` / ``edges``。

    保留的边须**同时**满足：
    1. 至少出现在一个 ``SubCycle.edges`` 中（不在任何子循环中的边剔除）；
    2. ``mass_flow`` 非 ``None`` 且在容差意义下非零（无流量边剔除）。

    返回的 ``nodes`` 为 ``Node`` 的浅拷贝字典，邻边槽中指向被剔除边的键置为 ``None``；
    ``edges`` 为仅含保留边的**新**字典，``Edge`` 对象与输入共享引用（只读使用）。
    """
    in_sub = _subcycle_edge_key_set(subcycles)
    kept_keys: set[str] = set()
    for ek in in_sub:
        if ek not in edges:
            continue
        if _edge_has_nonzero_mass_flow(edges[ek]):
            kept_keys.add(ek)

    filtered_edges = {k: edges[k] for k in sorted(kept_keys)}

    new_nodes: dict[int, Node] = {}
    for i, n in nodes.items():
        eu = n.edge_up if n.edge_up in filtered_edges else None
        ed = n.edge_down if n.edge_down in filtered_edges else None
        el = n.edge_left if n.edge_left in filtered_edges else None
        er = n.edge_right if n.edge_right in filtered_edges else None
        if (eu, ed, el, er) == (n.edge_up, n.edge_down, n.edge_left, n.edge_right):
            new_nodes[i] = n
        else:
            new_nodes[i] = replace(n, edge_up=eu, edge_down=ed, edge_left=el, edge_right=er)

    return new_nodes, filtered_edges


def _aggregate_chain_mass_flow(
    chain_edges: list[str],
    edges: dict[str, Edge],
    chain_label: str,
) -> float | None:
    """聚合链上一组原始边的 ``mass_flow``；全 ``None`` 返回 ``None``；存在非 ``None`` 但相互不一致抛 ``ValueError``。"""
    seen: float | None = None
    for ek in chain_edges:
        m = edges[ek].mass_flow
        if m is None:
            continue
        if seen is None:
            seen = float(m)
        elif not _mass_flow_close(seen, float(m)):
            raise ValueError(
                f"精简{chain_label}链 {chain_edges!r} 上 mass_flow 不一致：含 {seen!r} 与 {m!r}"
            )
    return seen


def _find_typed_chains(
    nodes: dict[int, Node],
    edges: dict[str, Edge],
    kind: Literal["mechanical", "heat"],
) -> list[dict]:
    """
    按 ``kind`` 找出该类型原始边形成的所有连通分量（每个分量是一条简单链，无分叉——PS 网格保证）。

    返回列表，每项为 ``{"nodes": ordered_list, "edges": ordered_list}``：
    ``nodes`` 按 PS 正向排序（机械：``P`` 升；换热：``S`` 升），``edges`` 为相邻节点对之间的边键，长度 ``= len(nodes) - 1``。
    """
    if kind == "mechanical":
        sort_key = lambda i: nodes[i].P
        slot_a = "edge_up"
        slot_b = "edge_down"
    else:
        sort_key = lambda i: nodes[i].S
        slot_a = "edge_left"
        slot_b = "edge_right"

    adj: dict[int, list[str]] = defaultdict(list)
    for i, n in nodes.items():
        for slot in (slot_a, slot_b):
            ek = getattr(n, slot)
            if ek is not None and edges[ek].kind == kind:
                adj[i].append(ek)

    visited: set[int] = set()
    chains: list[dict] = []
    for start in nodes:
        if start in visited or not adj[start]:
            continue
        component: set[int] = set()
        stack = [start]
        while stack:
            v = stack.pop()
            if v in component:
                continue
            component.add(v)
            for ek in adj[v]:
                e = edges[ek]
                other = e.head if e.tail == v else e.tail
                if other not in component:
                    stack.append(other)
        visited |= component
        ordered_nodes = sorted(component, key=sort_key)
        chain_edges_in_order: list[str] = []
        for k in range(len(ordered_nodes) - 1):
            a, b = ordered_nodes[k], ordered_nodes[k + 1]
            found = None
            for ek in adj[a]:
                e = edges[ek]
                if (e.tail == a and e.head == b) or (e.tail == b and e.head == a):
                    found = ek
                    break
            if found is None:
                raise RuntimeError(
                    f"{kind} 链上相邻节点 {a}/{b} 之间未找到边；PS 网格连通性被破坏"
                )
            chain_edges_in_order.append(found)
        chains.append({"nodes": ordered_nodes, "edges": chain_edges_in_order})
    return chains


def _simplify_chain(
    chain: dict,
    nodes: dict[int, Node],
    edges: dict[str, Edge],
    kind: Literal["mechanical", "heat"],
    is_mergeable: dict[int, bool],
    counter: int,
) -> tuple[list[tuple[str, SimplifiedEdge]], list[tuple[int, str]], int]:
    """
    把一条链按「切分点」切成若干精简边。切分点 = 链两端 + 链中间的非 mergeable 节点。
    相邻切分点之间的所有原始边合并为一条精简边，中间的 mergeable 节点登记到 ``merged_into``。
    """
    chain_nodes: list[int] = chain["nodes"]
    chain_edges: list[str] = chain["edges"]

    if len(chain_nodes) < 2:
        return [], [], counter

    cut_positions: list[int] = [0]
    for k in range(1, len(chain_nodes) - 1):
        if not is_mergeable[chain_nodes[k]]:
            cut_positions.append(k)
    cut_positions.append(len(chain_nodes) - 1)

    prefix = "SM" if kind == "mechanical" else "SH"
    chain_label = "机械" if kind == "mechanical" else "换热"

    simp_edges: list[tuple[str, SimplifiedEdge]] = []
    merged_entries: list[tuple[int, str]] = []

    for seg in range(len(cut_positions) - 1):
        i0, i1 = cut_positions[seg], cut_positions[seg + 1]
        ps_low_idx = chain_nodes[i0]
        ps_high_idx = chain_nodes[i1]
        seg_merged_in_ps_order = chain_nodes[i0 + 1 : i1]
        seg_edges_in_ps_order = chain_edges[i0:i1]

        aggregated_mf = _aggregate_chain_mass_flow(seg_edges_in_ps_order, edges, chain_label)

        if aggregated_mf is None or aggregated_mf == 0.0:
            tail = ps_low_idx
            head = ps_high_idx
            mf_out: float | None = aggregated_mf
            ordered_edges = tuple(seg_edges_in_ps_order)
            ordered_merged = tuple(seg_merged_in_ps_order)
        elif aggregated_mf > 0:
            tail = ps_low_idx
            head = ps_high_idx
            mf_out = aggregated_mf
            ordered_edges = tuple(seg_edges_in_ps_order)
            ordered_merged = tuple(seg_merged_in_ps_order)
        else:
            tail = ps_high_idx
            head = ps_low_idx
            mf_out = -aggregated_mf
            ordered_edges = tuple(reversed(seg_edges_in_ps_order))
            ordered_merged = tuple(reversed(seg_merged_in_ps_order))

        counter += 1
        simp_key = f"{prefix}{counter}"
        simp_edges.append(
            (
                simp_key,
                SimplifiedEdge(
                    kind=kind,
                    tail=tail,
                    head=head,
                    constituent_edges=ordered_edges,
                    merged_nodes=ordered_merged,
                    mass_flow=mf_out,
                ),
            )
        )
        for mn in ordered_merged:
            merged_entries.append((mn, simp_key))

    return simp_edges, merged_entries, counter


def build_simplified_topology(
    nodes: dict[int, Node],
    edges: dict[str, Edge],
    subcycles: Sequence[SubCycle],
) -> SimplifiedTopology:
    """
    在拓扑产物上构建精简图（**不修改入参**）。

    **步骤 0**：先 ``filter_topology_for_non_ideal(nodes, edges, subcycles)`` ——
    剔除不在任何子循环中的边，以及 ``mass_flow`` 为 ``None`` 或数值为零的边；并同步清空节点上指向被删边的邻边槽。

    **步骤 1**（链合并）：
    - 节点的 ``edge_up`` 与 ``edge_down`` 均为 ``None``（仅挂换热邻边） → 合并到换热链。
    - 节点的 ``edge_left`` 与 ``edge_right`` 均为 ``None``（仅挂机械邻边） → 合并到机械链。
    - 其余节点（含两类邻边都有、或都没有的孤立点）作为切分点或链端保留。
    - 精简边方向（``tail → head``）：聚合该精简边覆盖的所有原始 ``Edge.mass_flow``——非 ``None`` 值要求互相一致
      （``math.isclose``），否则抛 ``ValueError``；正值 → 沿 PS 正向，负值 → 反向并取正绝对值，全 ``None`` / 全 ``0`` → 沿 PS 正向且 ``mass_flow`` 保持原状。

    **步骤 2**：对过滤后 **四邻边槽全空** 且尚未出现在 ``merged_into`` 中的节点，追加
    ``(index, MERGED_ISOLATED_NODE_EDGE_KEY)``——不删 ``nodes_f`` 中的键，但使其不出现在 ``kept_nodes`` 中（与链上 ``merged_nodes`` 一并可经 ``merged_dict`` 查询）。

    **步骤 3**：``kept_nodes = { i ∈ nodes_f | i 不在 merged_into 的任一条目键中 }``。
    """
    nodes_f, edges_f = filter_topology_for_non_ideal(nodes, edges, subcycles)

    is_mech_only: dict[int, bool] = {}
    is_heat_only: dict[int, bool] = {}
    for i, n in nodes_f.items():
        has_mech = n.edge_up is not None or n.edge_down is not None
        has_heat = n.edge_left is not None or n.edge_right is not None
        is_mech_only[i] = has_mech and not has_heat
        is_heat_only[i] = has_heat and not has_mech

    mech_chains = _find_typed_chains(nodes_f, edges_f, "mechanical")
    heat_chains = _find_typed_chains(nodes_f, edges_f, "heat")

    simplified_edges: list[tuple[str, SimplifiedEdge]] = []
    merged_into: list[tuple[int, str]] = []

    mech_counter = 0
    for chain in mech_chains:
        segs, mlinks, mech_counter = _simplify_chain(
            chain,
            nodes_f,
            edges_f,
            kind="mechanical",
            is_mergeable=is_mech_only,
            counter=mech_counter,
        )
        simplified_edges.extend(segs)
        merged_into.extend(mlinks)

    heat_counter = 0
    for chain in heat_chains:
        segs, mlinks, heat_counter = _simplify_chain(
            chain,
            nodes_f,
            edges_f,
            kind="heat",
            is_mergeable=is_heat_only,
            counter=heat_counter,
        )
        simplified_edges.extend(segs)
        merged_into.extend(mlinks)

    merged_idx_so_far = {mn for mn, _ in merged_into}
    for i, n in nodes_f.items():
        if i in merged_idx_so_far:
            continue
        if (
            n.edge_up is None
            and n.edge_down is None
            and n.edge_left is None
            and n.edge_right is None
        ):
            merged_into.append((i, MERGED_ISOLATED_NODE_EDGE_KEY))

    merged_set: set[int] = {mn for mn, _ in merged_into}
    kept_nodes = frozenset(i for i in nodes_f if i not in merged_set)

    return SimplifiedTopology(
        kept_nodes=kept_nodes,
        simplified_edges=tuple(simplified_edges),
        merged_into=tuple(merged_into),
    )


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
    每次 ``analyze_topology()`` / ``commit_subcycle_mass_flows_to_topology()`` 末尾会自动调用 ``_rebuild_simplified()`` 同步 ``simplified``，并清空 ``non_ideal`` 快照容器；通过 ``ensure_non_ideal()`` 再生成与当前 ``simplified`` 同源的快照。

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
        self.simplified: SimplifiedTopology | None = None
        self.non_ideal: NonIdealClosedCycleLayer | None = None
        if auto_analyze:
            self.analyze_topology()

    def _invalidate_non_ideal(self) -> None:
        """理想层基准发生变化时清空非理想层挂载；外部不直接调用。"""
        self.non_ideal = None

    def _rebuild_simplified(self) -> None:
        """
        基于当前 ``nodes`` / ``edges`` / ``subcycles`` 重建 ``self.simplified``。

        子循环为空时不再调用 ``build_simplified_topology`` —— 直接置为空骨架
        ``SimplifiedTopology(kept_nodes=frozenset(), simplified_edges=(), merged_into=())``，
        并通过 ``warnings.warn`` 发出 ``RuntimeWarning`` 提示拓扑里没有子循环。
        """
        if not self.subcycles:
            warnings.warn(
                "ClosedCycleLayer: 当前没有子循环，simplified 置为空骨架（kept_nodes / simplified_edges / merged_into 均为空）。",
                RuntimeWarning,
                stacklevel=2,
            )
            self.simplified = SimplifiedTopology(
                kept_nodes=frozenset(),
                simplified_edges=(),
                merged_into=(),
            )
            return
        self.simplified = build_simplified_topology(self.nodes, self.edges, self.subcycles)

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
            self._rebuild_simplified()
            return
        self.quantize_subcycle_mass_flows()
        self.sync_subcycle_mass_flows_to_subcycles()
        self.assign_edge_mass_flows_from_subcycles()
        self._rebuild_simplified()

    def ensure_non_ideal(self) -> NonIdealClosedCycleLayer:
        """若尚无 ``non_ideal`` 则基于当前层构造并挂载；不修改理想层基准字段。须在 ``analyze_topology`` / ``commit`` 之后、理想层稳定时调用。"""
        from core.non_ideal_closed_cycle_layer import NonIdealClosedCycleLayer

        if self.non_ideal is None:
            self.non_ideal = NonIdealClosedCycleLayer.from_closed_cycle_layer(self)
        return self.non_ideal

    def analyze_topology(self) -> None:
        """构建拓扑与子循环；初始化 ``subcycle_mass_flows``（系数见根目录 ``config``）并同步，再赋边流量；末尾同步 ``simplified``。"""
        self._invalidate_non_ideal()
        self.nodes, self.edges, self.skipped_points = build_node_edge_topology(self.properties, self.input)
        self.subcycles = build_subcycles(self.nodes, self.edges)
        n = len(self.subcycles)
        if n == 0:
            self.subcycle_mass_flows = []
            self._rebuild_simplified()
            return
        mf = self.input.max_mass_flow
        q0 = cyges_config.SUBCYCLE_INITIAL_MASS_FLOW_FRACTION_OF_MAX * (0.0 if mf is None else float(mf))
        self.subcycle_mass_flows = [q0] * n
        self.sync_subcycle_mass_flows_to_subcycles()
        self.assign_edge_mass_flows_from_subcycles()
        self._rebuild_simplified()
