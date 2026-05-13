"""
非理想闭式循环层：依附于 ``ClosedCycleLayer``，使用冻结基准快照做后续分析；不修改理想层上的 ``nodes`` / ``edges`` / ``subcycles``。
父层在 ``analyze_topology`` 与 ``commit_subcycle_mass_flows_to_topology`` 时会清空 ``non_ideal``，须在理想层稳定后再 ``ensure_non_ideal``。

本层挂载：
- ``baseline``：理想层在 ensure_non_ideal 时刻的拓扑与流量基准快照（只读、值拷贝）。
- ``simplified``：在 baseline 拓扑上构建的"精简 PS 拓扑"——**先**按子循环与流量剔除无效边（见 ``filter_topology_for_non_ideal``），
  再按"另一类邻边槽全空"规则合并中间节点；每条精简边持有沿 mass_flow 的方向、覆盖的原始边链、统一效率（来自 ``config``，临时配置）。

简化规则与方向约定：
- ``edge_up`` 与 ``edge_down`` 均为 ``None`` 的节点（只挂换热边） → 合并到换热链。
- ``edge_left`` 与 ``edge_right`` 均为 ``None`` 的节点（只挂机械边） → 合并到机械链。
- 其余节点（同时挂两类邻边、或两类邻边都没有的孤立点）作为精简图的保留节点。
- 精简边方向：聚合链上所有原始 ``Edge.mass_flow``，沿 ``tail → head`` 方向使其恒非负；
  全 ``None`` / 全 ``0`` 时方向退化为 PS 正向（机械链按 ``P`` 升序，换热链按 ``S`` 升序）。
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, replace
from typing import Literal, Sequence

import config as cyges_config

from core.closed_cycle_layer import Edge, Node, SubCycle

# 过滤后四邻边槽全空、且未作为链上 ``merged_nodes`` 出现的节点，在 ``merged_into`` 中指向此占位键
# （无对应 ``SimplifiedEdge``；仅表示「不参与活跃子图」，与 ``kept_nodes`` 互斥）。
MERGED_ISOLATED_NODE_EDGE_KEY = "__CyGES_ISOLATED__"


@dataclass(frozen=True)
class BaselineEdgeRecord:
    """单条边的基准只读快照（不持有可变 ``Edge`` 引用）。"""

    key: str
    kind: Literal["mechanical", "heat"]
    tail: int
    head: int
    mass_flow: float | None


@dataclass(frozen=True)
class BaselineSubCycleRecord:
    """单个子循环的基准只读快照。"""

    nodes: tuple[int, int, int, int]
    edges: tuple[str, str, str, str]
    mass_flow: float | None


@dataclass(frozen=True)
class BaselineTopologySnapshot:
    """
    理想闭式循环层在某一时刻的拓扑与流量基准；字段均为不可变值，与之后对 ``ClosedCycleLayer`` 的原地修改脱钩。

    ``edge_records`` 与 ``subcycles`` / ``subcycle_mass_flows`` 为构造时拷贝；``node_items`` 引用 ``Node`` 对象（其自身 ``frozen``），
    若父层 ``analyze_topology`` 整表替换 ``nodes`` 字典，旧快照中的 ``Node`` 仍为历史状态视图。
    """

    node_items: tuple[tuple[int, Node], ...]
    edge_records: tuple[BaselineEdgeRecord, ...]
    subcycles: tuple[BaselineSubCycleRecord, ...]
    subcycle_mass_flows: tuple[float, ...]

    @classmethod
    def from_closed_cycle_layer(cls, layer) -> BaselineTopologySnapshot:
        ni = tuple(sorted(layer.nodes.items()))
        er = tuple(
            BaselineEdgeRecord(k, e.kind, e.tail, e.head, e.mass_flow)
            for k, e in sorted(layer.edges.items(), key=lambda kv: kv[0])
        )
        sc = tuple(
            BaselineSubCycleRecord(s.nodes, s.edges, s.mass_flow) for s in layer.subcycles
        )
        q = tuple(float(x) for x in layer.subcycle_mass_flows)
        return cls(ni, er, sc, q)


@dataclass(frozen=True)
class SimplifiedEdge:
    """
    精简边：一段连续同类型原始边链合并而成的有向边。

    - ``kind``：``"mechanical"`` 或 ``"heat"``。
    - ``tail`` / ``head``：精简边两个保留端点的原始 ``Node.index``；沿 ``mass_flow`` 方向 ``tail → head``。
    - ``constituent_edges``：覆盖的原始边键，按 ``tail → head`` 顺序。
    - ``merged_nodes``：被合并掉的中间节点 ``index``，按 ``tail → head`` 顺序；与 ``constituent_edges`` 长度相差 1。
    - ``mass_flow``：链上统一流量沿 ``tail → head`` 方向（恒非负）；若链上全部为 ``None`` 则保持 ``None``。
    - ``efficiency``：本边效率；机械边为等熵效率 η，换热边为压力保留比例 ``P_head / P_tail``。
    """

    kind: Literal["mechanical", "heat"]
    tail: int
    head: int
    constituent_edges: tuple[str, ...]
    merged_nodes: tuple[int, ...]
    mass_flow: float | None
    efficiency: float


@dataclass(frozen=True)
class SimplifiedTopology:
    """
    精简后的 PS 拓扑快照；与 ``BaselineTopologySnapshot`` 一样为冻结值拷贝。

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
    """边是否视为“有流量”：``mass_flow`` 非 ``None`` 且绝对值在容差意义下非零。"""
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
    efficiency: float,
    counter: int,
) -> tuple[list[tuple[str, SimplifiedEdge]], list[tuple[int, str]], int]:
    """
    把一条链按"切分点"切成若干精简边。切分点 = 链两端 + 链中间的非 mergeable 节点。
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
                    efficiency=efficiency,
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
    *,
    mech_eta: float,
    heat_eta: float,
) -> SimplifiedTopology:
    """
    在 baseline 拓扑产物上构建精简图（**不修改父层**）。

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
            efficiency=float(mech_eta),
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
            efficiency=float(heat_eta),
            counter=heat_counter,
        )
        simplified_edges.extend(segs)
        merged_into.extend(mlinks)

    # 过滤后四邻边槽全空的节点：不删字典项，记入 merged_into 占位键，使其不出现在 kept_nodes 中
    # （与链上 merged_nodes 互斥：已在 merged_into 中的索引不再追加）。
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

    # kept_nodes：未出现在 merged_into 中的节点（含链端点；不含链内 merged 点与孤立占位点）。
    merged_set: set[int] = {mn for mn, _ in merged_into}
    kept_nodes = frozenset(i for i in nodes_f if i not in merged_set)

    return SimplifiedTopology(
        kept_nodes=kept_nodes,
        simplified_edges=tuple(simplified_edges),
        merged_into=tuple(merged_into),
    )


@dataclass
class NonIdealClosedCycleLayer:
    """
    非理想分析容器；持有 ``BaselineTopologySnapshot`` 与 ``SimplifiedTopology``，
    后续偏移与子图等只写本对象字段，不回写父层。

    由 ``ClosedCycleLayer.ensure_non_ideal()`` 创建并赋给 ``parent.non_ideal``；父层 ``analyze_topology()`` 与 ``commit_subcycle_mass_flows_to_topology()`` 会将其清空。
    """

    baseline: BaselineTopologySnapshot
    simplified: SimplifiedTopology

    @classmethod
    def from_closed_cycle_layer(cls, layer) -> NonIdealClosedCycleLayer:
        baseline = BaselineTopologySnapshot.from_closed_cycle_layer(layer)
        simplified = build_simplified_topology(
            layer.nodes,
            layer.edges,
            layer.subcycles,
            mech_eta=float(cyges_config.NON_IDEAL_MECHANICAL_EFFICIENCY_DEFAULT),
            heat_eta=float(cyges_config.NON_IDEAL_HEAT_EFFICIENCY_DEFAULT),
        )
        return cls(baseline=baseline, simplified=simplified)
