"""非理想节点状态偏移（**旧实现**，待替换）。

当前换热 ``P ← P_ideal × σ^layer`` + 全表 ``PS`` 闭合，以及机械 ``PS→η→HP`` 支路 DFS，
自 ``non_ideal_closed_cycle_layer`` 拆出以便重写；**请勿**在新代码中依赖本模块语义。

算法细节见 ``docs/architecture.md`` 中换热/机械偏移章节（实现引用本文件）。
"""

from __future__ import annotations

import warnings
from collections import defaultdict
from dataclasses import replace
from typing import Mapping

import config as cyges_config

from core.closed_cycle_layer import Node, SimplifiedEdge
from core.fluid_property_solver import FluidPropertySolver
from core.non_ideal_closed_cycle_layer import (
    NonIdealClosedCycleLayer,
    SimplifiedDirectedGroup,
)


def _pressure_equal(P0: float, P1: float) -> bool:
    """两个压力是否在相对/绝对容差下视为相等。"""
    tol = max(1e-9, 1e-6 * max(abs(P0), abs(P1)))
    return abs(P1 - P0) <= tol


def _resolve_efficiency(
    param: float | None,
    stored: float | None,
    default: float,
    label: str,
) -> float:
    """参数 → 实例字段 → config 默认；校验 ``(0, 1]``。"""
    eta = param if param is not None else (stored if stored is not None else default)
    if not (0.0 < eta <= 1.0):
        raise ValueError(f"{label} 须在 (0, 1] 内，收到 {eta!r}")
    return eta


def _ensure_nodes(layer: NonIdealClosedCycleLayer) -> dict[int, Node]:
    if layer.nodes is None:
        layer.nodes = {i: replace(n) for i, n in layer.ideal_nodes.items()}
    return layer.nodes


def _pick_mechanical_anchor(
    group: SimplifiedDirectedGroup,
    ideal_nodes: Mapping[int, Node],
) -> int:
    """机械组基准节点选择（锚点整点状态不再修改）。

    决策树：一级 ∈ ``upstream_special_nodes`` → 该一级（index 最小）→
    ``min(upstream_special_nodes)`` → 任一一级 → ``min(layer==0)``。详见
    ``docs/architecture.md §机械锚点决策``。
    """
    group_nodes = group.nodes()
    primaries = sorted(v for v in group_nodes if ideal_nodes[v].parent is None)
    specials = group.upstream_special_nodes
    if primaries and (not specials or primaries[0] in specials):
        return primaries[0]
    if specials:
        return min(specials)
    if primaries:
        return primaries[0]
    layers = group.depth_dict()
    layer0 = sorted(v for v in group_nodes if layers[v] == 0)
    if layer0:
        return layer0[0]
    raise ValueError(f"机械有向组 {sorted(group.edge_keys)!r} 无法确定基准节点")


def _refresh_node_from_ps(
    nodes: dict[int, Node],
    properties: FluidPropertySolver,
    index: int,
) -> None:
    """以当前 ``(P, S)`` 通过 ``PS`` 闪蒸刷新节点 ``T,H``（``P,S`` 以求解器返回值为准）。"""
    n = nodes[index]
    st = properties.state("PS", n.P, n.S)
    nodes[index] = replace(n, T=st["T"], P=st["P"], H=st["H"], S=st["S"])


def _mechanical_step_known_to_unknown(
    *,
    nodes: dict[int, Node],
    properties: FluidPropertySolver,
    known: int,
    unknown: int,
    edge: SimplifiedEdge,
    eta_is: float,
) -> None:
    """沿精简机械边由 ``known`` 推 ``unknown`` 一步：``PS→η→HP``。

    公式细节见 ``docs/architecture.md §机械步公式``。
    """
    n_known = nodes[known]
    n_unknown = nodes[unknown]
    p_unknown = n_unknown.P
    h_known = n_known.H
    s_known = n_known.S

    if edge.tail == known and edge.head == unknown:
        p_tail, p_head = n_known.P, n_unknown.P
    elif edge.head == known and edge.tail == unknown:
        p_tail, p_head = n_unknown.P, n_known.P
    else:
        raise ValueError(
            f"机械边 {edge!r} 的端点 (tail={edge.tail}, head={edge.head}) 与"
            f"已知/待求 ({known}, {unknown}) 不匹配"
        )

    if _pressure_equal(p_tail, p_head):
        h_new = h_known
    else:
        h1 = properties.state("PS", p_unknown, s_known)["H"]
        if p_head > p_tail:
            h_new = (h1 - h_known) / eta_is + h_known
        else:
            h_new = (h1 - h_known) * eta_is + h_known

    st = properties.state("HP", h_new, p_unknown)
    nodes[unknown] = replace(
        n_unknown, T=st["T"], P=st["P"], H=st["H"], S=st["S"]
    )


def _walk_mechanical_branches(
    *,
    nodes: dict[int, Node],
    properties: FluidPropertySolver,
    group: SimplifiedDirectedGroup,
    edges_by_key: dict[str, SimplifiedEdge],
    anchor: int,
    eta_is: float,
) -> set[int]:
    """从锚点沿无向邻接 DFS，逐条支路单向推进（每条边只从已知端推到未知端）。"""
    adj: dict[int, list[tuple[int, str]]] = defaultdict(list)
    for ek in sorted(group.edge_keys):
        e = edges_by_key[ek]
        adj[e.tail].append((e.head, ek))
        adj[e.head].append((e.tail, ek))

    known: set[int] = {anchor}
    for start_nb, _ in sorted(adj.get(anchor, [])):
        if start_nb in known:
            continue
        stack: list[tuple[int, int]] = [(anchor, start_nb)]
        while stack:
            prev, cur = stack.pop()
            if cur in known:
                continue
            edge_key = next(ek for nb, ek in adj[prev] if nb == cur)
            edge = edges_by_key[edge_key]
            _mechanical_step_known_to_unknown(
                nodes=nodes,
                properties=properties,
                known=prev,
                unknown=cur,
                edge=edge,
                eta_is=eta_is,
            )
            known.add(cur)
            for nb, _ in sorted(adj[cur]):
                if nb not in known:
                    stack.append((cur, nb))
    return known


def apply_heat_pressure_offsets(
    layer: NonIdealClosedCycleLayer,
    heat_efficiency: float | None = None,
) -> NonIdealClosedCycleLayer:
    """（旧）按换热组层号修正 ``P`` 并全表 ``PS`` 闭合。详见 ``docs/architecture.md §换热偏移``。"""
    eta = _resolve_efficiency(
        heat_efficiency,
        layer.heat_efficiency,
        float(cyges_config.NON_IDEAL_HEAT_EFFICIENCY_DEFAULT),
        "换热总压恢复系数 σ",
    )

    nodes = _ensure_nodes(layer)

    for group in layer.heat_groups:
        layers = group.depth_dict()
        for v in group.nodes():
            p_ideal = layer.ideal_nodes[v].P
            p_new = p_ideal * (eta ** layers[v])
            nodes[v] = replace(nodes[v], P=p_new)

    for v in list(nodes.keys()):
        _refresh_node_from_ps(nodes, layer.properties, v)

    layer.heat_efficiency = eta
    return layer


def apply_mechanical_isentropic_offsets(
    layer: NonIdealClosedCycleLayer,
    mechanical_efficiency: float | None = None,
) -> NonIdealClosedCycleLayer:
    """（旧）机械组锚点 + DFS ``PS→η→HP``。须先调用 ``apply_heat_pressure_offsets``。"""
    eta = _resolve_efficiency(
        mechanical_efficiency,
        layer.mechanical_efficiency,
        float(cyges_config.NON_IDEAL_MECHANICAL_EFFICIENCY_DEFAULT),
        "机械等熵效率 η_is",
    )

    nodes = _ensure_nodes(layer)
    edges_by_key = layer.simplified.edges_dict()

    for group in layer.mechanical_groups:
        if not group.edge_keys:
            continue
        anchor = _pick_mechanical_anchor(group, layer.ideal_nodes)
        known = _walk_mechanical_branches(
            nodes=nodes,
            properties=layer.properties,
            group=group,
            edges_by_key=edges_by_key,
            anchor=anchor,
            eta_is=eta,
        )
        missing = group.nodes() - known
        if missing:
            warnings.warn(
                f"机械有向组 {sorted(group.edge_keys)!r} 中节点 {sorted(missing)} "
                f"与基准 {anchor} 在无向意义上不连通，已跳过",
                RuntimeWarning,
                stacklevel=2,
            )

    layer.mechanical_efficiency = eta
    return layer
