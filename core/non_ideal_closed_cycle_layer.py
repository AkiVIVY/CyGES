"""
非理想闭式循环层：依附于 ``ClosedCycleLayer``，使用冻结基准快照做后续分析；不修改理想层上的 ``nodes`` / ``edges`` / ``subcycles``。
父层在 ``analyze_topology`` 与 ``commit_subcycle_mass_flows_to_topology`` 时会清空 ``non_ideal``，须在理想层稳定后再 ``ensure_non_ideal``。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from core.closed_cycle_layer import Node


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


@dataclass
class NonIdealClosedCycleLayer:
    """
    非理想分析容器；持有 ``BaselineTopologySnapshot``，后续偏移与子图等只写本对象字段，不回写父层。

    由 ``ClosedCycleLayer.ensure_non_ideal()`` 创建并赋给 ``parent.non_ideal``；父层 ``analyze_topology()`` 与 ``commit_subcycle_mass_flows_to_topology()`` 会将其清空。
    """

    baseline: BaselineTopologySnapshot

    @classmethod
    def from_closed_cycle_layer(cls, layer) -> NonIdealClosedCycleLayer:
        return cls(baseline=BaselineTopologySnapshot.from_closed_cycle_layer(layer))
