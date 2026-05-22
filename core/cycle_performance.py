"""
闭式循环**性能统计**（只读第三层）：在精简拓扑上归类热/机械过程并汇总循环指标。

在 ``SimplifiedTopology`` 与节点状态（理想或非理想）上推导过程归类、边级记录与循环汇总；
不修改拓扑、不参与 ``analyze_topology`` / ``commit_*`` 失效链。

**分节顺序**：§1 数据模型 → §2 判据小工具 → §3 状态解析 → §4 统计计算。

公式与判据见 ``docs/architecture.md`` §10。
"""

from __future__ import annotations

import warnings
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING, Literal, Mapping

from core.closed_cycle_layer import Node, SimplifiedEdge, SimplifiedTopology

if TYPE_CHECKING:
    from core.closed_cycle_layer import ClosedCycleLayer
    from core.non_ideal_bias import NonIdealClosedCycleLayer


# ============================================================
# §1. 数据模型（过程类别 / 边记录 / 报告 / 输入视图）
# ============================================================


class ProcessCategory(str, Enum):
    """精简边 ``tail→head`` 上的热力学过程类别。"""

    COMPRESSION = "compression"
    """机械边且 ``P_head > P_tail``。"""
    EXPANSION = "expansion"
    """机械边且 ``P_head < P_tail``。"""
    HEAT_ABSORPTION = "heat_absorption"
    """换热边且 ``H_head > H_tail``（吸热）。"""
    HEAT_REJECTION = "heat_rejection"
    """换热边且 ``H_head < H_tail``（放热）。"""
    ISOBARIC_MECHANICAL = "isobaric_mechanical"
    """机械边且 ``P_head ≈ P_tail``（容差内等压）。"""


@dataclass(frozen=True)
class NodeStateSnapshot:
    """``kept_nodes`` 上某节点的状态快照（用于报告与端点关联）。"""

    index: int
    T: float
    P: float
    H: float
    S: float


@dataclass(frozen=True)
class ProcessRecord:
    """单条精简边（``SM*`` / ``SH*``）的性能记录。"""

    edge_key: str
    fluid: str
    """工质（CoolProp 流体名）。"""
    kind: Literal["mechanical", "heat"]
    category: ProcessCategory
    tail: int
    head: int
    mass_flow: float | None
    """沿 ``tail→head`` 的质量流量 [kg/s]；链上全 ``None`` 则为 ``None``。"""
    tail_state: NodeStateSnapshot
    head_state: NodeStateSnapshot
    delta_H: float
    """``H_head - H_tail`` [kJ/kg]。"""
    power_rate: float | None
    """``ṁ·ΔH`` [kW]；``mass_flow`` 为 ``None`` 时为 ``None``。"""


@dataclass(frozen=True)
class CycleTotals:
    """循环级汇总（精简边 ``power_rate`` 代数和；``mass_flow`` 为 ``None`` 的边不计入）。"""

    net_mechanical_power: float
    """所有机械边 ``ṁ·ΔH`` 之和 [kW]（压缩为正、膨胀为负）。"""
    net_heat_rate: float
    """所有换热边 ``ṁ·ΔH`` 之和 [kW]（吸热为正、放热为负）。"""
    compression_power: float
    """压缩过程 ``ṁ·ΔH`` 之和 [kW]。"""
    expansion_power: float
    """膨胀过程 ``ṁ·ΔH`` 之和 [kW]。"""
    heat_absorption_rate: float
    """吸热过程 ``ṁ·ΔH`` 之和 [kW]。"""
    heat_rejection_rate: float
    """放热过程 ``ṁ·ΔH`` 之和 [kW]。"""


@dataclass(frozen=True)
class CyclePerformanceReport:
    """一次性能统计的冻结结果。"""

    source: Literal["ideal", "non_ideal"]
    """节点状态来源：理想 ``layer.nodes`` 或非理想 ``NonIdealClosedCycleLayer.nodes``。"""
    by_edge: tuple[tuple[str, ProcessRecord], ...]
    """边键 → 单条过程记录（按边键排序）。"""
    by_category: tuple[tuple[ProcessCategory, tuple[ProcessRecord, ...]], ...]
    """按 ``ProcessCategory`` 分组的过程列表。"""
    nodes: tuple[tuple[int, NodeStateSnapshot], ...]
    """``kept_nodes`` 上的节点快照（按 index 排序）。"""
    cycle_totals: CycleTotals
    skipped_edges: tuple[str, ...]
    """预留：因异常跳过的边键（当前恒空）。"""


@dataclass(frozen=True)
class PerformanceContext:
    """``compute_cycle_performance`` 的输入视图。"""

    simplified: SimplifiedTopology
    nodes: Mapping[int, Node]
    source: Literal["ideal", "non_ideal"]
    fluid: str
    """工质（CoolProp 流体名）；由 ``resolve_performance_context`` 自动从层属性注入。"""


# ============================================================
# §2. 判据小工具（容差比较 / 快照 / 过程归类 / 功率）
# ============================================================


def _pressure_greater(P_tail: float, P_head: float) -> bool:
    """``P_head`` 是否显著大于 ``P_tail``。"""
    tol = max(1e-9, 1e-6 * max(abs(P_tail), abs(P_head)))
    return P_head - P_tail > tol


def _pressure_less(P_tail: float, P_head: float) -> bool:
    """``P_tail`` 是否显著大于 ``P_head``。"""
    tol = max(1e-9, 1e-6 * max(abs(P_tail), abs(P_head)))
    return P_tail - P_head > tol


def _enthalpy_greater(H_tail: float, H_head: float) -> bool:
    """``H_head`` 是否显著大于 ``H_tail``。"""
    tol = max(1e-12, 1e-9 * max(abs(H_tail), abs(H_head)))
    return H_head - H_tail > tol


def _enthalpy_less(H_tail: float, H_head: float) -> bool:
    """``H_tail`` 是否显著大于 ``H_head``。"""
    tol = max(1e-12, 1e-9 * max(abs(H_tail), abs(H_head)))
    return H_tail - H_head > tol


def _node_snapshot(node: Node) -> NodeStateSnapshot:
    """从 ``Node`` 提取不可变快照。"""
    return NodeStateSnapshot(index=node.index, T=node.T, P=node.P, H=node.H, S=node.S)


def _classify_edge(edge: SimplifiedEdge, tail_n: Node, head_n: Node) -> ProcessCategory:
    """沿 ``tail→head`` 判定过程类别；判据见 ``docs/architecture.md`` §10.2。"""
    if edge.kind == "mechanical":
        if _pressure_greater(tail_n.P, head_n.P):
            return ProcessCategory.COMPRESSION
        if _pressure_less(tail_n.P, head_n.P):
            return ProcessCategory.EXPANSION
        return ProcessCategory.ISOBARIC_MECHANICAL

    if _enthalpy_greater(tail_n.H, head_n.H):
        return ProcessCategory.HEAT_ABSORPTION
    if _enthalpy_less(tail_n.H, head_n.H):
        return ProcessCategory.HEAT_REJECTION
    warnings.warn(
        f"精简换热边 {edge.tail!r}→{edge.head!r} 的 ΔH 在容差内为零，归为吸热。",
        RuntimeWarning,
        stacklevel=3,
    )
    return ProcessCategory.HEAT_ABSORPTION


def _power_rate(mass_flow: float | None, delta_H: float) -> float | None:
    """``ṁ·ΔH`` [kW]；``mass_flow`` 为 ``None`` 时返回 ``None``。"""
    if mass_flow is None:
        return None
    return float(mass_flow) * delta_H


# ============================================================
# §3. 状态解析（理想 / 非理想节点表选取）
# ============================================================


def resolve_performance_context(
    layer: ClosedCycleLayer,
    *,
    non_ideal: NonIdealClosedCycleLayer | None = None,
) -> PerformanceContext:
    """从理想层（及可选非理想层）解析统计用 ``simplified`` 与节点表。

    - ``simplified``：``non_ideal.simplified``（若传入或挂载）否则 ``layer.simplified``。
    - ``nodes``：若 ``non_ideal.nodes`` 已写入（已 ``apply_combined_offsets``）→ 非理想；
      否则 ``layer.nodes``。
    - ``fluid``：从 ``layer.properties.fluid`` 获取。
    - 流量始终来自 ``SimplifiedEdge.mass_flow``（偏置不改变流量）。
    - 前置：``layer.simplified is not None``。

    返回 ``PerformanceContext``；不修改任何输入对象。
    """
    if layer.simplified is None:
        raise RuntimeError(
            "ClosedCycleLayer.simplified 尚未生成；须先 analyze_topology() 或 commit_subcycle_mass_flows_to_topology()。"
        )

    ni = non_ideal if non_ideal is not None else layer.non_ideal
    simplified = ni.simplified if ni is not None else layer.simplified

    if ni is not None and ni.nodes is not None:
        return PerformanceContext(
            simplified=simplified,
            nodes=ni.nodes,
            source="non_ideal",
            fluid=layer.properties.fluid,
        )

    return PerformanceContext(
        simplified=simplified,
        nodes=layer.nodes,
        source="ideal",
        fluid=layer.properties.fluid,
    )


# ============================================================
# §4. 统计计算（边级记录 + 分类汇总 + 循环总量）
# ============================================================


def compute_cycle_performance(ctx: PerformanceContext) -> CyclePerformanceReport:
    """在 ``ctx.simplified`` 的每条精简边上归类并计算性能量。

    仅 ``kept_nodes`` 参与；``merged_into`` 中的中间点不重复计功/热。
    边端点缺失于 ``ctx.nodes`` 时抛 ``KeyError``。
    """
    simp = ctx.simplified
    nodes_map = ctx.nodes

    node_snaps: dict[int, NodeStateSnapshot] = {}
    for idx in sorted(simp.kept_nodes):
        node_snaps[idx] = _node_snapshot(nodes_map[idx])

    by_edge: dict[str, ProcessRecord] = {}
    by_cat: dict[ProcessCategory, list[ProcessRecord]] = defaultdict(list)

    totals = {
        "net_mech": 0.0,
        "net_heat": 0.0,
        "compression": 0.0,
        "expansion": 0.0,
        "heat_abs": 0.0,
        "heat_rej": 0.0,
    }

    for edge_key, edge in simp.simplified_edges:
        tail_n = nodes_map[edge.tail]
        head_n = nodes_map[edge.head]
        category = _classify_edge(edge, tail_n, head_n)
        delta_H = head_n.H - tail_n.H
        pr = _power_rate(edge.mass_flow, delta_H)

        rec = ProcessRecord(
            edge_key=edge_key,
            fluid=ctx.fluid,
            kind=edge.kind,
            category=category,
            tail=edge.tail,
            head=edge.head,
            mass_flow=edge.mass_flow,
            tail_state=node_snaps[edge.tail],
            head_state=node_snaps[edge.head],
            delta_H=delta_H,
            power_rate=pr,
        )
        by_edge[edge_key] = rec
        by_cat[category].append(rec)

        if pr is not None:
            if edge.kind == "mechanical":
                totals["net_mech"] += pr
                if category == ProcessCategory.COMPRESSION:
                    totals["compression"] += pr
                elif category == ProcessCategory.EXPANSION:
                    totals["expansion"] += pr
            else:
                totals["net_heat"] += pr
                if category == ProcessCategory.HEAT_ABSORPTION:
                    totals["heat_abs"] += pr
                elif category == ProcessCategory.HEAT_REJECTION:
                    totals["heat_rej"] += pr

    cycle_totals = CycleTotals(
        net_mechanical_power=totals["net_mech"],
        net_heat_rate=totals["net_heat"],
        compression_power=totals["compression"],
        expansion_power=totals["expansion"],
        heat_absorption_rate=totals["heat_abs"],
        heat_rejection_rate=totals["heat_rej"],
    )

    by_category_sorted = tuple(
        (cat, tuple(by_cat[cat]))
        for cat in ProcessCategory
        if cat in by_cat
    )

    return CyclePerformanceReport(
        source=ctx.source,
        by_edge=tuple(sorted(by_edge.items())),
        by_category=by_category_sorted,
        nodes=tuple(sorted(node_snaps.items())),
        cycle_totals=cycle_totals,
        skipped_edges=(),
    )
