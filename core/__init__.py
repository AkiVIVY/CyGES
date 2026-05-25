"""CyGES 核心域模型与闭式循环相关逻辑。

拓扑边 ``kind`` 表示热力过程：``mechanical`` 为叶轮机械工作过程（非理想用等熵效率），
``heat`` 为换热过程（非理想用总压恢复系数）。详见 ``config``、``README`` 与非理想层模块说明。

为便于上层调用，将闭式循环层的常用类型与构建函数 re-export 至本包顶层；
原始全路径 ``from core.closed_cycle_layer import ...`` 等仍保持可用。
"""

from __future__ import annotations

from core.closed_cycle_layer import (
    ClosedCycleLayer,
    ClosedCycleTPInput,
    Edge,
    MERGED_ISOLATED_NODE_EDGE_KEY,
    Node,
    SimplifiedEdge,
    SimplifiedTopology,
    SkippedPoint,
    SubCycle,
    build_axis,
    build_node_edge_topology,
    build_simplified_topology,
    build_subcycles,
    filter_topology_for_non_ideal,
)
from core.fluid_property_solver import (
    CoolPropFluidPropertySolver,
    FluidPropertySolver,
    PropertyPair,
    ThermoLookup,
    ThermoStateTPHS,
)
from core.non_ideal_bias import (
    GroupDepthMetrics,
    NonIdealClosedCycleLayer,
    SimplifiedDirectedGroup,
    apply_combined_offsets,
    build_directed_groups,
    compute_group_downstream_depth,
    compute_group_downstream_reach,
    partition_simplified_edges_by_kind,
)
from core.cycle_performance import (
    CyclePerformanceReport,
    CycleTotals,
    NodeStateSnapshot,
    PerformanceContext,
    ProcessCategory,
    ProcessRecord,
    compute_cycle_performance,
    resolve_performance_context,
)
from core.postprocess import (
    EnthalpyLookup,
    HeatTQCurve,
    PinchResult,
    analyze_pinch,
    build_heat_tq_curves,
    compute_pinch,
)
from core.system import (
    CycleConfig,
    ExternalSourceInput,
    SystemInput,
    SystemPipeline,
    SystemResult,
    convert_sources,
)

__all__ = [
    "ClosedCycleLayer",
    "ClosedCycleTPInput",
    "Edge",
    "Node",
    "SkippedPoint",
    "SubCycle",
    "build_axis",
    "build_node_edge_topology",
    "build_subcycles",
    "CoolPropFluidPropertySolver",
    "FluidPropertySolver",
    "PropertyPair",
    "ThermoStateTPHS",
    "MERGED_ISOLATED_NODE_EDGE_KEY",
    "NonIdealClosedCycleLayer",
    "SimplifiedDirectedGroup",
    "build_directed_groups",
    "compute_group_downstream_depth",
    "compute_group_downstream_reach",
    "GroupDepthMetrics",
    "partition_simplified_edges_by_kind",
    "SimplifiedEdge",
    "SimplifiedTopology",
    "build_simplified_topology",
    "filter_topology_for_non_ideal",
    "apply_combined_offsets",
    "CyclePerformanceReport",
    "CycleTotals",
    "EnthalpyLookup",
    "HeatTQCurve",
    "NodeStateSnapshot",
    "PerformanceContext",
    "ProcessCategory",
    "ProcessRecord",
    "compute_cycle_performance",
    "resolve_performance_context",
    "PinchResult",
    "analyze_pinch",
    "build_heat_tq_curves",
    "compute_pinch",
    "CycleConfig",
    "ExternalSourceInput",
    "SystemInput",
    "SystemPipeline",
    "SystemResult",
    "convert_sources",
]
