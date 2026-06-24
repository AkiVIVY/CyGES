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
    InterpolatingHeliumSolver,
    PropertyPair,
    PropertyRegistry,
    ThermoStateTPHS,
)
from core.non_ideal_bias import (
    GroupDepthMetrics,
    NonIdealClosedCycleLayer,
    SimplifiedDirectedGroup,
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
    HeatTQCurve,
    PinchResult,
    PinchFixedResult,
    TQSegment,
    analyze_pinch,
    build_heat_tq_curves,
    compute_pinch,
    compute_pinch_fixed_alignment,
    split_tq_curve_to_records,
)
from core.heat_exchanger import (
    HXMatchResult,
    HXUnit,
    match_constructive,
    match_heat_exchanger_groups,
    match_heat_exchanger_staged,
)
from core.system import (
    CycleConfig,
    ExternalSourceInput,
    SystemInput,
    SystemPipeline,
    SystemResult,
    analyze_system_heat,
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
    "PropertyRegistry",
    "InterpolatingHeliumSolver",
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
    "CyclePerformanceReport",
    "CycleTotals",
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
    "TQSegment",
    "split_tq_curve_to_records",
    "CycleConfig",
    "ExternalSourceInput",
    "SystemInput",
    "SystemPipeline",
    "SystemResult",
    "HXUnit",
    "HXMatchResult",
    "match_heat_exchanger_groups",
    "match_constructive",
    "match_heat_exchanger_staged",
    "analyze_system_heat",
    "convert_sources",
]
