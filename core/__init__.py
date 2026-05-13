"""CyGES 核心域模型与闭式循环相关逻辑。

为便于上层调用，将闭式循环层的常用类型与构建函数 re-export 至本包顶层；
原始全路径 ``from core.closed_cycle_layer import ...`` 等仍保持可用。
"""

from core.closed_cycle_layer import (
    ClosedCycleLayer,
    ClosedCycleTPInput,
    Edge,
    Node,
    SkippedPoint,
    SubCycle,
    build_axis,
    build_node_edge_topology,
    build_subcycles,
)
from core.fluid_property_solver import (
    CoolPropFluidPropertySolver,
    FluidPropertySolver,
    PropertyPair,
    ThermoStateTPHS,
)
from core.non_ideal_closed_cycle_layer import (
    BaselineEdgeRecord,
    BaselineSubCycleRecord,
    BaselineTopologySnapshot,
    NonIdealClosedCycleLayer,
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
    "BaselineEdgeRecord",
    "BaselineSubCycleRecord",
    "BaselineTopologySnapshot",
    "NonIdealClosedCycleLayer",
]
