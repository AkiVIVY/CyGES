"""
系统层：整合外部冷热源与闭式循环，统一执行夹点分析与性能统计。

**分节**：§1 数据模型 → §2 冷热源转换 → §3 管线。
"""

from __future__ import annotations

from dataclasses import dataclass, field

from core.closed_cycle_layer import ClosedCycleLayer, ClosedCycleTPInput
from core.cycle_performance import (
    CyclePerformanceReport,
    NodeStateSnapshot,
    ProcessCategory,
    ProcessRecord,
)
from core.fluid_property_solver import PropertyRegistry, ThermoStateTPHS
from core.postprocess import (
    PinchResult,
    analyze_pinch,
)


# ============================================================
# §1. 数据模型
# ============================================================


@dataclass(frozen=True)
class ExternalSourceInput:
    """单个外部热源或冷源边界条件。

    构造时仅需 fluid、mass_flow 与进出口温度/压力；
    H / S 由 ``convert_sources`` 通过物性求解器补全。
    """

    fluid: str
    mass_flow: float
    T_in: float
    P_in: float
    T_out: float
    P_out: float


@dataclass
class CycleConfig:
    """单个闭式循环配置（可变）。"""

    input: ClosedCycleTPInput
    use_non_ideal: bool = False
    subcycle_mass_flows: list[float] = field(default_factory=list)
    internal_pinch_dT: float = 0.0
    """循环内部换热夹点温差 [K]；0 表示不做内部夹点。"""


@dataclass(frozen=True)
class SystemInput:
    """系统层输入：外部冷热源 + 闭式循环 + 系统级夹点温差。"""

    heat_sources: tuple[ExternalSourceInput, ...] = ()
    cold_sources: tuple[ExternalSourceInput, ...] = ()
    cycles: tuple[CycleConfig, ...] = ()
    delta_T_min: float = 0.0
    """系统级夹点最小温差 [K]；0 表示不执行系统级夹点。"""


@dataclass(frozen=True)
class SystemResult:
    """系统层一次运行输出的冻结结果。"""

    heat_source_records: tuple[ProcessRecord, ...]
    cold_source_records: tuple[ProcessRecord, ...]
    cycle_reports: tuple[CyclePerformanceReport, ...]
    cycle_pinch: PinchResult | None
    """循环内部夹点结果；``None`` = 未执行或无边。"""
    system_pinch: PinchResult | None
    """系统级夹点结果；``None`` = 未执行或无边。"""
    objective: float | None = None


# ============================================================
# §2. 冷热源转换
# ============================================================


def _source_to_record(
    src: ExternalSourceInput,
    edge_key: str,
    category: ProcessCategory,
    props: PropertyRegistry,
) -> ProcessRecord:
    """将单个外部源转换为 ``ProcessRecord``。

    通过 ``thermo_fn(fluid, "TP", T, P)`` 查询进出口完整状态（T/P/H/S）。
    index 用负值以区分拓扑节点。
    """
    in_state: ThermoStateTPHS = props(src.fluid, "TP", src.T_in, src.P_in)
    out_state: ThermoStateTPHS = props(src.fluid, "TP", src.T_out, src.P_out)

    tail_snap = NodeStateSnapshot(
        index=-1,
        T=in_state["T"],
        P=in_state["P"],
        H=in_state["H"],
        S=in_state["S"],
    )
    head_snap = NodeStateSnapshot(
        index=-2,
        T=out_state["T"],
        P=out_state["P"],
        H=out_state["H"],
        S=out_state["S"],
    )
    delta_H = head_snap.H - tail_snap.H

    return ProcessRecord(
        edge_key=edge_key,
        fluid=src.fluid,
        kind="heat",
        category=category,
        tail=-1,
        head=-2,
        mass_flow=src.mass_flow,
        tail_state=tail_snap,
        head_state=head_snap,
        delta_H=delta_H,
        power_rate=src.mass_flow * delta_H,
    )


def convert_sources(
    heat_sources: tuple[ExternalSourceInput, ...],
    cold_sources: tuple[ExternalSourceInput, ...],
    props: PropertyRegistry,
) -> tuple[tuple[ProcessRecord, ...], tuple[ProcessRecord, ...]]:
    """将外部冷/热源批量转换为 ``ProcessRecord`` 元组。

    热源 → ``HEAT_REJECTION``，冷源 → ``HEAT_ABSORPTION``。
    ``props(fluid, "TP", T, P)`` 须返回含 T/P/H/S 的完整状态。
    """
    hr: list[ProcessRecord] = []
    for i, src in enumerate(heat_sources):
        hr.append(_source_to_record(src, f"HS_{i}", ProcessCategory.HEAT_REJECTION, thermo_fn))

    cr: list[ProcessRecord] = []
    for i, src in enumerate(cold_sources):
        cr.append(_source_to_record(src, f"CS_{i}", ProcessCategory.HEAT_ABSORPTION, thermo_fn))

    return tuple(hr), tuple(cr)


# ============================================================
# §3. 管线
# ============================================================


class SystemPipeline:
    """系统管线：加载输入 → 运行循环 → 转换冷热源 → 夹点分析。

    一次 ``run()`` 返回冻结的 ``SystemResult``；可复用实例多次调用。
    """

    def __init__(self, system_input: SystemInput) -> None:
        self._input = system_input

    def run(self, thermo_fn: ThermoLookup) -> SystemResult:
        """执行系统级分析。

        :param props: 多工质物性注册表，用于 ``props(fluid, pair, x, y)`` 查询完整状态。
        """
        inp = self._input

        # 1. 转换外部冷热源
        heat_records, cold_records = convert_sources(
            inp.heat_sources, inp.cold_sources, thermo_fn,
        )

        # 2. 处理闭式循环
        cycle_reports: list[CyclePerformanceReport] = []
        cycle_all_rej: list[ProcessRecord] = []
        cycle_all_abs: list[ProcessRecord] = []

        for cfg in inp.cycles:
            layer = ClosedCycleLayer(cfg.input)
            if cfg.subcycle_mass_flows:
                layer.subcycle_mass_flows = cfg.subcycle_mass_flows
                layer.commit_subcycle_mass_flows_to_topology()

            if cfg.use_non_ideal:
                ni = layer.ensure_non_ideal()
                ni.apply_offsets()

            report = layer.performance_report()
            cycle_reports.append(report)

            for _, rec in report.by_edge:
                if rec.kind != "heat" or rec.power_rate is None:
                    continue
                if rec.category == ProcessCategory.HEAT_REJECTION:
                    cycle_all_rej.append(rec)
                elif rec.category == ProcessCategory.HEAT_ABSORPTION:
                    cycle_all_abs.append(rec)

        # 3. 循环内部夹点
        cycle_pinch: PinchResult | None = None
        cycle_unmatched_rej: list[ProcessRecord] = list(cycle_all_rej)
        cycle_unmatched_abs: list[ProcessRecord] = list(cycle_all_abs)

        if inp.cycles:
            cfg = inp.cycles[0]
            if cfg.internal_pinch_dT > 0.0 and cycle_all_rej and cycle_all_abs:
                cycle_pinch = analyze_pinch(
                    cycle_all_abs, cycle_all_rej, cfg.internal_pinch_dT, props,
                )

        # 4. 系统级夹点：循环换热 + 外部冷热源
        system_pinch: PinchResult | None = None
        if inp.delta_T_min > 0.0:
            all_abs = list(cold_records) + cycle_unmatched_abs
            all_rej = list(heat_records) + cycle_unmatched_rej
            if all_abs and all_rej:
                system_pinch = analyze_pinch(
                    all_abs, all_rej, inp.delta_T_min, props,
                )

        return SystemResult(
            heat_source_records=heat_records,
            cold_source_records=cold_records,
            cycle_reports=tuple(cycle_reports),
            cycle_pinch=cycle_pinch,
            system_pinch=system_pinch,
            objective=None,
        )
