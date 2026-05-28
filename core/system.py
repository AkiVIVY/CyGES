"""
系统层：整合外部冷热源与闭式循环，统一执行夹点分析与性能统计。
夹点分析独立于管线，由 ``analyze_system_heat()`` 完成。

**分节**：§1 数据模型 → §2 冷热源转换 → §3 管线 → §4 夹点分析。
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
from core.fluid_property_solver import FluidPropertySolver, PropertyRegistry, ThermoStateTPHS
from core.postprocess import (
    PinchResult,
    analyze_pinch,
    split_tq_curve_to_records,
)


# ============================================================
# §1. 数据模型
# ============================================================


@dataclass(frozen=True)
class ExternalSourceInput:
    fluid: str
    mass_flow: float
    T_in: float
    P_in: float
    T_out: float
    P_out: float


@dataclass
class CycleConfig:
    input: ClosedCycleTPInput
    use_non_ideal: bool = False
    subcycle_mass_flows: list[float] = field(default_factory=list)
    delta_T_min: float = 0.0
    heat_method: str | None = "pinch"
    mechanical_method: str | None = None


@dataclass(frozen=True)
class SystemInput:
    heat_sources: tuple[ExternalSourceInput, ...] = ()
    cold_sources: tuple[ExternalSourceInput, ...] = ()
    cycles: tuple[CycleConfig, ...] = ()
    delta_T_min: float = 0.0
    heat_method: str = "pinch"
    mechanical_method: str | None = None


@dataclass(frozen=True)
class SystemResult:
    """系统层一次运行输出的冻结结果。pinch 字段由 ``analyze_system_heat()`` 独立填充。"""

    heat_source_records: tuple[ProcessRecord, ...]
    cold_source_records: tuple[ProcessRecord, ...]
    cycle_reports: tuple[CyclePerformanceReport, ...]
    cycle_pinch: PinchResult | None = None
    system_pinch: PinchResult | None = None
    hot_match_pinch: PinchResult | None = None
    cold_match_pinch: PinchResult | None = None
    source_hot_pinch: PinchResult | None = None
    source_cold_pinch: PinchResult | None = None
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
    in_state: ThermoStateTPHS = props(src.fluid, "TP", src.T_in, src.P_in)
    out_state: ThermoStateTPHS = props(src.fluid, "TP", src.T_out, src.P_out)

    tail_snap = NodeStateSnapshot(index=-1, T=in_state["T"], P=in_state["P"], H=in_state["H"], S=in_state["S"])
    head_snap = NodeStateSnapshot(index=-2, T=out_state["T"], P=out_state["P"], H=out_state["H"], S=out_state["S"])
    delta_H = head_snap.H - tail_snap.H

    return ProcessRecord(
        edge_key=edge_key, fluid=src.fluid, kind="heat", category=category,
        tail=-1, head=-2, mass_flow=src.mass_flow,
        tail_state=tail_snap, head_state=head_snap,
        delta_H=delta_H, power_rate=src.mass_flow * delta_H,
    )


def convert_sources(
    heat_sources: tuple[ExternalSourceInput, ...],
    cold_sources: tuple[ExternalSourceInput, ...],
    props: PropertyRegistry,
) -> tuple[tuple[ProcessRecord, ...], tuple[ProcessRecord, ...]]:
    hr = tuple(_source_to_record(s, f"HS_{i}", ProcessCategory.HEAT_REJECTION, props) for i, s in enumerate(heat_sources))
    cr = tuple(_source_to_record(s, f"CS_{i}", ProcessCategory.HEAT_ABSORPTION, props) for i, s in enumerate(cold_sources))
    return hr, cr


# ============================================================
# §3. 管线
# ============================================================


class SystemPipeline:
    """系统管线：加载输入 → 运行循环 → 转换冷热源。夹点分析由 ``analyze_system_heat()`` 独立完成。"""

    def __init__(self, system_input: SystemInput) -> None:
        self._input = system_input

    def run(self, props: PropertyRegistry, *,
            cycle_properties: FluidPropertySolver | None = None) -> SystemResult:
        inp = self._input
        heat_records, cold_records = convert_sources(inp.heat_sources, inp.cold_sources, props)

        cycle_reports: list[CyclePerformanceReport] = []
        for cfg in inp.cycles:
            layer = ClosedCycleLayer(cfg.input, properties=cycle_properties)
            if cfg.subcycle_mass_flows:
                layer.subcycle_mass_flows = cfg.subcycle_mass_flows
                layer.commit_subcycle_mass_flows_to_topology()
            if cfg.use_non_ideal:
                ni = layer.ensure_non_ideal()
                ni.apply_offsets()
            cycle_reports.append(layer.performance_report())

        return SystemResult(
            heat_source_records=heat_records,
            cold_source_records=cold_records,
            cycle_reports=tuple(cycle_reports),
        )


# ============================================================
# §4. 夹点分析（独立于管线，按需求调用）
# ============================================================


def analyze_system_heat(
    raw_result: SystemResult,
    system_input: SystemInput,
    props: PropertyRegistry,
) -> SystemResult:
    """对系统执行夹点分析。不改变 ``raw_result``——返回新 ``SystemResult``。

    与 ``SystemPipeline.run()`` 解耦：可对同一 ``run()`` 输出调用不同模式。
    """
    heat_records = list(raw_result.heat_source_records)
    cold_records = list(raw_result.cold_source_records)

    # 提取循环换热边
    cycle_all_abs: list[ProcessRecord] = []
    cycle_all_rej: list[ProcessRecord] = []
    for report in raw_result.cycle_reports:
        for _, rec in report.by_edge:
            if rec.kind != "heat" or rec.power_rate is None:
                continue
            if rec.category == ProcessCategory.HEAT_ABSORPTION:
                cycle_all_abs.append(rec)
            elif rec.category == ProcessCategory.HEAT_REJECTION:
                cycle_all_rej.append(rec)

    # 循环内部夹点
    cycle_pinch: PinchResult | None = None
    cycle_unmatched_abs = list(cycle_all_abs)
    cycle_unmatched_rej = list(cycle_all_rej)

    if system_input.cycles:
        cfg = system_input.cycles[0]
        if cfg.heat_method == "pinch" and cycle_all_abs and cycle_all_rej:
            cycle_pinch = analyze_pinch(cycle_all_abs, cycle_all_rej, cfg.delta_T_min, props)
            cycle_unmatched_rej = []
            cycle_unmatched_abs = []
            if cycle_pinch.extra_rejection is not None:
                cycle_unmatched_rej = split_tq_curve_to_records(cycle_pinch.extra_rejection, props)
            if cycle_pinch.extra_absorption is not None:
                cycle_unmatched_abs = split_tq_curve_to_records(cycle_pinch.extra_absorption, props)

    # 系统级夹点
    system_pinch: PinchResult | None = None
    hot_match_pinch: PinchResult | None = None
    cold_match_pinch: PinchResult | None = None
    source_hot_pinch: PinchResult | None = None
    source_cold_pinch: PinchResult | None = None
    method = system_input.heat_method

    if method == "system_pinch":
        all_abs = cold_records + cycle_all_abs
        all_rej = heat_records + cycle_all_rej
        if all_abs and all_rej:
            system_pinch = analyze_pinch(all_abs, all_rej, system_input.delta_T_min, props)

    elif method == "pinch":
        all_abs = cold_records + cycle_unmatched_abs
        all_rej = heat_records + cycle_unmatched_rej
        if all_abs and all_rej:
            system_pinch = analyze_pinch(all_abs, all_rej, system_input.delta_T_min, props)

    elif method == "split_pinch":
        if cycle_unmatched_abs and heat_records:
            hot_match_pinch = analyze_pinch(cycle_unmatched_abs, heat_records, system_input.delta_T_min, props)
        if cold_records and cycle_unmatched_rej:
            cold_match_pinch = analyze_pinch(cold_records, cycle_unmatched_rej, system_input.delta_T_min, props)

    elif method == "source_pinch":
        residual_rej: list[ProcessRecord] = []
        residual_abs: list[ProcessRecord] = []
        if cycle_all_abs and heat_records:
            source_hot_pinch = analyze_pinch(cycle_all_abs, heat_records, system_input.delta_T_min, props)
            if source_hot_pinch.extra_rejection is not None:
                residual_rej.extend(split_tq_curve_to_records(source_hot_pinch.extra_rejection, props))
            if source_hot_pinch.extra_absorption is not None:
                residual_abs.extend(split_tq_curve_to_records(source_hot_pinch.extra_absorption, props))
        else:
            residual_rej = heat_records
            residual_abs = cycle_all_abs
        if cold_records and cycle_all_rej:
            source_cold_pinch = analyze_pinch(cold_records, cycle_all_rej, system_input.delta_T_min, props)
            if source_cold_pinch.extra_rejection is not None:
                residual_rej.extend(split_tq_curve_to_records(source_cold_pinch.extra_rejection, props))
            if source_cold_pinch.extra_absorption is not None:
                residual_abs.extend(split_tq_curve_to_records(source_cold_pinch.extra_absorption, props))
        if residual_abs and residual_rej:
            system_pinch = analyze_pinch(residual_abs, residual_rej, system_input.delta_T_min, props)

    return SystemResult(
        heat_source_records=raw_result.heat_source_records,
        cold_source_records=raw_result.cold_source_records,
        cycle_reports=raw_result.cycle_reports,
        cycle_pinch=cycle_pinch,
        system_pinch=system_pinch,
        hot_match_pinch=hot_match_pinch,
        cold_match_pinch=cold_match_pinch,
        source_hot_pinch=source_hot_pinch,
        source_cold_pinch=source_cold_pinch,
    )
