"""
闭式循环性能数据二次处理：T-Q 曲线构建与夹点分析。

T-Q 曲线在 ``ProcessRecord`` 上按温度节点离散累积；夹点分析通过平移吸热曲线
分离换热匹配与公用工程需求。

**分节**：§1 数据模型 → §2 T-Q 曲线构建 → §3 曲线插值与采样 → §4 夹点计算。
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import config as cyges_config
from core.cycle_performance import CyclePerformanceReport, ProcessCategory, ProcessRecord


EnthalpyLookup = Callable[[str, float, float], float]
"""(fluid: str, T[K], P[kPa]) -> H[kJ/kg]，用于 T-Q 曲线构建时按工质与温压查焓。"""


# ============================================================
# §1. 数据模型
# ============================================================


@dataclass(frozen=True)
class HeatTQCurve:
    """单条换热 T-Q 折线（吸热或放热）。"""

    category: ProcessCategory
    """``HEAT_ABSORPTION`` 或 ``HEAT_REJECTION``。"""
    q_points: tuple[float, ...]
    """累积热量 [kW]（从 0 开始逐段累加）。"""
    t_points: tuple[float, ...]
    """温度 [K]（与 ``q_points`` 一一对应）。"""


@dataclass(frozen=True)
class _PinchAnalysisResult:
    """单次夹点分析的冻结结果（底层，输入为 T-Q 曲线）。

    ``rejection`` / ``absorption`` 分别对应放热曲线 / 吸热曲线；
    夹点分析通过平移吸热曲线使两曲线最小温差为 ``delta_T_min``。
    """

    delta_T_min: float
    """给定夹点温差 [K]。"""
    delta_Q: float
    """吸热曲线平移量 [kW]：正值=吸热右移。"""
    pinch_T_hot: float
    """夹点处放热温度 [K]。"""
    pinch_T_cold: float
    """夹点处吸热温度 [K]（＝ ``pinch_T_hot - delta_T_min``）。"""
    pinch_Q_hot: float
    """夹点处放热 Q 坐标 [kW]。"""
    pinch_Q_cold: float
    """夹点处吸热 Q 坐标 [kW]（原坐标系）。"""

    rejection: HeatTQCurve
    """原放热曲线（未平移）。"""
    absorption_shifted: HeatTQCurve
    """平移后吸热曲线：``q_points[i] + delta_Q``，``t_points`` 不变。"""

    overlap_rejection: HeatTQCurve | None
    """重叠区放热段；无则为 ``None``。"""
    overlap_absorption: HeatTQCurve | None
    """重叠区吸热段（原坐标系）；无则为 ``None``。"""

    unmatched_rejection: HeatTQCurve | None
    """非重叠放热段（需冷公用）；无则为 ``None``。"""
    unmatched_absorption: HeatTQCurve | None
    """非重叠吸热段（需热公用）；无则为 ``None``。"""


@dataclass(frozen=True)
class PinchResult:
    """夹点分析结果：公用工程需求与匹配换热曲线。

    由 ``analyze_pinch`` 直接处理 ``ProcessRecord`` 列表生成；
    是 ``Compute_pinch`` 的上层封装。
    """

    delta_T_min: float
    """给定夹点温差 [K]。"""
    hot_utility_demand: float
    """需外部热源提供的热量 [kW]。"""
    cold_utility_demand: float
    """需外部冷源带走的热量 [kW]。"""
    matched_rejection: HeatTQCurve
    """匹配的放热段（内部换热）。"""
    matched_absorption: HeatTQCurve
    """匹配的吸热段（内部换热）。"""
    extra_absorption: HeatTQCurve | None
    """额外吸热曲线（需外部热源）；低于 ``config.PINCH_EXTRA_CURVE_FRACTION_THRESHOLD`` 时为 ``None``。"""
    extra_rejection: HeatTQCurve | None
    """额外放热曲线（需外部冷源）；低于阈值时为 ``None``。"""
    pinch_T_hot: float
    """夹点处放热温度 [K]。"""
    pinch_T_cold: float
    """夹点处吸热温度 [K]。"""
    delta_Q: float
    """吸热曲线平移量 [kW]。"""


# ============================================================
# §2. T-Q 曲线构建
# ============================================================


def _interp_pressure_by_temperature(rec: ProcessRecord, T: float) -> float:
    """在 ``rec.tail→head`` 的 T 空间线性插值 P，供区间端点 TP 查焓使用。"""
    t0 = float(rec.tail_state.T)
    t1 = float(rec.head_state.T)
    p0 = float(rec.tail_state.P)
    p1 = float(rec.head_state.P)
    if abs(t1 - t0) <= 1e-12:
        return p0
    alpha = (T - t0) / (t1 - t0)
    return p0 + alpha * (p1 - p0)


def _build_heat_tq_curve(
    records: list[ProcessRecord],
    category: ProcessCategory,
    enthalpy_fn: EnthalpyLookup,
) -> HeatTQCurve:
    """按温度节点离散并填充区间能量（TP 查焓），生成单条 T-Q 折线。

    对每条换热过程，在每个重叠温区的两个端点通过 TP 查焓并计算
    ``dq = |mass_flow| * |h(T_b, P_b) - h(T_a, P_a)|``，各区间累计后得到连续折线；
    Q 统一按绝对值，曲线位于 ``Q >= 0`` 右半轴。
    """
    selected = [rec for rec in records if rec.category == category]
    if not selected:
        return HeatTQCurve(category=category, q_points=(0.0,), t_points=(0.0,))

    temp_nodes = sorted(
        {float(rec.tail_state.T) for rec in selected}
        | {float(rec.head_state.T) for rec in selected}
    )
    if len(temp_nodes) == 1:
        return HeatTQCurve(
            category=category,
            q_points=(0.0,),
            t_points=(float(temp_nodes[0]),),
        )

    interval_q = [0.0] * (len(temp_nodes) - 1)
    for rec in selected:
        if rec.mass_flow is None:
            continue
        m_abs = abs(float(rec.mass_flow))
        t0 = float(rec.tail_state.T)
        t1 = float(rec.head_state.T)
        lo = min(t0, t1)
        hi = max(t0, t1)
        if hi - lo <= 1e-12:
            continue
        for i in range(len(temp_nodes) - 1):
            a = temp_nodes[i]
            b = temp_nodes[i + 1]
            overlap = max(0.0, min(hi, b) - max(lo, a))
            if overlap <= 0.0:
                continue
            Ta = max(lo, a)
            Tb = min(hi, b)
            Pa = _interp_pressure_by_temperature(rec, Ta)
            Pb = _interp_pressure_by_temperature(rec, Tb)
            h_a = float(enthalpy_fn(rec.fluid, Ta, Pa))
            h_b = float(enthalpy_fn(rec.fluid, Tb, Pb))
            interval_q[i] += m_abs * abs(h_b - h_a)

    q_points = [0.0]
    for dq in interval_q:
        q_points.append(q_points[-1] + dq)
    return HeatTQCurve(
        category=category,
        q_points=tuple(q_points),
        t_points=tuple(temp_nodes),
    )


def build_heat_tq_curves(
    report: CyclePerformanceReport,
    enthalpy_fn: EnthalpyLookup,
) -> tuple[HeatTQCurve, ...]:
    """从 ``CyclePerformanceReport`` 构建吸热和放热 T-Q 折线。

    :param report: 含 ``by_edge`` 的性能报告。
    :param enthalpy_fn: ``(fluid, T, P) → H``，通常为 ``lambda f, T, P: properties.state("TP", T, P)["H"]``。
    :returns: 一条或两条 ``HeatTQCurve``（吸热/放热）。
    """
    heat_records = [
        rec for _, rec in report.by_edge
        if rec.kind == "heat" and rec.power_rate is not None
    ]
    tq_curves: list[HeatTQCurve] = []
    for cat in (ProcessCategory.HEAT_ABSORPTION, ProcessCategory.HEAT_REJECTION):
        curve = _build_heat_tq_curve(heat_records, cat, enthalpy_fn)
        if len(curve.q_points) > 1 or curve.t_points[0] != 0.0:
            tq_curves.append(curve)
    return tuple(tq_curves)


# ============================================================
# §3. 曲线插值与采样工具
# ============================================================


def _interp_T_at_Q(curve: HeatTQCurve, Q: float) -> float:
    """在分段线性 T-Q 曲线上按 Q 插值 T；超出范围返回最近端点温度。"""
    qs = curve.q_points
    ts = curve.t_points
    n = len(qs)
    if n == 0:
        return 0.0
    if Q <= qs[0]:
        return float(ts[0])
    if Q >= qs[-1]:
        return float(ts[-1])
    for i in range(n - 1):
        if qs[i + 1] >= Q:
            denom = qs[i + 1] - qs[i]
            if denom <= 1e-15:
                return float(ts[i])
            alpha = (Q - qs[i]) / denom
            return float(ts[i] + alpha * (ts[i + 1] - ts[i]))
    return float(ts[-1])


def _interp_Q_at_T(curve: HeatTQCurve, T: float) -> float:
    """在分段线性 T-Q 曲线上按 T 插值 Q；超出范围返回最近端点 Q。"""
    qs = curve.q_points
    ts = curve.t_points
    n = len(ts)
    if n == 0:
        return 0.0
    if T <= ts[0]:
        return float(qs[0])
    if T >= ts[-1]:
        return float(qs[-1])
    for i in range(n - 1):
        if ts[i + 1] >= T:
            denom = ts[i + 1] - ts[i]
            if denom <= 1e-15:
                return float(qs[i])
            alpha = (T - ts[i]) / denom
            return float(qs[i] + alpha * (qs[i + 1] - qs[i]))
    return float(qs[-1])


def _sample_curve_segment(
    curve: HeatTQCurve,
    Q_start: float,
    Q_end: float,
) -> HeatTQCurve | None:
    """截取曲线 ``Q ∈ [Q_start, Q_end]`` 内的区段（含端点插值点）。

    ``Q_start`` / ``Q_end`` 在曲线坐标系内；
    截取区间与曲线无交集时返回 ``None``。
    """
    qs = curve.q_points
    ts = curve.t_points
    n = len(qs)
    if n == 0 or Q_end <= qs[0] or Q_start >= qs[-1] or Q_start >= Q_end:
        return None

    kept_q: list[float] = []
    kept_t: list[float] = []

    if Q_start > qs[0]:
        t_start = _interp_T_at_Q(curve, Q_start)
        kept_q.append(Q_start)
        kept_t.append(t_start)

    for i in range(n):
        if Q_start <= qs[i] <= Q_end:
            kept_q.append(float(qs[i]))
            kept_t.append(float(ts[i]))

    if Q_end < qs[-1]:
        t_end = _interp_T_at_Q(curve, Q_end)
        if kept_q and abs(kept_q[-1] - Q_end) < 1e-12:
            pass
        else:
            kept_q.append(Q_end)
            kept_t.append(t_end)

    if len(kept_q) < 2:
        return None
    return HeatTQCurve(
        category=curve.category,
        q_points=tuple(kept_q),
        t_points=tuple(kept_t),
    )


def _sample_curve_multi_segment(
    curve: HeatTQCurve,
    segments: list[tuple[float, float]],
) -> HeatTQCurve | None:
    """从曲线拼接多个 ``[(Q_start, Q_end), ...]`` 区段为一条新曲线。

    各区段间插入 ``(NaN, NaN)`` 标记断点（供绘图识别）。
    """
    qs = curve.q_points
    ts = curve.t_points
    n = len(qs)
    if n == 0 or not segments:
        return None

    kept_q: list[float] = []
    kept_t: list[float] = []
    first = True
    for seg_start, seg_end in segments:
        if seg_end <= qs[0] or seg_start >= qs[-1] or seg_start >= seg_end:
            continue
        if not first:
            kept_q.append(float("nan"))
            kept_t.append(float("nan"))
        first = False

        if seg_start > qs[0]:
            t_s = _interp_T_at_Q(curve, seg_start)
            kept_q.append(seg_start)
            kept_t.append(t_s)
        for i in range(n):
            if seg_start <= qs[i] <= seg_end:
                kept_q.append(float(qs[i]))
                kept_t.append(float(ts[i]))
        if seg_end < qs[-1]:
            t_e = _interp_T_at_Q(curve, seg_end)
            kept_q.append(seg_end)
            kept_t.append(t_e)

    if len(kept_q) < 2:
        return None
    return HeatTQCurve(
        category=curve.category,
        q_points=tuple(kept_q),
        t_points=tuple(kept_t),
    )


# ============================================================
# §4. 夹点计算
# ============================================================


def compute_pinch(
    rej_curve: HeatTQCurve,
    abs_curve: HeatTQCurve,
    delta_T_min: float = 0.0,
) -> _PinchAnalysisResult:
    """在放热/吸热 T-Q 曲线上执行夹点分析。

    - 将吸热曲线沿 Q 轴平移 ``delta_Q``，使 ``T_rej(Q) ≥ T_abs(Q - delta_Q) + delta_T_min`` 成立。
    - ``delta_Q`` 由两曲线顶点约束候选取最大值得出。
    - 夹点处 ``T_rej - T_abs_shifted = delta_T_min``。
    - 返回包含原曲线、平移后曲线、重叠/非重叠区段的冻结结果。

    前置：两曲线 ``category`` 应为 ``HEAT_REJECTION`` / ``HEAT_ABSORPTION``；
    ``q_points`` 从 0 递增，``t_points`` 升序排列。

    :param rej_curve: 放热曲线（热流）。
    :param abs_curve: 吸热曲线（冷流）。
    :param delta_T_min: 夹点最小温差 [K]，要求 ``≥ 0``。
    """
    if delta_T_min < 0.0:
        raise ValueError(f"delta_T_min 必须 ≥ 0，实际为 {delta_T_min}")

    q_h = rej_curve.q_points
    t_h = rej_curve.t_points
    q_c = abs_curve.q_points
    t_c = abs_curve.t_points

    n_h = len(q_h)
    n_c = len(q_c)

    candidates: list[tuple[float, float, float, float]] = []  # (δQ, pinch_T_hot, pinch_T_cold, pinch_Q_hot)

    for i in range(n_h):
        Q_h = float(q_h[i])
        T_h = float(t_h[i])
        T_target = T_h - delta_T_min
        Q_c_at_target = _interp_Q_at_T(abs_curve, T_target)
        dq = Q_h - Q_c_at_target
        candidates.append((dq, T_h, T_target, Q_h))

    for j in range(n_c):
        Q_c_j = float(q_c[j])
        T_c_j = float(t_c[j])
        T_target = T_c_j + delta_T_min
        Q_h_at_target = _interp_Q_at_T(rej_curve, T_target)
        dq = Q_h_at_target - Q_c_j
        candidates.append((dq, T_target, T_c_j, Q_h_at_target))

    delta_Q, pinch_T_hot, pinch_T_cold, pinch_Q_hot = max(candidates, key=lambda x: x[0])
    pinch_Q_cold = _interp_Q_at_T(abs_curve, pinch_T_cold)

    abs_shifted = HeatTQCurve(
        category=abs_curve.category,
        q_points=tuple(qc + delta_Q for qc in q_c),
        t_points=t_c,
    )

    Q_h_max = float(q_h[-1])
    Q_c_max = float(q_c[-1])
    overlap_start = max(0.0, delta_Q)
    overlap_end = min(Q_h_max, delta_Q + Q_c_max)

    if overlap_start < overlap_end:
        overlap_rej = _sample_curve_segment(rej_curve, overlap_start, overlap_end)
        overlap_abs = _sample_curve_segment(abs_curve, overlap_start - delta_Q, overlap_end - delta_Q)
    else:
        overlap_rej = None
        overlap_abs = None

    unmatched_rej_segs: list[tuple[float, float]] = []
    if delta_Q > 0.0 and delta_Q <= Q_h_max:
        unmatched_rej_segs.append((0.0, min(delta_Q, Q_h_max)))
    cold_end_shifted = delta_Q + Q_c_max
    if cold_end_shifted < Q_h_max:
        unmatched_rej_segs.append((cold_end_shifted, Q_h_max))
    unmatched_rejection = _sample_curve_multi_segment(rej_curve, unmatched_rej_segs)

    unmatched_abs_segs: list[tuple[float, float]] = []
    if delta_Q < 0.0:
        unmatched_abs_segs.append((0.0, min(-delta_Q, Q_c_max)))
    hot_end_shifted = Q_h_max - delta_Q
    if hot_end_shifted < Q_c_max:
        unmatched_abs_segs.append((hot_end_shifted, Q_c_max))
    unmatched_absorption = _sample_curve_multi_segment(abs_curve, unmatched_abs_segs)

    return _PinchAnalysisResult(
        delta_T_min=delta_T_min,
        delta_Q=delta_Q,
        pinch_T_hot=pinch_T_hot,
        pinch_T_cold=pinch_T_cold,
        pinch_Q_hot=pinch_Q_hot,
        pinch_Q_cold=pinch_Q_cold,
        rejection=rej_curve,
        absorption_shifted=abs_shifted,
        overlap_rejection=overlap_rej,
        overlap_absorption=overlap_abs,
        unmatched_rejection=unmatched_rejection,
        unmatched_absorption=unmatched_absorption,
    )


def analyze_pinch(
    abs_records: list[ProcessRecord],
    rej_records: list[ProcessRecord],
    delta_T_min: float,
    enthalpy_fn: EnthalpyLookup,
    *,
    threshold_fraction: float | None = None,
) -> PinchResult:
    """夹点分析入口：接收吸热/放热过程记录与夹点温差，输出公用工程需求与匹配/额外换热曲线。

    - 内部先构建吸热/放热 T-Q 曲线，再执行 ``compute_pinch`` 平移与区段分离。
    - 公用工程需求直接由非重叠区段 Q 跨度计算。
    - ``extra_absorption`` / ``extra_rejection`` 来自夹点平移后的非重叠区段；
      若其热量占该过程总换热的比例低于 ``threshold_fraction`` 则不输出（``None``）。

    :param abs_records: 吸热过程 ``ProcessRecord`` 列表。
    :param rej_records: 放热过程 ``ProcessRecord`` 列表。
    :param delta_T_min: 夹点最小温差 [K]。
    :param enthalpy_fn: ``(T, P) → H``，用于 T-Q 曲线构建时查焓。
    :param threshold_fraction: 额外曲线截断阈值比例 [0,1]；默认取 ``config.PINCH_EXTRA_CURVE_FRACTION_THRESHOLD``。
    """
    if threshold_fraction is None:
        threshold_fraction = cyges_config.PINCH_EXTRA_CURVE_FRACTION_THRESHOLD

    abs_curve = _build_heat_tq_curve(abs_records, ProcessCategory.HEAT_ABSORPTION, enthalpy_fn)
    rej_curve = _build_heat_tq_curve(rej_records, ProcessCategory.HEAT_REJECTION, enthalpy_fn)

    pa = compute_pinch(rej_curve, abs_curve, delta_T_min)

    Q_h_max = float(rej_curve.q_points[-1])
    Q_c_max = float(abs_curve.q_points[-1])
    dq = pa.delta_Q

    cold_utility = max(0.0, dq) + max(0.0, Q_h_max - dq - Q_c_max)
    hot_utility = max(0.0, -dq) + max(0.0, dq + Q_c_max - Q_h_max)

    # 额外曲线：低于阈值比例则不输出
    extra_abs = pa.unmatched_absorption
    if extra_abs is not None and Q_c_max > 1e-12:
        if abs(extra_abs.q_points[-1]) / Q_c_max < threshold_fraction:
            extra_abs = None

    extra_rej = pa.unmatched_rejection
    if extra_rej is not None and Q_h_max > 1e-12:
        if abs(extra_rej.q_points[-1]) / Q_h_max < threshold_fraction:
            extra_rej = None

    return PinchResult(
        delta_T_min=delta_T_min,
        hot_utility_demand=hot_utility,
        cold_utility_demand=cold_utility,
        matched_rejection=pa.overlap_rejection if pa.overlap_rejection is not None
        else HeatTQCurve(category=ProcessCategory.HEAT_REJECTION, q_points=(0.0,), t_points=(0.0,)),
        matched_absorption=pa.overlap_absorption if pa.overlap_absorption is not None
        else HeatTQCurve(category=ProcessCategory.HEAT_ABSORPTION, q_points=(0.0,), t_points=(0.0,)),
        extra_absorption=extra_abs,
        extra_rejection=extra_rej,
        pinch_T_hot=pa.pinch_T_hot,
        pinch_T_cold=pa.pinch_T_cold,
        delta_Q=dq,
    )
