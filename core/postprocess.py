"""
闭式循环性能数据二次处理：夹点分析（首项功能）。

在 ``HeatTQCurve`` 上做放热/吸热曲线平移、夹点定位与区段分离。
输入曲线来自 ``core.cycle_performance`` 的 ``CyclePerformanceReport.heat_tq_curves``。

**分节**：§1 数据模型 → §2 曲线插值与采样 → §3 夹点计算。
"""

from __future__ import annotations

from dataclasses import dataclass

from core.cycle_performance import HeatTQCurve, ProcessCategory


# ============================================================
# §1. 数据模型
# ============================================================


@dataclass(frozen=True)
class PinchAnalysisResult:
    """单次夹点分析的冻结结果。

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


# ============================================================
# §2. 曲线插值与采样工具
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

    # 左端点插值
    if Q_start > qs[0]:
        t_start = _interp_T_at_Q(curve, Q_start)
        kept_q.append(Q_start)
        kept_t.append(t_start)

    # 保留区间内原有点
    for i in range(n):
        if Q_start <= qs[i] <= Q_end:
            kept_q.append(float(qs[i]))
            kept_t.append(float(ts[i]))

    # 右端点插值
    if Q_end < qs[-1]:
        t_end = _interp_T_at_Q(curve, Q_end)
        if kept_q and abs(kept_q[-1] - Q_end) < 1e-12:
            pass  # 右端点与原有点重合，不重复添加
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
# §3. 夹点计算
# ============================================================


def compute_pinch(
    rej_curve: HeatTQCurve,
    abs_curve: HeatTQCurve,
    delta_T_min: float = 0.0,
) -> PinchAnalysisResult:
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

    # 收集所有候选 δQ
    candidates: list[tuple[float, float, float, float]] = []  # (δQ, pinch_T_hot, pinch_T_cold, pinch_Q_hot)

    # 候选 1: 遍历放热顶点
    for i in range(n_h):
        Q_h = float(q_h[i])
        T_h = float(t_h[i])
        T_target = T_h - delta_T_min
        Q_c_at_target = _interp_Q_at_T(abs_curve, T_target)
        dq = Q_h - Q_c_at_target
        candidates.append((dq, T_h, T_target, Q_h))

    # 候选 2: 遍历吸热顶点
    for j in range(n_c):
        Q_c_j = float(q_c[j])
        T_c_j = float(t_c[j])
        T_target = T_c_j + delta_T_min
        Q_h_at_target = _interp_Q_at_T(rej_curve, T_target)
        dq = Q_h_at_target - Q_c_j
        candidates.append((dq, T_target, T_c_j, Q_h_at_target))

    # 取最大 δQ；并列取首次出现
    delta_Q, pinch_T_hot, pinch_T_cold, pinch_Q_hot = max(candidates, key=lambda x: x[0])
    pinch_Q_cold = _interp_Q_at_T(abs_curve, pinch_T_cold)

    # 平移后吸热曲线
    abs_shifted = HeatTQCurve(
        category=abs_curve.category,
        q_points=tuple(qc + delta_Q for qc in q_c),
        t_points=t_c,
    )

    # 平移后重叠区：放热 [0, Q_h_max]，吸热 [delta_Q, delta_Q + Q_c_max]
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

    # 非重叠区段
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

    return PinchAnalysisResult(
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
