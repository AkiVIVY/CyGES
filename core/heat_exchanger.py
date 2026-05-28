"""
换热过程最优匹配：枚举星形换热器候选组（1热+k冷 或 k热+1冷），
逆流检验温差约束后，DP 最优集合打包，最大化匹配热量、最少换热器数。

**分节**：§1 数据模型 → §2 记录归一化 → §3 温度可行性 → §4 候选组枚举
         → §5 DP 打包 → §6 入口函数。

算法概述见 ``docs/architecture.md`` 换热匹配章节。
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

from core.cycle_performance import ProcessCategory, ProcessRecord


# ============================================================
# §1. 数据模型
# ============================================================


@dataclass(frozen=True)
class HXUnit:
    """一个逆流换热器单元：1热+k冷 或 k热+1冷。

    ``matched_heat`` = min(ΣQ_hot, ΣQ_cold) —— 实际交换热量 [kW]。
    ``residual`` = |ΣQ_hot - ΣQ_cold| —— 多出侧的未满足量 [kW]。
    ``internal_pinch`` = 组内所有检验点的最小 ΔT [K]。
    """

    hot_records: tuple[ProcessRecord, ...]
    cold_records: tuple[ProcessRecord, ...]
    matched_heat: float
    residual: float
    internal_pinch: float


@dataclass(frozen=True)
class HXMatchResult:
    """换热匹配完整结果。``total_unmatched`` 为优化目标（越小越好）。"""

    units: tuple[HXUnit, ...]
    unassigned_hots: tuple[ProcessRecord, ...]
    unassigned_colds: tuple[ProcessRecord, ...]
    total_matched: float
    total_unmatched: float
    num_units: int


# ============================================================
# §2. 记录归一化
# ============================================================


@dataclass
class _RecInfo:
    """内部轻量记录：含全局索引（用于 DP bitmask）。"""

    idx: int
    record: ProcessRecord
    is_hot: bool
    T_high: float
    T_low: float
    Q: float


def _normalize_records(
    hot_records: list[ProcessRecord],
    cold_records: list[ProcessRecord],
) -> tuple[list[_RecInfo], float, float]:
    """归一化所有输入记录，分配 idx=0..N-1，返回总热/冷量。"""
    recs: list[_RecInfo] = []
    total_hot = 0.0
    total_cold = 0.0

    for rec in hot_records:
        if rec.power_rate is None:
            continue
        Q = abs(float(rec.power_rate))
        if Q <= 1e-12:
            continue
        T_high = max(float(rec.tail_state.T), float(rec.head_state.T))
        T_low = min(float(rec.tail_state.T), float(rec.head_state.T))
        recs.append(_RecInfo(idx=len(recs), record=rec, is_hot=True,
                             T_high=T_high, T_low=T_low, Q=Q))
        total_hot += Q

    for rec in cold_records:
        if rec.power_rate is None:
            continue
        Q = abs(float(rec.power_rate))
        if Q <= 1e-12:
            continue
        T_high = max(float(rec.tail_state.T), float(rec.head_state.T))
        T_low = min(float(rec.tail_state.T), float(rec.head_state.T))
        recs.append(_RecInfo(idx=len(recs), record=rec, is_hot=False,
                             T_high=T_high, T_low=T_low, Q=Q))
        total_cold += Q

    return recs, total_hot, total_cold


# ============================================================
# §3. 温度可行性检验
# ============================================================


def _check_1hot_Ncold(h: _RecInfo, colds: list[_RecInfo],
                      dT_min: float) -> tuple[bool, float]:
    """检验 1 热 + k 冷的逆流链是否满足温差约束。

    冷按 T_high 降序排列（温度最优）；热从高温端逐段消耗。

    :returns: ``(feasible, min_pinch)``。
    """
    ordered = sorted(colds, key=lambda c: c.T_high, reverse=True)
    cum_q = 0.0
    min_pinch = float("inf")

    for c in ordered:
        if cum_q >= h.Q - 1e-12:
            break

        T_h_in = h.T_high - (h.T_high - h.T_low) * cum_q / h.Q
        remaining_hot = h.Q - cum_q
        if c.Q <= remaining_hot + 1e-12:
            T_h_out = h.T_high - (h.T_high - h.T_low) * (cum_q + c.Q) / h.Q
        else:
            T_h_out = h.T_low

        p1 = T_h_in - c.T_high
        p2 = T_h_out - c.T_low
        if p1 < dT_min - 1e-9 or p2 < dT_min - 1e-9:
            return False, min(p1, p2)
        min_pinch = min(min_pinch, p1, p2)
        cum_q += c.Q

    return True, min_pinch if min_pinch != float("inf") else 0.0


def _check_Nhot_1cold(hots: list[_RecInfo], c: _RecInfo,
                      dT_min: float) -> tuple[bool, float]:
    """检验 m 热 + 1 冷的逆流链是否满足温差约束。

    热按 T_high 降序排列；冷从高温端（Q_c 位置）向低温端消耗。

    :returns: ``(feasible, min_pinch)``。
    """
    ordered = sorted(hots, key=lambda h: h.T_high, reverse=True)
    cum_cold = 0.0
    min_pinch = float("inf")

    for h in ordered:
        if cum_cold >= c.Q - 1e-12:
            break

        pos_hot_end = c.Q - cum_cold
        remaining_cold = c.Q - cum_cold
        actual = min(h.Q, remaining_cold)
        pos_cold_end = pos_hot_end - actual

        T_c_hot = c.T_low + (c.T_high - c.T_low) * max(pos_hot_end, 0.0) / c.Q
        T_c_cold = c.T_low + (c.T_high - c.T_low) * max(pos_cold_end, 0.0) / c.Q

        p1 = h.T_high - T_c_hot
        p2 = h.T_low - T_c_cold
        if p1 < dT_min - 1e-9 or p2 < dT_min - 1e-9:
            return False, min(p1, p2)
        min_pinch = min(min_pinch, p1, p2)
        cum_cold += actual

    return True, min_pinch if min_pinch != float("inf") else 0.0


# ============================================================
# §4. 候选组枚举（含温度预过滤 + 组大小限制）
# ============================================================


def _enumerate_candidate_groups(
    hots: list[_RecInfo],
    colds: list[_RecInfo],
    dT_min: float,
    max_group_size: int | None = None,
) -> list:
    """枚举所有温度可行的星形候选组。

    温度预过滤：仅当逆流必要条件成立时才纳入子集枚举。
    ``max_group_size`` 限制组内非单一侧最大成员数（``None``=不限）。

    :returns: ``[(mask, matched, residual, pinch, hot_recs, cold_recs), ...]``
    """
    candidates: list = []

    # 温度预过滤矩阵：hot[i] 能否匹配 cold[j]（必要条件，仍需全检验）
    n_h, n_c = len(hots), len(colds)
    feasible_hc = [[False] * n_c for _ in range(n_h)]
    for i, h in enumerate(hots):
        for j, c in enumerate(colds):
            # 逆流必要条件：热入口 > 冷出口（首端约束；末端由全检验逐段判定）
            if h.T_high > c.T_high + dT_min - 1e-9:
                feasible_hc[i][j] = True

    # H→C 组：1 热 + k 冷
    for i, h in enumerate(hots):
        feasible_cold_indices = [j for j in range(n_c) if feasible_hc[i][j]]
        max_k = len(feasible_cold_indices)
        if max_group_size is not None:
            max_k = min(max_k, max_group_size)
        for k in range(1, max_k + 1):
            for combo in combinations(feasible_cold_indices, k):
                subset = [colds[j] for j in combo]
                ok, pinch = _check_1hot_Ncold(h, subset, dT_min)
                if not ok:
                    continue
                sum_c = sum(c.Q for c in subset)
                matched = min(h.Q, sum_c)
                residual = abs(h.Q - sum_c)
                mask = (1 << h.idx) | sum((1 << c.idx) for c in subset)
                candidates.append((mask, matched, residual, pinch,
                                   [h], list(subset)))

    # C→H 组：m 热 + 1 冷 (m ≥ 2；m=1 已在 H→C 中覆盖)
    for j, c in enumerate(colds):
        feasible_hot_indices = [i for i in range(n_h) if feasible_hc[i][j]]
        max_k = len(feasible_hot_indices)
        if max_group_size is not None:
            max_k = min(max_k, max_group_size)
        for k in range(2, max_k + 1):
            for combo in combinations(feasible_hot_indices, k):
                subset = [hots[i] for i in combo]
                ok, pinch = _check_Nhot_1cold(subset, c, dT_min)
                if not ok:
                    continue
                sum_h = sum(h.Q for h in subset)
                matched = min(sum_h, c.Q)
                residual = abs(sum_h - c.Q)
                mask = (1 << c.idx) | sum((1 << h.idx) for h in subset)
                candidates.append((mask, matched, residual, pinch,
                                   list(subset), [c]))

    return candidates


# ============================================================
# §5. 多起点随机贪心打包（O(K log K + R·K)，R=10）
# ============================================================


def _greedy_select(candidates: list, order: list[int]) -> tuple[list[int], int, float]:
    """按指定顺序贪心选取不重叠的候选组。

    :returns: ``(selected_indices, assigned_mask, residual_sum)``
    """
    used_mask = 0
    selected: list[int] = []
    residual_sum = 0.0
    for ci in order:
        mask_g = candidates[ci][0]
        if (used_mask & mask_g) != 0:
            continue
        selected.append(ci)
        used_mask |= mask_g
        residual_sum += candidates[ci][2]
    return selected, used_mask, residual_sum


def _solve_optimal_packing(
    N: int,
    candidates: list,
    n_restarts: int = 10,
) -> tuple[list[int], int, float]:
    """多起点随机贪心：按 matched 降序分桶、桶内随机打乱，取最佳匹配热。

    O(K log K + R·K)，R 次 restart 逼近 DP 质量，适合嵌入每轮优化评估。

    :returns: ``(selected_candidate_indices, assigned_mask, residual_sum)``
    """
    import random as _random

    if not candidates:
        return [], 0, 0.0

    # 分桶：按 matched_heat 四舍五入到 0.1kW 分组（同值归入一桶）
    indexed = list(range(len(candidates)))
    indexed.sort(key=lambda ci: candidates[ci][1], reverse=True)  # matched DESC

    buckets: list[list[int]] = []
    cur_bucket: list[int] = []
    cur_matched = None
    for ci in indexed:
        m = round(candidates[ci][1], 1)
        if cur_matched is None:
            cur_matched = m
            cur_bucket = [ci]
        elif abs(m - cur_matched) < 1e-6:
            cur_bucket.append(ci)
        else:
            buckets.append(cur_bucket)
            cur_matched = m
            cur_bucket = [ci]
    if cur_bucket:
        buckets.append(cur_bucket)

    best_masked = 0.0
    best_selected: list[int] = []
    best_mask = 0
    best_resid = 0.0

    rng = _random.Random(42)

    for _ in range(n_restarts):
        order: list[int] = []
        for b in buckets:
            b_copy = list(b)
            rng.shuffle(b_copy)
            order.extend(b_copy)

        sel, mask, resid = _greedy_select(candidates, order)
        total_masked = sum(candidates[ci][1] for ci in sel)
        if total_masked > best_masked + 1e-9:
            best_masked = total_masked
            best_selected = sel
            best_mask = mask
            best_resid = resid

    return best_selected, best_mask, best_resid


# ============================================================
# §6. 入口函数
# ============================================================


def match_heat_exchanger_groups(
    hot_records: list[ProcessRecord],
    cold_records: list[ProcessRecord],
    dT_min: float = 10.0,
    max_group_size: int | None = None,
) -> HXMatchResult:
    """主入口：将所有换热过程组建为星形逆流换热器，最大化匹配热。

    :param hot_records: 放热过程 (``HEAT_REJECTION``) 列表。
    :param cold_records: 吸热过程 (``HEAT_ABSORPTION``) 列表。
    :param dT_min: 最小温差 [K]（默认 10K）。
    :param max_group_size: 组内非单一侧最大成员数（``None``=不限）。
    :returns: 含换热器单元列表、未分配记录与总不匹配量的冻结结果。
    """
    recs, total_hot, total_cold = _normalize_records(hot_records, cold_records)
    N = len(recs)

    if N == 0:
        return HXMatchResult(
            units=(), unassigned_hots=(), unassigned_colds=(),
            total_matched=0.0, total_unmatched=0.0, num_units=0,
        )

    hots = [r for r in recs if r.is_hot]
    colds = [r for r in recs if not r.is_hot]

    candidates = _enumerate_candidate_groups(hots, colds, dT_min, max_group_size)

    if not candidates:
        ua_hots = tuple(r.record for r in hots)
        ua_colds = tuple(r.record for r in colds)
        return HXMatchResult(
            units=(), unassigned_hots=ua_hots, unassigned_colds=ua_colds,
            total_matched=0.0, total_unmatched=total_hot + total_cold,
            num_units=0,
        )

    selected_indices, assigned_mask, base_unmatched = _solve_optimal_packing(N, candidates)

    # 计算未分配记录的 Q
    for i in range(N):
        if not (assigned_mask & (1 << i)):
            base_unmatched += recs[i].Q

    # 构建 HXUnit 列表
    units: list[HXUnit] = []
    for ci in selected_indices:
        _, matched, residual, pinch, h_recs, c_recs = candidates[ci]
        units.append(HXUnit(
            hot_records=tuple(r.record for r in h_recs),
            cold_records=tuple(r.record for r in c_recs),
            matched_heat=matched,
            residual=residual,
            internal_pinch=pinch,
        ))

    # 收集未分配记录
    ua_hot: list[ProcessRecord] = []
    ua_cold: list[ProcessRecord] = []
    for r in hots:
        if not (assigned_mask & (1 << r.idx)):
            ua_hot.append(r.record)
    for r in colds:
        if not (assigned_mask & (1 << r.idx)):
            ua_cold.append(r.record)

    total_matched = sum(u.matched_heat for u in units)

    return HXMatchResult(
        units=tuple(units),
        unassigned_hots=tuple(ua_hot),
        unassigned_colds=tuple(ua_cold),
        total_matched=total_matched,
        total_unmatched=base_unmatched,
        num_units=len(units),
    )
