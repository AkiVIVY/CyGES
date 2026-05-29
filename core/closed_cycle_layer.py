"""
闭式循环**理想层**：离散 TP/PS 拓扑、子循环流量、活跃子图精简，并为非理想层提供基准。

``ClosedCycleLayer`` 持有 ``nodes`` / ``edges`` / ``subcycles`` 与 ``simplified``；不实现非理想效率
或状态偏移，后者由 :mod:`core.non_ideal_bias` 在 ``ensure_non_ideal()`` 时快照
``simplified`` 后承接。

失效语义：``analyze_topology()`` 与 ``commit_subcycle_mass_flows_to_topology()`` 会重建
``simplified`` 并清空 ``non_ideal``；非理想分析须在理想层稳定后再 ``ensure_non_ideal()``。

单位约定（字段名不缀单位）：T[K]、P[kPa]、H[kJ/kg]、S[kJ/(kg·K)]。分位序列由调用方保证
互不重复；本模块不对 ``ClosedCycleTPInput`` 做防御性校验。

算法细节（PS 平面约定、拓扑流水线、子循环模板、链合并规则、流量汇聚）见
``docs/architecture.md``。

**分节顺序**与 ``analyze_topology`` 流水线一致：数据模型 → 小工具 → 基础拓扑构建
（``build_node_edge_topology``）→ 子循环（``build_subcycles``）→ 精简拓扑
（``build_simplified_topology``）→ ``ClosedCycleLayer``。
"""

from __future__ import annotations

import math
import warnings
from collections import defaultdict
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING, Literal, Sequence

if TYPE_CHECKING:
    # 仅类型检查导入，避免与 non_ideal_bias（反向 import 本模块）循环依赖；
    # 运行时 ``ensure_non_ideal()`` 内再延迟 import 并构造实例。
    from core.non_ideal_bias import NonIdealClosedCycleLayer

import config as cyges_config

from core.fluid_property_solver import CoolPropFluidPropertySolver, FluidPropertySolver


# ============================================================
# §1. 数据模型（输入 / 节点 / 边 / 子循环 / 精简拓扑）
# ============================================================

MERGED_ISOLATED_NODE_EDGE_KEY = "__CyGES_ISOLATED__"
"""精简拓扑中四邻边槽全空、且未上链的节点占位键：在 ``SimplifiedTopology.merged_into`` 中
作为「不参与活跃子图」的标记，无对应 ``SimplifiedEdge``。"""


@dataclass(frozen=True)
class ClosedCycleTPInput:
    """闭式循环层拓扑输入：含工质，可单独作为创建该层的完整信息。"""

    fluid: str
    """工质标识，与物性后端（如 CoolProp）流体名一致。"""
    t_min: float
    t_max: float
    p_min: float
    p_max: float
    t_quantiles: Sequence[float] = field(default_factory=tuple)
    """[0,1] 分位，在 [t_min, t_max] 内线性插值；须互不重复（与端点亦不重复）。"""
    p_quantiles: Sequence[float] = field(default_factory=tuple)
    """[0,1] 分位，在 [p_min, p_max] 内线性插值；须互不重复（与端点亦不重复）。"""
    subcycle_mass_flow_initial: float | None = None
    """子循环初始质量流量 [kg/s]；``None`` 时取 ``config.SUBCYCLE_INITIAL_MASS_FLOW_DEFAULT``。"""


@dataclass(frozen=True)
class Node:
    """
    拓扑状态节点；一级（TP 网格）与二级（等熵延伸）共用本结构，用 ``parent`` 区分。

    - ``index``：全局唯一序号，亦作 ``nodes`` 字典键。
    - ``T, P, H, S``：状态量（项目约定单位）。
    - ``parent``：一级为 ``None``；二级为所属一级节点的 ``index``。
    - ``edge_up`` / ``edge_down`` / ``edge_left`` / ``edge_right``：PS 平面上邻接边在 ``edges`` 中的键；拓扑生成末段由 ``tail``/``head`` 与边 ``kind`` 填入，初值为 ``None``。
    """

    index: int
    T: float
    P: float
    H: float
    S: float
    parent: int | None = None
    edge_up: str | None = None
    edge_down: str | None = None
    edge_left: str | None = None
    edge_right: str | None = None


@dataclass
class Edge:
    """
    拓扑快照中的有向边（非物性计算边）。

    - ``kind``：过程类型，与 PS 离散方向对应——
      ``"mechanical"``：**叶轮机械工作过程**（压缩/膨胀等），在 PS 平面上沿压力轴相邻离散；
      ``"heat"``：**换热过程**（加热/冷却等），在 PS 平面上沿等压线（熵 ``S`` 方向）相邻离散。
    - ``tail`` / ``head``：端点节点序号，方向为尾 → 头（工质沿此方向流经该过程）。
    - ``mass_flow``：质量流 [kg/s] 等，初始 ``None``；经 ``ClosedCycleLayer.assign_edge_mass_flows_from_subcycles`` 后为沿 ``tail → head`` 的代数和（未被子循环覆盖的边保持 ``None``）。

    非理想修正所用效率见 ``config``（机械边：等熵效率；换热边：总压恢复系数）及非理想层。
    """

    kind: Literal["mechanical", "heat"]
    tail: int
    head: int
    mass_flow: float | None = None


@dataclass
class SubCycle:
    """
    最小子循环（4 节点、4 边）。非 ``frozen``，便于后续写入流量等。

    - ``nodes``：``(左下, 左上, 右上, 右下)`` 顺时针，各元素为 ``Node.index``。
    - ``edges``：``(左, 上, 右, 下)`` 各为 ``edges`` 字典键；左/右为**机械边**（叶轮机械过程），上/下为**换热边**（换热过程）（由 ``build_subcycles`` 走边保证）。
    - ``mass_flow``：子循环质量流 [kg/s] 等标量，可为负；由 ``ClosedCycleLayer`` 的 ``subcycle_mass_flows`` 同步写入（见 ``sync_subcycle_mass_flows_to_subcycles``）。
    """

    nodes: tuple[int, int, int, int]
    edges: tuple[str, str, str, str]
    mass_flow: float | None = None


@dataclass(frozen=True)
class SkippedPoint:
    """
    拓扑构建中被物性求解器异常静默跳过的候选点；仅用于诊断，不参与拓扑/流量计算。

    - ``stage``：``"primary_TP"`` 为一级网格点构建阶段（按 ``(T, P)`` 求 ``H, S``）；
      ``"secondary_PS"`` 为沿一级节点等熵延伸阶段（按 ``(P, S)`` 求 ``T, H``）。
    - ``T``：``stage == "primary_TP"`` 时为请求的 T；``stage == "secondary_PS"`` 时为 ``None``（T 由 PS 解出，失败时不存在）。
    - ``P``：请求时的压力。
    - ``S``：``stage == "secondary_PS"`` 时为父节点 S；``stage == "primary_TP"`` 时为 ``None``。
    - ``reason``：触发跳过的异常字符串表示，便于复盘是否落在包线/相区外等。
    """

    stage: Literal["primary_TP", "secondary_PS"]
    T: float | None
    P: float
    S: float | None
    reason: str


@dataclass(frozen=True)
class SimplifiedEdge:
    """
    精简边：过滤后子图中，一段**同类型**原始边链合并成的单条有向边。

    - ``kind``：``"mechanical"``（叶轮机械过程链，沿 P）或 ``"heat"``（换热过程链，沿 S）。
    - ``tail`` / ``head``：保留端点 ``Node.index``；有向为 ``tail → head``。
      由链上 ``mass_flow`` 符号规范化：正值沿 PS 正向（机械 P 升、换热 S 升），负值则反向并取绝对值。
    - ``constituent_edges``：被合并的原始 ``M*`` / ``H*`` 键，按 ``tail → head`` 顺序。
    - ``merged_nodes``：链内被吞并的中间点 index，顺序与 ``constituent_edges`` 一致；``len = len(edges) - 1``。
    - ``mass_flow``：链上聚合流量，**恒非负**（方向已体现在 tail/head）；链上全 ``None`` 则为 ``None``。

    非理想参数不在此结构存储：机械链用等熵效率、换热链用总压恢复系数，见 ``config.NON_IDEAL_*`` 与非理想层。
    """

    kind: Literal["mechanical", "heat"]
    tail: int
    head: int
    constituent_edges: tuple[str, ...]
    merged_nodes: tuple[int, ...]
    mass_flow: float | None


@dataclass(frozen=True)
class SimplifiedTopology:
    """
    活跃子图上的精简拓扑快照（``frozen``，可安全共享引用）。

    - ``kept_nodes``：仍参与精简图的节点 index（链端点等；**不**含 ``merged_into`` 中的中间点与孤立占位点）。
    - ``simplified_edges``：``("SM1", SimplifiedEdge(...))`` 等；机械/换热分开编号。
    - ``merged_into``：被合并进某条精简边的中间点 → 对应 ``SM*``/``SH*``；或
      ``MERGED_ISOLATED_NODE_EDGE_KEY``（过滤后四邻边槽皆空、且未出现在任何 ``merged_nodes`` 中的点）。

    由 ``ClosedCycleLayer._rebuild_simplified()`` 在每次 analyze/commit 后写入；``ensure_non_ideal()`` 读取此刻快照。
    """

    kept_nodes: frozenset[int]
    simplified_edges: tuple[tuple[str, SimplifiedEdge], ...]
    merged_into: tuple[tuple[int, str], ...]

    def edges_dict(self) -> dict[str, SimplifiedEdge]:
        """便捷访问：把精简边以 ``dict`` 形式返回（构造时为 ``tuple`` 以维持冻结性）。"""
        return dict(self.simplified_edges)

    def merged_dict(self) -> dict[int, str]:
        """便捷访问：节点 ``index`` → 所在精简边键，或 ``MERGED_ISOLATED_NODE_EDGE_KEY``（过滤后四邻边全空且未上链）。"""
        return dict(self.merged_into)


# ============================================================
# §2. 流量与边小工具（容差比较、子循环边集、链上 mass_flow 聚合）
# ============================================================


def _mass_flow_close(a: float | None, b: float | None) -> bool:
    """链上 ``mass_flow`` 一致性比较：``None`` 视为缺测，与任意值兼容；两个非 ``None`` 用
    ``math.isclose``。"""
    if a is None or b is None:
        return True
    return math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12)


def _subcycle_edge_key_set(subcycles: Sequence[SubCycle]) -> frozenset[str]:
    """所有子循环 ``SubCycle.edges`` 中出现过的原始边键并集。"""
    out: set[str] = set()
    for sc in subcycles:
        for ek in sc.edges:
            out.add(ek)
    return frozenset(out)


def _edge_has_nonzero_mass_flow(e: Edge) -> bool:
    """边是否视为「有流量」：``mass_flow`` 非 ``None`` 且绝对值在容差意义下非零。"""
    if e.mass_flow is None:
        return False
    return not math.isclose(float(e.mass_flow), 0.0, rel_tol=0.0, abs_tol=1e-12)


# ============================================================
# §3. 基础拓扑构建（轴、节点、定向边、邻槽回写）
# ============================================================


def _attach_edges_to_nodes_ps(nodes: dict[int, Node], edges: dict[str, Edge]) -> dict[int, Node]:
    """将边键写回各节点 PS 四邻槽。

    机械边沿 P 向上：``tail.edge_up`` / ``head.edge_down``；换热边沿 S 向右：
    ``tail.edge_right`` / ``head.edge_left``。同槽重复时抛 ``ValueError``。
    """
    slot: dict[int, dict[str, str | None]] = {
        i: {"up": None, "down": None, "left": None, "right": None} for i in nodes
    }

    def put(ni: int, side: str, ek: str) -> None:
        cur = slot[ni][side]
        if cur is not None:
            raise ValueError(f"node {ni} duplicate {side} edge: {cur!r} and {ek!r}")
        slot[ni][side] = ek

    for ek, e in edges.items():
        # 机械：tail 在 head 下方（P 较低）→ tail.edge_up / head.edge_down
        if e.kind == "mechanical":
            put(e.tail, "up", ek)
            put(e.head, "down", ek)
        else:
            # 换热：tail 在 head 左侧（S 较低）→ tail.edge_right / head.edge_left
            put(e.tail, "right", ek)
            put(e.head, "left", ek)

    return {
        i: replace(
            n,
            edge_up=slot[i]["up"],
            edge_down=slot[i]["down"],
            edge_left=slot[i]["left"],
            edge_right=slot[i]["right"],
        )
        for i, n in nodes.items()
    }


def build_axis(min_v: float, max_v: float, quantiles: Sequence[float]) -> list[float]:
    """由 ``[min_v, max_v]`` 及分位线性插值生成排序后的采样点序列。

    调用方须保证 ``quantiles`` 与端点互不重复；本函数不做去重。
    """
    span = max_v - min_v
    pts = [min_v, max_v, *[min_v + q * span for q in quantiles]]
    return sorted(pts)


def build_node_edge_topology(
    solver: FluidPropertySolver,
    inp: ClosedCycleTPInput,
) -> tuple[dict[int, Node], dict[str, Edge], list[SkippedPoint]]:
    """一次性生成基础拓扑：一级 TP 网格 → 等熵二级点 → 机械边 ``M*`` → 等压换热边 ``H*`` → 邻槽回写。

    返回 ``(nodes, edges, skipped_points)``；``skipped_points`` 仅供诊断。流水线细节见
    ``docs/architecture.md §拓扑构建流水线``。
    """

    def oriented_edge(
        kind: Literal["mechanical", "heat"],
        na: Node,
        nb: Node,
        *,
        rtol_p: float = 1e-9,
        rtol_s: float = 1e-12,
    ) -> Edge:
        """在容差下判定 tail→head，使 P、S 沿边非减（PS 平面向上且向右）。"""
        eps_p = rtol_p * max(1.0, abs(na.P), abs(nb.P))
        eps_s = rtol_s * max(1.0, abs(na.S), abs(nb.S))
        a_sw = na.P <= nb.P + eps_p and na.S <= nb.S + eps_s
        b_sw = nb.P <= na.P + eps_p and nb.S <= na.S + eps_s
        if a_sw and not b_sw:
            return Edge(kind=kind, tail=na.index, head=nb.index, mass_flow=None)
        if b_sw and not a_sw:
            return Edge(kind=kind, tail=nb.index, head=na.index, mass_flow=None)
        raise RuntimeError(
            f"PS 方向不可判定：节点 {na.index!r}/{nb.index!r} 互不严格小于；"
            "上游算法（机械边按 P 升序成链、换热边按 S 升序分桶）不应抵达此分支"
        )

    skipped_points: list[SkippedPoint] = []

    # —— 一级 TP 网格：每个 (T,P) 闪蒸得 H,S，失败点记入 skipped_points 不建节点 ——
    t_list = build_axis(inp.t_min, inp.t_max, inp.t_quantiles)
    p_list = build_axis(inp.p_min, inp.p_max, inp.p_quantiles)
    grid: list[Node] = []
    idx = 0
    for Tk in t_list:
        for Pk in p_list:
            try:
                st = solver.state("TP", Tk, Pk)
            except Exception as exc:
                skipped_points.append(
                    SkippedPoint(stage="primary_TP", T=Tk, P=Pk, S=None, reason=repr(exc))
                )
                continue
            grid.append(Node(index=idx, T=Tk, P=Pk, H=st["H"], S=st["S"], parent=None))
            idx += 1

    # —— 熵序警告：(t_max, p_max) 的熵不应低于 (t_min, p_min)，否则 PS 平面子循环矩形可能不连通 ——
    def _node_matches(n: Node, T: float, P: float) -> bool:
        return abs(n.T - T) < 1e-9 * max(1.0, abs(T)) and abs(n.P - P) < 1e-9 * max(1.0, abs(P))

    hot_hi: Node | None = None
    cold_lo: Node | None = None
    for node in grid:
        if _node_matches(node, inp.t_max, inp.p_max):
            hot_hi = node
        if _node_matches(node, inp.t_min, inp.p_min):
            cold_lo = node
        if hot_hi is not None and cold_lo is not None:
            break
    if hot_hi is not None and cold_lo is not None and hot_hi.S < cold_lo.S:
        warnings.warn(
            f"S(t_max={inp.t_max:.1f}K, p_max={inp.p_max:.1f}kPa)={hot_hi.S:.3f} < "
            f"S(t_min={inp.t_min:.1f}K, p_min={inp.p_min:.1f}kPa)={cold_lo.S:.3f}: "
            "熵序反转可能导致子循环无法连通",
            RuntimeWarning,
        )

    # —— 二级等熵 + 机械边：沿每个一级点的 S 在压力轴上延伸，同 parent 的二级点按 P 排序连成 M* ——
    p_axis = build_axis(inp.p_min, inp.p_max, inp.p_quantiles)
    secondary: list[Node] = []
    mechanical: dict[str, Edge] = {}
    m = 0
    idx = len(grid)

    for prim in grid:
        batch: list[Node] = []
        for p2 in p_axis:
            if abs(p2 - prim.P) < 1e-9 * max(1.0, abs(p2)):
                continue
            try:
                st = solver.state("PS", p2, prim.S)
                t2 = st["T"]
            except Exception as exc:
                skipped_points.append(
                    SkippedPoint(stage="secondary_PS", T=None, P=p2, S=prim.S, reason=repr(exc))
                )
                continue
            if not (inp.t_min <= t2 <= inp.t_max):
                continue
            node = Node(
                index=idx,
                T=t2,
                P=p2,
                H=st["H"],
                S=prim.S,
                parent=prim.index,
            )
            secondary.append(node)
            batch.append(node)
            idx += 1

        # 一级 + 其二级子点按 P 排序，相邻对生成一条机械边（压缩或膨胀由 P 大小决定方向）
        chain = sorted([prim, *batch], key=lambda n: (n.P, n.index))
        for k in range(len(chain) - 1):
            a, b = chain[k], chain[k + 1]
            m += 1
            mechanical[f"M{m}"] = oriented_edge("mechanical", a, b)

    nodes: dict[int, Node] = {}
    for n in grid:
        nodes[n.index] = n
    for n in secondary:
        nodes[n.index] = n

    # —— 换热边：等压线上按 S 排序，相邻节点连成 H*（沿熵增方向） ——
    by_p: dict[float, list[Node]] = defaultdict(list)
    for n in nodes.values():
        by_p[n.P].append(n)
    heat: dict[str, Edge] = {}
    h = 0
    for lst in by_p.values():
        lst.sort(key=lambda n: (n.S, n.index))
        for k in range(len(lst) - 1):
            a, b = lst[k], lst[k + 1]
            h += 1
            heat[f"H{h}"] = oriented_edge("heat", a, b)

    all_edges = {**mechanical, **heat}
    # 将边键写回各节点 edge_up/down/left/right，供 build_subcycles 沿邻槽走环
    nodes = _attach_edges_to_nodes_ps(nodes, all_edges)
    return nodes, all_edges, skipped_points


# ============================================================
# §4. 子循环枚举（最小 4 节点环 / 顺时针模板）
# ============================================================


def build_subcycles(nodes: dict[int, Node], edges: dict[str, Edge]) -> list[SubCycle]:
    """枚举所有最小 4 节点子循环；同物理单元至多输出一次。

    走法：``edge_up`` → ``edge_right`` → ``edge_down``（取该边 ``tail``） → ``edge_left``（取
    ``tail``）须回到起点。模板细节见 ``docs/architecture.md §子循环模板``。
    """
    seen_cells: set[frozenset[int]] = set()
    out: list[SubCycle] = []

    for n0 in nodes.values():
        # 固定走法：左下 n0 → e_left(M↑) → 左上 n1 → e_top(H→) → 右上 n2
        #           → e_right(M↓) → 右下 n3 → e_bottom(H←) → 回到 n0
        ku = n0.edge_up
        if ku is None:
            continue
        e_left = edges[ku]
        if e_left.kind != "mechanical":
            continue
        n1 = nodes[e_left.head]
        kr = n1.edge_right
        if kr is None:
            continue
        e_top = edges[kr]
        if e_top.kind != "heat":
            continue
        n2 = nodes[e_top.head]
        kd = n2.edge_down
        if kd is None:
            continue
        e_right = edges[kd]
        if e_right.kind != "mechanical":
            continue
        n3 = nodes[e_right.tail]
        kl = n3.edge_left
        if kl is None:
            continue
        e_bottom = edges[kl]
        if e_bottom.kind != "heat":
            continue
        if e_bottom.tail != n0.index:
            continue

        idxs = (n0.index, n1.index, n2.index, n3.index)
        if len(set(idxs)) != 4:
            continue

        # 四角集合去重：同一物理单元只输出一个 SubCycle
        cell = frozenset(idxs)
        if cell in seen_cells:
            continue
        seen_cells.add(cell)

        out.append(
            SubCycle(
                nodes=(n0.index, n1.index, n2.index, n3.index),
                edges=(ku, kr, kd, kl),
            )
        )

    return out


# ============================================================
# §5. 精简拓扑算法（过滤 + 同类型链合并 + 占位）
# ----------------------------------------------------------------
# 详见 docs/architecture.md §精简拓扑：过滤规则、切分点、方向规范化与 MERGED_ISOLATED_NODE_EDGE_KEY。
# ============================================================


def filter_topology_for_non_ideal(
    nodes: dict[int, Node],
    edges: dict[str, Edge],
    subcycles: Sequence[SubCycle],
) -> tuple[dict[int, Node], dict[str, Edge]]:
    """剔除不在任何子循环或 ``mass_flow`` 为 ``None`` / 零的边，返回拷贝（**不修改入参**）。

    返回 ``(nodes, edges)``：``nodes`` 为 ``Node`` 浅拷贝且邻边槽指向被剔除边时置 ``None``；
    ``edges`` 为仅含保留边的新字典，``Edge`` 对象与输入共享引用。
    """
    # 子循环边键并集：不在任何 4 元环上的边不参与非理想精简
    in_sub = _subcycle_edge_key_set(subcycles)
    kept_keys: set[str] = set()
    for ek in in_sub:
        if ek not in edges:
            continue
        # 仅保留 mass_flow 有定义且非零的边（零流量视为该过程未激活）
        if _edge_has_nonzero_mass_flow(edges[ek]):
            kept_keys.add(ek)

    filtered_edges = {k: edges[k] for k in sorted(kept_keys)}

    # 节点浅拷贝：邻边槽若指向已剔除边则置 None，避免后续链遍历误入死边
    new_nodes: dict[int, Node] = {}
    for i, n in nodes.items():
        eu = n.edge_up if n.edge_up in filtered_edges else None
        ed = n.edge_down if n.edge_down in filtered_edges else None
        el = n.edge_left if n.edge_left in filtered_edges else None
        er = n.edge_right if n.edge_right in filtered_edges else None
        if (eu, ed, el, er) == (n.edge_up, n.edge_down, n.edge_left, n.edge_right):
            new_nodes[i] = n
        else:
            new_nodes[i] = replace(n, edge_up=eu, edge_down=ed, edge_left=el, edge_right=er)

    return new_nodes, filtered_edges


def _aggregate_chain_mass_flow(
    chain_edges: list[str],
    edges: dict[str, Edge],
    chain_label: str,
) -> float | None:
    """
    聚合链上一组原始边的 ``mass_flow``，供精简边方向规范化使用。

    链上全 ``None`` 则返回 ``None``；出现多个非 ``None`` 且数值不一致时抛 ``ValueError``
    （同一精简段内原始边流量应已由子循环汇聚一致）。
    """
    seen: float | None = None
    for ek in chain_edges:
        m = edges[ek].mass_flow
        if m is None:
            continue
        if seen is None:
            seen = float(m)
        elif not _mass_flow_close(seen, float(m)):
            raise ValueError(
                f"精简{chain_label}链 {chain_edges!r} 上 mass_flow 不一致：含 {seen!r} 与 {m!r}"
            )
    return seen


def _chain_neighbor_slots(kind: Literal["mechanical", "heat"]) -> tuple[str, str]:
    """该 ``kind`` 在 ``Node`` 上用于无向链遍历的一对邻边槽名。"""
    if kind == "mechanical":
        return "edge_up", "edge_down"
    return "edge_left", "edge_right"


def _chain_sort_coordinate(nodes: dict[int, Node], kind: Literal["mechanical", "heat"], index: int) -> float:
    """链上节点按 PS 正向排序用的标量（机械 ``P``、换热 ``S``）。"""
    return nodes[index].P if kind == "mechanical" else nodes[index].S


def _find_typed_chains(
    nodes: dict[int, Node],
    edges: dict[str, Edge],
    kind: Literal["mechanical", "heat"],
) -> list[dict]:
    """返回 ``kind`` 类型边的全部简单链，每项 ``{"nodes": [...], "edges": [...]}``。

    ``nodes`` 按 PS 正向（机械 ``P`` 升，换热 ``S`` 升）排序；``edges`` 为相邻节点对之间的边
    键，长度 ``= len(nodes) - 1``。PS 网格保证连通分量均为无分叉简单链；找不到相邻边时抛
    ``RuntimeError``。
    """
    slot_a, slot_b = _chain_neighbor_slots(kind)
    sort_key = lambda i: _chain_sort_coordinate(nodes, kind, i)

    # 仅用该 kind 的邻边槽建无向邻接（每条原始边在两端各出现一次）
    adj: dict[int, list[str]] = defaultdict(list)
    for i, n in nodes.items():
        for slot in (slot_a, slot_b):
            ek = getattr(n, slot)
            if ek is not None and edges[ek].kind == kind:
                adj[i].append(ek)

    visited: set[int] = set()
    chains: list[dict] = []
    for start in nodes:
        if start in visited or not adj[start]:
            continue
        # 无向 DFS 收集一个连通分量（PS 网格保证为简单链、无分叉）
        component: set[int] = set()
        stack = [start]
        while stack:
            v = stack.pop()
            if v in component:
                continue
            component.add(v)
            for ek in adj[v]:
                e = edges[ek]
                other = e.head if e.tail == v else e.tail
                if other not in component:
                    stack.append(other)
        visited |= component
        # 按 PS 正向排序后，相邻节点对之间应恰有一条该 kind 的边
        ordered_nodes = sorted(component, key=sort_key)
        chain_edges_in_order: list[str] = []
        for k in range(len(ordered_nodes) - 1):
            a, b = ordered_nodes[k], ordered_nodes[k + 1]
            found = None
            for ek in adj[a]:
                e = edges[ek]
                if (e.tail == a and e.head == b) or (e.tail == b and e.head == a):
                    found = ek
                    break
            if found is None:
                raise RuntimeError(
                    f"{kind} 链上相邻节点 {a}/{b} 之间未找到边；PS 网格连通性被破坏"
                )
            chain_edges_in_order.append(found)
        chains.append({"nodes": ordered_nodes, "edges": chain_edges_in_order})
    return chains


def _orient_chain_segment(
    ps_low_idx: int,
    ps_high_idx: int,
    aggregated_mf: float | None,
    seg_edges_in_ps_order: list[str],
    seg_merged_in_ps_order: list[int],
) -> tuple[int, int, float | None, tuple[str, ...], tuple[int, ...]]:
    """按聚合 ``mass_flow`` 符号确定精简段 ``tail→head`` 与元组顺序（沿 PS 正向或反向）。"""
    if aggregated_mf is None or aggregated_mf == 0.0 or aggregated_mf > 0:
        return (
            ps_low_idx,
            ps_high_idx,
            aggregated_mf,
            tuple(seg_edges_in_ps_order),
            tuple(seg_merged_in_ps_order),
        )
    return (
        ps_high_idx,
        ps_low_idx,
        -aggregated_mf,
        tuple(reversed(seg_edges_in_ps_order)),
        tuple(reversed(seg_merged_in_ps_order)),
    )


def _merge_chains_of_kind(
    chains: list[dict],
    *,
    nodes_f: dict[int, Node],
    edges_f: dict[str, Edge],
    kind: Literal["mechanical", "heat"],
    is_mergeable: dict[int, bool],
    counter: int,
) -> tuple[list[tuple[str, SimplifiedEdge]], list[tuple[int, str]], int]:
    """将同 ``kind`` 的全部简单链合并为精简边并累加 ``merged_into`` 条目。"""
    simplified_edges: list[tuple[str, SimplifiedEdge]] = []
    merged_into: list[tuple[int, str]] = []
    for chain in chains:
        segs, mlinks, counter = _simplify_chain(
            chain,
            nodes_f,
            edges_f,
            kind=kind,
            is_mergeable=is_mergeable,
            counter=counter,
        )
        simplified_edges.extend(segs)
        merged_into.extend(mlinks)
    return simplified_edges, merged_into, counter


def _simplify_chain(
    chain: dict,
    nodes: dict[int, Node],
    edges: dict[str, Edge],
    kind: Literal["mechanical", "heat"],
    is_mergeable: dict[int, bool],
    counter: int,
) -> tuple[list[tuple[str, SimplifiedEdge]], list[tuple[int, str]], int]:
    """将一条简单链按「不可合并节点」切分为若干 ``SimplifiedEdge``，并返回 ``merged_into`` 条目。

    切分点 = 链两端 + 链上 ``is_mergeable[v] is False`` 的节点；切分点之间的原始边串合并为
    一条精简边，方向由聚合 ``mass_flow`` 规范化（详见 ``docs/architecture.md §精简拓扑``）。
    """
    chain_nodes: list[int] = chain["nodes"]
    chain_edges: list[str] = chain["edges"]

    if len(chain_nodes) < 2:
        return [], [], counter

    # 切分点 = 链端点 + 链内不可合并节点（仍挂两类邻边或孤立，须保留为独立状态点）
    cut_positions: list[int] = [0]
    for k in range(1, len(chain_nodes) - 1):
        if not is_mergeable[chain_nodes[k]]:
            cut_positions.append(k)
    cut_positions.append(len(chain_nodes) - 1)

    prefix = "SM" if kind == "mechanical" else "SH"
    chain_label = "机械" if kind == "mechanical" else "换热"

    simp_edges: list[tuple[str, SimplifiedEdge]] = []
    merged_entries: list[tuple[int, str]] = []

    for seg in range(len(cut_positions) - 1):
        i0, i1 = cut_positions[seg], cut_positions[seg + 1]
        ps_low_idx = chain_nodes[i0]
        ps_high_idx = chain_nodes[i1]
        seg_merged_in_ps_order = chain_nodes[i0 + 1 : i1]
        seg_edges_in_ps_order = chain_edges[i0:i1]

        aggregated_mf = _aggregate_chain_mass_flow(seg_edges_in_ps_order, edges, chain_label)
        tail, head, mf_out, ordered_edges, ordered_merged = _orient_chain_segment(
            ps_low_idx,
            ps_high_idx,
            aggregated_mf,
            seg_edges_in_ps_order,
            seg_merged_in_ps_order,
        )

        counter += 1
        simp_key = f"{prefix}{counter}"
        simp_edges.append(
            (
                simp_key,
                SimplifiedEdge(
                    kind=kind,
                    tail=tail,
                    head=head,
                    constituent_edges=ordered_edges,
                    merged_nodes=ordered_merged,
                    mass_flow=mf_out,
                ),
            )
        )
        for mn in ordered_merged:
            merged_entries.append((mn, simp_key))

    return simp_edges, merged_entries, counter


def build_simplified_topology(
    nodes: dict[int, Node],
    edges: dict[str, Edge],
    subcycles: Sequence[SubCycle],
) -> SimplifiedTopology:
    """在过滤后的子图上做同类型链合并，返回 ``SimplifiedTopology``（**不修改入参**）。

    流程：``filter_topology_for_non_ideal`` → 标记仅挂一类邻边的可合并中间点 →
    ``_find_typed_chains`` 拆链 → ``_simplify_chain`` 切段并规范方向 → 四邻全空且未上链的节点
    记入 ``MERGED_ISOLATED_NODE_EDGE_KEY`` 占位。``kept_nodes`` 为 ``nodes_f`` 中不在
    ``merged_into`` 键集合里的节点。

    精简边方向、流量聚合、占位语义详见 ``docs/architecture.md §精简拓扑``。同精简段内原始边
    ``mass_flow`` 不一致时抛 ``ValueError``。
    """
    nodes_f, edges_f = filter_topology_for_non_ideal(nodes, edges, subcycles)

    # 仅挂一类邻边的中间点可吞并进精简边；两类都有或四邻皆空的点作切分/占位
    is_mech_only: dict[int, bool] = {}
    is_heat_only: dict[int, bool] = {}
    for i, n in nodes_f.items():
        has_mech = n.edge_up is not None or n.edge_down is not None
        has_heat = n.edge_left is not None or n.edge_right is not None
        is_mech_only[i] = has_mech and not has_heat
        is_heat_only[i] = has_heat and not has_mech

    simplified_edges: list[tuple[str, SimplifiedEdge]] = []
    merged_into: list[tuple[int, str]] = []

    mech_segs, mech_mlinks, _ = _merge_chains_of_kind(
        _find_typed_chains(nodes_f, edges_f, "mechanical"),
        nodes_f=nodes_f,
        edges_f=edges_f,
        kind="mechanical",
        is_mergeable=is_mech_only,
        counter=0,
    )
    simplified_edges.extend(mech_segs)
    merged_into.extend(mech_mlinks)

    heat_segs, heat_mlinks, _ = _merge_chains_of_kind(
        _find_typed_chains(nodes_f, edges_f, "heat"),
        nodes_f=nodes_f,
        edges_f=edges_f,
        kind="heat",
        is_mergeable=is_heat_only,
        counter=0,
    )
    simplified_edges.extend(heat_segs)
    merged_into.extend(heat_mlinks)

    # 过滤后四邻全空、且未上链的节点：记入 merged_into 占位，不出现在 kept_nodes
    merged_idx_so_far = {mn for mn, _ in merged_into}
    for i, n in nodes_f.items():
        if i in merged_idx_so_far:
            continue
        if (
            n.edge_up is None
            and n.edge_down is None
            and n.edge_left is None
            and n.edge_right is None
        ):
            merged_into.append((i, MERGED_ISOLATED_NODE_EDGE_KEY))

    merged_set: set[int] = {mn for mn, _ in merged_into}
    kept_nodes = frozenset(i for i in nodes_f if i not in merged_set)

    return SimplifiedTopology(
        kept_nodes=kept_nodes,
        simplified_edges=tuple(simplified_edges),
        merged_into=tuple(merged_into),
    )


# ============================================================
# §6. 入口类 ClosedCycleLayer
# ============================================================


class ClosedCycleLayer:
    """闭式循环**理想层**入口：拓扑构建、子循环流量、精简图与非理想快照挂载点。

    字段
    ----
    - ``nodes`` / ``edges``：完整基准拓扑（``M*`` / ``H*``）。
    - ``subcycles`` / ``subcycle_mass_flows``：最小 4 元环及其流量向量（优化器优先改后者）。
    - ``skipped_points``：物性失败点，仅诊断。
    - ``simplified``：最近一次 analyze/commit 后的 ``SimplifiedTopology``；无子循环时为空骨架。
    - ``non_ideal``：``ensure_non_ideal()`` 挂载的快照层；analyze/commit 时置 ``None``。

    构造仅需 ``ClosedCycleTPInput``；默认 ``properties=CoolPropFluidPropertySolver(fluid)``，
    ``auto_analyze=True`` 时构造末尾执行一次 ``analyze_topology()``。
    """

    def __init__(
        self,
        inp: ClosedCycleTPInput,
        *,
        properties: FluidPropertySolver | None = None,
        auto_analyze: bool = True,
    ) -> None:
        self.input = inp
        self.properties = properties if properties is not None else CoolPropFluidPropertySolver(inp.fluid)
        self.nodes: dict[int, Node] = {}
        self.edges: dict[str, Edge] = {}
        self.subcycles: list[SubCycle] = []
        self.subcycle_mass_flows: list[float] = []
        self.skipped_points: list[SkippedPoint] = []
        self.simplified: SimplifiedTopology | None = None
        self.non_ideal: NonIdealClosedCycleLayer | None = None
        if auto_analyze:
            self.analyze_topology()

    def _invalidate_non_ideal(self) -> None:
        """理想层基准发生变化时清空非理想层挂载；外部不直接调用。"""
        self.non_ideal = None

    def _rebuild_simplified(self) -> None:
        """基于当前 ``nodes`` / ``edges`` / ``subcycles`` 重建 ``self.simplified``。

        清空既有 ``non_ideal`` 快照；``subcycles`` 为空时直接置空骨架并发出 ``RuntimeWarning``。
        """
        self._invalidate_non_ideal()
        if not self.subcycles:
            warnings.warn(
                "ClosedCycleLayer: 当前没有子循环，simplified 置为空骨架（kept_nodes / simplified_edges / merged_into 均为空）。",
                RuntimeWarning,
                stacklevel=2,
            )
            self.simplified = SimplifiedTopology(
                kept_nodes=frozenset(),
                simplified_edges=(),
                merged_into=(),
            )
            return
        self.simplified = build_simplified_topology(self.nodes, self.edges, self.subcycles)
        if not self.simplified.simplified_edges:
            warnings.warn(
                "ClosedCycleLayer: 子循环存在但精简后无边（可能所有边 mass_flow 为零），"
                "simplified 置为空骨架。",
                RuntimeWarning,
                stacklevel=2,
            )

    def sync_subcycle_mass_flows_to_subcycles(self) -> None:
        """将 ``subcycle_mass_flows[i]`` 写入 ``subcycles[i].mass_flow``；长度须一致。"""
        if len(self.subcycle_mass_flows) != len(self.subcycles):
            raise ValueError(
                f"subcycle_mass_flows 长度 {len(self.subcycle_mass_flows)} "
                f"与 subcycles 长度 {len(self.subcycles)} 不一致"
            )
        for i, sc in enumerate(self.subcycles):
            sc.mass_flow = self.subcycle_mass_flows[i]

    def assign_edge_mass_flows_from_subcycles(self) -> None:
        """按顺时针约定将子循环环量代数汇聚到各边 ``mass_flow``。

        ``mass_flow > 0`` 沿顺时针回路、``< 0`` 沿逆时针、``None`` 按 ``0`` 计；先清空所有边
        再仅对至少出现一次的边写入累加结果。子循环边方向与顺时针段不一致时抛 ``ValueError``。
        汇聚规则详见 ``docs/architecture.md §边流量汇聚``。
        """
        for e in self.edges.values():
            e.mass_flow = None
        totals: dict[str, float] = defaultdict(float)

        for sc in self.subcycles:
            q = 0.0 if sc.mass_flow is None else float(sc.mass_flow)
            n0, n1, n2, n3 = sc.nodes
            # 顺时针四段：(左下→左上)、(左上→右上)、(右上→右下)、(右下→左下)
            segment_edge: tuple[tuple[tuple[int, int], str], ...] = (
                ((n0, n1), sc.edges[0]),
                ((n1, n2), sc.edges[1]),
                ((n2, n3), sc.edges[2]),
                ((n3, n0), sc.edges[3]),
            )
            for (u, v), ek in segment_edge:
                edge = self.edges[ek]
                # 段方向与 edge tail→head 同向则 +q，反向则 -q（子循环环量叠加）
                if edge.tail == u and edge.head == v:
                    totals[ek] += q
                elif edge.tail == v and edge.head == u:
                    totals[ek] -= q
                else:
                    raise ValueError(
                        f"子循环边 {ek!r} 的 tail/head {(edge.tail, edge.head)!r} 与顺时针段 {(u, v)!r} 不一致"
                    )

        for ek, tot in totals.items():
            self.edges[ek].mass_flow = tot

    def commit_subcycle_mass_flows_to_topology(self) -> None:
        """写回 ``SubCycle`` → 汇聚到边 → 重建 ``simplified`` 一气呵成，并清空 ``non_ideal``。

        ``len(subcycle_mass_flows)`` 与 ``len(subcycles)`` 不一致时抛 ``ValueError``。无子循环
        时仅清空各边 ``mass_flow`` 并重建空骨架 ``simplified``。
        """
        self._invalidate_non_ideal()
        n_sc = len(self.subcycles)
        n_q = len(self.subcycle_mass_flows)
        if n_sc != n_q:
            raise ValueError(f"subcycle_mass_flows 长度 {n_q} 与 subcycles 长度 {n_sc} 不一致")
        if n_sc == 0:
            self.assign_edge_mass_flows_from_subcycles()
            self._rebuild_simplified()
            return
        # 写回 SubCycle → 汇聚到 Edge → 按新流量重建 simplified
        self.sync_subcycle_mass_flows_to_subcycles()
        self.assign_edge_mass_flows_from_subcycles()
        self._rebuild_simplified()

    def ensure_non_ideal(self) -> NonIdealClosedCycleLayer:
        """创建或返回 ``non_ideal`` 快照层；不修改理想层字段。

        要求 ``self.simplified`` 已生成。
        """
        from core.non_ideal_bias import NonIdealClosedCycleLayer

        if self.non_ideal is None:
            self.non_ideal = NonIdealClosedCycleLayer.from_closed_cycle_layer(self)
        return self.non_ideal

    def performance_report(self):
        """精简拓扑性能统计；非理想偏置已写入时采用非理想节点状态。

        委托 :mod:`core.cycle_performance`；不缓存结果。
        """
        from core.cycle_performance import (
            compute_cycle_performance,
            resolve_performance_context,
        )

        return compute_cycle_performance(resolve_performance_context(self))

    def analyze_topology(self) -> None:
        """重建完整拓扑并初始化子循环流量。

        ``build_node_edge_topology`` → ``build_subcycles`` → 初值 ``subcycle_mass_flows``
        （``subcycle_mass_flow_initial`` 优先，否则取 ``config.SUBCYCLE_INITIAL_MASS_FLOW_DEFAULT``）
        → 同步与汇聚 → ``_rebuild_simplified()``。覆盖既有拓扑并清空 ``non_ideal``；无子循环时
        ``simplified`` 为空骨架并 ``warn``。
        """
        self._invalidate_non_ideal()
        self.nodes, self.edges, self.skipped_points = build_node_edge_topology(self.properties, self.input)
        self.subcycles = build_subcycles(self.nodes, self.edges)
        n = len(self.subcycles)
        if n == 0:
            self.subcycle_mass_flows = []
            self._rebuild_simplified()
            return
        # 子循环统一初值（优先用输入直给，否则取 config 默认值）
        q0 = self.input.subcycle_mass_flow_initial
        if q0 is None:
            q0 = float(cyges_config.SUBCYCLE_INITIAL_MASS_FLOW_DEFAULT)
        self.subcycle_mass_flows = [q0] * n
        self.sync_subcycle_mass_flows_to_subcycles()
        self.assign_edge_mass_flows_from_subcycles()
        self._rebuild_simplified()
