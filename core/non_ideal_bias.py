"""闭式循环**非理想偏置**（``core.non_ideal_bias``）：理想 ``simplified`` 快照、有向组、组内深度，以及**单步**节点偏置。

由 ``ClosedCycleLayer.ensure_non_ideal()`` 构造，**不修改**父层 ``nodes`` / ``edges``。父层
``analyze_topology()`` / ``commit_*`` 时 ``layer.non_ideal`` 被置 ``None``，已生成的实例仍指向
ensure 时刻的快照。

``ideal_nodes`` 为 ensure 时刻对父层 ``nodes`` 的只读引用；首次调用 ``apply_offsets`` /
:func:`_apply_combined_offsets` 时整表拷贝到 ``nodes`` 再写修正量（理想层不变）。

模块分节：§1 分组 → §2 有向组数据类 → §3 深度 → §4 有向组构造 → §5 快照容器
``NonIdealClosedCycleLayer`` → §6 节点偏置 ``apply_offsets``。

``reach`` / ``layer`` / 换热 σ / 机械 η_is 公式见 ``docs/architecture.md``。
"""

from __future__ import annotations

import warnings
from collections import defaultdict, deque
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Literal, Mapping

import config as cyges_config

from core.closed_cycle_layer import Node, SimplifiedEdge, SimplifiedTopology
from core.fluid_property_solver import FluidPropertySolver

if TYPE_CHECKING:
    from core.closed_cycle_layer import ClosedCycleLayer


# ============================================================
# §1. 分组与有向邻接（按 kind 拆分 + 并查集求无向连通分量）
# ============================================================


def partition_simplified_edges_by_kind(
    topology: SimplifiedTopology,
) -> tuple[tuple[frozenset[str], ...], tuple[frozenset[str], ...]]:
    """按 ``kind`` 拆分精简边，再各自求**无向**连通分量。

    返回 ``(mechanical_edge_groups, heat_edge_groups)``，元素为边键 ``frozenset``；组顺序按
    组内最小边键字符串排序稳定化。不改变边的有向语义。
    """
    mech: list[tuple[str, int, int]] = []
    heat: list[tuple[str, int, int]] = []
    for ek, se in topology.simplified_edges:
        if se.kind == "mechanical":
            mech.append((ek, se.tail, se.head))
        else:
            heat.append((ek, se.tail, se.head))

    def _components(items: list[tuple[str, int, int]]) -> tuple[frozenset[str], ...]:
        """并查集：同一 kind 下 tail/head 无向连通的精简边归为同一组。"""
        if not items:
            return ()
        parent: dict[int, int] = {}

        def find(x: int) -> int:
            if x not in parent:
                parent[x] = x
                return x
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]

        def union(a: int, b: int) -> None:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra

        for _, tail, head in items:
            union(tail, head)

        root_to_keys: dict[int, list[str]] = defaultdict(list)
        for ek, tail, head in items:
            root_to_keys[find(tail)].append(ek)

        groups = [frozenset(sorted(keys)) for keys in root_to_keys.values()]
        groups.sort(key=lambda g: min(g) if g else "")
        return tuple(groups)

    return _components(mech), _components(heat)

# ============================================================
# §2. 有向组数据类 `SimplifiedDirectedGroup`
# ============================================================


@dataclass(frozen=True)
class SimplifiedDirectedGroup:
    """同 ``kind`` 精简边的一个**无向连通分量**与组内有向层号、最上游特殊节点。

    字段
    ----
    - ``kind`` / ``edge_keys``：组类型与本组 ``SM*`` / ``SH*`` 边键集合。
    - ``node_depth``：``(node_index, layer)`` 元组，按 index 升序；最上游 ``0``，最下游
      ``max_depth``。
    - ``upstream_special_nodes``：``reach`` 最大的节点（``frozenset``，允许并列）。

    另一 ``kind`` 的组对同一 ``Node.index`` 可能有不同层号；查询请用 ``depth_dict()``。
    """

    kind: Literal["mechanical", "heat"]
    edge_keys: frozenset[str]
    node_depth: tuple[tuple[int, int], ...]
    upstream_special_nodes: frozenset[int]

    def depth_dict(self) -> dict[int, int]:
        """本组内 ``node_index → 层号``（勿与另一 kind 的组混用）。"""
        return dict(self.node_depth)

    def nodes(self) -> frozenset[int]:
        """组内作为精简边端点出现的所有 ``Node.index``。"""
        return frozenset(i for i, _ in self.node_depth)

    @property
    def max_depth(self) -> int:
        """组内层号的最大值（最下游）；``node_depth`` 为空时为 ``0``。"""
        if not self.node_depth:
            return 0
        return max(d for _, d in self.node_depth)
# ============================================================
# §3. 组内深度（reach 下游最长路径 / layer 主脊分层）
# ----------------------------------------------------------------
# 定义与示例详见 docs/architecture.md §有向组与深度。
# ============================================================


@dataclass(frozen=True)
class GroupDepthMetrics:
    """组内有向深度的中间结果。

    - ``reach``：从该点沿 ``tail→head`` 的最长下游路径边数（用于 ``upstream_special_nodes``）。
    - ``layer``：主脊分层刻度，见 ``docs/architecture.md §7.3``（换热 ``σ^layer`` 等偏置使用）。
    """

    reach: dict[int, int]
    layer: dict[int, int]


def _group_adjacency(
    edges_by_key: dict[str, SimplifiedEdge],
    edge_keys: frozenset[str],
) -> tuple[dict[int, list[int]], set[int], list[tuple[int, int]]]:
    """从一组精简边键构建有向邻接表 ``tail→head``、节点集与有向边列表。"""
    adj: dict[int, list[int]] = defaultdict(list)
    nodes: set[int] = set()
    directed_edges: list[tuple[int, int]] = []
    for ek in edge_keys:
        se = edges_by_key[ek]
        adj[se.tail].append(se.head)
        nodes.add(se.tail)
        nodes.add(se.head)
        directed_edges.append((se.tail, se.head))
    return adj, nodes, directed_edges


def _reach_from_adj(
    adj: dict[int, list[int]],
    nodes: set[int],
    edge_keys: frozenset[str],
) -> dict[int, int]:
    """在已有有向邻接表上计算 ``reach``（供 ``_compute_group_depth_metrics`` 复用邻接）。"""
    memo: dict[int, int] = {}
    visiting: set[int] = set()

    def reach_dfs(v: int) -> int:
        if v in memo:
            return memo[v]
        if v in visiting:
            raise ValueError(
                f"精简边组 {sorted(edge_keys)!r} 内存在有向环（涉及节点 {v}），无法定义下游深度"
            )
        visiting.add(v)
        best = 0
        for u in adj.get(v, ()):
            best = max(best, 1 + reach_dfs(u))
        visiting.remove(v)
        memo[v] = best
        return best

    return {v: reach_dfs(v) for v in nodes}


def compute_group_downstream_reach(
    edges_by_key: dict[str, SimplifiedEdge],
    edge_keys: frozenset[str],
) -> dict[int, int]:
    """计算组内各端点的**下游最长路径边数** ``reach``。

    ``reach(v) = 0`` 表示 ``v`` 无出边。组内存在有向环时抛 ``ValueError``。
    """
    if not edge_keys:
        return {}
    adj, nodes, _ = _group_adjacency(edges_by_key, edge_keys)
    return _reach_from_adj(adj, nodes, edge_keys)


_LAYER_UNSET = -1


def _compute_layer_by_spine(
    adj: dict[int, list[int]],
    nodes: set[int],
    directed_edges: list[tuple[int, int]],
) -> dict[int, int]:
    """主脊分层：最长路径定标（并列时最小起点）→ 脊上 ``0…L`` → 支流/汇流对齐。

    详见 ``docs/architecture.md §7.3``。
    """
    if not nodes:
        return {}

    # ── 阶段 1：Kahn 拓扑排序 → 各节点自上而下的最长距离 dist ──
    in_degree: dict[int, int] = {v: 0 for v in nodes}
    for u, w in directed_edges:
        in_degree[w] += 1

    dist: dict[int, int] = {v: 0 for v in nodes}
    q: deque[int] = deque(sorted(v for v, d in in_degree.items() if d == 0))
    remaining = dict(in_degree)
    while q:
        u = q.popleft()
        for w in adj.get(u, ()):
            if dist[u] + 1 > dist[w]:
                dist[w] = dist[u] + 1
            remaining[w] -= 1
            if remaining[w] == 0:
                q.append(w)

    max_d = max(dist.values())

    # ── 阶段 2：反向邻接 + 回溯 → 找距 max_d 最远的最小起点 ──
    rev: dict[int, list[int]] = defaultdict(list)
    for u, w in directed_edges:
        if dist[u] + 1 == dist[w]:
            rev[w].append(u)

    def starts_at(v: int, d: int) -> list[int]:
        if d == 0:
            return [v]
        out: list[int] = []
        for u in rev.get(v, ()):
            if dist[u] == d - 1:
                out.extend(starts_at(u, d - 1))
        return out

    best_start: int | None = None
    for end in nodes:
        if dist[end] != max_d:
            continue
        for s in starts_at(end, max_d):
            if best_start is None or s < best_start:
                best_start = s
    if best_start is None:
        best_start = min(nodes)

    # 「能延续到 dist == max_d」的节点集合：仅这些节点是合法的脊路径中继点
    on_longest: set[int] = set()
    rev_stack: list[int] = [v for v in nodes if dist[v] == max_d]
    while rev_stack:
        v = rev_stack.pop()
        if v in on_longest:
            continue
        on_longest.add(v)
        for u in rev.get(v, ()):
            if u not in on_longest:
                rev_stack.append(u)

    # ── 阶段 3：从 best_start 沿 dist 递增走主脊（后继必须 on_longest） ──
    spine: list[int] = [best_start]
    cur = best_start
    while dist[cur] < max_d:
        # 必须挑既递增一层、又能继续走到 max_d 的后继（避免选到岔出的支流）
        nxt = min(
            w
            for w in adj[cur]
            if dist[w] == dist[cur] + 1 and w in on_longest
        )
        spine.append(nxt)
        cur = nxt

    layer: dict[int, int] = {v: _LAYER_UNSET for v in nodes}
    spine_set = frozenset(spine)
    for i, v in enumerate(spine):
        layer[v] = i

    # ── 阶段 4：迭代正反向松弛 → 支流/汇流对齐，主脊节点不受影响 ──
    for _ in range(len(nodes) + 1):
        changed = False
        for u, w in directed_edges:
            if layer[u] >= 0:
                nw = layer[u] + 1
                if layer[w] < 0:
                    layer[w] = nw
                    changed = True
                elif w not in spine_set and layer[w] < nw:
                    layer[w] = nw
                    changed = True
            if layer[w] >= 0 and u not in spine_set:
                nu = layer[w] - 1
                if layer[u] < 0:
                    layer[u] = nu
                    changed = True
                elif layer[u] > nu:
                    layer[u] = nu
                    changed = True
        if not changed:
            break

    # 收敛后仍悬空的孤立/不连通节点 → 归零
    for v in nodes:
        if layer[v] < 0:
            layer[v] = 0
    return layer


def _compute_group_depth_metrics(
    edges_by_key: dict[str, SimplifiedEdge],
    edge_keys: frozenset[str],
) -> GroupDepthMetrics:
    """同时计算组内 ``reach`` 与 ``layer``（单次构建有向邻接）。"""
    if not edge_keys:
        return GroupDepthMetrics(reach={}, layer={})

    adj, nodes, directed_edges = _group_adjacency(edges_by_key, edge_keys)
    reach = _reach_from_adj(adj, nodes, edge_keys)
    if not reach:
        return GroupDepthMetrics(reach={}, layer={})

    layer = _compute_layer_by_spine(adj, nodes, directed_edges)
    return GroupDepthMetrics(reach=reach, layer=layer)


def compute_group_downstream_depth(
    edges_by_key: dict[str, SimplifiedEdge],
    edge_keys: frozenset[str],
) -> dict[int, int]:
    """计算组内各端点 ``layer``（``node_depth`` 公开语义）。

    定义与示例详见 ``docs/architecture.md §有向组与深度``。
    """
    return _compute_group_depth_metrics(edges_by_key, edge_keys).layer


# ============================================================
# §4. 有向组构造（build_directed_groups*）
# ============================================================

def _directed_groups_from_edge_keys(
    kind: Literal["mechanical", "heat"],
    edge_key_groups: tuple[frozenset[str], ...],
    edges_by_key: dict[str, SimplifiedEdge],
) -> tuple[SimplifiedDirectedGroup, ...]:
    """由已分好的边键组构建 ``SimplifiedDirectedGroup``（避免重复 partition）。"""
    out: list[SimplifiedDirectedGroup] = []
    for keys in edge_key_groups:
        metrics = _compute_group_depth_metrics(edges_by_key, keys)
        if metrics.reach:
            r_max = max(metrics.reach.values())
            special = frozenset(v for v, r in metrics.reach.items() if r == r_max)
            node_depth = tuple(sorted(metrics.layer.items()))
        else:
            special = frozenset()
            node_depth = ()
        out.append(
            SimplifiedDirectedGroup(
                kind=kind,
                edge_keys=keys,
                node_depth=node_depth,
                upstream_special_nodes=special,
            )
        )
    return tuple(out)


def build_directed_groups(
    topology: SimplifiedTopology,
    kind: Literal["mechanical", "heat"],
) -> tuple[SimplifiedDirectedGroup, ...]:
    """按 ``kind`` 切分 → 逐组算层号与 ``upstream_special_nodes``。

    若同时需要机械与换热两组，请用 ``build_directed_groups_both`` 或自行只 partition 一次。
    """
    mech_keys, heat_keys = partition_simplified_edges_by_kind(topology)
    edge_key_groups = mech_keys if kind == "mechanical" else heat_keys
    return _directed_groups_from_edge_keys(kind, edge_key_groups, topology.edges_dict())


def build_directed_groups_both(
    topology: SimplifiedTopology,
) -> tuple[tuple[SimplifiedDirectedGroup, ...], tuple[SimplifiedDirectedGroup, ...]]:
    """一次 partition 后同时构建机械、换热有向组（``ensure_non_ideal`` 使用）。"""
    mech_keys, heat_keys = partition_simplified_edges_by_kind(topology)
    edges_by_key = topology.edges_dict()
    return (
        _directed_groups_from_edge_keys("mechanical", mech_keys, edges_by_key),
        _directed_groups_from_edge_keys("heat", heat_keys, edges_by_key),
    )
# ============================================================
# §5. 非理想分析容器
# ----------------------------------------------------------------
# 仅持有快照与分组；状态偏移由 §6 ``apply_offsets()`` 完成。
# ============================================================


@dataclass
class NonIdealClosedCycleLayer:
    """非理想分析容器（可变 dataclass，便于挂载偏移结果）。

    字段
    ----
    - ``simplified``：ensure 时刻 ``layer.simplified`` 的同一引用（frozen，不可变）。
    - ``mechanical_groups`` / ``heat_groups``：按 ``kind`` 划分的 ``SimplifiedDirectedGroup`` 元组。
    - ``ideal_nodes``：ensure 时刻 ``layer.nodes`` 的只读映射。
    - ``nodes``：首次 ``apply_offsets()`` 时从 ``ideal_nodes`` 拷贝；写非理想状态。
    - ``heat_efficiency`` / ``mechanical_efficiency``：最近一次偏置采用的 ``σ`` / ``η_is``；
      未偏置时为 ``None``。
    - ``properties``：与父层相同的物性求解器。

     由 ``from_closed_cycle_layer`` / ``ClosedCycleLayer.ensure_non_ideal()`` 创建；
      偏置请调用 ``apply_offsets()``。
    """

    simplified: SimplifiedTopology
    mechanical_groups: tuple[SimplifiedDirectedGroup, ...]
    heat_groups: tuple[SimplifiedDirectedGroup, ...]
    ideal_nodes: Mapping[int, Node]
    properties: FluidPropertySolver
    nodes: dict[int, Node] | None = None
    heat_efficiency: float | None = None
    mechanical_efficiency: float | None = None

    @classmethod
    def from_closed_cycle_layer(cls, layer: ClosedCycleLayer) -> NonIdealClosedCycleLayer:
        """读取理想层 ``simplified`` / ``nodes``，构建机械/换热有向组快照。

        前置：``layer.simplified is not None``（已 ``analyze_topology`` 或 ``commit_*``）。
        """
        if layer.simplified is None:
            raise RuntimeError(
                "ClosedCycleLayer.simplified 尚未生成；须先 analyze_topology() 或 commit_subcycle_mass_flows_to_topology()。"
            )
        simp = layer.simplified
        mech_groups, heat_groups = build_directed_groups_both(simp)
        return cls(
            simplified=simp,
            mechanical_groups=mech_groups,
            heat_groups=heat_groups,
            ideal_nodes=layer.nodes,
            properties=layer.properties,
        )

    def apply_offsets(
        self,
        sigma: float | None = None,
        eta_is: float | None = None,
    ) -> None:
        """单步应用换热 ``σ`` 与机械 ``η_is`` 偏置，写入 ``self.nodes``。

        详见 ``docs/architecture.md`` §8。
        """
        _apply_combined_offsets(self, sigma=sigma, eta_is=eta_is)

# ============================================================
# §6. 节点偏置（单步：换热 σ → 机械组 PS 重置 → DFS ``PS→η→HP``）
# ----------------------------------------------------------------
# 公式与 anchor 决策见 ``docs/architecture.md`` §8。
# ============================================================


def _pressure_equal(P0: float, P1: float) -> bool:
    """两个压力是否在相对/绝对容差下视为相等。"""
    tol = max(1e-9, 1e-6 * max(abs(P0), abs(P1)))
    return abs(P1 - P0) <= tol


def _resolve_efficiency(
    param: float | None,
    stored: float | None,
    default: float,
    label: str,
) -> float:
    """参数 → 实例字段 → config 默认；校验 ``(0, 1]``。"""
    eta = param if param is not None else (stored if stored is not None else default)
    if not (0.0 < eta <= 1.0):
        raise ValueError(f"{label} 须在 (0, 1] 内，收到 {eta!r}")
    return eta


def _ensure_nodes(layer: NonIdealClosedCycleLayer) -> dict[int, Node]:
    """首次调用时由 ``ideal_nodes`` 整表拷贝出可写 ``nodes`` 表。"""
    if layer.nodes is None:
        layer.nodes = {i: replace(n) for i, n in layer.ideal_nodes.items()}
    return layer.nodes


def _apply_heat_pressure(layer: NonIdealClosedCycleLayer, sigma: float) -> None:
    """**步骤 1**：换热组内 ``P_new = P_ideal × σ ** layer``；仅写 P。

    始终以 ``ideal_nodes[v].P`` 为底，避免重复偏移；不调用 PS 闪蒸，``T/H/S`` 保留理想值，
    供步骤 2 取一级节点理想 ``T``。
    """
    nodes = _ensure_nodes(layer)
    for group in layer.heat_groups:
        layers = group.depth_dict()
        for v in group.nodes():
            p_ideal = layer.ideal_nodes[v].P
            p_new = p_ideal * (sigma ** layers[v])
            nodes[v] = replace(nodes[v], P=p_new)


def _pick_anchor_and_sbase(
    layer: NonIdealClosedCycleLayer,
    group: SimplifiedDirectedGroup,
) -> tuple[int, float]:
    """含一级 → 步骤 2；否则 → 步骤 3。返回 ``(anchor, S_base)``。"""
    nodes = _ensure_nodes(layer)
    group_nodes = group.nodes()
    primaries = sorted(v for v in group_nodes if layer.ideal_nodes[v].parent is None)
    if primaries:
        anchor = primaries[0]
        s_base = layer.properties.state(
            "TP", layer.ideal_nodes[anchor].T, nodes[anchor].P
        )["S"]
        return anchor, s_base

    depths = group.depth_dict()
    zeros = sorted(v for v in group_nodes if depths.get(v, -1) == 0)
    if zeros:
        anchor = zeros[0]
        return anchor, nodes[anchor].S

    raise ValueError(
        f"机械有向组 {sorted(group.edge_keys)!r} 既无一级节点也无 layer==0 节点，无法确定 anchor"
    )


def _reset_group_to_base_entropy(
    layer: NonIdealClosedCycleLayer,
    group: SimplifiedDirectedGroup,
    s_base: float,
) -> None:
    """**步骤 2/3 共用**：组内每个节点 ``state("PS", P_v, S_base)`` 重置 ``T,H,S``。"""
    nodes = _ensure_nodes(layer)
    for v in group.nodes():
        st = layer.properties.state("PS", nodes[v].P, s_base)
        nodes[v] = replace(nodes[v], T=st["T"], P=st["P"], H=st["H"], S=st["S"])


def _mechanical_step_known_to_unknown(
    *,
    nodes: dict[int, Node],
    properties: FluidPropertySolver,
    known: int,
    unknown: int,
    edge: SimplifiedEdge,
    eta_is: float,
) -> None:
    """沿一条精简机械边由 ``known`` 推 ``unknown``：``PS→η→HP``（见 architecture §8）。"""
    n_known = nodes[known]
    n_unknown = nodes[unknown]
    p_unknown = n_unknown.P
    h_known = n_known.H
    s_known = n_known.S

    if edge.tail == known and edge.head == unknown:
        p_tail, p_head = n_known.P, n_unknown.P
    elif edge.head == known and edge.tail == unknown:
        p_tail, p_head = n_unknown.P, n_known.P
    else:
        raise ValueError(
            f"机械边 {edge!r} 的端点 (tail={edge.tail}, head={edge.head}) 与"
            f"已知/待求 ({known}, {unknown}) 不匹配"
        )

    if _pressure_equal(p_tail, p_head):
        h_new = h_known
    else:
        h1 = properties.state("PS", p_unknown, s_known)["H"]
        if p_head > p_tail:
            h_new = (h1 - h_known) / eta_is + h_known
        else:
            h_new = (h1 - h_known) * eta_is + h_known

    st = properties.state("HP", h_new, p_unknown)
    nodes[unknown] = replace(
        n_unknown, T=st["T"], P=st["P"], H=st["H"], S=st["S"]
    )


def _walk_mechanical_branches(
    *,
    layer: NonIdealClosedCycleLayer,
    group: SimplifiedDirectedGroup,
    edges_by_key: dict[str, SimplifiedEdge],
    anchor: int,
    eta_is: float,
) -> set[int]:
    """从 anchor 沿组内无向邻接 DFS，每条精简边由已知端推未知端。"""
    nodes = _ensure_nodes(layer)
    adj: dict[int, list[tuple[int, str]]] = defaultdict(list)
    for ek in sorted(group.edge_keys):
        e = edges_by_key[ek]
        adj[e.tail].append((e.head, ek))
        adj[e.head].append((e.tail, ek))

    known: set[int] = {anchor}
    for start_nb, _ in sorted(adj.get(anchor, [])):
        if start_nb in known:
            continue
        stack: list[tuple[int, int]] = [(anchor, start_nb)]
        while stack:
            prev, cur = stack.pop()
            if cur in known:
                continue
            edge_key = next(ek for nb, ek in adj[prev] if nb == cur)
            edge = edges_by_key[edge_key]
            _mechanical_step_known_to_unknown(
                nodes=nodes,
                properties=layer.properties,
                known=prev,
                unknown=cur,
                edge=edge,
                eta_is=eta_is,
            )
            known.add(cur)
            for nb, _ in sorted(adj[cur]):
                if nb not in known:
                    stack.append((cur, nb))
    return known


def _apply_mechanical_group(
    layer: NonIdealClosedCycleLayer,
    group: SimplifiedDirectedGroup,
    edges_by_key: dict[str, SimplifiedEdge],
    eta_is: float,
) -> None:
    """单机械组：选 anchor + ``S_base`` → 组内 PS 重置 → 从 anchor DFS 非理想机械步。"""
    if not group.edge_keys:
        return
    anchor, s_base = _pick_anchor_and_sbase(layer, group)
    _reset_group_to_base_entropy(layer, group, s_base)
    known = _walk_mechanical_branches(
        layer=layer,
        group=group,
        edges_by_key=edges_by_key,
        anchor=anchor,
        eta_is=eta_is,
    )
    missing = group.nodes() - known
    if missing:
        warnings.warn(
            f"机械有向组 {sorted(group.edge_keys)!r} 中节点 {sorted(missing)} "
            f"与 anchor {anchor} 在无向意义上不连通，已跳过非理想机械步",
            RuntimeWarning,
            stacklevel=3,
        )


def _apply_combined_offsets(
    layer: NonIdealClosedCycleLayer,
    *,
    sigma: float | None = None,
    eta_is: float | None = None,
) -> NonIdealClosedCycleLayer:
    """单步骤非理想节点偏置：换热 σ 改 P → 机械组 PS 重置 → 从 anchor DFS 非理想机械步。

    ``sigma``：换热总压恢复系数；``None`` → ``layer.heat_efficiency`` →
    ``config.NON_IDEAL_HEAT_EFFICIENCY_DEFAULT``。须在 ``(0, 1]``。
    ``eta_is``：机械等熵效率；``None`` → ``layer.mechanical_efficiency`` →
    ``config.NON_IDEAL_MECHANICAL_EFFICIENCY_DEFAULT``。须在 ``(0, 1]``。

    返回 ``layer`` 自身（``nodes`` 已写入新状态；``heat_efficiency`` /
    ``mechanical_efficiency`` 记录本次 ``σ`` / ``η_is``）。不修改 ``ideal_nodes`` 与父层。
    """
    sigma_v = _resolve_efficiency(
        sigma,
        layer.heat_efficiency,
        float(cyges_config.NON_IDEAL_HEAT_EFFICIENCY_DEFAULT),
        "换热总压恢复系数 σ",
    )
    eta_v = _resolve_efficiency(
        eta_is,
        layer.mechanical_efficiency,
        float(cyges_config.NON_IDEAL_MECHANICAL_EFFICIENCY_DEFAULT),
        "机械等熵效率 η_is",
    )

    _ensure_nodes(layer)
    _apply_heat_pressure(layer, sigma_v)

    edges_by_key = layer.simplified.edges_dict()
    for group in layer.mechanical_groups:
        _apply_mechanical_group(layer, group, edges_by_key, eta_v)

    layer.heat_efficiency = sigma_v
    layer.mechanical_efficiency = eta_v
    return layer
