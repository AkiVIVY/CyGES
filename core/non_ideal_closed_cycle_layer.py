"""闭式循环**非理想层**：理想 ``simplified`` 快照 + 按 ``kind`` 划分的有向连通组与组内深度。

由 ``ClosedCycleLayer.ensure_non_ideal()`` 构造，**不修改**父层 ``nodes`` / ``edges``。父层
``analyze_topology()`` / ``commit_*`` 时 ``layer.non_ideal`` 被置 ``None``，已生成的实例仍指向
ensure 时刻的快照。

``ideal_nodes`` 为 ensure 时刻对父层 ``nodes`` 的只读引用；``nodes`` / 偏移系数字段预留给
后续新偏移实现。旧版换热/机械节点偏置已迁至
``core.non_ideal_node_offsets_legacy``（待重写，勿作长期 API）。

``reach`` / ``layer`` 定义见 ``docs/architecture.md``。
"""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass
from typing import Literal, Mapping

from core.closed_cycle_layer import Node, SimplifiedEdge, SimplifiedTopology
from core.fluid_property_solver import FluidPropertySolver


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
# §2. 组内深度（reach 下游最长路径 / layer 主脊分层）
# ----------------------------------------------------------------
# 定义与示例详见 docs/architecture.md §有向组与深度。
# ============================================================


@dataclass(frozen=True)
class GroupDepthMetrics:
    """组内有向深度的中间结果。

    - ``reach``：从该点沿 ``tail→head`` 的最长下游路径边数（用于 ``upstream_special_nodes``）。
    - ``layer``：主脊分层刻度，见 ``docs/architecture.md §7.3``（新偏移实现可复用）。
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

    spine: list[int] = [best_start]
    cur = best_start
    while dist[cur] < max_d:
        nxt = min(w for w in adj[cur] if dist[w] == dist[cur] + 1)
        spine.append(nxt)
        cur = nxt

    layer: dict[int, int] = {v: _LAYER_UNSET for v in nodes}
    spine_set = frozenset(spine)
    for i, v in enumerate(spine):
        layer[v] = i

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
# §3. 有向组数据类与构造
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
# §4. 非理想分析容器（节点偏置实现见 non_ideal_node_offsets_legacy）
# ============================================================


@dataclass
class NonIdealClosedCycleLayer:
    """非理想分析容器（分组/深度快照；节点偏置待新实现）。

    字段
    ----
    - ``simplified``：ensure 时刻 ``layer.simplified`` 的同一引用（frozen，不可变）。
    - ``mechanical_groups`` / ``heat_groups``：按 ``kind`` 划分的 ``SimplifiedDirectedGroup`` 元组。
    - ``ideal_nodes``：ensure 时刻 ``layer.nodes`` 的只读映射。
    - ``nodes`` / ``heat_efficiency`` / ``mechanical_efficiency``：预留给新偏移管线。
    - ``properties``：与父层相同的物性求解器。

    由 ``from_closed_cycle_layer`` / ``ClosedCycleLayer.ensure_non_ideal()`` 创建。
    临时沿用旧偏移请 import ``core.non_ideal_node_offsets_legacy``。
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
    def from_closed_cycle_layer(cls, layer) -> NonIdealClosedCycleLayer:
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
