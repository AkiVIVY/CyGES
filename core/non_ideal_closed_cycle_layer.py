"""闭式循环**非理想层**：理想 ``simplified`` 快照 + 按 ``kind`` 划分的有向连通组与组内深度。

由 ``ClosedCycleLayer.ensure_non_ideal()`` 构造，**不修改**父层 ``nodes`` / ``edges``。父层
``analyze_topology()`` / ``commit_*`` 时 ``layer.non_ideal`` 被置 ``None``，已生成的实例仍指向
ensure 时刻的快照。

节点状态：``ideal_nodes`` 为 ensure 时刻对父层 ``nodes`` 的只读引用；``nodes`` 在首次偏移
方法调用时整表拷贝，再写修正量，理想层不变。同一 ``Node.index`` 可同时属于一个机械组与一个
换热组，深度按组查询（``group.depth_dict()``）。

算法细节（``reach`` / ``layer`` 定义、换热 σ 偏移、机械锚点决策、机械步公式）见
``docs/architecture.md``。
"""

from __future__ import annotations

import warnings
from collections import defaultdict, deque
from dataclasses import dataclass, replace
from typing import Literal, Mapping

import config as cyges_config

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
# §2. 组内深度（reach 下游最长路径 / layer 拓扑序 DP 最长上游路径）
# ----------------------------------------------------------------
# 定义与示例详见 docs/architecture.md §有向组与深度。
# ============================================================


@dataclass(frozen=True)
class GroupDepthMetrics:
    """组内有向深度的中间结果。

    - ``reach``：从该点沿 ``tail→head`` 的最长下游路径边数（用于 ``upstream_special_nodes``）。
    - ``layer``：从任一入度 0 源点到该点的最长上游路径边数（用于换热 ``σ^layer`` 叠乘）。
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
        # 记忆化 DFS：沿 tail→head 取最长下游路径；回边检测用于发现组内有向环
        best = 0
        for u in adj.get(v, ()):
            best = max(best, 1 + reach_dfs(u))
        visiting.remove(v)
        memo[v] = best
        return best

    return {v: reach_dfs(v) for v in nodes}


def _compute_group_depth_metrics(
    edges_by_key: dict[str, SimplifiedEdge],
    edge_keys: frozenset[str],
) -> GroupDepthMetrics:
    """同时计算组内 ``reach`` 与 ``layer``。

    ``layer`` 用拓扑序 DP（Kahn）：源点 ``layer=0``，每条 ``u→v`` 满足
    ``layer[v] = max(layer[v], layer[u] + 1)``。``reach`` 的计算已保证组内无有向环。
    """
    reach = compute_group_downstream_reach(edges_by_key, edge_keys)
    if not reach:
        return GroupDepthMetrics(reach={}, layer={})

    adj, _, directed_edges = _group_adjacency(edges_by_key, edge_keys)
    # 入度 0 的节点 = 组内源点（无上游精简边指向），层号从 0 起算
    in_degree: dict[int, int] = {v: 0 for v in reach}
    for _, w in directed_edges:
        in_degree[w] += 1

    layer: dict[int, int] = {v: 0 for v in reach}
    q: deque[int] = deque(sorted(v for v, d in in_degree.items() if d == 0))
    remaining = dict(in_degree)
    while q:
        u = q.popleft()
        for w in adj.get(u, ()):
            # 多源汇入时取最长上游路径：layer[w] = max(经各 u 的 layer[u]+1)
            if layer[u] + 1 > layer[w]:
                layer[w] = layer[u] + 1
            remaining[w] -= 1
            if remaining[w] == 0:
                q.append(w)

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


def build_directed_groups(
    topology: SimplifiedTopology,
    kind: Literal["mechanical", "heat"],
) -> tuple[SimplifiedDirectedGroup, ...]:
    """按 ``kind`` 切分 → 逐组算层号与 ``upstream_special_nodes``，输出完整 ``SimplifiedDirectedGroup``。

    ``partition_simplified_edges_by_kind`` 只分边键；本函数另带层号与特殊节点，供
    ``NonIdealClosedCycleLayer`` 持有。
    """
    mech_keys, heat_keys = partition_simplified_edges_by_kind(topology)
    edge_key_groups = mech_keys if kind == "mechanical" else heat_keys
    edges_by_key = topology.edges_dict()

    out: list[SimplifiedDirectedGroup] = []
    for keys in edge_key_groups:
        metrics = _compute_group_depth_metrics(edges_by_key, keys)
        if metrics.reach:
            # upstream_special：下游延伸 reach 最大者（机械锚点、绘图高亮用）
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


# ============================================================
# §4. 偏移实现（机械锚点 / 机械步公式 / 支路 DFS / 换热 σ + PS 闭合）
# ----------------------------------------------------------------
# 详见 docs/architecture.md §换热偏移、§机械偏移。
# ============================================================


def _pressure_equal(P0: float, P1: float) -> bool:
    """两个压力是否在相对/绝对容差下视为相等。"""
    tol = max(1e-9, 1e-6 * max(abs(P0), abs(P1)))
    return abs(P1 - P0) <= tol


def _pick_mechanical_anchor(
    group: SimplifiedDirectedGroup,
    ideal_nodes: Mapping[int, Node],
) -> int:
    """机械组基准节点选择（锚点整点状态不再修改）。

    决策树：一级 ∈ ``upstream_special_nodes`` → 该一级（index 最小）→
    ``min(upstream_special_nodes)`` → 任一一级 → ``min(layer==0)``。详见
    ``docs/architecture.md §机械锚点决策``。
    """
    group_nodes = group.nodes()
    primaries = sorted(v for v in group_nodes if ideal_nodes[v].parent is None)
    specials = group.upstream_special_nodes
    # 一级且在 reach 最大集合内 → 熵不动点与网格原点一致
    if primaries and (not specials or primaries[0] in specials):
        return primaries[0]
    # 分叉组：一级在侧枝末端时，改用下游延伸最长的 special 作锚
    if specials:
        return min(specials)
    if primaries:
        return primaries[0]
    layers = group.depth_dict()
    layer0 = sorted(v for v in group_nodes if layers[v] == 0)
    if layer0:
        return layer0[0]
    raise ValueError(f"机械有向组 {sorted(group.edge_keys)!r} 无法确定基准节点")


def _refresh_node_from_ps(
    nodes: dict[int, Node],
    properties: FluidPropertySolver,
    index: int,
) -> None:
    """以当前 ``(P, S)`` 通过 ``PS`` 闪蒸刷新节点 ``T,H``（``P,S`` 以求解器返回值为准）。"""
    n = nodes[index]
    st = properties.state("PS", n.P, n.S)
    nodes[index] = replace(n, T=st["T"], P=st["P"], H=st["H"], S=st["S"])


def _mechanical_step_known_to_unknown(
    *,
    nodes: dict[int, Node],
    properties: FluidPropertySolver,
    known: int,
    unknown: int,
    edge: SimplifiedEdge,
    eta_is: float,
) -> None:
    """沿精简机械边由 ``known`` 推 ``unknown`` 一步：``PS→η→HP``。

    ``H1 = state("PS", P_unknown, S_known)["H"]``；按边方向（压缩 ``/η_is``、膨胀 ``×η_is``、
    等压保持）算 ``H2``；``state("HP", H2, P_unknown)`` 写回 ``T,P,H,S``。仅修改
    ``nodes[unknown]``。公式细节见 ``docs/architecture.md §机械步公式``。
    """
    n_known = nodes[known]
    n_unknown = nodes[unknown]
    p_unknown = n_unknown.P
    h_known = n_known.H
    s_known = n_known.S

    # 统一得到边的 P_tail、P_head（与精简边 tail/head 几何一致），用于判压缩/膨胀
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
        # 等熵假想：待求端压力下沿用已知端熵（H1）；再按 η_is 偏离到实际焓 H2
        h1 = properties.state("PS", p_unknown, s_known)["H"]
        if p_head > p_tail:
            h_new = (h1 - h_known) / eta_is + h_known
        else:
            h_new = (h1 - h_known) * eta_is + h_known

    # HP 闭合：非理想态 S 一般 ≠ S_known
    st = properties.state("HP", h_new, p_unknown)
    nodes[unknown] = replace(
        n_unknown, T=st["T"], P=st["P"], H=st["H"], S=st["S"]
    )


def _walk_mechanical_branches(
    *,
    nodes: dict[int, Node],
    properties: FluidPropertySolver,
    group: SimplifiedDirectedGroup,
    edges_by_key: dict[str, SimplifiedEdge],
    anchor: int,
    eta_is: float,
) -> set[int]:
    """从锚点沿无向邻接 DFS，逐条支路单向推进（每条边只从已知端推到未知端）。

    返回处理完毕的节点集合（含锚点）。与锚点无向不连通的节点不会被更新，由调用方决定如何告警。
    """
    # 无向邻接：(邻居节点, 连接边键)，便于从锚点向各支路 DFS
    adj: dict[int, list[tuple[int, str]]] = defaultdict(list)
    for ek in sorted(group.edge_keys):
        e = edges_by_key[ek]
        adj[e.tail].append((e.head, ek))
        adj[e.head].append((e.tail, ek))

    known: set[int] = {anchor}
    for start_nb, _ in sorted(adj.get(anchor, [])):
        if start_nb in known:
            continue
        # 每条从锚点出发的支路：栈中 (上一已知点, 待更新点)，只沿远离锚点方向推进
        stack: list[tuple[int, int]] = [(anchor, start_nb)]
        while stack:
            prev, cur = stack.pop()
            if cur in known:
                continue
            edge_key = next(ek for nb, ek in adj[prev] if nb == cur)
            edge = edges_by_key[edge_key]
            _mechanical_step_known_to_unknown(
                nodes=nodes,
                properties=properties,
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


@dataclass
class NonIdealClosedCycleLayer:
    """非理想分析容器（可变 dataclass，便于挂载偏移结果）。

    字段
    ----
    - ``simplified``：ensure 时刻 ``layer.simplified`` 的同一引用（frozen，不可变）。
    - ``mechanical_groups`` / ``heat_groups``：按 ``kind`` 划分的 ``SimplifiedDirectedGroup`` 元组。
    - ``ideal_nodes``：ensure 时刻 ``layer.nodes`` 的只读映射。
    - ``nodes``：非理想节点表，首次偏移调用时从 ``ideal_nodes`` 拷贝。
    - ``properties``：与父层相同的物性求解器。
    - ``heat_efficiency`` / ``mechanical_efficiency``：最近一次偏移采用的系数；未偏移时为
      ``None``。

    仅由 ``from_closed_cycle_layer`` / ``ClosedCycleLayer.ensure_non_ideal()`` 创建。偏移须分别
    调用 ``apply_heat_pressure_offsets()`` / ``apply_mechanical_isentropic_offsets()``。
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
        return cls(
            simplified=simp,
            mechanical_groups=build_directed_groups(simp, "mechanical"),
            heat_groups=build_directed_groups(simp, "heat"),
            ideal_nodes=layer.nodes,
            properties=layer.properties,
        )

    def _ensure_nodes(self) -> dict[int, Node]:
        if self.nodes is None:
            self.nodes = {i: replace(n) for i, n in self.ideal_nodes.items()}
        return self.nodes

    def apply_heat_pressure_offsets(
        self,
        heat_efficiency: float | None = None,
    ) -> NonIdealClosedCycleLayer:
        """按换热组层号修正 ``P``：``P ← P_ideal × σ^layer``，再对全部节点 ``PS(P,S)`` 闭合 ``T,H``。

        ``σ`` 取参数 → ``self.heat_efficiency`` → ``config.NON_IDEAL_HEAT_EFFICIENCY_DEFAULT``；
        须在 ``(0, 1]`` 内。始终以理想 ``P`` 为底，避免重复叠加；不修改 ``ideal_nodes`` 与父层。
        细节见 ``docs/architecture.md §换热偏移``。
        """
        if heat_efficiency is not None:
            eta = heat_efficiency
        elif self.heat_efficiency is not None:
            eta = self.heat_efficiency
        else:
            eta = float(cyges_config.NON_IDEAL_HEAT_EFFICIENCY_DEFAULT)
        if not (0.0 < eta <= 1.0):
            raise ValueError(f"换热总压恢复系数 σ 须在 (0, 1] 内，收到 {eta!r}")

        nodes = self._ensure_nodes()

        # 仅换热组成员按层号叠乘 σ；始终以 ideal P 为底，避免重复偏移
        for group in self.heat_groups:
            layers = group.depth_dict()
            for v in group.nodes():
                p_ideal = self.ideal_nodes[v].P
                p_new = p_ideal * (eta ** layers[v])
                nodes[v] = replace(nodes[v], P=p_new)

        # 全表 PS 闭合：P 已变、S 仍理想拷贝，刷新 T,H 供后续机械步使用
        for v in list(nodes.keys()):
            _refresh_node_from_ps(nodes, self.properties, v)

        self.heat_efficiency = eta
        return self

    def apply_mechanical_isentropic_offsets(
        self,
        mechanical_efficiency: float | None = None,
    ) -> NonIdealClosedCycleLayer:
        """在每个机械组从锚点 DFS 单向推进 ``PS→η→HP``；锚点整点状态不变。

        ``η_is`` 取参数 → ``self.mechanical_efficiency`` →
        ``config.NON_IDEAL_MECHANICAL_EFFICIENCY_DEFAULT``；须在 ``(0, 1]`` 内。**须先**调用
        ``apply_heat_pressure_offsets()`` 让 ``P,H`` 反映换热后状态。不修改 ``ideal_nodes`` 与
        父层。机械组中与锚点无向不连通的节点 → ``RuntimeWarning``。

        锚点决策与机械步公式详见 ``docs/architecture.md §机械锚点决策、§机械步公式``。
        """
        if mechanical_efficiency is not None:
            eta = mechanical_efficiency
        elif self.mechanical_efficiency is not None:
            eta = self.mechanical_efficiency
        else:
            eta = float(cyges_config.NON_IDEAL_MECHANICAL_EFFICIENCY_DEFAULT)
        if not (0.0 < eta <= 1.0):
            raise ValueError(f"机械等熵效率 η_is 须在 (0, 1] 内，收到 {eta!r}")

        nodes = self._ensure_nodes()
        edges_by_key = self.simplified.edges_dict()

        for group in self.mechanical_groups:
            if not group.edge_keys:
                continue
            anchor = _pick_mechanical_anchor(group, self.ideal_nodes)
            known = _walk_mechanical_branches(
                nodes=nodes,
                properties=self.properties,
                group=group,
                edges_by_key=edges_by_key,
                anchor=anchor,
                eta_is=eta,
            )
            missing = group.nodes() - known
            if missing:
                warnings.warn(
                    f"机械有向组 {sorted(group.edge_keys)!r} 中节点 {sorted(missing)} "
                    f"与基准 {anchor} 在无向意义上不连通，已跳过",
                    RuntimeWarning,
                    stacklevel=2,
                )

        self.mechanical_efficiency = eta
        return self
