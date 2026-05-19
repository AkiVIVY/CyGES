"""
非理想闭式循环层：理想层 ``simplified`` 的**只读快照** + 按过程类型划分的**有向连通组**索引。

边类型的物理含义
----------------
精简边 ``kind`` 与理想层 ``Edge`` 一致，表示不同的热力过程单元：

- ``mechanical``（``SM*``）：**叶轮机械工作过程**（压缩/膨胀等）。非理想参数为**等熵效率**
  ``η_is``（``config.NON_IDEAL_MECHANICAL_EFFICIENCY_DEFAULT``）；``apply_mechanical_isentropic_offsets``
  沿组内边修正焓并以 ``HP`` 闭合 ``T,P,H,S``。
- ``heat``（``SH*``）：**换热过程**（加热/冷却等）。非理想参数为**总压恢复系数**
  ``σ``（``config.NON_IDEAL_HEAT_EFFICIENCY_DEFAULT``），表示沿过程流向的总压保留程度；当前已实现
  按组内层号对节点压力 ``P`` 的偏移 ``P ← P_ideal × σ^layer``，并在所有节点上以 ``PS(P_new, S)``
  闭合 ``T,H``（``apply_heat_pressure_offsets``）。

与理想层的分工
--------------
- **理想层**（``closed_cycle_layer``）：生成/更新 ``SimplifiedTopology``（过滤、链合并、``SM*``/``SH*``）。
- **本模块**：在 ``ClosedCycleLayer.ensure_non_ideal()`` 时拷贝引用并**派生**组内结构，供后续
  非理想效率、节点状态偏移、约束装配使用；**不**回写父层 ``nodes`` / ``edges``。

挂载与失效
----------
- ``from_closed_cycle_layer`` 要求 ``layer.simplified`` 已存在。
- 父层 ``analyze_topology()`` / ``commit_subcycle_mass_flows_to_topology()`` 将 ``layer.non_ideal = None``；
  旧 ``NonIdealClosedCycleLayer`` 实例及其 ``simplified`` 引用仍指向 ensure 时刻，与父层新 ``simplified`` 解耦。

组内深度（重要）
----------------
在每个 ``SimplifiedDirectedGroup`` 内，**仅使用该组、该 kind** 的精简边（``tail→head``）建有向图：

- **层号**（``node_depth`` 第二分量）：最上游为 ``0``，沿流向递增；由 ``reach``、锚点顺流
  BFS 与逆流补全得到相对深度 ``rel`` 后，``layer(v) = max(rel) - rel(v)``。
- ``upstream_special_nodes`` = ``reach`` 最大的节点（下游延伸最长，``frozenset`` 可并列）。

同一 ``Node.index`` 可同时出现在一个机械组与一个换热组中，此时有两套独立的
深度，数值可以不同。偏移应通过 ``group.depth_dict()`` 按组查询，
勿用全局 ``node_index → 深度`` 单表混用。详见项目 ``AGENTS.md`` §5。

节点状态
--------
- ``ideal_nodes``：ensure 时刻对 ``layer.nodes`` 的只读引用（理想 ``T,P,H,S``）。
- ``nodes``：首次偏移方法调用时整表拷贝，再写修正量；理想层不变。
- ``properties``：与理想层相同的物性求解器，机械偏移用 ``HP`` 闭合状态。
"""

from __future__ import annotations

import warnings
from collections import defaultdict, deque
from dataclasses import dataclass, replace
from typing import Literal, Mapping

import config as cyges_config

from core.closed_cycle_layer import Node, SimplifiedEdge, SimplifiedTopology
from core.fluid_property_solver import FluidPropertySolver


def partition_simplified_edges_by_kind(
    topology: SimplifiedTopology,
) -> tuple[tuple[frozenset[str], ...], tuple[frozenset[str], ...]]:
    """
    将 ``topology.simplified_edges`` 按 ``SimplifiedEdge.kind`` 拆成机械、换热两套边表，再各自求**无向**连通分量。

    无向连通：每条精简边连接 ``tail`` 与 ``head``，忽略 ``tail→head`` 方向；**不**改变边的有向语义。

    返回
    ----
    ``(mechanical_edge_groups, heat_edge_groups)``，元素为边键 ``frozenset``。
    组顺序稳定：按组内最小边键字符串排序。

    后续 ``build_directed_groups`` 对每个边键集合再算组内有向层号。
    """
    mech: list[tuple[str, int, int]] = []
    heat: list[tuple[str, int, int]] = []
    for ek, se in topology.simplified_edges:
        if se.kind == "mechanical":
            mech.append((ek, se.tail, se.head))
        else:
            heat.append((ek, se.tail, se.head))

    def _components(items: list[tuple[str, int, int]]) -> tuple[frozenset[str], ...]:
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


@dataclass(frozen=True)
class GroupDepthMetrics:
    """组内有向深度中间量：``reach``（下游延伸）与公开 ``layer``（层号）。"""

    reach: dict[int, int]
    layer: dict[int, int]


def _group_adjacency(
    edges_by_key: dict[str, SimplifiedEdge],
    edge_keys: frozenset[str],
) -> tuple[dict[int, list[int]], set[int], list[tuple[int, int]]]:
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
    """
    在**单一有向连通组**上计算各端点的**下游可达边数** ``reach``。

    ``reach(v)`` = 从 ``v`` 出发沿 ``tail→head`` 能延伸的最长路径边数（无出边为 ``0``）。
    组内有向环时抛出 ``ValueError``。
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
    """
    ``reach``：下游最长有向路径（边数）；
    ``layer``：组内**任一**入度 0 源点到 ``v`` 的有向**最长路径**（边数），即 ``v`` 上游
    沿组内边累计经过的过程步数。

    在多源 DAG 上等价于「拓扑序 DP」：源点 ``layer=0``，对每条 ``u→v`` 有
    ``layer[v] = max(layer[v], layer[u] + 1)``。``compute_group_downstream_reach``
    已经保证组内无有向环。
    """
    reach = compute_group_downstream_reach(edges_by_key, edge_keys)
    if not reach:
        return GroupDepthMetrics(reach={}, layer={})

    adj, _, directed_edges = _group_adjacency(edges_by_key, edge_keys)
    in_degree: dict[int, int] = {v: 0 for v in reach}
    for _, w in directed_edges:
        in_degree[w] += 1

    layer: dict[int, int] = {v: 0 for v in reach}
    q: deque[int] = deque(sorted(v for v, d in in_degree.items() if d == 0))
    remaining = dict(in_degree)
    while q:
        u = q.popleft()
        for w in adj.get(u, ()):
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
    """
    在**单一有向连通组**上计算各端点**层号** ``layer``（``node_depth`` 公开语义）。

    ``layer(v)`` = 组内任一入度 0 源点到 ``v`` 的**有向最长路径**长度（边数）；
    源点本身 ``layer = 0``。基于拓扑序 DP：对每条 ``u → v`` 有
    ``layer[v] = max(layer[v], layer[u] + 1)``。

    示例：``A→B``、``B→C``、``A→D``、``E→D`` 时层号为 A=0,B=1,C=2,D=1,E=0
    （源点 A 与 E 同为 0，``D`` 取入边来源中最大的 layer 后 +1）。
    """
    return _compute_group_depth_metrics(edges_by_key, edge_keys).layer


@dataclass(frozen=True)
class SimplifiedDirectedGroup:
    """
    同 ``kind`` 精简边的一个**无向连通分量** + 组内**有向**层号与最上游特殊节点。

    字段
    ----
    - ``kind``：``"mechanical"``（叶轮机械过程组）或 ``"heat"``（换热过程组），与 ``SimplifiedEdge.kind`` 一致。
    - ``edge_keys``：本组 ``SM*`` / ``SH*`` 边键集合。
    - ``node_depth``：``(node_index, 层号)``，按 index 升序；最上游为 ``0``，最下游为 ``max_depth``。
    - ``upstream_special_nodes``：``reach`` 最大的节点（``frozenset``，允许并列）。

    查询
    ----
    - ``depth_dict()``：``index → 层号``，仅限本组语义。
    - ``nodes()`` / ``max_depth``：组内节点集与最大层号。

    注意：另一 ``kind`` 的组对同一 ``index`` 可能有不同层号，须分表使用。
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
    """
    对 ``topology`` 中指定 ``kind`` 的全部精简边：先无向连通划分，再逐组换算层号与 ``upstream_special_nodes``。

    与 ``partition_simplified_edges_by_kind`` 的关系：后者只分边键；本函数产出完整的
    ``SimplifiedDirectedGroup``（含层号与特殊节点），供 ``NonIdealClosedCycleLayer`` 持有。
    """
    mech_keys, heat_keys = partition_simplified_edges_by_kind(topology)
    edge_key_groups = mech_keys if kind == "mechanical" else heat_keys
    edges_by_key = topology.edges_dict()

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


def _pressure_equal(P0: float, P1: float) -> bool:
    tol = max(1e-9, 1e-6 * max(abs(P0), abs(P1)))
    return abs(P1 - P0) <= tol


def _pick_mechanical_anchor(
    group: SimplifiedDirectedGroup,
    ideal_nodes: Mapping[int, Node],
) -> int:
    """
    机械组基准节点（整点状态不再修改）：

    1. 若一级节点（``parent is None``）同时位于 ``upstream_special_nodes``，取该一级节点
       （并列时 index 最小）——一级节点本就是熵不动点。
    2. 否则使用 ``min(upstream_special_nodes)``——分叉组中一级位于侧枝时的熵不动点退化。
    3. 否则使用任一一级节点；再否则取层号 0（index 最小）。
    """
    group_nodes = group.nodes()
    primaries = sorted(v for v in group_nodes if ideal_nodes[v].parent is None)
    specials = group.upstream_special_nodes
    if primaries and (not specials or primaries[0] in specials):
        return primaries[0]
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
    """以当前 ``P,S`` 通过 ``PS`` 闪蒸刷新节点的 ``T,H``（``P,S`` 回写以求解器返回值为准）。"""
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
    """
    沿精简机械边 ``tail→head`` 对**待求端点**做一步非理想求解：

    1. 取已知端熵 ``S_known``，与待求端压力 ``P_unknown`` 一起做 ``PS`` 闪蒸，得**等熵焓** ``H1``；
    2. 按等熵效率 ``η_is`` 与边方向（``P_head`` vs ``P_tail``）取真实焓：

       - 压缩（``P_head > P_tail``）：``H2 = (H1 - H_known) / η_is + H_known``
       - 膨胀（``P_head < P_tail``）：``H2 = (H1 - H_known) × η_is + H_known``
       - 等压：``H2 = H_known``（等熵假设下 ``H1 ≡ H_known``）；

    3. 用 ``HP(H2, P_unknown)`` 写回待求端的完整 ``T,P,H,S``（``S`` 一般不再等于 ``S_known``）。

    本函数仅修改 ``nodes[unknown]``。
    """
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
    nodes: dict[int, Node],
    properties: FluidPropertySolver,
    group: SimplifiedDirectedGroup,
    edges_by_key: dict[str, SimplifiedEdge],
    anchor: int,
    eta_is: float,
) -> set[int]:
    """
    从锚点出发，沿无向邻接逐条支路单向 DFS 推进；每条边只从已知端推到未知端。

    锚点处的每个邻居各开启一条支路；在每个新已知节点上若仍有未访问邻居则继续沿链
    （子分叉同理 DFS）。返回处理完毕的节点集合（含锚点）。
    """
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
    """
    非理想分析容器（可变 dataclass，便于在本对象上挂载偏移结果）。

    字段
    ----
    - ``simplified``：与 ``ensure_non_ideal()`` 时刻 ``layer.simplified`` **同一引用**
      （``SimplifiedTopology`` 为 frozen，内容不可变）。
    - ``mechanical_groups`` / ``heat_groups``：按 kind 划分的 ``SimplifiedDirectedGroup`` 元组
      （含 ``node_depth``、``max_depth``、``upstream_special_nodes``）。
    - ``ideal_nodes``：ensure 时刻 ``layer.nodes`` 的只读映射（理想状态）。
    - ``nodes``：非理想节点表；``apply_heat_pressure_offsets()`` 首次调用时从 ``ideal_nodes`` 拷贝。
    - ``heat_efficiency``：最近一次换热压力偏移使用的总压恢复系数 ``σ``；未偏移时为 ``None``。
    - ``properties``：物性求解器（与 ``ClosedCycleLayer.properties`` 相同引用）。
    - ``mechanical_efficiency``：最近一次机械等熵偏移使用的 ``η_is``；未偏移时为 ``None``。

    构造
    ----
    仅通过 ``from_closed_cycle_layer`` 或 ``ClosedCycleLayer.ensure_non_ideal()`` 创建。
    换热 / 机械偏移须分别调用 ``apply_heat_pressure_offsets()``、``apply_mechanical_isentropic_offsets()``。
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
        """
        读取理想层当前 ``simplified`` 与 ``nodes``，构建机械/换热有向组快照。

        前置条件：``layer.simplified is not None``（已 ``analyze_topology`` 或 ``commit_*``）。
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
        """
        按**换热过程**有向组层号修正节点压力，并对全部节点用 ``PS`` 重新闭合 ``T,H``。

        换热边表示换热过程；``σ``（``heat_efficiency`` 参数）为**总压恢复系数**，
        ``P_non_ideal(v) = P_ideal(v) × σ ** layer(v)``（``layer`` 见 ``heat_groups`` 的 ``depth_dict``）。

        ``P`` 写完后，对 ``self.nodes`` 中**所有**节点用当前 ``(P, S)`` 调用 ``state("PS", P, S)``
        得到一致的 ``T,H``；这样后续机械步可在新 ``P`` 下直接使用各节点的 ``S``。

        - ``σ`` 默认 ``config.NON_IDEAL_HEAT_EFFICIENCY_DEFAULT``；可用参数或 ``self.heat_efficiency`` 覆盖。
        - 以各节点 **理想** ``P`` 为底，避免重复叠加。
        - 不修改 ``ideal_nodes`` 与父层 ``ClosedCycleLayer.nodes``。
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

        for group in self.heat_groups:
            layers = group.depth_dict()
            for v in group.nodes():
                p_ideal = self.ideal_nodes[v].P
                p_new = p_ideal * (eta ** layers[v])
                nodes[v] = replace(nodes[v], P=p_new)

        for v in list(nodes.keys()):
            _refresh_node_from_ps(nodes, self.properties, v)

        self.heat_efficiency = eta
        return self

    def apply_mechanical_isentropic_offsets(
        self,
        mechanical_efficiency: float | None = None,
    ) -> NonIdealClosedCycleLayer:
        """
        按**叶轮机械过程**有向组，从基准点出发沿每条支路**单向**逐边推进非理想偏置。

        基准节点（整点状态不再修改）：见 :func:`_pick_mechanical_anchor`，简言之优先一级节点；
        若一级位于侧枝（不在 ``upstream_special_nodes``），退化为 ``min(upstream_special_nodes)``。

        每一步沿一条精简机械边由已知端 ``k`` 推未知端 ``u``：

        1. 取 ``S_k``，在 ``P_u`` 下做 ``PS`` 闪蒸得**等熵焓** ``H1``；
        2. 按边方向（``P_head`` vs ``P_tail``）以 ``η_is`` 得真实焓 ``H2``
           （压缩 ``/η_is``、膨胀 ``×η_is``、等压保持）；
        3. ``HP(H2, P_u)`` 写回 ``T,P,H,S``。

        组内传播采用「从基准向每个邻居各开一条支路，沿无向链 DFS 远离基准单向推进」的方式，
        每条精简边只被使用一次。

        **须先** ``apply_heat_pressure_offsets()`` 再调用本方法（``P`` 与 ``H`` 用换热后值）。
        不修改 ``ideal_nodes`` 与父层 ``nodes``。
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
