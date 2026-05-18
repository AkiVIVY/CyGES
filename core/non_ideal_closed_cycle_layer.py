"""
非理想闭式循环层：理想层 ``simplified`` 的**只读快照** + 按机械/换热划分的**有向连通组**索引。

与理想层的分工
--------------
- **理想层**（``closed_cycle_layer``）：生成/更新 ``SimplifiedTopology``（过滤、链合并、``SM*``/``SH*``）。
- **本模块**：在 ``ClosedCycleLayer.ensure_non_ideal()`` 时拷贝引用并**派生**组内结构，供后续
  非理想效率、节点状态偏移、约束装配使用；**不**回写 ``nodes`` / ``edges``。

挂载与失效
----------
- ``from_closed_cycle_layer`` 要求 ``layer.simplified`` 已存在。
- 父层 ``analyze_topology()`` / ``commit_subcycle_mass_flows_to_topology()`` 将 ``layer.non_ideal = None``；
  旧 ``NonIdealClosedCycleLayer`` 实例及其 ``simplified`` 引用仍指向 ensure 时刻，与父层新 ``simplified`` 解耦。

下游深度（重要）
----------------
在每个 ``SimplifiedDirectedGroup`` 内，**仅使用该组、该 kind** 的精简边（``tail→head``）建有向图：

- **深度** = 从节点 ``v`` 出发能走的最长有向路径**边数**（无出边为 ``0``）；存于 ``node_depth`` 的第二个分量。
- ``upstream_special_nodes`` = 组内深度取 **最大值** 的节点集合（**可并列**）。

同一 ``Node.index`` 可同时出现在一个机械组与一个换热组中，此时有两套独立的
深度（机械深度 / 换热深度），数值可以不同。偏移与约束应通过 ``group.depth_dict()`` 按组查询，
勿用全局 ``node_index → 深度`` 单表混用。详见项目 ``AGENTS.md`` §5。
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Literal

from core.closed_cycle_layer import SimplifiedEdge, SimplifiedTopology


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

    后续 ``build_directed_groups`` 对每个边键集合再算组内有向深度。
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


def compute_group_downstream_depth(
    edges_by_key: dict[str, SimplifiedEdge],
    edge_keys: frozenset[str],
) -> dict[int, int]:
    """
    在**单一有向连通组**（同一 ``kind``、无向连通的一组精简边）上计算各端点的下游深度。

    定义
    ----
    ``head`` 为边的下游端点。**深度** = 从 ``v`` 出发沿 ``tail→head`` 能延伸的最长路径边数。

    示例（机械组内）：边 ``A→B``、``B→C``、``A→D``、``E→D`` 时，
    A 的深度为 2（经 B 到 C），B 与 E 为 1，C 与 D 为 0；深度最大者为 ``{A}``。

    异常
    ----
    组内存在有向环时无法定义深度，抛出 ``ValueError``。

    返回
    ----
    仅包含出现在 ``edge_keys`` 中的 ``tail``/``head`` 节点 index。
    """
    if not edge_keys:
        return {}

    adj: dict[int, list[int]] = defaultdict(list)
    nodes: set[int] = set()
    for ek in edge_keys:
        se = edges_by_key[ek]
        adj[se.tail].append(se.head)
        nodes.add(se.tail)
        nodes.add(se.head)

    memo: dict[int, int] = {}
    visiting: set[int] = set()

    def depth(v: int) -> int:
        if v in memo:
            return memo[v]
        if v in visiting:
            raise ValueError(
                f"精简边组 {sorted(edge_keys)!r} 内存在有向环（涉及节点 {v}），无法定义下游深度"
            )
        visiting.add(v)
        best = 0
        for u in adj.get(v, ()):
            best = max(best, 1 + depth(u))
        visiting.remove(v)
        memo[v] = best
        return best

    return {v: depth(v) for v in nodes}


@dataclass(frozen=True)
class SimplifiedDirectedGroup:
    """
    同 ``kind`` 精简边的一个**无向连通分量** + 组内**有向**下游深度与最上游特殊节点。

    字段
    ----
    - ``kind``：``"mechanical"`` 或 ``"heat"``，与 ``SimplifiedEdge.kind`` 一致。
    - ``edge_keys``：本组 ``SM*`` / ``SH*`` 边键集合。
    - ``node_depth``：``(node_index, 深度)`` 元组，按 index 升序；仅含与本组边相连的端点。
    - ``upstream_special_nodes``：深度等于 ``max_depth`` 的节点（``frozenset``，允许并列）。

    查询
    ----
    - ``depth_dict()``：``index → 深度``，仅限本组语义。
    - ``nodes()`` / ``max_depth``：组内节点集与最大深度。

    注意：另一 ``kind`` 的组对同一 ``index`` 可能有不同深度，须分表使用。
    """

    kind: Literal["mechanical", "heat"]
    edge_keys: frozenset[str]
    node_depth: tuple[tuple[int, int], ...]
    upstream_special_nodes: frozenset[int]

    def depth_dict(self) -> dict[int, int]:
        """本组内 ``node_index → 深度``（勿与另一 kind 的组混用）。"""
        return dict(self.node_depth)

    def nodes(self) -> frozenset[int]:
        """组内作为精简边端点出现的所有 ``Node.index``。"""
        return frozenset(i for i, _ in self.node_depth)

    @property
    def max_depth(self) -> int:
        """组内深度的最大值；``node_depth`` 为空时为 ``0``。"""
        if not self.node_depth:
            return 0
        return max(d for _, d in self.node_depth)


def build_directed_groups(
    topology: SimplifiedTopology,
    kind: Literal["mechanical", "heat"],
) -> tuple[SimplifiedDirectedGroup, ...]:
    """
    对 ``topology`` 中指定 ``kind`` 的全部精简边：先无向连通划分，再逐组计算深度与 ``upstream_special_nodes``。

    与 ``partition_simplified_edges_by_kind`` 的关系：后者只分边键；本函数产出完整的
    ``SimplifiedDirectedGroup``（含深度与特殊节点），供 ``NonIdealClosedCycleLayer`` 持有。
    """
    mech_keys, heat_keys = partition_simplified_edges_by_kind(topology)
    edge_key_groups = mech_keys if kind == "mechanical" else heat_keys
    edges_by_key = topology.edges_dict()

    out: list[SimplifiedDirectedGroup] = []
    for keys in edge_key_groups:
        depths = compute_group_downstream_depth(edges_by_key, keys)
        if depths:
            d_max = max(depths.values())
            special = frozenset(v for v, d in depths.items() if d == d_max)
        else:
            special = frozenset()
        node_depth = tuple(sorted(depths.items()))
        out.append(
            SimplifiedDirectedGroup(
                kind=kind,
                edge_keys=keys,
                node_depth=node_depth,
                upstream_special_nodes=special,
            )
        )
    return tuple(out)


@dataclass
class NonIdealClosedCycleLayer:
    """
    非理想分析容器（可变 dataclass，便于后续在本对象上挂载偏移结果等字段）。

    字段
    ----
    - ``simplified``：与 ``ensure_non_ideal()`` 时刻 ``layer.simplified`` **同一引用**
      （``SimplifiedTopology`` 为 frozen，内容不可变）。
    - ``mechanical_groups`` / ``heat_groups``：按 kind 划分的 ``SimplifiedDirectedGroup`` 元组
      （含 ``node_depth``、``max_depth``、``upstream_special_nodes``），组序与 ``partition_simplified_edges_by_kind`` 一致。

    构造
    ----
    仅通过 ``from_closed_cycle_layer`` 或 ``ClosedCycleLayer.ensure_non_ideal()`` 创建。
    """

    simplified: SimplifiedTopology
    mechanical_groups: tuple[SimplifiedDirectedGroup, ...]
    heat_groups: tuple[SimplifiedDirectedGroup, ...]

    @classmethod
    def from_closed_cycle_layer(cls, layer) -> NonIdealClosedCycleLayer:
        """
        读取理想层当前 ``simplified``，构建机械/换热有向组快照。

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
        )
