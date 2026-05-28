# `core/` 已知问题备忘

> ✅ = 已修复

## `non_ideal_bias.py`

### 高严重

### 1. `_mechanical_step_known_to_unknown` — 反向行走效率逆推是近似
**`line 581-585`**

```python
# 正向（tail→head）：精确
# 压缩: h_out = h_in + (h_isentropic_out - h_in) / η_is       ✓
# 膨胀: h_out = h_in + η_is * (h_isentropic_out - h_in)       ✓

# 反向（head→tail，出口推入口）：近似
# 膨胀反推真解: h_in = (h_out - η_is * h_isentropic_out) / (1 - η_is)
# 代码用: h_in = h_out + η_is * (h_isentropic_in - h_out)     ← 不等价
```

原因：DFS 行走时不携带 `S_inlet`，无法真算 `h_isentropic_out`。η_is 远离 1.0 时误差累积加剧。

**建议**：评估典型工况下误差量级；如可接受，在代码中加注释说明近似边界条件。

---

## 中严重

### 2. `_pick_anchor_and_sbase` — 无 primary 节点时直接用理想 S
**`line 529-533`**

```python
depths = group.depth_dict()
zeros = sorted(v for v in group_nodes if depths.get(v, -1) == 0)
if zeros:
    anchor = zeros[0]
    return anchor, nodes[anchor].S    # ← 此时仍是理想 S
```

有 primary 节点时，`s_base = state("TP", T_ideal, P_modified).S`，同时吸收了 σ 对 S 的间接影响。

无 primary 时，`s_base = nodes[anchor].S`（阶段 1 只改 P，S 仍是理想值）。此组的机械偏置完全基于理想熵启动。

**建议**：无 primary 时也按 `state("TP", T_ideal, P_modified)` 取 S，与 primary 路径一致。

---

### 3. `_pick_anchor_and_sbase` — 多 primary 节点选点无物理准则
**`line 522-523`**

```python
primaries = sorted(v for v in group_nodes if layer.ideal_nodes[v].parent is None)
if primaries:
    anchor = primaries[0]   # 确定性但任取最小 index
```

不同 anchor 给出不同 `s_base`，影响后续整个机械组的非理想熵偏移。应基于物理准则（如压力最高 / 温度最高 / 最长路径可达节点最多）选择。

**建议**：如架构文档未指定 tie-breaker，至少选 reach 最大的 primary 节点（最上游），因为 anchor 应尽可能"源"端。

---

### 4. `_walk_mechanical_branches` — 机械组有环时 DFS 静默丢弃边
**`line 613-631`**

```python
known: set[int] = {anchor}
stack = [(anchor, start_nb)]
while stack:
    prev, cur = stack.pop()
    if cur in known:       # 防重访 → 有环时回边被跳过
        continue
```

无向邻接的 DFS 以 `known` 集合防重访。如果机械组存在简化拓扑形成的环，环上某条边永远不会被访问，对应的非理想贡献丢失。当前没有任何环检测和 warning。

**建议**：加入环检测 —— DFS 完成后若 `group.nodes() - known` 非空，说明有连通的未被访问边（可能是环），发出 warning。

---

## 低严重

### 5. `_walk_mechanical_branches` — 外层 for 循环大部分迭代是 no-op
**`line 610`**

```python
for start_nb, _ in sorted(adj.get(anchor, [])):
    if start_nb in known:   # 首个邻居的 DFS 已遍历全连通分量，后续全是 continue
        continue
    stack = [(anchor, start_nb)]
    while stack:
        ...
```

`known` 集合跨循环共享，首个 DFS 就遍历完整个连通分量。后续邻居全是 no-op。功能正确但结构误导（暗示有多个独立 DFS 根）。

**建议**：如果图保证连通，外层循环可只取第一个邻居；否则留下也没问题。

---

### 6. `_apply_heat_pressure` — 中间状态 T/H/S 与 P 不一致
**`line 510`**

```python
# 只写 P，T/H/S 保留理想值
nodes[v] = replace(nodes[v], P=p_new)
```

阶段 1 后的 `nodes` 中 T/H/S 是理想值、P 是非理想值 —— 此状态物理上不一致。注释已标明"供步骤 2 取一级节点理想 T"，但若未来有代码在阶段 1 和 2 之间读取 node 状态，会读到不一致数据。

**建议**：加断言或标记位保护，或文档化此为严格不变式。

---

### 7. `_compute_layer_by_spine` — `min()` 在异常数据上可能抛空序列错误
**`line 276`**

```python
nxt = min(
    w for w in adj[cur]
    if dist[w] == dist[cur] + 1 and w in on_longest
)
```

如果 `on_longest` 计算在某条路径上出 bug，导致 `cur` 没有满足条件的后继，`min()` 抛 `ValueError: min() arg is an empty sequence`。提示信息不明。

**建议**：加保护检查 + 有意义的错误消息，或确保 `on_longest` 不变式成立。

---

### 8. `_resolve_efficiency` — 中文错误消息
**`line 488`**

```python
raise ValueError(f"{label} 须在 (0, 1] 内，收到 {eta!r}")
```

其余代码库异常用英文。

---

### 9. `_apply_combined_offsets` — `_ensure_nodes` 重复调用
**`line 692`**

```python
_ensure_nodes(layer)          # 外层调用
_apply_heat_pressure(layer)   # 内部再次 _ensure_nodes → 第二次是 no-op
```

第二次调用内 `layer.nodes is not None` → 直接返回，无实际副作用，仅冗余。

---

### 10. `build_directed_groups` — 同次 partition 算了两组 key 但只用一个
**`line 386-388`**

```python
mech_keys, heat_keys = partition_simplified_edges_by_kind(topology)
edge_key_groups = mech_keys if kind == "mechanical" else heat_keys
```

同时算了两组，丢一组。调用方若两种都需要，用 `build_directed_groups_both` 更高效。

---

# `core/closed_cycle_layer.py` 低优问题备忘

### ✅ C0a. `_rebuild_simplified` 不调用 `_invalidate_non_ideal()`
**`line 842`** — 已修复：`b669576` 开头新增 `self._invalidate_non_ideal()`。

### ✅ C0b. 子循环存在但重组后拓扑为空时静默
**`line 860`** — 已修复：`b669576` `build_simplified_topology` 返回空边时发出 `RuntimeWarning`。

### ✅ C0c. `sync_subcycle_mass_flows_to_subcycles` 无长度校验
**`line 863`** — 已修复：`b669576` 开头新增长度校验抛出 `ValueError`。

---

### C1. `_merge_chains_of_kind` 返回的 counter 被调用方丢弃
**`line 635, 751, 762`**

```python
def _merge_chains_of_kind(..., counter: int) -> tuple[...]:
    ...
    return simplified_edges, merged_into, counter    # 返回值含 counter

# 调用方
mech_segs, mech_mlinks, _ = _merge_chains_of_kind(...)   # 丢弃 counter
heat_segs, heat_mlinks, _ = _merge_chains_of_kind(...)   # 同样丢弃
```

因为 mech/heat 用不同前缀 `SM` vs `SH` 不冲突，所以每次从 0 重新编号没问题。但函数签名暗示 counter 有意义却无人使用。

**建议**：移除 counter 参数和返回，改为函数内部局部变量。或至少加注释说明。

---

### C2. `assign_edge_mass_flows_from_subcycles` 每次全量清空边流量
**`line 870-871`**

```python
for e in self.edges.values():
    e.mass_flow = None
```

O(E) 操作，调用频繁时可用差分更新（只重算变化的子循环）。当前全量重建是设计意图（简单可靠），留作低优。

---

### C3. `_find_typed_chains` adj 构建每节点遍历两槽
**`line 556-561`**

每次调用遍历所有节点、每节点 2 个槽做 `getattr`。节点数 O(100) 量级，实际无性能影响，仅结构上略为重复（`_group_adjacency` 在 non_ideal_bias 中有类似逻辑可复用）。

---

### C4. `_attach_edges_to_nodes_ps` 预分配所有节点 4 槽
**`line 234`**

```python
slot = {i: {"up": None, "down": None, "left": None, "right": None} for i in nodes}
```

每个节点预初始化全部 4 个槽为 None，但多数节点只有 2-3 个邻边。节点数少时无实际浪费。

---

### C5. `subcycle_mass_flows` 与 `subcycles[i].mass_flow` 双重存储
**`line 828` vs `SubCycle.mass_flow`**

同一数据存两处，需通过 `sync_subcycle_mass_flows_to_subcycles()` 手动同步。赋值 → commit 的模式已固化，但直接修改 `subcycles[i].mass_flow` 而不调 sync 会导致不一致。

**建议**：长期考虑统一单一数据源，或将 `SubCycle.mass_flow` 改为 property 代理到 `subcycle_mass_flows`。

---

### C6. `build_subcycles` 迭代 `nodes.values()` 顺序无物理含义
**`line 407`**

```python
for n0 in nodes.values():
```

Python 3.7+ dict 保序，顺序来自 `build_node_edge_topology` 的插入顺序（一级 TP 网格 → 二级等熵点）。这个顺序无物理意义。如果代码依赖 `subcycles[0]` 位置表示特定子循环（左下角等），这不安全。

**建议**：如未来有子循环排序需求，应在 `build_subcycles` 出结果后按 PS 坐标显式排序。
