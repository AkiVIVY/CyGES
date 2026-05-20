# CyGES — Agent 交接说明

面向**接续开发的 Agent / 开发者**。

- 用户向使用说明：[`README.md`](README.md)。
- 算法细节单一信息源：[`docs/architecture.md`](docs/architecture.md)。**所有 PS 约定、链合并规则、有向组深度定义、机械锚点决策、机械步公式、换热 σ 计算只在 architecture 文档中维护一份**，本文件与代码 docstring 仅引用。

---

## 1. 当前实现范围

- 在给定工质与温压包线上构建 PS 平面离散拓扑、最小 4 节点子循环、子循环质量流向量、精简拓扑。
- 非理想偏移：换热 `P ← P_ideal × σ^layer` + 全表 `PS(P, S)` 闭合 `T, H`；机械 `PS → η_is → HP` 沿支路单向推进。
- 物性走 CoolProp 统一 `state(pair, x, y)` 接口。

**未实现**（请勿在文档中描述为已完成）：

- 非理想方程/约束装配、多目标优化器。
- 换热网络（HEN）与多热源/多冷源边界耦合。
- 按边 / 按组分别赋值 `σ`、`η_is`（当前共用 `config` 默认）。

---

## 2. 失效与快照不变量

| 触发 | 行为 |
|------|------|
| `ClosedCycleLayer.analyze_topology()` | 重建 `nodes` / `edges` / `subcycles` / `subcycle_mass_flows`，重建 `simplified`，**置 `non_ideal = None`** |
| `ClosedCycleLayer.commit_subcycle_mass_flows_to_topology()` | 量化 → 同步 → 汇聚 → 重建 `simplified`，**置 `non_ideal = None`** |
| `ClosedCycleLayer.ensure_non_ideal()` | 创建 `NonIdealClosedCycleLayer` 快照；要求 `layer.simplified is not None` |
| 父层 analyze / commit | 已生成的 `NonIdealClosedCycleLayer` 仍指向 ensure 时刻的 `simplified`，与父层新数据解耦 |

`closed_cycle_layer` 对 `NonIdealClosedCycleLayer` 使用 `TYPE_CHECKING` + `ensure_non_ideal()` 内延迟 import，避免循环依赖。

**空子循环**：`len(subcycles) == 0` 时 `simplified` 为空骨架 + `RuntimeWarning`，不调用 `build_simplified_topology`。

---

## 3. 设计取舍：同节点两套深度

同一 `Node.index` 可同时位于一个机械组与一个换热组中，在两个不同有向子图上各自得到 `layer`，**数值可以不同**——这是预期行为。

- 深度挂在 `SimplifiedDirectedGroup.node_depth`（已含 `kind`）。
- 偏移/约束**必须**按 `kind + 组` 通过 `group.depth_dict()` 查询，不要构造全局 `node_index → 深度` 单表。
- 若下游确实需要按节点汇总，应显式维护 `depth_mech` / `depth_heat` 两套字典。

定义与示例见 [`docs/architecture.md §7`](docs/architecture.md)。

---

## 4. 待办（优先级供参考）

1. **换热 `layer` 改为主脊分层**（已知偏差，见 [`docs/architecture.md §7.5`](docs/architecture.md)）：现实现为「多源入度 0 均为 0」；应为「最长路径定标 + 支流/汇流对齐」，例如 `A→B→C, A→D, E→C` 时 `E` 应为 1 而非 0。
2. 非理想方程/约束装配（边动量 / 能量平衡、约束闭包）。
3. 按边 / 按组分别赋值 `σ`、`η_is`（替代当前 config 默认值的"全局共用"）。
4. API 收紧：若确认每组只需一个特殊节点，将 `upstream_special_nodes: frozenset` 改为 `int` 并统一 tie-break。
5. 熵单调性硬约束（当前算法仅在大多数物理情况下保证沿机械边 `ΔS ≥ 0`）；如需强约束可在偏移后校验或换闭合策略。
6. 多目标优化器与 HEN 边界耦合。

---

## 5. 协作约定

- 不要自动 `git commit` / `push`，除非用户明确要求。
- 风格：`dataclass`、`frozen` 快照、中文 docstring、小步 diff。代码 docstring 只写**契约**（输入/输出/前后置/抛错），算法细节一律放到 [`docs/architecture.md`](docs/architecture.md)。
- 修改算法主干时同步更新 architecture 文档，README / AGENTS 仅在结构性变更时改动。
- 仓库其它目录 [`inputs/`](inputs/)、[`solvers/`](solvers/)、[`optimize/`](optimize/)、[`oldFile/`](oldFile/) 为历史或占位，**不作为当前 API 依据**。

---

*最后更新：文档分层（architecture 为算法 SoT）、`layer` 主脊语义已知问题见 architecture §7.5。*
