# CyGES — Agent 交接说明

面向**接续开发的 Agent / 开发者**。用户向说明见 [`README.md`](README.md)；算法与实现细节**只在** [`docs/architecture.md`](docs/architecture.md) 维护一份，本文件与代码 docstring **仅引用、不重复公式**。

---

## 1. 文档分工

| 文件 | 读者 | 写什么 |
|------|------|--------|
| [`README.md`](README.md) | 用户 | 环境、快速开始、config、测试命令 |
| **本文件** | Agent | 范围、不变量、模块地图、协作约定 |
| [`docs/architecture.md`](docs/architecture.md) | 实现者 | PS 约定、流水线、精简、深度、偏置公式 |

改算法主干 → 先改 `architecture.md`，再改代码；README / AGENTS 仅在结构性变更时动。

---

## 2. 当前实现范围

**已实现**

- 理想层：PS 离散拓扑、最小 4 节点子循环、子循环流量、活跃子图精简（[`closed_cycle_layer.py`](core/closed_cycle_layer.py)）。
- 非理想偏置：有向组、组内 `reach` / `layer`（主脊分层）、单步 [`NonIdealClosedCycleLayer.apply_offsets()`](core/non_ideal_bias.py)（换热 `σ` → 机械组 PS 重置 → DFS `PS→η→HP`）。
- 性能统计：精简边过程归类与循环汇总（[`cycle_performance.py`](core/cycle_performance.py)），纯计算、零外部依赖；[`ClosedCycleLayer.performance_report()`](core/closed_cycle_layer.py) 为薄封装。
- 夹点分析：T-Q 曲线构建（`build_heat_tq_curves`）、夹点平移（`compute_pinch`）、公用工程需求（`analyze_pinch`），全部位于 [`postprocess.py`](core/postprocess.py)。
- 系统集成：外部冷热源转换、`SystemPipeline.run(props)` 四种夹点模式（`system_pinch` / `pinch` / `split_pinch` / `source_pinch`）位于 [`system.py`](core/system.py)。
- 物性：`PropertyRegistry` 多工质注册表 + CoolProp，`props(fluid, pair, x, y)` / `props.enthalpy(fluid, T, P)`。

**未实现**（勿写进 README 为已完成）

- 非理想方程/约束装配、多目标优化器。
- HEN / 多热源多冷源边界耦合。
- 按边 / 按组分别赋值 `σ`、`η_is`（当前共用 `config` 默认）。

---

## 3. 模块地图与分节

### [`core/closed_cycle_layer.py`](core/closed_cycle_layer.py)（理想层）

§1 数据模型 → §2 小工具 → §3 基础拓扑 → §4 子循环 → §5 精简 → §6 `ClosedCycleLayer`。

典型流水线：`analyze_topology` = `build_node_edge_topology` → `build_subcycles` → 初值流量 → 汇聚 → `_rebuild_simplified`。

### [`core/non_ideal_bias.py`](core/non_ideal_bias.py)（非理想偏置）

§1 分组 → §2 `SimplifiedDirectedGroup` → §3 深度 → §4 `build_directed_groups*` → §5 `NonIdealClosedCycleLayer` → §6 `apply_offsets`。

### [`core/cycle_performance.py`](core/cycle_performance.py)（性能统计）

§1 数据模型 → §2 判据小工具 → §3 状态解析 → §4 统计计算。只读；不参与 analyze/commit 失效链。公式见 architecture §10。

### [`core/postprocess.py`](core/postprocess.py)（二次处理）

§1 数据模型 → §2 T-Q 曲线构建 → §3 曲线插值与采样 → §4 夹点计算。T-Q 在 ``ProcessRecord`` 上构建（``build_heat_tq_curves``），夹点分析含底层 ``compute_pinch`` 与高层 ``analyze_pinch``。公式见 architecture §11。

### [`core/system.py`](core/system.py)（系统层）

§1 数据模型 → §2 冷热源转换 → §3 管线。``SystemPipeline.run(props)`` 四种夹点模式。公式见 architecture §9。

### [`optimize/`](optimize/)（优化器）

``Optimizer``: CMA-ES (BIPOP 重启) + DE，高斯基编码（维度固定），分位/流量/T 边界离散化，持久 ``ProcessPoolExecutor`` 多进程并行评估，``@register`` 目标函数注册表。

---

## 4. 失效与快照不变量

| 触发 | 行为 |
|------|------|
| `ClosedCycleLayer.analyze_topology()` | 重建 `nodes` / `edges` / `subcycles` / 初值流量 / `simplified`；**`non_ideal = None`** |
| `ClosedCycleLayer.commit_subcycle_mass_flows_to_topology()` | 量化 → 同步 → 汇聚 → 重建 `simplified`；**`non_ideal = None`** |
| `ClosedCycleLayer.ensure_non_ideal()` | 创建 `NonIdealClosedCycleLayer` 快照；要求 `layer.simplified is not None` |
| 父层 analyze / commit 之后 | 已存在的 `NonIdealClosedCycleLayer` 仍指向 ensure 时刻快照，与父层新数据**解耦** |

- `closed_cycle_layer` 对 `NonIdealClosedCycleLayer`：`TYPE_CHECKING` + `ensure_non_ideal()` 内延迟 import，避免循环依赖。
- `len(subcycles) == 0`：`simplified` 为空骨架 + `RuntimeWarning`，不调 `build_simplified_topology`。

---

## 5. 设计约束（改代码前必读）

**同节点两套深度**：同一 `Node.index` 可同时在一个机械组与一个换热组；`layer` 数值可不同。偏移/约束**必须**按 `kind + 组` 用 `group.depth_dict()` 查询，勿建全局 `node_index → 深度` 单表。见 architecture §7.4。

**非理想偏置顺序**：必须先 `ensure_non_ideal()`，再 `ni.apply_offsets()`；步骤 1 只改换热组 `P`，步骤 2/3 在机械组内 PS 重置后 DFS。见 architecture §8。

**主脊 `layer`**：最长路径定标 + 并列最小起点；脊延伸时后继须在「能到达 `max_d`」的集合内。见 architecture §7.3。

---

## 6. 测试（仓库内仅保留）

| 文件 | 用途 |
|------|------|
| [`tests/test_system_full.py`](tests/test_system_full.py) | 全管线可视化测试：冷热源 TS/PS、循环 TS/PS、循环夹点、循环性能、系统夹点（5 图） |
| [`tests/test_optimize.py`](tests/test_optimize.py) | 系统优化测试：CMA-ES 优化 + result 图 + Excel 累积日志 |

勿恢复已删的旧测试文件，除非用户明确要求。

---

## 7. 协作约定

- 不要自动 `git commit` / `push`，除非用户明确要求。
- 风格：`dataclass`、`frozen` 快照、中文 docstring、小步 diff；代码 docstring 只写**契约**（输入/输出/前后置/抛错）。
- 仓库 [`inputs/`](inputs/)、[`solvers/`](solvers/)、[`optimize/`](optimize/)、[`oldFile/`](oldFile/) 为历史或占位，**不作为当前 API 依据**。

---

## 8. 待办（优先级供参考）

1. 非理想方程/约束装配（边动量 / 能量平衡、约束闭包）。
2. 按边 / 按组分别赋值 `σ`、`η_is`。
3. `upstream_special_nodes: frozenset` → 单 `int` + 统一 tie-break（若 API 收紧）。
4. 机械边 `ΔS ≥ 0` 硬约束（可选）。
5. 多目标优化器与 HEN 边界耦合。
