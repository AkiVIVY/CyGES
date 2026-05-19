# CyGES

闭式循环能源系统相关原型代码。当前主干在 **固定温压包线与分位轴** 上构建 **TP 离散拓扑**：一级 TP 网格、沿压力轴的 **等熵二级节点**、**机械边**（`M*`，叶轮机械工作过程）与 **换热边**（`H*`，换热过程），在 PS 约定下枚举 **最小子循环**（4 节点、4 边），并支持 **子循环质量流向量** 与 **量化步长**。物性由 **CoolProp** 经统一 `state` 接口提供。

拓扑上的「机械 / 换热」是**过程类型**的抽象，而非几何类别本身：非理想修正时，机械边对应**等熵效率** `η_is`，换热边对应**总压恢复系数** `σ`（见 [`config.py`](config.py) 与非理想层）。

---

## 运行环境与依赖

| 依赖 | 说明 |
|------|------|
| Python | 3.12+（当前验证环境） |
| [CoolProp](https://www.coolprop.org/) | `import CoolProp.CoolProp`，物性计算必需 |
| pytest | 单元测试 |
| matplotlib | 仅 `tests/` 下绘图用例需要 |

```bash
pip install CoolProp matplotlib pytest
```

**重要**：从仓库根运行代码或测试时，须将**项目根**加入 `PYTHONPATH`，否则 `import core` 与 `import config` 会失败（[`core`](core/) 为包，[`config.py`](config.py) 在根目录）。

```bash
cd /path/to/CyGES
set PYTHONPATH=.          # Windows CMD
# 或
$env:PYTHONPATH="."       # PowerShell
python -m pytest tests/ -q
```

单独运行绘图用例：

```bash
# 理想循环
python -m pytest tests/test_tp_topology.py::test_helium_topology_overview_plot -v
# 非理想（σ、η_is 偏移后）
python -m pytest tests/test_non_ideal_offsets_plot.py -v
```

输出 PNG（本地 pytest 生成，未纳入 git）：

| 文件 | 说明 |
|------|------|
| [`tests/ts_topology_he.png`](tests/ts_topology_he.png) | 理想层：T–S / P–S，节点、边、子循环 |
| [`tests/non_ideal_offsets_he.png`](tests/non_ideal_offsets_he.png) | 2×2：理想 vs 非理想（T–S、P–S） |
| [`tests/non_ideal_offsets_nodes_9_14_35_44.png`](tests/non_ideal_offsets_nodes_9_14_35_44.png) | 指定节点局部放大 |

---

## 仓库结构

### 当前主干（README 以这些为准）

| 路径 | 说明 |
|------|------|
| [`config.py`](config.py) | 全局可调数值：子循环初值系数、量化步长、非理想精简边统一效率默认值（见下表） |
| [`core/closed_cycle_layer.py`](core/closed_cycle_layer.py) | 闭式循环层：`ClosedCycleTPInput`、`Node`、`Edge`、`SubCycle`、`ClosedCycleLayer`、`SimplifiedEdge` / `SimplifiedTopology`；`build_axis`、`build_node_edge_topology`、`build_subcycles`、`filter_topology_for_non_ideal`、`build_simplified_topology`（每次 `analyze` / `commit` 后由 `ClosedCycleLayer._rebuild_simplified()` 自动同步） |
| [`core/non_ideal_closed_cycle_layer.py`](core/non_ideal_closed_cycle_layer.py) | 非理想层：`NonIdealClosedCycleLayer`（`ideal_nodes`/`nodes`、`apply_heat_pressure_offsets`）、`SimplifiedDirectedGroup`（层号含锚点重算）；`build_directed_groups`、`compute_group_downstream_depth` / `compute_group_downstream_reach` |
| [`core/fluid_property_solver.py`](core/fluid_property_solver.py) | `FluidPropertySolver` 协议与 `CoolPropFluidPropertySolver`（`state(pair,x,y)` → `T,P,H,S`） |
| [`tests/test_tp_topology.py`](tests/test_tp_topology.py) | 理想 He 循环绘图 → `ts_topology_he.png` |
| [`tests/test_non_ideal_offsets_plot.py`](tests/test_non_ideal_offsets_plot.py) | 非理想（换热 σ + 机械 η_is）偏移前后对比绘图 |
| [`core/__init__.py`](core/__init__.py) | `core` 包说明 |

### 其他目录（扩展/历史）

| 路径 | 说明 |
|------|------|
| [`inputs/`](inputs/)、[`solvers/`](solvers/)、[`optimize/`](optimize/) | 输入模式、求解器封装、优化占位等，与当前拓扑单测无强耦合 |
| [`oldFile/`](oldFile/) | 归档与旧实现，**不作为**当前 API 依据 |

---

## [`config.py`](config.py) 可调参数

| 常量 | 默认值 | 含义 |
|------|--------|------|
| `SUBCYCLE_INITIAL_MASS_FLOW_FRACTION_OF_MAX` | `0.1` | `analyze_topology` 后子循环初值：`每项 = 该系数 × max_mass_flow`；`max_mass_flow` 为 `None` 时按 `0` 参与乘法 |
| `SUBCYCLE_MASS_FLOW_STEP_FRACTION_DEFAULT` | `0.01` | `ClosedCycleTPInput` **未显式传入** `subcycle_mass_flow_step_fraction` 时采用；量化步长 `step = subcycle_mass_flow_step_fraction × max_mass_flow` |
| `NON_IDEAL_MECHANICAL_EFFICIENCY_DEFAULT` | `0.85` | 精简**机械边**（`SM*`，叶轮机械过程）统一**等熵效率** `η_is`；`apply_mechanical_isentropic_offsets` 默认读取：以已知端 `S` 在待求端 `P` 下 `PS` 得等熵焓 `H1`，再按 `η_is` 取真实焓 `H2`，最后 `HP` 闭合 |
| `NON_IDEAL_HEAT_EFFICIENCY_DEFAULT` | `0.99` | 精简**换热边**（`SH*`，换热过程）统一**总压恢复系数** `σ`；`apply_heat_pressure_offsets` 默认读取，`P ← P_ideal × σ^layer`，再对全部节点用 `PS(P,S)` 闭合 `T,H` |

### 边类型与非理想参数（物理含义）

| 拓扑 `kind` | 代表过程 | PS 离散方向 | 非理想参数 | 含义（概要） |
|-------------|----------|-------------|------------|----------------|
| `mechanical` / `M*`、`SM*` | 叶轮机械工作（压缩、膨胀等） | 沿压力轴（`edge_up` / `edge_down`） | `η_is` | **等熵效率**：实际过程相对等熵过程的效率 |
| `heat` / `H*`、`SH*` | 换热（加热、冷却等） | 沿等压线、熵方向（`edge_left` / `edge_right`） | `σ` | **总压恢复系数**：流经该换热过程后工质总压的保留程度 |

当前已实现：换热组内按层号用 `σ` 修正 `P`，并对**所有节点**用 `PS(P,S)` 闭合 `T,H`；机械组在各 `SimplifiedDirectedGroup` 内**从基准点**（优先一级节点，分叉时退化为 `min(upstream_special_nodes)`）**沿无向支路逐边单向**推进，每步由已知端的 `S` 在待求端 `P` 下 `PS` 得等熵 `H1`，再按等熵效率得真实 `H2` 并 `HP` 闭合。建议先换热、后机械。

---

## 单位约定

字段名不缀单位：

| 量 | 单位 | 典型位置 |
|----|------|----------|
| T | K | `Node.T`、`t_min` / `t_max` |
| P | kPa | `Node.P`、`p_min` / `p_max` |
| H | kJ/kg | `Node.H`、物性字典 `"H"` |
| S | kJ/(kg·K) | `Node.S`、物性字典 `"S"` |

CoolProp 内部为 SI；[`CoolPropFluidPropertySolver`](core/fluid_property_solver.py) 在边界与上述约定对齐。

---

## 物性 API

[`CoolPropFluidPropertySolver.state(pair, x, y)`](core/fluid_property_solver.py) 返回 **`ThermoStateTPHS`**（`T,P,H,S`）。

| `pair` | `x` | `y` |
|--------|-----|-----|
| `"TP"` | T [K] | P [kPa] |
| `"PS"` | P [kPa] | S [kJ/(kg·K)] |
| `"HP"` | H [kJ/kg] | P [kPa] |
| `"HS"` | H [kJ/kg] | S [kJ/(kg·K)] |

实现复用 `AbstractState("HEOS", fluid)`，适合拓扑生成中的高频调用。

---

## 闭式循环层：拓扑与子循环

### 推荐调用顺序

1. 构造 [`ClosedCycleTPInput`](core/closed_cycle_layer.py)：`fluid`、温压上下限、`t_quantiles` / `p_quantiles`（\[0,1\] 分位，调用方保证与端点互不重复）、可选 `max_mass_flow`、可选 `subcycle_mass_flow_step_fraction`（未传则用 [`config.SUBCYCLE_MASS_FLOW_STEP_FRACTION_DEFAULT`](config.py)）。
2. 构造 [`ClosedCycleLayer(inp, properties=..., auto_analyze=True)`](core/closed_cycle_layer.py)；默认 `properties` 为 `CoolPropFluidPropertySolver(inp.fluid)`；**默认 `auto_analyze=True`** 时构造末尾自动调用一次 **`analyze_topology()`**（`build_node_edge_topology` → `build_subcycles` → 初始化 **`subcycle_mass_flows`** 并 **`sync_subcycle_mass_flows_to_subcycles()`** → **`assign_edge_mass_flows_from_subcycles()`**）。若 `auto_analyze=False`，须随后自行调用 **`analyze_topology()`**。
3. 需要**重建**拓扑与子循环初值时，再次调用 **`analyze_topology()`**（会覆盖 `nodes` / `edges` / `subcycles` / `subcycle_mass_flows` 等）。

### PS 平面与有向边

**P** 自下而上增大，**S** 自左而右增大。边上方向 **`tail → head`**，且 **`P`、`S` 在容差意义下不减**。**机械边**（叶轮机械过程）对应节点 **`edge_up` / `edge_down`**；**换热边**（换热过程）对应 **`edge_right` / `edge_left`**（值为 `edges` 字典键）。

### `ClosedCycleLayer` 结果字段

- **`nodes: dict[int, Node]`**：键为全局 `index`；`parent is None` 为一级 TP 点，`parent` 为一级 `index` 为二级等熵点。
- **`edges: dict[str, Edge]`**：键 `M1,M2,…`（机械 / 叶轮机械过程）、`H1,H2,…`（换热 / 换热过程）；完成 `analyze_topology`（含构造时默认自动分析）后，`Edge.mass_flow` 由子循环汇聚写入；未分析前为 `None`。
- **`subcycles: list[SubCycle]`**：最小 4 节点 4 边环；`nodes` 顺序为左下→左上→右上→右下；`edges` 为左、上、右、下；`mass_flow` 由层同步。
- **`subcycle_mass_flows: list[float]`**：与 `subcycles[i]` 同索引；优化场景下优先改此列表，再 **`commit_subcycle_mass_flows_to_topology()`**（内部：量化 → `sync` → 赋边 → `_rebuild_simplified`），或分步 **`quantize_subcycle_mass_flows()`**、**`sync_subcycle_mass_flows_to_subcycles()`**、**`assign_edge_mass_flows_from_subcycles()`**。无子循环时为空列表。再次 **`analyze_topology()`** 会重置拓扑与本列表。**`analyze_topology()`** 与 **`commit_subcycle_mass_flows_to_topology()`** 之后均会重建 **`simplified`** 并清空 **`non_ideal`** 快照容器，须在理想层稳定后再 **`ensure_non_ideal()`**。
- **`skipped_points: list[SkippedPoint]`**：**仅诊断用**，记录 `build_node_edge_topology` 中被 CoolProp 异常静默跳过的候选点（一级 `TP` 或二级 `PS` 阶段），含 `stage / T / P / S / reason`；不参与拓扑与流量计算。每次 **`analyze_topology()`** 重置。
- **`simplified: SimplifiedTopology | None`**：精简 PS 拓扑（保留节点 `kept_nodes`、`SimplifiedEdge` 序列、合并占位映射 `merged_into`）。由 **`_rebuild_simplified()`** 在每次 `analyze_topology()` / `commit_subcycle_mass_flows_to_topology()` 末尾自动构建——先 **`filter_topology_for_non_ideal(nodes, edges, subcycles)`** 剔除不在任何子循环中的边及 `mass_flow` 为 `None` / 0 的边并同步清空节点邻边槽（**不修改**父层），再做同类型链合并；过滤后四邻边槽全空的节点记入 `simplified.merged_into`，值为占位常量 **`MERGED_ISOLATED_NODE_EDGE_KEY`**（无对应 `SimplifiedEdge`），且不出现在 `kept_nodes` 中。`len(subcycles) == 0` 时退化为空骨架并发出 **`RuntimeWarning`**。
- **`ensure_non_ideal()`** → [`NonIdealClosedCycleLayer`](core/non_ideal_closed_cycle_layer.py)：持有 `simplified`、`ideal_nodes`、`properties`，并构建机械/换热有向组。另调 **`apply_heat_pressure_offsets()`**（`σ`、修正 `P`，再 `PS(P,S)` 全表闭合 `T,H`）与 **`apply_mechanical_isentropic_offsets()`**（`η_is`、从基准沿支路单向 `PS→HP` 闭合）；均写入 `nodes` 拷贝，不改父层。父层 `analyze` / `commit` 之后会清空 `non_ideal`。

### 子循环流量（B 模式）

- **输入校验在外部**：本层不对 `ClosedCycleTPInput` 任何字段（`t/p` 上下限、分位序列、`max_mass_flow`、`subcycle_mass_flow_step_fraction` 等）做防御性校验，由调用方在业务侧自行保证合法性后再传入；本层只在算法不变量被破坏时抛错。
- **初值**：存在子循环时，`subcycle_mass_flows` 每项为 **`SUBCYCLE_INITIAL_MASS_FLOW_FRACTION_OF_MAX × max_mass_flow`**（`max_mass_flow` 为 `None` 时按 `0`）。
- **量化与提交**：`quantize_subcycle_mass_flows` 就地舍入，使用 `step = subcycle_mass_flow_step_fraction × max_mass_flow`；调用方须确保 `max_mass_flow` 与 `subcycle_mass_flow_step_fraction` 在量化前为有效正数，本层不做额外检查。**`commit_subcycle_mass_flows_to_topology()`** 在后续分析前一次性完成「量化 → 同步子循环 → 赋边」，保证列表与子循环、边一致。
- **长度一致性**：**`commit_subcycle_mass_flows_to_topology()`** 要求 ``len(subcycle_mass_flows)==len(subcycles)``，否则抛出 ``ValueError``（这是本层为数不多的运行时强校验）。

### 仅拓扑、不经 `ClosedCycleLayer`

可调用 **`build_node_edge_topology(solver, inp)`**、**`build_subcycles(nodes, edges)`**、**`build_axis(min, max, quantiles)`**；子循环初值向量与同步需自行处理。

---

## 测试一览

当前 `tests/` 仅保留**绘图**用例（需 CoolProp + matplotlib）：

| 用例 | 输出 |
|------|------|
| `test_helium_topology_overview_plot` | `ts_topology_he.png`（理想循环 T–S / P–S） |
| `test_helium_ideal_and_non_ideal_offsets_topology_plot` | `non_ideal_offsets_he.png`（理想 vs 偏移后 2×2） |
| `test_helium_focus_nodes_9_14_35_44_offsets_plot` | `non_ideal_offsets_nodes_9_14_35_44.png`（局部节点） |

---

## 开发与扩展

- 改网格/边：`build_node_edge_topology`；子循环枚举：`build_subcycles`。
- 改全局默认系数或步长比例：[`config.py`](config.py)。
- 换物性后端：实现 `FluidPropertySolver`，`ClosedCycleLayer(..., properties=...)` 注入。

若本文与源码不一致，以仓库内实际文件为准。
