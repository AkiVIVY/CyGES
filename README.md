# CyGES

闭式循环能源系统相关原型代码。当前主干在 **固定温压包线与分位轴** 上构建 **TP 离散拓扑**：一级 TP 网格、沿压力轴的 **等熵二级节点**、**机械边**（`M*`）与 **换热边**（`H*`），在 PS 约定下枚举 **最小子循环**（4 节点、4 边），并支持 **子循环质量流向量** 与 **量化步长**。物性由 **CoolProp** 经统一 `state` 接口提供。

---

## 运行环境与依赖

| 依赖 | 说明 |
|------|------|
| Python | 3.12+（当前验证环境） |
| [CoolProp](https://www.coolprop.org/) | `import CoolProp.CoolProp`，物性计算必需 |
| pytest | 单元测试 |
| matplotlib | 仅 [`tests/test_tp_topology.py`](tests/test_tp_topology.py) 中绘图用例需要 |

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

单独运行氦气拓扑综合图用例：

```bash
python -m pytest tests/test_tp_topology.py::test_helium_topology_overview_plot -v
```

成功后生成 [`tests/ts_topology_he.png`](tests/ts_topology_he.png)（双子图 **T–S** 与 **P–S**：节点、机械/换热边、子循环多边形及编号）。

---

## 仓库结构

### 当前主干（README 以这些为准）

| 路径 | 说明 |
|------|------|
| [`config.py`](config.py) | 全局可调数值：子循环初值系数、量化步长比例默认值（见下表） |
| [`core/closed_cycle_layer.py`](core/closed_cycle_layer.py) | 闭式循环层：`ClosedCycleTPInput`、`Node`、`Edge`、`SubCycle`、`ClosedCycleLayer`；`build_axis`、`build_node_edge_topology`、`build_subcycles` |
| [`core/fluid_property_solver.py`](core/fluid_property_solver.py) | `FluidPropertySolver` 协议与 `CoolPropFluidPropertySolver`（`state(pair,x,y)` → `T,P,H,S`） |
| [`tests/test_tp_topology.py`](tests/test_tp_topology.py) | 拓扑与子循环流量相关测试及 PNG 输出 |
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

**P** 自下而上增大，**S** 自左而右增大。边上方向 **`tail → head`**，且 **`P`、`S` 在容差意义下不减**。机械边对应节点 **`edge_up` / `edge_down`**；换热边对应 **`edge_right` / `edge_left`**（值为 `edges` 字典键）。

### `ClosedCycleLayer` 结果字段

- **`nodes: dict[int, Node]`**：键为全局 `index`；`parent is None` 为一级 TP 点，`parent` 为一级 `index` 为二级等熵点。
- **`edges: dict[str, Edge]`**：键 `M1,M2,…`（机械）、`H1,H2,…`（换热）；完成 `analyze_topology`（含构造时默认自动分析）后，`Edge.mass_flow` 由子循环汇聚写入；未分析前为 `None`。
- **`subcycles: list[SubCycle]`**：最小 4 节点 4 边环；`nodes` 顺序为左下→左上→右上→右下；`edges` 为左、上、右、下；`mass_flow` 由层同步。
- **`subcycle_mass_flows: list[float]`**：与 `subcycles[i]` 同索引；优化场景下优先改此列表，再 **`commit_subcycle_mass_flows_to_topology()`**（内部：量化 → `sync` → 赋边），或分步 **`quantize_subcycle_mass_flows()`**、**`sync_subcycle_mass_flows_to_subcycles()`**、**`assign_edge_mass_flows_from_subcycles()`**。无子循环时为空列表。再次 **`analyze_topology()`** 会重置拓扑与本列表。**`analyze_topology()`** 与 **`commit_subcycle_mass_flows_to_topology()`** 之后均会清空 **`non_ideal`**，须在理想层稳定后再 **`ensure_non_ideal()`**。
- **`skipped_points: list[SkippedPoint]`**：**仅诊断用**，记录 `build_node_edge_topology` 中被 CoolProp 异常静默跳过的候选点（一级 `TP` 或二级 `PS` 阶段），含 `stage / T / P / S / reason`；不参与拓扑与流量计算。每次 **`analyze_topology()`** 重置。

### 子循环流量（B 模式）

- **输入校验在外部**：本层不对 `ClosedCycleTPInput` 任何字段（`t/p` 上下限、分位序列、`max_mass_flow`、`subcycle_mass_flow_step_fraction` 等）做防御性校验，由调用方在业务侧自行保证合法性后再传入；本层只在算法不变量被破坏时抛错。
- **初值**：存在子循环时，`subcycle_mass_flows` 每项为 **`SUBCYCLE_INITIAL_MASS_FLOW_FRACTION_OF_MAX × max_mass_flow`**（`max_mass_flow` 为 `None` 时按 `0`）。
- **量化与提交**：`quantize_subcycle_mass_flows` 就地舍入，使用 `step = subcycle_mass_flow_step_fraction × max_mass_flow`；调用方须确保 `max_mass_flow` 与 `subcycle_mass_flow_step_fraction` 在量化前为有效正数，本层不做额外检查。**`commit_subcycle_mass_flows_to_topology()`** 在后续分析前一次性完成「量化 → 同步子循环 → 赋边」，保证列表与子循环、边一致。
- **长度一致性**：**`commit_subcycle_mass_flows_to_topology()`** 要求 ``len(subcycle_mass_flows)==len(subcycles)``，否则抛出 ``ValueError``（这是本层为数不多的运行时强校验）。

### 仅拓扑、不经 `ClosedCycleLayer`

可调用 **`build_node_edge_topology(solver, inp)`**、**`build_subcycles(nodes, edges)`**、**`build_axis(min, max, quantiles)`**；子循环初值向量与同步需自行处理。

---

## 测试一览

| 用例 | 作用 |
|------|------|
| `test_subcycle_mass_flow_defaults_zero_when_max_mass_flow_none` | `max_mass_flow` 为 `None` 时初值全 0 |
| `test_subcycle_mass_flow_step_fraction_defaults_from_config` | 未传 `subcycle_mass_flow_step_fraction` 时等于 `config` 默认 |
| `test_auto_analyze_false_then_manual_analyze` | `auto_analyze=False` 后手动 `analyze_topology` |
| `test_non_ideal_cleared_on_analyze` | `analyze_topology` 后清空 `non_ideal` |
| `test_baseline_snapshot_edge_mass_flow_detached` | 基准边快照与父层 `Edge` 解耦 |
| `test_non_ideal_cleared_after_commit` | `commit` 后清空 `non_ideal` |
| `test_commit_subcycle_mass_flows_len_mismatch_raises` | `commit_subcycle_mass_flows_to_topology` 长度不一致抛错 |
| `test_helium_topology_overview_plot` | He 宽网格、双子图 PNG |

---

## 开发与扩展

- 改网格/边：`build_node_edge_topology`；子循环枚举：`build_subcycles`。
- 改全局默认系数或步长比例：[`config.py`](config.py)。
- 换物性后端：实现 `FluidPropertySolver`，`ClosedCycleLayer(..., properties=...)` 注入。

若本文与源码不一致，以仓库内实际文件为准。
