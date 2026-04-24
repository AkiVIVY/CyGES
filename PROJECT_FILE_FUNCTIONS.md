# CyGES 项目文件与函数说明

本文档用于说明当前项目中各文件的作用，以及主要类/函数的功能边界，方便后续重构与协作。

---

## 1. 项目整体结构

当前项目已切换为新框架主线（`core/` + `SystemModel`），旧实现文件已从仓库工作区清理。

---

## 2. 根目录文件说明

### 已清理的旧文件

以下旧架构文件和历史测试文件已移除，不再作为工作入口：

- `propertySolver.py`
- `createModel.py`
- `test.py`
- `CycleModel.py`
- `creatCycle.py`

### `config.py`

**作用**：项目默认参数配置。  
**主要常量**：

- `DEFAULT_SUBCYCLE_M_DOT`
- `DEFAULT_SUBCYCLE_M_DOT_MIN`
- `DEFAULT_SUBCYCLE_M_DOT_MAX`
- `DEFAULT_EDGE_EFF`
- `DEFAULT_TOLERANCE`

---

### `.gitignore`

**作用**：忽略 Python 缓存文件。  
当前规则：

- `__pycache__/`
- `*.pyc`

---
---------------------------------------------------------------------------------------------------
## 3. 新架构目录说明

### `core/`（核心对象层）

#### `core/system_model.py`

**类**：`SystemModel`（系统总入口）

- `from_spec(spec)`
  - 从输入字典构造系统对象（热源、冷源、闭式循环）
- `build()`
  - 触发每个闭式循环生成拓扑
- `solve()`
  - 汇总各循环功热结果
  - 输出系统总功、总吸热、总放热

---

#### `core/closed_cycle.py`

**类**：`ClosedCycle`（闭式循环对象）

关键字段（双态设计）：

- `nodes_raw / nodes_working`
- `edges_raw / edges_working`
- `node_groups_in_p / node_groups_in_s`
- `subcycles`
- `topology_diagnostics`

主要方法：

- `generate_topology()`
  - 按边界与分位点生成拓扑：
    1. 输入校验
    2. TP 一级节点
    3. PS 二级节点
    4. 去重
    5. 分组
    6. 全局边生成
    7. 子循环提取
  - 分位点输入仅要求至少1个，函数内部自动补齐 `0` 与 `1`
  - 同时输出拓扑诊断信息（A/B节点数、去重点数量、无效子循环数、分组数量）
  - 每个 `ClosedCycle` 仅创建一次 `CoolPropSolver`，并在后续复用

- `calc_cycle_balance()`
  - 汇总全部子循环的功率与换热

- `_snapshot_topology_diagnostics()`
  - 记录节点/边/子循环数量诊断
  - 包含：
    - `n_nodes_a`
    - `n_nodes_b`
    - `n_nodes_dedup_removed`
    - `n_invalid_subcycles`
    - `n_pressure_groups`
    - `n_entropy_groups`

- `_normalized_key(value, tolerance)`
  - 浮点归一化辅助
  - 小数位精度与 `tolerance` 相关，不再写死固定小数位

- `_complete_levels(levels, level_name)`
  - 校验分位点输入并自动补齐 `0` 和 `1`

- `_group_nodes(node_list, key_attr, tolerance)`
  - 按 `P` 或 `S` 分组节点

- `_append_group_edges(grouped_nodes, sort_attr, family)`
  - 在分组中排序后创建边

- `_extract_subcycles(group_in_p)`
  - 从网格提取 `SubCycle` 对象
  - 返回 `subcycles` 与 `invalid_subcycle_count`
  - 子循环提取仅依赖节点分布，不依赖边信息

---

#### `core/subcycle.py`

**类**：`SubCycle`（可优化基本单元）

字段：

- `nodes_raw / nodes_working`
- `m_dot`, `m_dot_min`, `m_dot_max`
- `metadata`

方法：

- `calc_power()`
  - 直接由子循环左右两侧节点焓差计算功率指标
- `calc_heat()`
  - 直接由子循环上下两侧节点焓差计算换热指标

---

#### `core/process_edge.py`

**类**：`ProcessEdge`

字段：

- `edge_id`, `edge_type`
- `upstream`, `downstream`
- `role`（`left/right/top/bottom/unknown`）
- `eff`, `m_dot`, `constraints_ref`

方法：

- `delta_h()`：焓差
- `qw_dot()`：热/功率流率（`m_dot * delta_h`）

---

#### `core/state_node.py`

**类**：`StateNode`

字段：

- `node_id`, `fluid`, `T`, `P`, `H`, `S`
- `source_tag`

方法：

- `clone(node_id=None)`
  - 复制节点（用于 raw/working 双轨）

---

#### `core/external_stream.py`

**类**：`ExternalStream`（热源/冷源统一抽象）

- 字段：`stream_id`, `stream_type`, `fluid`, `m_dot`, `segments`
- 方法：`to_tq_segments()`（返回 TQ 段数据）

---

#### `core/units.py`

**作用**：单位转换工具（工程单位约定）。

- `pa_to_kpa()`
- `kpa_to_pa()`
- `jpkg_to_kjpkg()`
- `kjpkg_to_jpkg()`

---

#### `core/__init__.py`

**作用**：对外导出核心类入口。

---

### `inputs/`（输入校验层）

#### `inputs/input_schema.py`

**函数**：`validate_system_spec(spec)`

校验内容：

- 顶层键是否齐全：
  - `hot_streams`
  - `cold_streams`
  - `closed_cycles`
- 流体流信息必填键：
  - `stream_id`, `fluid`, `m_dot`
- 闭式循环必填键：
  - `cycle_id`, `fluid`, `boundary`, `levels`
- 分位点必填键：
  - `TLevel`, `PLevel`（每项至少1个值）

#### `inputs/__init__.py`

导出 `validate_system_spec`

---

### `solvers/`（求解器适配层）

#### `solvers/coolprop_solver.py`

**类**：`CoolPropSolver`（新架构原生物性求解器）

- `__init__(substance_list)`
  - 初始化并缓存 CoolProp `AbstractState`
  - 校验工质名是否合法
- `solve(pair, substance, input1, input2)`
  - 支持 `HP/TP/HS/PS`
  - 完成工程单位与 CoolProp SI 单位转换
  - 返回标准物性字典 `T/P/H/S/substance`
- `_to_coolprop_inputs(pair, input1, input2)`
  - 不同输入对下的单位/顺序转换辅助函数

#### `solvers/__init__.py`

导出 `CoolPropSolver`

---

### `optimize/`

#### `optimize/__init__.py`

当前为占位，后续用于放置优化变量、目标函数、约束接口。

---

### `tests/`

#### `tests/smoke_build.py`

**作用**：新架构最小冒烟测试。  
流程：

1. 生成 `sample_spec()`（示例分位点仅输入中间值，`0/1` 在闭式循环内自动补齐）
2. `validate_system_spec()`
3. `SystemModel.from_spec()`
4. `build()`
5. `solve()`
6. 绘制闭式循环 `T-S` 图（节点、边、子循环）
7. 输出构建与功率信息

---

## 4. 图像与资料文件说明

旧架构生成的图片与历史资料文件已从当前工作区清理。

---

## 5. 当前推荐入口

- **验证新架构**：`tests/smoke_build.py`
- **系统总入口**：`core/system_model.py`
- **物性求解器（新）**：`solvers/coolprop_solver.py`

---

## 6. 后续重构建议（简版）

1. 将 `createModel.py` 的旧逻辑逐步迁到 `core/closed_cycle.py`
2. 新增 `TQ` 数据模型，接入 `ExternalStream` 与 `SubCycle`
3. 在 `optimize/` 中定义 `decision variables / objective / constraints`
4. 为 `SubCycle` 增加效率扰动下的 `raw -> working` 迁移方法

