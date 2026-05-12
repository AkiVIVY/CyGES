# CyGES — 闭式循环层与物性（当前实现摘要）

本文档供后续 Agent / 开发者快速接续工作，描述仓库中与 **闭式循环层拓扑**、**物性**、**测试** 相关的约定与入口。

## 环境

- **Python 3.12+**（项目当前验证环境）
- **CoolProp**：物性计算必需（`import CoolProp.CoolProp`）
- **matplotlib**：仅 `tests/test_tp_topology.py::test_helium_tp_topology_ts_plot` 作图需要

安装示例：

```bash
pip install CoolProp matplotlib pytest
```

运行测试时请将项目根加入 `PYTHONPATH`：

```bash
cd /path/to/CyGES
set PYTHONPATH=.   # Windows CMD
# 或 $env:PYTHONPATH="."  # PowerShell
python -m pytest tests/test_tp_topology.py -q
```

仅跑氦气 T–S 图用例：

```bash
python -m pytest tests/test_tp_topology.py::test_helium_tp_topology_ts_plot -v
```

图输出路径：`tests/ts_topology_he.png`。

## 目录与核心文件

| 路径 | 作用 |
|------|------|
| [`core/closed_cycle_layer.py`](core/closed_cycle_layer.py) | 闭式循环层：`ClosedCycleLayer`、`ClosedCycleTPInput`、`Node`、`Edge`；`analyze_topology()` 填充 `nodes` 与快照边 |
| [`core/fluid_property_solver.py`](core/fluid_property_solver.py) | CoolProp `AbstractState("HEOS", fluid)`；统一 `state(pair, x, y)` 返回 `T,P,H,S` 字典 |
| [`tests/test_tp_topology.py`](tests/test_tp_topology.py) | 拓扑与物性相关单元测试及氦气 T–S 可视化 |

## 单位约定（全项目一致）

| 量 | 单位 | 说明 |
|----|------|------|
| T | K | `Node.T`、输入 `t_min` / `t_max` |
| P | kPa | `Node.P`、输入 `p_min` / `p_max` |
| H | kJ/kg | `Node.H`、物性字典 `"H"` |
| S | kJ/(kg·K) | `Node.S`、物性字典 `"S"` |

物性求解器内部用 SI，在边界换算。

## 物性 API（`CoolPropFluidPropertySolver`）

- **`state(pair, x, y) -> dict`**，键为 **`"T","P","H","S"`**（与 `Node` 字段对应）。
- **`pair`**：`"TP"` | `"PS"` | `"HP"` | `"HS"`，`(x,y)` 含义见该文件中文 docstring。

## 闭式循环层（`ClosedCycleLayer`）

### 输入

- **`ClosedCycleTPInput`**：含 **`fluid`**、温压上下限、**`t_quantiles` / `p_quantiles`**（\[0,1\] 分位，调用方保证不重复）、可选 **`max_mass_flow`**。

### 拓扑分析

- **`analyze_topology() -> None`**：仅修改层实例，**无返回值**。
- 写 **`self.nodes: dict[int, Node]`**（键 = 全局 **`index`**，不要求顺序）。
- **`parent is None`**：一级 TP 网格点；**`parent == 某一级 index`**：二级等熵延伸点。
- **`self.mechanical_edges`**：等熵离散链上按压力排序的相邻有向边（`tail` → `head`）。
- **`self.heat_edges`**：全节点按等 **P** 分桶、桶内按 **T** 排序后的相邻换热边。
- **`Edge.mass_flow`**：快照阶段为 **`None`**，后续业务再赋值；**`kind`** 为 `"mechanical"` / `"heat"`。

以上为 **拓扑快照**；后续可在副本或新结构上扩展「用于实际计算的边」。

## 扩展时建议

1. 改网格/边逻辑：主要动 [`core/closed_cycle_layer.py`](core/closed_cycle_layer.py) 中 `build_grid_nodes`、`expand_isentropic_nodes`、`build_heat_edges`。
2. 改物性后端：实现与 **`FluidPropertySolver`** 协议兼容的类，在 **`ClosedCycleLayer(..., properties=...)`** 注入。
3. 新测试：沿用 **`PYTHONPATH=.`**；大图可只跑氦气用例以节省时间。

## 版本说明

本 README 反映截至编写时的代码结构；若与源码不一致，以仓库内实际文件为准。
