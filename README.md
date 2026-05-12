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
| [`core/closed_cycle_layer.py`](core/closed_cycle_layer.py) | 闭式循环层：`ClosedCycleLayer`、`ClosedCycleTPInput`、`Node`、`Edge`；`analyze_topology()` 内部调用 ``build_node_edge_topology`` 填充 `nodes` 与 `edges` |
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
- **`self.nodes: dict[int, Node]`**：键 = 全局 **`index`**，不要求顺序。**`parent is None`** 为一级 TP 网格点；**`parent == 某一级 index`** 为二级等熵延伸点。每个节点另有 **`edge_up` / `edge_down` / `edge_left` / `edge_right`**（值为 **`self.edges`** 中的键或 **`None`**）：机械边写入 **上/下**，换热边写入 **右/左**（与 PS 约定一致）。
- **`self.edges: dict[str, Edge]`**：拓扑快照，**机械边**键为 **`M1`**, **`M2`**, …，**换热边**键为 **`H1`**, **`H2`**, …，均在同一字典中；**`Edge.kind`** 为 `"mechanical"` / `"heat"`，**`Edge.mass_flow`** 快照阶段为 **`None`**。
- **PS 平面约定**（与本模块边方向一致）：**P** 由小到大为自下而上，**S** 由小到大为自左而右；有向边 **`tail → head`** 满足 **`P`、`S` 均不减**（机械边沿等熵链 **P** 增大；换热边等 **P**、桶内按 **S** 相邻）。

以上为 **拓扑快照**；后续可在副本或新结构上扩展「用于实际计算的边」。

## 扩展时建议

1. 改网格/边逻辑：主要改 [`core/closed_cycle_layer.py`](core/closed_cycle_layer.py) 中的 ``build_node_edge_topology``（及轴采样 ``build_axis`` 若需调整分位策略）。
2. 改物性后端：实现与 **`FluidPropertySolver`** 协议兼容的类，在 **`ClosedCycleLayer(..., properties=...)`** 注入。
3. 新测试：沿用 **`PYTHONPATH=.`**；大图可只跑氦气用例以节省时间。

## 版本说明

本 README 反映截至编写时的代码结构；若与源码不一致，以仓库内实际文件为准。
