# CyGES

闭式循环能源系统相关原型代码：当前主干围绕 **TP 离散拓扑**、**CoolProp 物性** 与 **最小子循环** 展开，便于在固定温压包线与分位轴上生成状态节点、有向边及 PS 意义下的子循环单元。

---

## 运行环境

| 依赖 | 说明 |
|------|------|
| Python | 3.12+（当前验证版本） |
| [CoolProp](https://www.coolprop.org/) | 物性计算（`import CoolProp.CoolProp`） |
| pytest | 单元测试 |
| matplotlib | 仅拓扑可视化测试需要 |

安装示例：

```bash
pip install CoolProp matplotlib pytest
```

从仓库根目录运行测试时，请将项目根加入 `PYTHONPATH`，否则 `import core` 会失败：

```bash
cd /path/to/CyGES
set PYTHONPATH=.          # Windows CMD
# 或
$env:PYTHONPATH="."      # PowerShell
python -m pytest tests/test_tp_topology.py -q
```

仅运行氦气拓扑综合图用例：

```bash
python -m pytest tests/test_tp_topology.py::test_helium_topology_overview_plot -v
```

成功后会生成图像文件 [`tests/ts_topology_he.png`](tests/ts_topology_he.png)（双子图：T–S 与 P–S，含节点、机械/换热边、最小子循环多边形及编号）。

---

## 仓库结构（主干）

| 路径 | 作用 |
|------|------|
| [`core/closed_cycle_layer.py`](core/closed_cycle_layer.py) | 闭式循环层：`ClosedCycleTPInput`、`Node`、`Edge`、`SubCycle`、`ClosedCycleLayer`；公开函数 `build_axis`、`build_node_edge_topology`、`build_subcycles` |
| [`core/fluid_property_solver.py`](core/fluid_property_solver.py) | 物性协议 `FluidPropertySolver` 与 CoolProp 实现 `CoolPropFluidPropertySolver`（`state(pair, x, y)` → `T,P,H,S`） |
| [`tests/test_tp_topology.py`](tests/test_tp_topology.py) | 氦气拓扑单测与上述 PNG 输出 |
| [`core/__init__.py`](core/__init__.py) | 包说明 |

仓库中另含 `oldFile/`、`solvers/` 等历史或扩展目录；**当前 README 描述以 `core/closed_cycle_layer.py` 与 `core/fluid_property_solver.py` 为准**。

---

## 单位约定

字段名不缀单位，全项目一致：

| 量 | 单位 | 典型出现位置 |
|----|------|----------------|
| T | K | `Node.T`、`ClosedCycleTPInput.t_min` / `t_max` |
| P | kPa | `Node.P`、`p_min` / `p_max` |
| H | kJ/kg | `Node.H`、物性字典 `"H"` |
| S | kJ/(kg·K) | `Node.S`、物性字典 `"S"` |

CoolProp 内部为 SI；[`CoolPropFluidPropertySolver`](core/fluid_property_solver.py) 在边界做 kPa、kJ/kg 等与上述约定对齐。

---

## 物性 API

[`CoolPropFluidPropertySolver.state(pair, x, y)`](core/fluid_property_solver.py) 返回 **`ThermoStateTPHS`** 字典，键为 **`"T"`, `"P"`, `"H"`, `"S"`**。

| `pair` | `x` | `y` |
|--------|-----|-----|
| `"TP"` | T [K] | P [kPa] |
| `"PS"` | P [kPa] | S [kJ/(kg·K)] |
| `"HP"` | H [kJ/kg] | P [kPa] |
| `"HS"` | H [kJ/kg] | S [kJ/(kg·K)] |

实现类长期复用 `AbstractState("HEOS", fluid)`，适合在拓扑生成中高频调用。

---

## 闭式循环层拓扑

### 入口与流程

1. 构造 [`ClosedCycleTPInput`](core/closed_cycle_layer.py)（工质名、温压上下限、`t_quantiles` / `p_quantiles` 为 \[0,1\] 分位且互不重复、可选 `max_mass_flow`）。
2. 构造 [`ClosedCycleLayer(inp, properties=...)`](core/closed_cycle_layer.py)；默认物性后端为 `CoolPropFluidPropertySolver`。
3. 调用 **`analyze_topology()`**：内部依次执行 **`build_node_edge_topology`**（一级 TP 网格 → 等熵二级与机械边 → 换热边 → 将边键写回节点）与 **`build_subcycles`**（最小 4 节点 4 边子循环枚举）。

### PS 平面与有向边

在 PS 平面上约定：**P** 由小到大为自下而上，**S** 由小到大为自左而右。拓扑边上的方向为 **`tail → head`**，且满足 **`P`、`S` 在容差意义下不减**（机械边沿等熵链压力升高；换热边在等压桶内沿熵增大方向）。

### `ClosedCycleLayer` 上的结果字段

- **`nodes: dict[int, Node]`**：键为全局 **`index`**。`parent is None` 为一级 TP 网格点；`parent` 为某一级 `index` 时为二级等熵延伸点。节点上 **`edge_up` / `edge_down` / `edge_left` / `edge_right`** 存邻边在 **`edges`** 中的键（无邻边为 `None`）：机械边对应上/下，换热边对应右/左（与上述 PS 约定一致）。
- **`edges: dict[str, Edge]`**：机械边键 **`M1`, `M2`, …**；换热边键 **`H1`, `H2`, …**。`Edge.kind` 为 `"mechanical"` 或 `"heat"`；`Edge.mass_flow` 快照阶段为 `None`，供后续业务填写。
- **`subcycles: list[SubCycle]`**：最小子循环。每个 [`SubCycle`](core/closed_cycle_layer.py) 含 **`nodes`**（左下、左上、右上、右下，顺时针）、**`edges`**（左、上、右、下，与 `edges` 字典键对应）、**`mass_flow`**（子循环标量质量流，初值 `None`，可为负，类非 `frozen` 可就地修改）。枚举规则为：从某节点沿 `edge_up` → `edge_right` → `edge_down`（取该边 `tail`）→ `edge_left`（取该边 `tail`）回到起点，四角互异，并用 `frozenset` 去重；**不对 P/S 几何形状做额外校验**。

若需单独复现拓扑而不通过 `ClosedCycleLayer`，可直接调用模块级函数 **`build_node_edge_topology(solver, inp)`** 与 **`build_subcycles(nodes, edges)`**；轴采样使用 **`build_axis(min, max, quantiles)`**。

---

## 开发与扩展

- 修改网格或边生成逻辑：主要编辑 [`core/closed_cycle_layer.py`](core/closed_cycle_layer.py) 中的 `build_node_edge_topology`；子循环策略编辑 `build_subcycles`。
- 替换物性后端：实现与 **`FluidPropertySolver`** 协议兼容的类型，并在 **`ClosedCycleLayer(..., properties=...)`** 注入。
- 若源码与本文不一致，以仓库内实际文件为准。
