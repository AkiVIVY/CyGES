# CyGES

闭式循环能源系统原型代码。在固定温压包线与分位轴上构建 **PS 离散拓扑**（一级 TP 网格、沿压力轴的等熵二级节点、机械边 `M*` 与换热边 `H*`）、枚举**最小子循环**（4 节点 4 边），并支持**子循环质量流向量**与精简拓扑。物性由 **CoolProp** 经统一 `state` 接口提供。

非理想修正中，机械边对应**等熵效率** `η_is`，换热边对应**总压恢复系数** `σ`。算法细节统一见 [`docs/architecture.md`](docs/architecture.md)。

---

## 运行环境与依赖

| 依赖 | 说明 |
|------|------|
| Python | 3.12+（当前验证环境） |
| [CoolProp](https://www.coolprop.org/) | `import CoolProp.CoolProp`，物性计算必需 |
| pytest | 跑测试用例 |
| matplotlib | 仅 `tests/` 下绘图用例需要 |

```bash
pip install CoolProp matplotlib pytest
```

**重要**：从仓库根运行代码或测试时须将项目根加入 `PYTHONPATH`，否则 `import core` 与 `import config` 会失败。

```bash
cd /path/to/CyGES
set PYTHONPATH=.          # Windows CMD
# 或
$env:PYTHONPATH="."       # PowerShell
python -m pytest tests/ -q
```

---

## 单位约定

字段名不缀单位。

| 量 | 单位 | 典型位置 |
|----|------|----------|
| T | K | `Node.T`、`t_min` / `t_max` |
| P | kPa | `Node.P`、`p_min` / `p_max` |
| H | kJ/kg | `Node.H`、物性字典 `"H"` |
| S | kJ/(kg·K) | `Node.S`、物性字典 `"S"` |

CoolProp 内部为 SI，[`CoolPropFluidPropertySolver`](core/fluid_property_solver.py) 在边界与上述约定对齐。

---

## 物性 API

[`CoolPropFluidPropertySolver.state(pair, x, y)`](core/fluid_property_solver.py) 返回 `ThermoStateTPHS`（`T,P,H,S`）。

| `pair` | `x` | `y` |
|--------|-----|-----|
| `"TP"` | T [K] | P [kPa] |
| `"PS"` | P [kPa] | S [kJ/(kg·K)] |
| `"HP"` | H [kJ/kg] | P [kPa] |
| `"HS"` | H [kJ/kg] | S [kJ/(kg·K)] |

实现复用 `AbstractState("HEOS", fluid)`，适合拓扑生成中的高频调用。

---

## 输入：`ClosedCycleTPInput`

| 字段 | 含义 |
|------|------|
| `fluid` | 工质标识（与 CoolProp 流体名一致） |
| `t_min` / `t_max` | 温度上下限 [K] |
| `p_min` / `p_max` | 压力上下限 [kPa] |
| `t_quantiles` / `p_quantiles` | `[0, 1]` 分位，在对应范围内线性插值；调用方保证与端点互不重复 |
| `max_mass_flow` | 全层质量流上界（可选），用于子循环初值与量化步长基准 |
| `subcycle_mass_flow_step_fraction` | 量化步长占 `max_mass_flow` 的比例；未传时取 `config.SUBCYCLE_MASS_FLOW_STEP_FRACTION_DEFAULT` |

**输入校验在调用方**：本层不对 `ClosedCycleTPInput` 做防御性检查，仅在算法不变量被破坏时抛错。

---

## 推荐调用顺序

```python
from core import ClosedCycleLayer, ClosedCycleTPInput

inp = ClosedCycleTPInput(
    fluid="Helium",
    t_min=..., t_max=..., p_min=..., p_max=...,
    t_quantiles=(...,), p_quantiles=(...,),
    max_mass_flow=...,
)
layer = ClosedCycleLayer(inp)              # 构造时自动 analyze_topology()
layer.subcycle_mass_flows = [...]          # 优化器主控向量；可任意浮点
layer.commit_subcycle_mass_flows_to_topology()   # 量化 → 同步 → 汇聚 → 重建 simplified
non_ideal = layer.ensure_non_ideal()       # 创建非理想快照层
non_ideal.apply_heat_pressure_offsets()    # σ：P ← P × σ^layer + 全表 PS 闭合
non_ideal.apply_mechanical_isentropic_offsets()  # η_is：PS → η → HP（须在换热之后）
```

`auto_analyze=False` 时须随后自行 `layer.analyze_topology()`。

`ClosedCycleLayer` / `NonIdealClosedCycleLayer` 的字段定义见对应模块的类 docstring；算法细节（PS 平面、链合并、有向组深度、机械锚点、机械步公式、换热 σ）统一见 [`docs/architecture.md`](docs/architecture.md)。

---

## [`config.py`](config.py) 可调常量

| 常量 | 默认值 | 含义 |
|------|--------|------|
| `SUBCYCLE_INITIAL_MASS_FLOW_FRACTION_OF_MAX` | `0.1` | `analyze_topology` 后每个子循环初值 = 该系数 × `max_mass_flow`；`max_mass_flow=None` 时按 `0` |
| `SUBCYCLE_MASS_FLOW_STEP_FRACTION_DEFAULT` | `0.01` | `ClosedCycleTPInput` 未显式提供 `subcycle_mass_flow_step_fraction` 时采用 |
| `NON_IDEAL_MECHANICAL_EFFICIENCY_DEFAULT` | `0.85` | 机械精简边统一**等熵效率** `η_is` |
| `NON_IDEAL_HEAT_EFFICIENCY_DEFAULT` | `0.99` | 换热精简边统一**总压恢复系数** `σ` |

---

## 仓库结构

| 路径 | 说明 |
|------|------|
| [`core/closed_cycle_layer.py`](core/closed_cycle_layer.py) | 理想层：`ClosedCycleLayer`、数据模型、拓扑构建、子循环、精简 |
| [`core/non_ideal_closed_cycle_layer.py`](core/non_ideal_closed_cycle_layer.py) | 非理想层：`NonIdealClosedCycleLayer`、有向组、`σ` / `η_is` 偏移 |
| [`core/fluid_property_solver.py`](core/fluid_property_solver.py) | `FluidPropertySolver` 协议与 CoolProp 实现 |
| [`config.py`](config.py) | 全局可调常量 |
| [`docs/architecture.md`](docs/architecture.md) | 算法细节单一信息源 |
| [`tests/`](tests/) | 绘图用例 |
| [`AGENTS.md`](AGENTS.md) | 接续 agent / 开发者交接说明 |
| `inputs/` / `solvers/` / `optimize/` / `oldFile/` | 历史或占位，**非当前 API 依据** |

---

## 测试

`tests/` 当前仅保留**绘图**用例（需 CoolProp + matplotlib）。生成的 PNG 不纳入 git。

| 用例 | 输出 |
|------|------|
| [`tests/test_tp_topology.py`](tests/test_tp_topology.py)::`test_helium_topology_overview_plot` | `ts_topology_he.png`（理想 T–S / P–S） |
| [`tests/test_non_ideal_offsets_plot.py`](tests/test_non_ideal_offsets_plot.py)::`test_helium_ideal_and_non_ideal_offsets_topology_plot` | `non_ideal_offsets_he.png`（理想 vs 偏移后 2×2） |
| [`tests/test_non_ideal_offsets_plot.py`](tests/test_non_ideal_offsets_plot.py)::`test_helium_focus_nodes_9_14_35_44_offsets_plot` | `non_ideal_offsets_nodes_9_14_35_44.png`（局部节点） |

---

## 开发与扩展

- 改网格/边：[`build_node_edge_topology`](core/closed_cycle_layer.py)；子循环枚举：[`build_subcycles`](core/closed_cycle_layer.py)。
- 改全局默认系数或步长比例：[`config.py`](config.py)。
- 换物性后端：实现 [`FluidPropertySolver`](core/fluid_property_solver.py) 协议，`ClosedCycleLayer(..., properties=...)` 注入。

若本文与源码不一致，以仓库内实际文件为准。
