# CyGES

闭式循环能源系统原型：在温压包线与分位轴上生成 **PS 平面离散拓扑**、**最小子循环**与**精简拓扑**，并支持单步**非理想节点偏置**（换热 `σ`、机械 `η_is`）。物性由 **CoolProp** 提供。

算法与实现细节见 [`docs/architecture.md`](docs/architecture.md)；接续开发见 [`AGENTS.md`](AGENTS.md)。

---

## 环境

| 依赖 | 用途 |
|------|------|
| Python 3.12+ | 当前验证环境 |
| CoolProp | 物性（必需） |
| pytest | 测试 |
| matplotlib | 绘图测试（可选） |

```bash
pip install CoolProp matplotlib pytest
```

从仓库根目录运行代码或测试时，须将项目根加入 `PYTHONPATH`：

```bash
cd /path/to/CyGES
set PYTHONPATH=.          # Windows CMD
# 或
$env:PYTHONPATH="."       # PowerShell
python -m pytest tests/ -q
```

---

## 单位

| 量 | 单位 |
|----|------|
| T | K |
| P | kPa |
| H | kJ/kg |
| S | kJ/(kg·K) |

字段名不缀单位；CoolProp 内部为 SI，[`CoolPropFluidPropertySolver`](core/fluid_property_solver.py) 在边界对齐上述约定。

---

## 快速开始

```python
from core import ClosedCycleLayer, ClosedCycleTPInput
from core.non_ideal_bias import apply_combined_offsets

inp = ClosedCycleTPInput(
    fluid="He",
    t_min=100.0,
    t_max=900.0,
    p_min=1000.0,
    p_max=9000.0,
    t_quantiles=(0.3, 0.7),
    p_quantiles=(0.3, 0.7),
    max_mass_flow=10.0,
)
layer = ClosedCycleLayer(inp)                    # 默认自动 analyze_topology()
layer.subcycle_mass_flows = [...]                # 优化器主控向量
layer.commit_subcycle_mass_flows_to_topology()   # 量化 → 写回 → 汇聚 → 重建 simplified

non_ideal = layer.ensure_non_ideal()
apply_combined_offsets(non_ideal)                # 单步：σ 改 P → 机械组 PS 重置 → η_is 步进

report = layer.performance_report()            # 精简边过程归类与循环汇总（非理想已偏置时用非理想状态）
```

- `auto_analyze=False` 时须自行调用 `layer.analyze_topology()`。
- `ClosedCycleLayer` / `NonIdealClosedCycleLayer` 字段见类 docstring。
- 非理想偏置步骤与 anchor 规则见 [`docs/architecture.md` §8](docs/architecture.md)。

---

## 输入：`ClosedCycleTPInput`

| 字段 | 含义 |
|------|------|
| `fluid` | 工质（CoolProp 流体名） |
| `t_min` / `t_max` | 温度范围 [K] |
| `p_min` / `p_max` | 压力范围 [kPa] |
| `t_quantiles` / `p_quantiles` | `[0,1]` 分位，在范围内线性插值 |
| `max_mass_flow` | 质量流上界（可选）；用于子循环初值与量化步长 |
| `subcycle_mass_flow_step_fraction` | 量化步长占 `max_mass_flow` 的比例；默认见 `config` |

分位须与端点互不重复；本层不对输入做防御性校验。

---

## 可调常量（[`config.py`](config.py)）

| 常量 | 默认 | 含义 |
|------|------|------|
| `SUBCYCLE_INITIAL_MASS_FLOW_FRACTION_OF_MAX` | `0.1` | 子循环初值 = 系数 × `max_mass_flow` |
| `SUBCYCLE_MASS_FLOW_STEP_FRACTION_DEFAULT` | `0.01` | 流量量化步长比例 |
| `NON_IDEAL_MECHANICAL_EFFICIENCY_DEFAULT` | `0.9` | 机械等熵效率 `η_is` |
| `NON_IDEAL_HEAT_EFFICIENCY_DEFAULT` | `0.95` | 换热总压恢复系数 `σ` |

---

## 仓库结构

| 路径 | 说明 |
|------|------|
| [`core/closed_cycle_layer.py`](core/closed_cycle_layer.py) | 理想层：拓扑、子循环、精简、`ClosedCycleLayer` |
| [`core/non_ideal_bias.py`](core/non_ideal_bias.py) | 非理想偏置：有向组、深度、`apply_combined_offsets` |
| [`core/cycle_performance.py`](core/cycle_performance.py) | 性能统计：过程归类、`compute_cycle_performance` |
| [`core/fluid_property_solver.py`](core/fluid_property_solver.py) | 物性求解器 |
| [`config.py`](config.py) | 全局常量 |
| [`docs/architecture.md`](docs/architecture.md) | 算法与实现（单一信息源） |
| [`tests/`](tests/) | 理想 / 非理想两工况绘图测试 |
| `inputs/`、`solvers/`、`optimize/`、`oldFile/` | 历史或占位，**非当前 API** |

---

## 测试

需 CoolProp；绘图用例另需 matplotlib。生成的 PNG **不纳入 git**。

| 用例 | 命令 | 输出 |
|------|------|------|
| 理想 He 拓扑 | `pytest tests/test_tp_topology.py -q` | `tests/ts_topology_he.png` |
| 非理想 Case A / B | `pytest tests/test_non_ideal_two_cases.py -q` | `tests/non_ideal_offsets_case_a_he_wide.png`、`tests/non_ideal_offsets_case_b_he_mid_fine.png` |

全量：`pytest tests/ -q`（共 7 个用例）。也可直接运行：

```bash
python tests/test_non_ideal_two_cases.py
```

---

## 扩展入口

- 改网格与边：[`build_node_edge_topology`](core/closed_cycle_layer.py)
- 改默认 `σ` / `η_is`：[`config.py`](config.py)
- 换物性后端：实现 [`FluidPropertySolver`](core/fluid_property_solver.py)，构造时 `ClosedCycleLayer(..., properties=...)`

若本文与源码不一致，以仓库内实际文件为准。
