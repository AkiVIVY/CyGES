# CyGES

闭式循环能源系统原型：在温压包线与分位轴上生成 **PS 平面离散拓扑**、**最小子循环**与**精简拓扑**，并支持单步**非理想节点偏置**（换热 `σ`、机械 `η_is`）、**T-Q 换热曲线**与**夹点分析**。物性由 **CoolProp** 提供。

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

inp = ClosedCycleTPInput(
    fluid="He",
    t_min=100.0,
    t_max=900.0,
    p_min=1000.0,
    p_max=9000.0,
    t_quantiles=(0.3, 0.7),
    p_quantiles=(0.3, 0.7),
    mass_flow_min=-5.0,
    mass_flow_max=10.0,
)
layer = ClosedCycleLayer(inp)                    # 默认自动 analyze_topology()
layer.subcycle_mass_flows = [...]                # 优化器主控向量
layer.commit_subcycle_mass_flows_to_topology()   # 量化 → 写回 → 汇聚 → 重建 simplified

non_ideal = layer.ensure_non_ideal()
non_ideal.apply_offsets()                       # 单步：σ 改 P → 机械组 PS 重置 → η_is 步进

report = layer.performance_report()            # 精简边过程归类与循环汇总（非理想已偏置时用非理想状态）

# ── T-Q 曲线与夹点分析 ──
from core import PropertyRegistry, build_heat_tq_curves, analyze_pinch, ProcessCategory
props = PropertyRegistry()
tq = build_heat_tq_curves(report, props)              # 从报告构建吸热/放热 T-Q 曲线
result = analyze_pinch(                                # 夹点分析：输入 ProcessRecord 列表
    report.by_category[ProcessCategory.HEAT_ABSORPTION],
    report.by_category[ProcessCategory.HEAT_REJECTION],
    delta_T_min=10.0, props=props,
)
print(result.hot_utility_demand, result.cold_utility_demand)  # 公用工程需求 [kW]
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
| `mass_flow_min` / `mass_flow_max` | 子循环质量流量边界 [kg/s]；默认见 `config.MASS_FLOW_MIN_DEFAULT` / `MASS_FLOW_MAX_DEFAULT` |
| `subcycle_mass_flow_step_fraction` | 量化步长占 `mass_flow_max - mass_flow_min` 的比例；默认见 `config` |
| `subcycle_mass_flow_initial` | 子循环初始质量流量 [kg/s]；`None` 则用 config 公式计算 |

分位须与端点互不重复；本层不对输入做防御性校验。

---

## 可调常量（[`config.py`](config.py)）

| 常量 | 默认 | 含义 |
|------|------|------|
| `SUBCYCLE_INITIAL_MASS_FLOW_FRACTION_OF_MAX` | `0.1` | 子循环初值 = `mass_flow_min + 系数 × (mass_flow_max - mass_flow_min)` |
| `SUBCYCLE_MASS_FLOW_STEP_FRACTION_DEFAULT` | `0.01` | 量化步长占流量边界的比例 |
| `MASS_FLOW_MIN_DEFAULT` | `-5.0` | 子循环质量流量下限默认值 [kg/s] |
| `MASS_FLOW_MAX_DEFAULT` | `10.0` | 子循环质量流量上限默认值 [kg/s] |
| `NON_IDEAL_MECHANICAL_EFFICIENCY_DEFAULT` | `0.9` | 机械等熵效率 `η_is` |
| `NON_IDEAL_HEAT_EFFICIENCY_DEFAULT` | `0.95` | 换热总压恢复系数 `σ` |
| `PINCH_EXTRA_CURVE_FRACTION_THRESHOLD` | `0.01` | 额外换热曲线截断阈值（低于此比例不输出） |

---

## 仓库结构

| 路径 | 说明 |
|------|------|
| [`core/closed_cycle_layer.py`](core/closed_cycle_layer.py) | 理想层：拓扑、子循环、精简、`ClosedCycleLayer` |
| [`core/non_ideal_bias.py`](core/non_ideal_bias.py) | 非理想偏置：有向组、深度、`NonIdealClosedCycleLayer.apply_offsets()` |
| [`core/cycle_performance.py`](core/cycle_performance.py) | 性能统计：过程归类、循环汇总（纯计算，零外部依赖） |
| [`core/postprocess.py`](core/postprocess.py) | 二次处理：T-Q 曲线构建、夹点分析（`analyze_pinch`） |
| [`core/system.py`](core/system.py) | 系统层：冷热源整合、`SystemPipeline.run(props)`（三种夹点模式） |
| [`core/fluid_property_solver.py`](core/fluid_property_solver.py) | 物性求解器 + `PropertyRegistry`（多工质注册表） |
| [`config.py`](config.py) | 全局常量 |
| [`docs/architecture.md`](docs/architecture.md) | 算法与实现（单一信息源） |
| [`tests/test_system_full.py`](tests/test_system_full.py) | 全管线可视化测试（5 图） |

---

## 测试

```bash
pytest tests/test_system_full.py -q
```

输出 5 张 PNG：冷热源 TS/PS、循环 TS/PS、循环夹点、循环性能、系统夹点。

---

## 扩展入口

- 改网格与边：[`build_node_edge_topology`](core/closed_cycle_layer.py)
- 改默认 `σ` / `η_is`：[`config.py`](config.py)
- 换物性后端：实现 [`FluidPropertySolver`](core/fluid_property_solver.py)，构造时 `ClosedCycleLayer(properties=...)`
- 多工质场景：[`PropertyRegistry`](core/fluid_property_solver.py)，`SystemPipeline.run(props)` 传入

若本文与源码不一致，以仓库内实际文件为准。
