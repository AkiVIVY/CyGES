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
| scipy | L-BFGS-B 内层优化（可选） |

```bash
pip install CoolProp matplotlib pytest scipy
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
| `s_quantiles` | S 分位节点位置（`[0,1]`）— 在 primary nodes 的 S 范围内创建 PS 延伸节点 |
| `subcycle_mass_flow_initial` | 子循环初始质量流量 [kg/s]；`None` 则取 `config.SUBCYCLE_INITIAL_MASS_FLOW_DEFAULT` |

分位须与端点互不重复；本层不对输入做防御性校验。

---

## 外部源：`ExternalSourceInput`

| 字段 | 含义 |
|------|------|
| `fluid` | 冷/热源工质（任意 CoolProp 流体名，如 `"Air"`, `"Hydrogen"`, `"Methane"` 等） |
| `mass_flow` | 质量流量 [kg/s] |
| `T_in` / `T_out` | 进出口温度 [K] |
| `P_in` / `P_out` | 进出口压力 [kPa] |

---

## 可调常量（[`config.py`](config.py)）

| 常量 | 默认 | 含义 |
|------|------|------|
| `SUBCYCLE_INITIAL_MASS_FLOW_DEFAULT` | `0.0` | 子循环初始流量默认值 [kg/s] |
| `NON_IDEAL_MECHANICAL_EFFICIENCY_DEFAULT` | `0.9` | 机械等熵效率 `η_is` |
| `NON_IDEAL_HEAT_EFFICIENCY_DEFAULT` | `0.98` | 换热总压恢复系数 `σ` |

---

## 仓库结构

| 路径 | 说明 |
|------|------|
| [`core/closed_cycle_layer.py`](core/closed_cycle_layer.py) | 理想层：拓扑、子循环、精简、`ClosedCycleLayer` |
| [`core/non_ideal_bias.py`](core/non_ideal_bias.py) | 非理想偏置：有向组、深度、`NonIdealClosedCycleLayer.apply_offsets()` |
| [`core/cycle_performance.py`](core/cycle_performance.py) | 性能统计：过程归类、循环汇总（纯计算，零外部依赖） |
| [`core/postprocess.py`](core/postprocess.py) | 二次处理：T-Q 曲线构建、夹点分析（`analyze_pinch`） |
| [`core/heat_exchanger.py`](core/heat_exchanger.py) | HX 匹配：星型分组 + 串联夹点（`match_series_pinch`） |
| [`core/system.py`](core/system.py) | 系统层：冷热源整合、`SystemPipeline.run(props)`（四种夹点模式） |
| [`core/fluid_property_solver.py`](core/fluid_property_solver.py) | 物性求解器 + `PropertyRegistry`（多工质注册表，含 memo 缓存） |
| [`optimize/`](optimize/) | 优化器：CMA-ES/DE、高斯基编码、目标注册表、多进程并行 |
| [`config.py`](config.py) | 全局常量 |
| [`docs/architecture.md`](docs/architecture.md) | 算法与实现（单一信息源） |
| [`tests/test_system_full.py`](tests/test_system_full.py) | 全管线可视化测试（5 图） |
| [`tests/test_optimize.py`](tests/test_optimize.py) | 优化测试：CMA-ES + 高斯基编码 + Excel 日志 |
| [`tests/test_layered_opt.py`](tests/test_layered_opt.py) | 分层优化：LHS+DE 外层 + CMA/L-BFGS-B 内层，多模式 |
| [`tests/test_inner_compare.py`](tests/test_inner_compare.py) | CMA-ES vs L-BFGS-B 内层方法对比 |
| [`tests/optimization_example.py`](tests/optimization_example.py) | CH₄ 冷源效率优化示例 |

---

## 测试

```bash
# 全管线测试（5 图）
pytest tests/test_system_full.py -q

# 分层优化测试
pytest -s tests/test_layered_opt.py::test_layered_0p0t0s_series

# CMA-ES 优化测试
pytest -s tests/test_optimize.py

# 内层方法对比（CMA vs L-BFGS-B）
pytest -s tests/test_inner_compare.py::test_compare_inner

# CH₄ 冷源效率优化
pytest -s tests/optimization_example.py::test_h2_scan_efficiency
```

输出图表位于 `tests/run_layered_fast/`、`tests/inner_compare/`、`tests/optimization_example/`。

---

## 扩展入口

- 改网格与边：[`build_node_edge_topology`](core/closed_cycle_layer.py)
- 改默认 `σ` / `η_is`：[`config.py`](config.py)
- 换物性后端：实现 [`FluidPropertySolver`](core/fluid_property_solver.py)，构造时 `ClosedCycleLayer(properties=...)`
- 多工质场景：[`PropertyRegistry`](core/fluid_property_solver.py)，`SystemPipeline.run(props)` 传入；冷源工质通过 `cold_fluid` 参数切换
- 优化器配置：[`_DEFAULT_HP`](tests/test_layered_opt.py) 集中所有可调参数（内层方法、L-BFGS-B 多起点、冷源工质等）
- 性能加速：`PropertyRegistry` 自动 memo 缓存；`InterpolatingHeliumSolver` 可通过 `use_interp_he=True` 启用

若本文与源码不一致，以仓库内实际文件为准。
