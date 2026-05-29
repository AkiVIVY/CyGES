# CyGES — Agent 交接说明

面向**接续开发的 Agent / 开发者**。用户向说明见 [`README.md`](README.md)；算法与实现细节**只在** [`docs/architecture.md`](docs/architecture.md) 维护一份，本文件与代码 docstring **仅引用、不重复公式**。

---

## 1. 文档分工

| 文件 | 读者 | 写什么 |
|------|------|--------|
| [`README.md`](README.md) | 用户 | 环境、快速开始、config、测试命令 |
| **本文件** | Agent | 范围、不变量、模块地图、架构决策、协作约定 |
| [`docs/architecture.md`](docs/architecture.md) | 实现者 | PS 约定、流水线、精简、深度、偏置公式 |

改算法主干 → 先改 `architecture.md`，再改代码；README / AGENTS 仅在结构性变更时动。

---

## 2. 模块地图

```
CyGES/
├── core/
│   ├── closed_cycle_layer.py    # 理想层：拓扑构建（含S分位）、子循环枚举、简化
│   ├── non_ideal_bias.py        # 非理想偏置：有向组、深度、ni.apply_offsets()
│   ├── cycle_performance.py     # 性能统计：ProcessRecord, CyclePerformanceReport (纯计算)
│   ├── postprocess.py           # T-Q 构建 + 夹点分析 (build_tq / analyze_pinch)
│   ├── heat_exchanger.py        # HX 匹配：星型拓扑 + 串联夹点(match_series_pinch)
│   ├── system.py                # 系统层：Pipeline + analyze_system_heat (解耦)
│   └── fluid_property_solver.py # 物性：PropertyRegistry + InterpolatingHeliumSolver
├── optimize/
│   ├── solver.py                # Optimizer: CMA-ES(BIPOP)+DE, 高斯基编码, 并行
│   ├── objective.py             # @register 目标函数注册表（含 hx_unmatched）
│   └── types.py                 # OptimizationResult
├── config.py
├── tests/
│   ├── test_system_full.py     # 全管线可视化 5 图
│   ├── test_optimize.py        # 优化测试 + Excel 日志 + run_XXX 输出
│   ├── test_optimize_hx.py     # HX 匹配优化: 4 图 + dT_min 参数化
│   ├── test_hx_match.py        # HX 匹配单元测试
│   └── test_layered_opt.py     # 分层优化: LHS+DE外层 + CMA内层, TS/PS/HX图
├── docs/architecture.md
│   └── core_bugs.md            # 已知问题备忘
└── AGENTS.md                   # 本文件
```

---

## 3. 当前实现范围

**已实现**

- 理想层：PS 离散拓扑（含 S 分位节点）、最小子循环、简化边。
- S 分位节点：`ClosedCycleTPInput.s_quantiles` — 在 primary nodes 的 S_min/S_max 范围内按分位值创建额外 PS 延伸节点（`parent="S"`），扩展子循环拓扑。
- 非理想偏置：`ni.apply_offsets()` 单步 σ 改 P → PS 重置 → DFS 机械步。
- 性能统计：`ProcessRecord`/`CyclePerformanceReport`（纯计算，零外部依赖）。
- T-Q 构建 + 夹点分析：`build_heat_tq_curves`、`compute_pinch`（底层）、`analyze_pinch`（高层）。
- 系统集成：`SystemPipeline.run(props, layers=...)` 产纯数据，支持预建 layer 复用拓扑；`analyze_system_heat(raw, inp, props)` 独立做夹点。
- 多工质：`PropertyRegistry`（自动缓存 CoolProp 求解器）。
- 物性加速：`InterpolatingHeliumSolver`（通过 `use_interp_he` 开关启用，默认关闭）。
- **HX 串联夹点匹配**：`match_series_pinch()` — 所有放热/吸热过程按温度降序在 Q 上串联为单一逆流单元，在每条过程边界检查 dT_min 违规。与星型匹配并行共存，通过 `hx_series=True` 开关切换。
  - 目标函数：`obj = unmatched_ratio + 1e-2 * pinch_violation`（串联模式）
  - 目标函数：`obj = unmatched_ratio + 1e-3 * num_unmatched`（星型模式）
- 星型 HX 匹配：`match_heat_exchanger_groups()` / `match_constructive()` / `match_heat_exchanger_staged()`。
- 分层优化 `test_layered_opt.py`：外层 LHS+DE 搜索拓扑参数（含 s_q），内层 CMA 复用固定拓扑优化流量 + H2 + h2_T_out。支持 `h2_T_out` 固定、`h2_mf` 固定等模式。`_DEFAULT_HP` 集中所有可调参数。
- 可视化：理想/非理想 TS/PS（子循环多边形叠加）、HX T-Q（概览 + 每单元逆流）、DE 收敛曲线。

**HX 匹配模式**

| 模式 | 调用方式 | 说明 |
|------|------|------|
| 星型分组 | `hx_series=False` (默认) | 1H+kC / mH+1C，枚举或构造法 |
| 串联夹点 | `hx_series=True` | 全局串联为单一逆流单元，边界夹点检查 |

**未实现**

- 非理想方程/约束装配、多目标优化器。
- HEN / 多热源多冷源边界耦合。
- 按边/按组分别赋值 σ、η_is。

---

## 4. 关键设计决策

**夹点分析与管线解耦**：`SystemPipeline.run()` 产纯数据，`analyze_system_heat()` 独立做夹点。

**高斯基编码**：子循环流量由 PS 空间上的高斯核加权和生成，维度固定，与子循环数无关。

**离散化在优化器层**：流量/分位/T 边界步长全在 `Optimizer` 中控制。

**S 分位节点**：`ClosedCycleTPInput.s_quantiles` 与 t_quantiles/p_quantiles 平级。S 分位节点由 primary nodes 的 S_min/S_max 线性插值得到，对每个 P 轴 做 PS 延伸（`solver.state("PS", P, S_q)`），`parent="S"`。新节点通过同 parent 的 P 排序链形成机械边，并参与换热边构建。S 分位点的离散化步长与 T/P 分位共用 `qstep`。

**串联夹点匹配**：所有过程按 T_high 降序排列，在累计 Q 上构成两条链（hots + colds 均降序）。单一单元检查每个过程边界的逆流温差。目标函数 = `unmatched_ratio + 1e-2*pinch_violation`，直接暴露系统级的功率不平衡与温差违规。

**分层优化**：外层 LHS/DE 搜索拓扑参数（t_min/t_max/p_min/p_max/t_q/p_q/s_q），内层 CMA 在固定拓扑下优化流量 + H2 + h2_T_out。支持 `h2_T_fixed`（`h2_T_out_lo==h2_T_out_hi`）和 `h2_fixed`（`h2_mf_lo==h2_mf_hi`）模式，自动调整 CMA 维度。`LayerResult` 数据类含 `s_q` 字段。

**物性加速**：`InterpolatingHeliumSolver` 默认关闭（`_DEFAULT_HP.use_interp_he=False`），所有测试用纯 CoolProp。

**n_sc==0 短路**：`_inner_cma_fast` 中 n_sc==0 时直接返回 obj=1.0，跳过 CMA 迭代。

**熵序警告**：`build_node_edge_topology` 中 S(t_max,p_max) < S(t_min,p_min) 时发出 RuntimeWarning。

---

## 5. 测试

| 文件 | 用途 | 运行 |
|------|------|------|
| `test_layered_opt.py` | 分层优化，多个预设配置 | `pytest -s tests/test_layered_opt.py::<test_name>` |
| `test_optimize.py` | CMA-ES + DE 优化测试 | `pytest -s tests/test_optimize.py` |

可选测试函数（`test_layered_opt.py`）：

| 函数 | Topo | HX | H2 | 说明 |
|------|------|----|----|------|
| `test_layered_0p0t_hx1` | 1P+0T+1S | series | [2,5.5] | 串联夹点 |
| `test_layered_0p0t0s_series` | 0P0T0S | series | [2,5.5] | 最简串联 |
| `test_layered_0p0t0s_h2_4` | 0P0T0S | series | 4固定 | H2=4 串联 |
| `test_layered_0p1t0s_h2_4` | 0T+1P+0S | series | 4固定 | +P分位 |
| `test_layered_0p0t_sq1` | 0P0T+1S | group HX2 | [3,6] | S分位测试 |
| `test_layered_1p1t_ideal` | 1T+1P | group HX2 | [3,6] | 1P1T 基准 |

输出：`tests/run_layered_fast/` 下图表。

---

## 6. 失效与快照不变量

| 触发 | 行为 |
|------|------|
| `ClosedCycleLayer.analyze_topology()` | 重建 topology + `simplified`，清空 `non_ideal` |
| `ClosedCycleLayer.commit_subcycle_mass_flows_to_topology()` | 同步子循环流量 → 汇聚到边 → 重建 `simplified`，清空 `non_ideal` |
| `ClosedCycleLayer._rebuild_simplified()` | 清空 `non_ideal` + 重建 `simplified`（子循环为空或全零流量时 warn） |
| `ClosedCycleLayer.ensure_non_ideal()` | 创建 `NonIdealClosedCycleLayer` 快照；要求 `simplified` 存在 |
| 父层变动后 | 已存在 `NonIdealClosedCycleLayer` 仍指向 ensure 时刻快照，与父层**解耦** |
| `_inner_cma_fast` n_sc==0 | 直接返回 obj=1.0, n_evals=0，跳过所有 CMA 迭代 |

---

## 7. 设计约束

- **同节点两套深度**：Node.index 可同时在一机械组与一换热组；偏移必须按 `kind+组` 用 `group.depth_dict()` 查询。
- **非理想偏置顺序**：先 `ensure_non_ideal()`，再 `ni.apply_offsets()`。
- **Node.parent 类型**：`int | str | None` — 支持 `"S"` 标记 S 分位节点。
- **无自动化 `git commit/push`**。

---

## 8. 当前最优参数 & 测试结果

### 来源

`pytest -s tests/test_layered_opt.py::test_layered_0p0t0s_series`

**配置**: 0P0T0S, 理想, 串联夹点, H2∈[2,5.5], h2_T_out∈[400,900], p_min=2000kPa 固定

**优化器**: LHS=30, DE popsize=10 × maxiter=60, CMA restarts=3 × maxiter=15

| 参数 | 值 |
|------|:--|
| t_max | 938K |
| t_min | 349K |
| p_max | 8609kPa |
| p_min | 2000kPa |
| n_sc | 1 |
| flow | 29.7 kg/s |
| H2 | 5.50 kg/s |
| h2_T_out | 805K |
| obj | **0.00000** |
| 匹配 | 103.6MW |
| 未匹配 | **4kW** |
| HX 单元 | 1 (串联) |
| 耗时 | 14.9s, 610 evals |

### H2=4 固定测试 (0P0T0S)

| 参数 | 值 |
|------|:--|
| t_max | 1100K |
| t_min | 176K |
| p_max | 12000kPa |
| n_sc | 1 |
| H2 | 4.00 kg/s |
| h2_T_out | 404K |
| obj | 0.04095 |
| 匹配 | 120.7MW |
| 未匹配 | 11.1MW |

### H2=4 固定 + 1P 分位 (1P0T0S, 100gen)

| 参数 | 值 |
|------|:--|
| t_max | 978K |
| t_min | 50K |
| p_max | 8831kPa |
| p_q | 0.213 |
| n_sc | 4 |
| flows | [31.1, 40.2, 1.7, 18.2] kg/s |
| H2 | 4.00 kg/s |
| h2_T_out | 832K |
| obj | 0.02168 |
| 匹配 | 180.4MW |
| 未匹配 | 8.0MW |
| 耗时 | 126.9s, 1515 evals |

### 关键结论

- **H2 流量是核心瓶颈**：H2=5.5 串联夹点可实现近乎精确解（4kW 未匹配），H2=4 固定时未匹配跳至 11MW。
- **串联夹点模式有效**：直接暴露系统级功率不平衡和温差违规，单一逆流单元检查所有过程边界。
- **P 分位增加拓扑自由度**（n_sc 1→4），在 H2 受约束时可提升总匹配功率（120→180MW），但多加的换热过程会放大功率差。
