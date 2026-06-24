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
│   └── fluid_property_solver.py # 物性：PropertyRegistry(含memo缓存) + InterpolatingHeliumSolver
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
│   ├── test_layered_opt.py     # 分层优化: LHS+DE外层 + CMA/L-BFGS-B内层, TS/PS/HX图
│   ├── test_inner_compare.py   # 内层方法对比: CMA-ES vs L-BFGS-B (多起点LHS)
│   ├── test_outer_quantile_opt.py # 外层分位点优化: 网格扫描 + L-BFGS-B中心验证
│   └── optimization_example.py # 效率优化示例: CH4冷源, η最大化
│   └── test_feasibility_de.py # 新框架: h2_T_out 外层纳入 + 内层 flows-only
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
- 多工质：`PropertyRegistry` 含 memo 缓存 `(fluid, pair, x, y)` 去重；冷源支持任意 CoolProp 流体（`cold_fluid` 参数）。
- 物性加速：`InterpolatingHeliumSolver`（通过 `use_interp_he` 开关启用，默认关闭）。
- **HX 串联夹点匹配**：`match_series_pinch()` — 所有放热/吸热过程按温度降序在 Q 上串联为单一逆流单元，在每条过程边界检查 dT_min 违规。
- 星型 HX 匹配：`match_heat_exchanger_groups()` / `match_constructive()` / `match_heat_exchanger_staged()`。
- **分层优化** `test_layered_opt.py`：外层 LHS+DE 搜索拓扑参数（含 s_q），内层支持三种方法（CMA-ES / L-BFGS-B / hybrid 两阶段），`_DEFAULT_HP` 集中所有可调参数。
- **内层 L-BFGS-B 多起点优化**：LHS 采样初始点 + best-track 记录最小值 + scipy L-BFGS-B（`ftol=1e-8`）+ 多进程并行（`w=8` 确定性验证通过）。
- **Smooth softplus 惩罚**：`log(1+exp(k*(util-tol)))/k` 替代不可微 ReLU，使 L-BFGS-B 梯度法可用；`penalty_k=10`。
- **目标函数模式**：`_eval_fast` 支持 4 种 `obj_mode`。
- **h2_T_out 外层新框架** `test_feasibility_de.py`：将 h2_T_out 从内层 L-BFGS-B 变量提升为外层 DE 变量，内层仅优化子循环流量（dim=n_sc）。```outer(DE 7D: [拓扑 + h2_T_out]) → inner(L-BFGS-B S48 w=8, dim=n_sc flows only)```。h2_T_out 搜索边界 [500,1000]K。
- **0P1S / 1P0S 拓扑对比**：0P1S (S分位, n_sc=2) 远优于 1P0S — η 从 0.530 → 0.621。S 分位在 PS空间创造更高效换热匹配。0P1S 内层 spread=0（确定性收敛），1P0S spread=0.006（S64 最佳）。
- 可视化：理想/非理想 TS/PS（子循环多边形叠加）、HX T-Q（概览 + 每单元逆流）、DE 收敛曲线、夹点复合曲线。

**目标函数模式**

| obj_mode | 公式 | 用途 |
|------|------|------|
| `group` | `unmatched_ratio + 1e-3 * num_unmatched` | 星型 HX 匹配（默认） |
| `series` | `unmatched_ratio + 1e-2 * pinch_violation` | 串联夹点 HX 匹配 |
| `pinch` | `(hot_util + cold_util) / q_source` | 直接最小化公用工程量（**内层推荐**） |
| `eff_pinch` | `-η + penalty_w * softplus(util,tol,k) / q_source` | 最大效率 + 公用工程容差 |
| `eff_pinch_split` | `-η + penalty_w × (softplus(hot) + softplus(cold)) / q_source` | 分拆冷热惩罚（实验性） |
| `pinch_aligned` | `pen_w × softplus(util−tol,k)/q_source − pen_dt × min_dT / 100` | 夹点最大化 + 公用工程惩罚（**内层新推荐**） |

**`pinch_aligned` 模式详解**

| 特性 | 说明 |
|:--|:--|
| 夹点分析 | `compute_pinch_fixed_alignment(hot_curve, cold_curve)` — Q=0 左对齐，扫重叠区取 min ΔT |
| 目标 | 公用工程越小越好 + 夹点温差越大越好 |
| 关键优势 | min_dT 项在 util=0 时**仍提供非零梯度** → 消除 util=0 面上 obj 全等的平坦问题 |
| 引入位置 | `core/postprocess.py` 新增 `PinchFixedResult` + `compute_pinch_fixed_alignment` |
| 验证结果 | 0P1S S32 w=16: 9 seeds **全部收敛到同一 flows** [10.9, 22.4], spread=0.000, pinch dT=24.2K |

> 原 `pinch` 模式在 util=0 时 obj 全等、不同 seed 收敛到不同 flows 组合；`pinch_aligned` 模式通过 min_dT 奖励创建唯一最优解，消除多盆地。

**内层优化方法**

| inner_method | 说明 |
|------|------|
| `lbfgsb` | L-BFGS-B + LHS 多起点 + best-track + 并行（唯一方案） |

**内层默认配置（2024-06-30 最终）**

| 参数 | 值 | 说明 |
|------|:--|------|
| starts | 32 | S32 balance spread vs cost; S64 for 1P0S(spread 0.006) |
| maxiter | 20 | ftol 先于 maxiter 收敛 |
| workers | 16 | 内层并行 |
| ftol | 1e-8 | 不可降 — best 退化 3~5× |
| eps | default | 不可改 — 差分精度最优 |
| obj_mode | pinch_aligned | util penalty + pinch reward (pen_dt=0.5) |
| sampler | LHS | Sobol 反效果 |

**内层优化方法对比 (1P0S pinch_aligned pen0.5)**

| 方法 | best | spread | 评价 |
|:--|:--|:--|:--|
| **L-BFGS-B S64** | **-0.044** | **0.006** | 首选 |
| L-BFGS-B S256 | -0.047 | 0.004 | 成本 8× |
| DE pop30×35 | -0.043 | 0.042 | spread 差 |
| Adam | +57→0.78 | 583 | 不可用 |

**L-BFGS-B starts 收敛:** S32(0.014) → S64(0.006) → S128(0.005) → S256(0.004) — 递减改善。

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

**高斯基编码**：子循环流量由 PS 空间上的高斯核加权和生成，维度固定，与子循环数无关。**优化器层已去除量化步长**（flow_step/h2_step/qstep），直接使用连续值。

**离散化在优化器层**：流量/分位/T 边界步长全在 `Optimizer` 中控制。

**S 分位节点**：`ClosedCycleTPInput.s_quantiles` 与 t_quantiles/p_quantiles 平级。S 分位节点由 primary nodes 的 S_min/S_max 线性插值得到，对每个 P 轴做 PS 延伸（`solver.state("PS", P, S_q)`），`parent="S"`。新节点通过同 parent 的 P 排序链形成机械边，并参与换热边构建。

**串联夹点匹配**：所有过程按 T_high 降序排列，在累计 Q 上构成两条链（hots + colds 均降序）。单一单元检查每个过程边界的逆流温差。目标函数 = `unmatched_ratio + 1e-2*pinch_violation`，直接暴露系统级的功率不平衡与温差违规。

**分层优化**：外层 LHS/DE 搜索拓扑参数（t_min/t_max/p_min/p_max/t_q/p_q/s_q），内层在固定拓扑下优化流量 + H2 + h2_T_out。支持 `h2_T_fixed` 和 `h2_fixed` 模式。`_DEFAULT_HP` 中 `inner_method` 选内层方法，`lbfgsb_starts`/`lbfgsb_maxiter`/`lbfgsb_workers` 控制 L-BFGS-B。

**多工质冷源**：`_make_system_input()`/`_make_h2_source()` 支持 `cold_fluid` 参数切换冷源工质（H₂/CH₄ 等任意 CoolProp 流体）。

**Smooth softplus 惩罚**：`_eval_fast` 在 `obj_mode='eff_pinch'` 时使用 `softplus(k*(util-tol))/k` 的光滑惩罚，数值安全分段处理（excess>50→线性, excess<-10→0）。梯度连续可微，L-BFGS-B 可用。

**PropertyRegistry memo 缓存**：`(fluid, pair, x, y)` 键 → round-4 精度值，避免重复 CoolProp 查询，累计提速 11.5%。

**n_sc==0 短路**：`_inner_cma_fast`/`_inner_lbfgsb_fast`/`_inner_hybrid_fast` 中 n_sc==0 时直接返回 obj=1.0，跳过迭代。

**内层 obj_mode 选择**：当 h2_T_out 固定（外层搜索或新框架），η 为常数 → `eff_pinch` 的 `-η + penalty_w × softplus / q_source` 等价于 `const + C × util`，其中 penalty_w 和 k 造成无效比例尺混乱。改用 `obj_mode='pinch'`（纯 `util/q_source`），spread 从 0.29 → 0.0014（200× 改善），evals 减少 30%。新框架 `test_feasibility_de.py` 内层应用 `obj_mode='pinch'`。

**`pinch_aligned` 模式**：`pen_w × softplus(util−tol) / q_source − pen_dt × min_dT / 100`。`compute_pinch_fixed_alignment`（Q=0 左对齐，扫重叠区 min ΔT）替代原 δQ 优化对齐。min_dT 项在 util=0 时仍提供非零梯度 → 消除 util=0 面上 obj 全等的平坦问题。0P1S S32 验证 9 seeds 全部收敛到同一 flows [10.9, 22.4], spread=0, pinch dT=24.2K。

**熵序警告**：`build_node_edge_topology` 中 S(t_max,p_max) < S(t_min,p_min) 时发出 RuntimeWarning。

**best-track 必需**：记录优化过程中的 obj 最小值，非最终迭代值。

---

## 5. 测试

| 文件 | 用途 | 运行 |
|------|------|------|
| `test_layered_opt.py` | 分层优化，多个预设配置 | `pytest -s tests/test_layered_opt.py::<test_name>` |
| `test_optimize.py` | CMA-ES + DE 优化测试 | `pytest -s tests/test_optimize.py` |
| `test_inner_compare.py` | CMA-ES vs L-BFGS-B 内层对比 | `pytest -s tests/test_inner_compare.py::test_compare_inner` |
| `optimization_example.py` | CH₄ 效率优化示例 | `pytest -s tests/optimization_example.py::test_h2_scan_efficiency` |
| `test_feasibility_de.py` | h2_T_out 外层纳入新框架 | `pytest -s tests/test_feasibility_de.py::test_1p0s_baseline` |

可选测试函数（`test_layered_opt.py`）：

| 函数 | Topo | HX | H2 | 说明 |
|------|------|----|----|------|
| `test_layered_0p0t_hx1` | 1P+0T+1S | series | [2,5.5] | 串联夹点 |
| `test_layered_0p0t0s_series` | 0P0T0S | series | [2,5.5] | 最简串联 |
| `test_layered_0p0t0s_h2_4` | 0P0T0S | series | 4固定 | H2=4 串联 |
| `test_layered_0p1t0s_h2_4` | 0T+1P+0S | series | 4固定 | +P分位 |
| `test_layered_0p0t_sq1` | 0P0T+1S | group HX2 | [3,6] | S分位测试 |
| `test_layered_1p1t_ideal` | 1T+1P | group HX2 | [3,6] | 1P1T 基准 |
| `test_layered_0p1s_h2_4_lbfgsb` | 0P1S | eff_pinch | 4固定 | L-BFGS-B 内层 |
| `test_layered_1p1s_h2_4_lbfgsb` | 1P1S | eff_pinch | 4固定 | L-BFGS-B 内层 |
| `test_inner_method_compare` | 1P1S | eff_pinch | 3.5固定 | CMA/L-BFGS-B/hybrid 三法对比 |
| `test_inner_global_search_0p1s` | 0P1S | eff_pinch | 3.5固定 | 多seed全局搜索 |
| `test_h2_scan_1p1s` | 1P1S | eff_pinch | 扫描 | H₂流量扫描 |
| `test_h2_scan_1p1s_pinch` | 1P1S | pinch | 扫描 | pinch模式扫描 |

输出：`tests/run_layered_fast/`、`tests/inner_compare/`、`tests/optimization_example/` 下图表。

---

## 6. 失效与快照不变量

| 触发 | 行为 |
|------|------|
| `ClosedCycleLayer.analyze_topology()` | 重建 topology + `simplified`，清空 `non_ideal` |
| `ClosedCycleLayer.commit_subcycle_mass_flows_to_topology()` | 同步子循环流量 → 汇聚到边 → 重建 `simplified`，清空 `non_ideal` |
| `ClosedCycleLayer._rebuild_simplified()` | 清空 `non_ideal` + 重建 `simplified`（子循环为空或全零流量时 warn） |
| `ClosedCycleLayer.ensure_non_ideal()` | 创建 `NonIdealClosedCycleLayer` 快照；要求 `simplified` 存在 |
| 父层变动后 | 已存在 `NonIdealClosedCycleLayer` 仍指向 ensure 时刻快照，与父层**解耦** |
| 内层 n_sc==0 | 直接返回 obj=1.0, n_evals=0，跳过所有迭代 |

---

## 7. 设计约束

- **同节点两套深度**：Node.index 可同时在一机械组与一换热组；偏移必须按 `kind+组` 用 `group.depth_dict()` 查询。
- **非理想偏置顺序**：先 `ensure_non_ideal()`，再 `ni.apply_offsets()`。
- **Node.parent 类型**：`int | str | None` — 支持 `"S"` 标记 S 分位节点。
- **无自动化 `git commit/push`**。
- **内层方法选择**：`_DEFAULT_HP.inner_method` ∈ `{"cma","lbfgsb","hybrid"}`；`lbfgsb` 为推荐默认。
- **并行确定性**：`lbfgsb_workers=8` 不引入非确定性（w=1 vs w=8 结果完全一致）。
- **冷源工质**：`_DEFAULT_HP.cold_fluid` 支持任意 CoolProp 流体名。

---

## 8. 当前最优参数 & 测试结果

### 来源 1: 最简串联夹点 (0P0T0S, series)

`pytest -s tests/test_layered_opt.py::test_layered_0p0t0s_series`

**配置**: 0P0T0S, 理想, 串联夹点, H2∈[2,5.5], h2_T_out∈[400,900], p_min=2000kPa 固定
**优化器**: LHS=30, DE popsize=10 × maxiter=60, CMA restarts=3 × maxiter=15

| 参数 | 值 |
|------|:--|
| t_max | 938K |
| t_min | 349K |
| p_max | 8609kPa |
| n_sc | 1 |
| flow | 29.7 kg/s |
| H2 | 5.50 kg/s |
| h2_T_out | 805K |
| obj | **0.00000** |
| 匹配 | 103.6MW |
| 未匹配 | **4kW** |
| 耗时 | 14.9s, 610 evals |

### 来源 2: L-BFGS-B 内层对比 (1P1S, eff_pinch)

`pytest -s tests/test_inner_compare.py::test_compare_inner`

**配置**: 1P1S, p_q=s_q=0.5, H2=3.5kg/s 固定, obj_mode=eff_pinch, 理想

| 方法 | η | evals | 时间 |
|------|:--|:--|:--|
| CMA-b100 | -0.069 | 584 | 0.94s |
| CMA-b250 | 0.334 | 1,387 | 0.85s |
| CMA-b500 | 0.483 | 2,750 | 1.69s |
| CMA-b1000 | 0.547 | 5,129 | 3.11s |
| **L-BFGS-B×16** | **0.679** | 16,432 | 2.45s |

> L-BFGS-B 远超 CMA-b1000（2.45s vs 3.11s, η 0.679 vs 0.547）。

### 来源 3: CH₄ 冷源效率优化 (1P1S, eff_pinch)

`pytest -s tests/optimization_example.py::test_h2_scan_efficiency`

**配置**: 1P1S, CH₄ 120K→900K, Air 1000K→500K, p_min=2000kPa 固定

| 参数 | 值 |
|------|:--|
| CH₄ | 11.0 kg/s |
| t_max | 932.6K |
| t_min | 150K |
| p_max | 12000kPa |
| n_sc | 4 |
| flows | 37.7 kg/s 总和 |
| **η** | **0.431** |
| hot_util | 0.9 kW |
| cold_util | 0 kW |

### 来源 4: 0P1S 分拆惩罚确定性收敛

`pytest -s tests/test_layered_opt.py::test_inner_global_search_0p1s`

**配置**: 0P1S, H2=3.5kg/s 固定, split penalty (hot/cold 分开), L-BFGS-B S32

| 参数 | 值 |
|------|:--|
| **η** | **0.627** |
| spread | 0.0002 (所有 seed 一致) |
| flows | [13.8, 26.3] kg/s |
| h2_T_out | 641K |

### 来源 5: 1P0S 天花板

**配置**: 1P0S, H2=3.5kg/s 固定, L-BFGS-B S96

| 参数 | 值 |
|------|:--|
| **η** | **0.521** |
| spread | 0.024 |
| mean_η | 0.501~0.513 |

> P 分位在此工况下非活跃变量（η 0.627→0.521）。

### 来源 6: 1P0S S256 多 seed 更大规模对比

**配置**: 1P0S, H2=3.5kg/s 固定, L-BFGS-B S256 w=8 maxiter=80, 9 seeds

| 参数 | 值 |
|------|:--|
| **best η** | **0.523** |
| spread (9 seeds) | 0.050 |
| spread (除 seed=800 outlier) | **0.020** |
| mean evals/seed | ~180k |
| 总耗时 | 110.8s |
| outlier (seed=800) | η=0.473, h2_T_out=891K (脱靶) |

> S256 将内聚度提升到 8/9 seeds η∈[0.503,0.523]，但仍有 1 个 outlier 因 h2_T_out 脱出 40K 峰区。

### 来源 7: h2_T_out 扫描（固定拓扑）

**配置**: 1P0S, t_max=1000, p_max=10000, p_q=0.5, H2=3.5

| h2_T_out [K] | η | 备注 |
|:--|:--|:--|
| 500~775 | 负 (obj>0) | 不可行区 |
| 800 | ≈0 | 临界点 |
| **825** | **0.514** | 峰顶 |
| 850~1000 | 0.50→0.40 | 衰减区 |

> h2_T_out 是单峰光滑变量，峰宽~40K，500K 范围中有效区仅占~8%。验证了 h2_T_out 应从内层 LHS 随机搜索中移除。

### 来源 8: LHS vs Sobol (1P0S S96)

| 采样 | best η | mean η | spread | evals |
|:--|:--|:--|:--|:--|
| LHS | 0.5194 | 0.5056 | **0.030** | 632k |
| Sobol | 0.5216 | 0.5029 | 0.064 | 605k |

> Sobol 在 n=96(非2^k) 时退化，且 seed=500 出现 η=0.458 的严重离群值。

### 来源 9: 变量分块测试 (h2_T_out=812 固定, dim=4, 9 seeds)

| 指标 | dim=5 (含 h2_T_out) | dim=4 (固定 h2_T_out=812K) |
|:--|:--|:--|
| best η | 0.519 | 0.522 |
| spread | **0.030** | **0.068** |

> 固定 h2_T_out 后 6/9 seeds 收敛到同一最优值（η=0.522），但 3 个 seed 陷入更差盆地。变量分块无效 — h2_T_out 与 flows 交互不可忽略。

### 来源 10: 惩罚参数扫描 (1P0S, h2_T_out=750K 固定, dim=4)

**配置**: S48 w=8 maxiter=80, 9 seeds, 对比三种目标

| obj_mode | best | mean | spread | evals/seed |
|:--|:--|:--|:--|:--|
| eff_pinch w=1000 k=10 | 2.810 | 2.899 | **0.290** | ~26k |
| **pinch (direct)** | **0.0034** | **0.0038** | **0.0014** | **~18k** |
| eff_pinch w=2000 k=5 | 6.253 | 7.221 | 3.540 | ~27k |

> `pinch` 模式的 spread 比 `eff_pinch` 低 200×，evals 少 30%。h2_T_out 固定后 η 为常数，`eff_pinch` 的 penalty_w×softplus 等价于 `const + C × util`，引入了无效比例尺。`pinch` 直接 `util/q_source` 给出干净梯度信号。

### 来源 11: flows 精度验证

3 seed 打印 6 位小数：flows 收敛到不同连续值（差异 0.07~0.38 kg/s），非截断所致。1 位小数的显示（原 `{v:.1f}` 格式）是之前"看起来一致"的原因。

### 来源 12: pinch_aligned 验证 (0P1S, h2tout=800K)

**配置**: 0P1S, h2tout=800K 固定, S32 w=16, 9 seeds, pinch vs pinch_aligned

| obj_mode | best | spread | flows |
|:--|:--|:--|:--|
| pinch | 0.00000 | 0.000 | [13.6,21.6], [8.7,23.1]... 多个不同 |
| **pinch_aligned** | **-0.24193** | **0.000** | **[10.9,22.4] 全部相同** |

> pinch_aligned 通过 min_dT 奖励 (24.2K) 创建唯一最优解。9 seeds 全部收敛到同一 flows，彻底消除 pinch 模式在 util=0 面上的多盆地问题。

### 关键结论

- **H2 流量是核心瓶颈**：H2=5.5 串联夹点可实现近乎精确解（4kW 未匹配）。
- **L-BFGS-B 远超 CMA-ES**：固定拓扑下梯度法 + smooth softplus 收敛远优于演化算法。
- **分拆惩罚仅低维有效**：0P1S (d=3) 确定性收敛，1P0S (d=5) 反效果。
- **局部极值随维度增长**：spread d=3(0.000)→d=5(0.024)→d=7(0.12)，每增 2 维扩大 5-6×。
- **量化步长阻塞梯度**：flow_step/h2_step/qstep 已全部移除，L-BFGS-B 直接连续优化。

**h2_T_out 变量归属**：h2_T_out 经扫描和 LHS 多 seed 验证为~40K 宽单峰变量（峰区 800~840K，边界外 η 崩盘）。不应放入内层随机 LHS 多起点搜索（dim=5 时峰区覆盖率仅 ~11%），应作为外层系统级变量或以确定性 Brent 线搜索处理。`test_feasibility_de.py` 将 h2_T_out 纳入外层 DE 7D 向量。

**LHS vs Sobol 采样**：Sobol 在 n≠2^k 时平衡性退化，d=5 时对多盆地地貌反效果（spread 0.030→0.064）。LHS 仍是当前推荐起点生成方法。

**变量分块不适用**：尝试固定 h2_T_out 将 dim 从 5 降为 4，spread 反而从 0.030 恶化到 0.068 — 因为 flows 盆地无法通过调节 h2_T_out 补偿。h2_T_out 地貌分层：`h2_T_out`（主导，窄峰）>> `flows`（次级，多盆地）。

**内层起点数与收敛**：S96→S128→S256 spread 不归零（0.030→0.033→0.050，除 outlier 外 0.020）。ftol=1e-8 已超越 maxiter=80，提高代数无改善。瓶颈在地貌而非收敛精度。

**内层 obj_mode**：h2_T_out 固定后 η 为常数，`pinch`（`util/q_source`）优于 `eff_pinch`（spread 0.0014 vs 0.29，200× 改善）。`pinch_aligned`（`pen_w × softplus(util) / q_source − pen_dt × min_dT / 100`）更优——min_dT 项在 util=0 时仍提供非零梯度 → 9 seeds 全部收敛到同一 flows。新框架 `test_feasibility_de.py` 内层应用 `obj_mode='pinch'`。

**核心架构问题（2024-06-24 确认）**：所有问题的核心出在内层计算的不确定性上。当内层为优化结果时，如果内层存在局部解，外层的边界可能导致内层计算结果不连续，最终使得整个计算失效。0P1S (n_sc=2) 内层完全确定 (spread=0)，1P0S (n_sc=4) 仍有 0.006~0.014 的不确定性——这在嵌套优化架构下会传播到外层目标函数的噪声中。下一版本考虑用核函数方法等替代当前纯 LHS 多起点+L-BFGS-B 的内层搜索策略，以消除剩余的地貌不确定性。

**变量分块不适用**：尝试固定 h2_T_out 将 dim 从 5 降为 4，spread 反而从 0.030 恶化到 0.068 — 因为 flows 盆地无法通过调节 h2_T_out 补偿。h2_T_out 地貌分层：`h2_T_out`（主导，窄峰）>> `flows`（次级，多盆地）。

**内层起点数与收敛**：S96→S128→S256 spread 不归零（0.030→0.033→0.050，除 outlier 外 0.020）。ftol=1e-8 已超越 maxiter=80，提高代数无改善。瓶颈在地貌而非收敛精度。

**内层 obj_mode**：h2_T_out 固定后 η 为常数，`pinch`（`util/q_source`）优于 `eff_pinch`（spread 0.0014 vs 0.29，200× 改善）。`pinch_aligned`（`pen_w × softplus(util) / q_source − pen_dt × min_dT / 100`）更优——min_dT 项在 util=0 时仍提供非零梯度 → 9 seeds 全部收敛到同一 flows。新框架 `test_feasibility_de.py` 内层应用 `obj_mode='pinch'`。
