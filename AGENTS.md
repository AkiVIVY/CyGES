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
│   ├── closed_cycle_layer.py    # 理想层：拓扑构建、子循环枚举、简化
│   ├── non_ideal_bias.py        # 非理想偏置：有向组、深度、ni.apply_offsets()
│   ├── cycle_performance.py     # 性能统计：ProcessRecord, LoopReport (纯计算)
│   ├── postprocess.py           # T-Q 构建 + 夹点分析 (build_tq / analyze_pinch)
│   ├── heat_exchanger.py        # HX 匹配：星型拓扑, 贪婪打包, dT_min 约束
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
│   └── test_layered_opt.py     # 分层优化: LHS+DE外层 + CMA内层, TS/PS/HX图表
└── docs/architecture.md
    └── core_bugs.md             # 已知问题备忘
```

---

## 3. 当前实现范围

**已实现**

- 理想层：PS 离散拓扑、最小子循环、简化边。
- 非理想偏置：`ni.apply_offsets()` 单步 σ 改 P → PS 重置 → DFS 机械步。
- 性能统计：`ProcessRecord`/`CyclePerformanceReport`（纯计算，零外部依赖）。
- T-Q 构建 + 夹点分析：`build_heat_tq_curves`、`compute_pinch`（底层）、`analyze_pinch`（高层）。
- 系统集成：`SystemPipeline.run(props, layers=...)` 产纯数据，支持预建 layer 复用拓扑（只变流量时 25× 加速）；`analyze_system_heat(raw, inp, props)` 独立做夹点（五种模式）。
- 多工质：`PropertyRegistry`（自动缓存 CoolProp 求解器）。
- 物性加速：`InterpolatingHeliumSolver`（预计算 (T,P) 网格，双线性插值 + 二分查找替代 CoolProp 迭代；通过 `use_interp_he` 开关启用）。InterpHe 拓扑断裂时自动回退 CoolProp。
- HX 匹配：星型拓扑（1 热 + k 冷 / m 热 + 1 冷）；两种候选生成方式：
  - **枚举法** `match_heat_exchanger_groups()` — `itertools.combinations` + 贪心多起点打包
  - **构造法** `match_constructive()` — Q 优先排序 + 残差最小化 + 双向遍历
  - **三阶段分步** `match_heat_exchanger_staged()` — 热源↔吸热→冷源↔放热→内部消纳
- 优化器：CMA-ES (BIPOP 重启) + DE；分层优化 `test_layered_opt.py`：外层 LHS+DE 搜索拓扑参数，内层 CMA 复用固定拓扑优化流量 + H2 + h2_T_out。
- 可视化：理想/非理想 TS/PS（子循环多边形叠加）、HX T-Q（概览 + 每单元逆流）、DE 收敛曲线。

**夹点模式结论**

| 模式 | 状态 | 说明 |
|------|:--:|------|
| `system_pinch` | ✅ **最优** | 所有换热统一匹配，H₂ 最低、功率最高 |
| `pinch`（循环+系统） | ✅ 可用 | 两级匹配，略差于 system_pinch |
| `split_pinch` | ⚠️ 未验证 | 热/冷源独立匹配 |
| `staged_pinch` | 🆕 新增 | 三阶段：外部热源↔循环吸热 → 外部冷源↔循环放热 → 残差消纳 |
| `source_pinch` | ❌ 无效 | 分步匹配破坏全局热平衡 |

**未实现**

- 非理想方程/约束装配、多目标优化器。
- HEN / 多热源多冷源边界耦合。
- 按边/按组分别赋值 σ、η_is。

---

## 4. 关键设计决策

**夹点分析与管线解耦**：`SystemPipeline.run()` 产纯数据，`analyze_system_heat()` 独立做夹点。可对同一 `run()` 输出调用不同模式。

**高斯基编码**：子循环流量由 PS 空间上的高斯核加权和生成，维度固定（3×3=9），与子循环数无关。比直接优化 mf[] 维度稳定，收敛更快。

**离散化在优化器层**：流量/分位/T 边界步长、分位合并比全在 `Optimizer` 中控制。模型中不保留任何离散化逻辑。

**HX 匹配与管线解耦**：HX 匹配在目标函数 `hx_unmatched` 内调用匹配函数，管线本身只产 `ProcessRecord` 数据。三种匹配方式：枚举法 `match_heat_exchanger_groups`（默认 `max_group_size=3`）、构造法 `match_constructive`（`n_restarts=50`）、三阶段分步 `match_heat_exchanger_staged`（外部热源优先消纳）。

**物性加速（interp He）**：`InterpolatingHeliumSolver` 预计算 (T,P) 网格，用双线性插值 + 等压线二分查找替代 CoolProp 迭代求根。通过 `Optimizer(use_interp_he=True)` 启用，每个 worker 缓存一份网格实例。精度（H 误差 <0.03%，S 误差 <0.4%）对优化影响极小，但拓扑可能因 S 累计误差产生少量差异，最终解宜用 CoolProp 验证。

**分层优化**：外层 LHS/DE 搜索拓扑参数（t_min/t_max/p_min/p_max/分位），内层 CMA 在固定拓扑下优化流量 + H2 + h2_T_out。拓扑复用实现 25× 内层加速。`test_layered_opt.py` 是全参数可调的完整实现，`_DEFAULT_HP` 集中所有参数。

**构造式 HX 匹配**：`match_constructive` 用 Q 优先排序 + 残差最小化 + 双向遍历替代 `combinations` 枚举，O(R·N·avg_k) 无组合爆炸。`match_heat_exchanger_staged` 三阶段分步匹配（外部源优先消纳）。

---

## 5. 测试

| 文件 | 用途 | 运行 |
|------|------|------|
| `test_optimize.py` | 优化测试：CMA-ES + 结果图 + Excel 日志 | `pytest -s tests/test_optimize.py` |
| `test_system_full.py` | 全管线可视化 5 图 | `pytest tests/test_system_full.py` |
| `test_optimize_hx.py` | HX 匹配优化：4 图 + dT_min/group_size 参数化 | `pytest -s tests/test_optimize_hx.py` |
| `test_hx_match.py` | HX 匹配单元测试（合成数据） | `pytest tests/test_hx_match.py` |
| `test_layered_opt.py` | 分层优化：LHS+DE外层 + CMA内层, TS/PS/HX图 | `pytest -s tests/test_layered_opt.py::test_layered_0p0t_wide` |

输出：`tests/run_XXX/` 下 6 张图 + `tests/opt_log.xlsx` 累积日志。

---

## 6. 失效与快照不变量

| 触发 | 行为 |
|------|------|
| `ClosedCycleLayer.analyze_topology()` | 重建 topology + `simplified`，清空 `non_ideal` |
| `ClosedCycleLayer.commit_subcycle_mass_flows_to_topology()` | 同步子循环流量 → 汇聚到边 → 重建 `simplified`，清空 `non_ideal` |
| `ClosedCycleLayer._rebuild_simplified()` | 清空 `non_ideal` + 重建 `simplified`（子循环为空或全零流量时 warn） |
| `ClosedCycleLayer.ensure_non_ideal()` | 创建 `NonIdealClosedCycleLayer` 快照；要求 `simplified` 存在 |
| 父层变动后 | 已存在 `NonIdealClosedCycleLayer` 仍指向 ensure 时刻快照，与父层**解耦** |

---

## 7. 设计约束

- **同节点两套深度**：Node.index 可同时在一机械组与一换热组；偏移必须按 `kind+组` 用 `group.depth_dict()` 查询。
- **非理想偏置顺序**：先 `ensure_non_ideal()`，再 `ni.apply_offsets()`。
- **主脊 layer**：最长路径定标 + 并列最小起点。
- **无自动化 `git commit/push`**。

---

## 8. 当前最优参数

### 0P0T_wide（h2_T_out 内层优化，无非理想，HX dT_min=10K）

| 参数 | 值 |
|------|:--|
| He t∈[50,500]K, t_max∈[700,1100]K, p∈[1000,4000]kPa, p_max∈[6000,15000]kPa |
| 分位：0p+0t, step=0.001 |
| H2流量∈[2,6], h2_T_out∈[400,900]K |
| 内层 CMA: dim=3, sigma0=15, maxiter=10, restarts=2 |
| 外层 DE: popsize=15, maxiter=40, 6 workers |
| obj=0.00067, t_max≈1063K, t_min≈210K, p_max≈12730kPa, p_min≈3909kPa |
| flows≈21.5, H2≈5.86kg/s, h2_T_out≈590K, HX: 2单元 |
