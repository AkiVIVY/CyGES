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
│   ├── system.py                # 系统层：Pipeline + analyze_system_heat (解耦)
│   └── fluid_property_solver.py # 物性：PropertyRegistry (多工质)
├── optimize/
│   ├── solver.py                # Optimizer: CMA-ES(BIPOP)+DE, 高斯基编码, 并行
│   ├── objective.py             # @register 目标函数注册表
│   └── types.py                 # OptimizationResult
├── config.py
├── tests/
│   ├── test_system_full.py     # 全管线可视化 5 图
│   └── test_optimize.py        # 优化测试 + Excel 日志 + run_XXX 输出
└── docs/architecture.md
```

---

## 3. 当前实现范围

**已实现**

- 理想层：PS 离散拓扑、最小子循环、简化边。
- 非理想偏置：`ni.apply_offsets()` 单步 σ 改 P → PS 重置 → DFS 机械步。
- 性能统计：`ProcessRecord`/`CyclePerformanceReport`（纯计算，零外部依赖）。
- T-Q 构建 + 夹点分析：`build_heat_tq_curves`、`compute_pinch`（底层）、`analyze_pinch`（高层）。
- 系统集成：`SystemPipeline.run(props)` 产纯数据，`analyze_system_heat(raw, inp, props)` 独立做夹点（四种模式）。
- 多工质：`PropertyRegistry`（自动缓存 CoolProp 求解器）。
- 优化器：CMA-ES (BIPOP 5 重启) + DE，高斯基编码（维度固定），分位/流量/T 边界离散化，持久 `ProcessPoolExecutor` 6 进程并行。

**夹点模式结论**

| 模式 | 状态 | 说明 |
|------|:--:|------|
| `system_pinch` | ✅ **最优** | 所有换热统一匹配，H₂ 最低、功率最高 |
| `pinch`（循环+系统） | ✅ 可用 | 两级匹配，略差于 system_pinch |
| `split_pinch` | ⚠️ 未验证 | 热/冷源独立匹配 |
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

---

## 5. 测试

| 文件 | 用途 | 运行 |
|------|------|------|
| `test_optimize.py` | 优化测试：CMA-ES + 结果图 + Excel 日志 | `pytest -s tests/test_optimize.py` |
| `test_system_full.py` | 全管线可视化 5 图 | `pytest tests/test_system_full.py` |

输出：`tests/run_XXX/` 下 6 张图 + `tests/opt_log.xlsx` 累积日志。

---

## 6. 失效与快照不变量

| 触发 | 行为 |
|------|------|
| `ClosedCycleLayer.analyze_topology()` | 重建 topology + `simplified`，清空 `non_ideal` |
| `ClosedCycleLayer.commit_subcycle_mass_flows_to_topology()` | 同步子循环流量 → 汇聚到边 → 重建 `simplified`，清空 `non_ideal` |
| `ClosedCycleLayer.ensure_non_ideal()` | 创建 `NonIdealClosedCycleLayer` 快照；要求 `simplified` 存在 |
| 父层变动后 | 已存在 `NonIdealClosedCycleLayer` 仍指向 ensure 时刻快照，与父层**解耦** |

---

## 7. 设计约束

- **同节点两套深度**：Node.index 可同时在一机械组与一换热组；偏移必须按 `kind+组` 用 `group.depth_dict()` 查询。
- **非理想偏置顺序**：先 `ensure_non_ideal()`，再 `ni.apply_offsets()`。
- **主脊 layer**：最长路径定标 + 并列最小起点。
- **无自动化 `git commit/push`**。

---

## 8. 当前最优参数（system_pinch 模式, H₂=3.5）

| 参数 | 值 |
|------|:--|
| σ / η_is | 0.98 / 0.9 |
| He t∈[40, 1000]K, p∈[2, 10] MPa |
| 分位：1t+2p, step=0.01 |
| 基函数：3×3, mf∈[0, 50], step=0.01 |
| CMA-ES, popsize auto, maxiter=300, 6 workers |
| obj=0.0000, t_max≈1090K, net_mech≈-38.9MW |
