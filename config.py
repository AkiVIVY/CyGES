"""
CyGES 项目级可调数值参数（不含密钥）。

运行拓扑等代码时请保证项目根在 ``PYTHONPATH`` 中，以便 ``import config``。
"""

# 子循环质量流在 ``ClosedCycleLayer.analyze_topology`` 中的初值：
# 每项 = mass_flow_min + SUBCYCLE_INITIAL_MASS_FLOW_FRACTION_OF_MAX * (mass_flow_max - mass_flow_min)
# （``mass_flow_min`` / ``mass_flow_max`` 见 ``ClosedCycleTPInput``）
SUBCYCLE_INITIAL_MASS_FLOW_FRACTION_OF_MAX: float = 0.1

# 子循环质量流量化：步长 = subcycle_mass_flow_step_fraction * (mass_flow_max - mass_flow_min)；
# ``ClosedCycleTPInput`` 未显式传入 ``subcycle_mass_flow_step_fraction`` 时使用本默认值。
SUBCYCLE_MASS_FLOW_STEP_FRACTION_DEFAULT: float = 0.001

# 全层子循环质量流量边界默认值 [kg/s]；``ClosedCycleTPInput`` 未显式传入时使用。
MASS_FLOW_MIN_DEFAULT: float = -5.0
MASS_FLOW_MAX_DEFAULT: float = 10.0

# 非理想精简「机械边」统一等熵效率 η_is ∈ (0, 1]。
# 拓扑上 mechanical 边表示叶轮机械工作过程（压缩/膨胀等）；非理想层后续按该 η_is
# 修正状态（等熵效率定义：实际功与等熵功之比，或等价的焓升/降相对理想过程的折减）。
# 临时配置：组内所有精简机械边 SM* 共用该值，后续会按边分别赋值。
NON_IDEAL_MECHANICAL_EFFICIENCY_DEFAULT: float = 0.9

# 非理想精简「换热边」统一总压恢复系数 σ ∈ (0, 1]（沿 tail→head 的压力保留比）。
# 拓扑上 heat 边表示换热过程（加热/冷却）；σ 表征流经该过程后工质总压相对理想
# 换热路径的保留程度（``non_ideal_bias._apply_combined_offsets`` 步骤 1：P ← P_ideal × σ^layer）。
# 临时配置：组内所有精简换热边 SH* 共用该值，后续会按边分别赋值。
NON_IDEAL_HEAT_EFFICIENCY_DEFAULT: float = 0.95
