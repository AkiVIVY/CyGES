"""
CyGES 项目级可调数值参数（不含密钥）。

运行拓扑等代码时请保证项目根在 ``PYTHONPATH`` 中，以便 ``import config``。
"""

# 子循环质量流在 ``ClosedCycleLayer.analyze_topology`` 中的初值：
# 每项 = SUBCYCLE_INITIAL_MASS_FLOW_FRACTION_OF_MAX * max_mass_flow
# （``max_mass_flow`` 见 ``ClosedCycleTPInput``；若为 ``None`` 则按 0 参与乘法）
SUBCYCLE_INITIAL_MASS_FLOW_FRACTION_OF_MAX: float = 0.1

# 子循环质量流量化：步长 = subcycle_mass_flow_step_fraction * max_mass_flow；
# ``ClosedCycleTPInput`` 未显式传入 ``subcycle_mass_flow_step_fraction`` 时使用本默认值。
SUBCYCLE_MASS_FLOW_STEP_FRACTION_DEFAULT: float = 0.01
