"""
非理想闭式循环层：仅作为理想层 ``ClosedCycleLayer.simplified`` 的快照容器。

简化拓扑的生成（``filter_topology_for_non_ideal`` / ``build_simplified_topology`` 等）已下沉到
``core.closed_cycle_layer``，由理想层在 ``analyze_topology()`` 与 ``commit_subcycle_mass_flows_to_topology()``
末尾通过 ``_rebuild_simplified()`` 自动同步至 ``ClosedCycleLayer.simplified``。

本层只持有：

- ``simplified``：理想层「ensure 时刻」的 ``SimplifiedTopology`` 引用。该结构是 ``frozen=True``
  且所有字段为 ``tuple`` / ``frozenset``，本身不可变；与父层后续 ``analyze`` / ``commit`` 之后生成的新
  ``SimplifiedTopology`` 对象解耦（旧引用仍指向 ensure 当时的值快照）。
"""

from __future__ import annotations

from dataclasses import dataclass

from core.closed_cycle_layer import SimplifiedTopology


@dataclass
class NonIdealClosedCycleLayer:
    """
    非理想分析容器；仅持有理想层 ``simplified`` 的快照引用。

    由 ``ClosedCycleLayer.ensure_non_ideal()`` 创建并赋给 ``parent.non_ideal``；
    父层 ``analyze_topology()`` 与 ``commit_subcycle_mass_flows_to_topology()`` 会清空 ``parent.non_ideal``，
    保证「ensure 时刻」的快照不被静默替换——若需重新对齐当前理想层，应重新 ``ensure_non_ideal()``。
    """

    simplified: SimplifiedTopology

    @classmethod
    def from_closed_cycle_layer(cls, layer) -> NonIdealClosedCycleLayer:
        """以理想层当前 ``simplified`` 构造快照容器；``layer.simplified`` 须已生成。"""
        if layer.simplified is None:
            raise RuntimeError(
                "ClosedCycleLayer.simplified 尚未生成；须先 analyze_topology() 或 commit_subcycle_mass_flows_to_topology()。"
            )
        return cls(simplified=layer.simplified)
