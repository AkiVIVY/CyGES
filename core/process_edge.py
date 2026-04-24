from dataclasses import dataclass
from typing import Literal

from .state_node import StateNode


EdgeType = Literal["compress", "expand", "heat_in", "heat_out", "unknown"]


@dataclass
class ProcessEdge:
    edge_id: str
    edge_type: EdgeType
    upstream: StateNode
    downstream: StateNode
    role: Literal["left", "right", "top", "bottom", "unknown"] = "unknown"
    eff: float = 1.0
    m_dot: float = 0.0
    constraints_ref: dict | None = None

    def delta_h(self) -> float:
        return self.downstream.H - self.upstream.H

    def qw_dot(self) -> float:
        return self.m_dot * self.delta_h()

