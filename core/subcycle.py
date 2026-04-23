from dataclasses import dataclass, field

from .process_edge import ProcessEdge
from .state_node import StateNode


@dataclass
class SubCycle:
    subcycle_id: str
    nodes_raw: dict[str, StateNode]
    nodes_working: dict[str, StateNode]
    edges_raw: dict[str, ProcessEdge]
    edges_working: dict[str, ProcessEdge]
    m_dot: float = 1.0
    m_dot_min: float = 0.0
    m_dot_max: float = 100.0
    metadata: dict = field(default_factory=dict)

    def _edge(self, role: str) -> ProcessEdge | None:
        return self.edges_working.get(role)

    def calc_power(self) -> dict[str, float]:
        left = self._edge("left")
        right = self._edge("right")
        w_left = left.qw_dot() if left else 0.0
        w_right = right.qw_dot() if right else 0.0
        return {"left_power": w_left, "right_power": w_right, "net_power": w_right - w_left}

    def calc_heat(self) -> dict[str, float]:
        top = self._edge("top")
        bottom = self._edge("bottom")
        q_top = top.qw_dot() if top else 0.0
        q_bottom = bottom.qw_dot() if bottom else 0.0
        return {"heat_in_top": q_top, "heat_out_bottom": q_bottom}

