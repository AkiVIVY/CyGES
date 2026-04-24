from dataclasses import dataclass, field

from .state_node import StateNode


@dataclass
class SubCycle:
    subcycle_id: str
    nodes_raw: dict[str, StateNode]
    nodes_working: dict[str, StateNode]
    m_dot: float = 1.0
    m_dot_min: float = 0.0
    m_dot_max: float = 100.0
    metadata: dict = field(default_factory=dict)

    def calc_power(self) -> dict[str, float]:
        n = self.nodes_working
        left_dh = n["left_top"].H - n["left_bottom"].H
        right_dh = n["right_top"].H - n["right_bottom"].H
        w_left = self.m_dot * left_dh
        w_right = self.m_dot * right_dh
        return {"left_power": w_left, "right_power": w_right, "net_power": w_right - w_left}

    def calc_heat(self) -> dict[str, float]:
        n = self.nodes_working
        top_dh = n["right_top"].H - n["left_top"].H
        bottom_dh = n["right_bottom"].H - n["left_bottom"].H
        q_top = self.m_dot * top_dh
        q_bottom = self.m_dot * bottom_dh
        return {"heat_in_top": q_top, "heat_out_bottom": q_bottom}

