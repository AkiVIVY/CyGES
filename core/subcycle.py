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
    edge_flow_contributions: dict = field(default_factory=dict)

    def _ts_frame_nodes(self) -> dict[str, StateNode]:
        """
        按 T-S 图中的几何位置重建四角：
        - 左侧两点：S较小的两点，再按T分上下
        - 右侧两点：S较大的两点，再按T分上下
        """
        node_list = list(self.nodes_working.values())
        if len(node_list) != 4:
            raise ValueError(f"{self.subcycle_id}: subcycle must contain exactly 4 nodes")

        sorted_by_s = sorted(node_list, key=lambda n: (n.S, n.T))
        left_pair = sorted(sorted_by_s[:2], key=lambda n: n.T)
        right_pair = sorted(sorted_by_s[2:], key=lambda n: n.T)
        return {
            "left_bottom": left_pair[0],
            "left_top": left_pair[1],
            "right_bottom": right_pair[0],
            "right_top": right_pair[1],
        }

    def oriented_edges_for_flow(self) -> list[dict]:
        """
        根据子循环流量方向，给出四条边的有向流信息。
        约定：
        - m_dot > 0：在TS图上顺时针
        - m_dot < 0：在TS图上逆时针
        """
        n = self._ts_frame_nodes()
        flow_abs = abs(self.m_dot)
        clockwise_edges = [
            ("left", "left_bottom", "left_top"),
            ("top", "left_top", "right_top"),
            ("right", "right_top", "right_bottom"),
            ("bottom", "right_bottom", "left_bottom"),
        ]
        edges = clockwise_edges if self.m_dot >= 0 else list(reversed([
            (role, end_key, start_key) for role, start_key, end_key in clockwise_edges
        ]))

        result = []
        for role, start_key, end_key in edges:
            start_node = n[start_key]
            end_node = n[end_key]
            result.append({
                "role": role,
                "start_node_id": start_node.node_id,
                "end_node_id": end_node.node_id,
                "start_node": start_node,
                "end_node": end_node,
                "flow": flow_abs,
                "subcycle_id": self.subcycle_id,
            })
        self.edge_flow_contributions = {item["role"]: item for item in result}
        return result

    def calc_power(self) -> dict[str, float]:
        n = self._ts_frame_nodes()
        left_dh = n["left_top"].H - n["left_bottom"].H
        right_dh = n["right_top"].H - n["right_bottom"].H
        w_left = self.m_dot * left_dh
        w_right = self.m_dot * right_dh
        return {"left_power": w_left, "right_power": w_right, "net_power": w_right - w_left}

    def calc_heat(self) -> dict[str, float]:
        n = self._ts_frame_nodes()
        top_dh = n["right_top"].H - n["left_top"].H
        bottom_dh = n["right_bottom"].H - n["left_bottom"].H
        q_top = self.m_dot * top_dh
        q_bottom = self.m_dot * bottom_dh
        return {"heat_in_top": q_top, "heat_out_bottom": q_bottom}

