from dataclasses import dataclass, field

from config import DEFAULT_SUBCYCLE_M_DOT, DEFAULT_SUBCYCLE_M_DOT_MAX, DEFAULT_SUBCYCLE_M_DOT_MIN, DEFAULT_TOLERANCE

from .process_edge import ProcessEdge
from .state_node import StateNode
from .subcycle import SubCycle
from solvers import CoolPropSolver


@dataclass
class ClosedCycle:
    cycle_id: str
    boundary: dict
    levels: dict
    fluid: str
    subcycles: dict[str, SubCycle] = field(default_factory=dict)
    nodes_raw: dict[str, StateNode] = field(default_factory=dict)
    nodes_working: dict[str, StateNode] = field(default_factory=dict)
    edges_raw: dict[str, ProcessEdge] = field(default_factory=dict)
    edges_working: dict[str, ProcessEdge] = field(default_factory=dict)
    node_groups_in_p: dict[float, list[StateNode]] = field(default_factory=dict)
    node_groups_in_s: dict[float, list[StateNode]] = field(default_factory=dict)
    topology_diagnostics: dict = field(default_factory=dict)

    def generate_topology(self) -> None:
        """
        根据边界与分位点生成闭式循环拓扑（节点、边、子循环）。
        """
        t_range = self.boundary.get("TRange")
        p_range = self.boundary.get("PRange")
        t_level = self.levels.get("TLevel")
        p_level = self.levels.get("PLevel")
        tolerance = self.boundary.get("tolerance", DEFAULT_TOLERANCE)

        if not isinstance(t_range, (list, tuple)) or len(t_range) != 2:
            raise ValueError(f"{self.cycle_id}: boundary.TRange must be [Tmin, Tmax]")
        if not isinstance(p_range, (list, tuple)) or len(p_range) != 2:
            raise ValueError(f"{self.cycle_id}: boundary.PRange must be [Pmin, Pmax]")
        if not isinstance(t_level, (list, tuple)) or len(t_level) < 2:
            raise ValueError(f"{self.cycle_id}: levels.TLevel must have at least 2 values")
        if not isinstance(p_level, (list, tuple)) or len(p_level) < 2:
            raise ValueError(f"{self.cycle_id}: levels.PLevel must have at least 2 values")

        t_min, t_max = t_range
        p_min, p_max = p_range
        if t_min >= t_max or p_min >= p_max:
            raise ValueError(f"{self.cycle_id}: invalid boundary ranges")

        solver = CoolPropSolver([self.fluid])

        p_list = [p_min + (p_max - p_min) * x for x in p_level]
        t_list = [t_min + (t_max - t_min) * x for x in t_level]

        # Step 1: generate A nodes from TP
        a_nodes: list[StateNode] = []
        for p_idx, px in enumerate(p_list):
            for t_idx, tx in enumerate(t_list):
                prop = solver.solve("TP", self.fluid, tx, px)
                a_nodes.append(
                    StateNode(
                        node_id=f"N_A_{p_idx}_{t_idx}",
                        fluid=self.fluid,
                        T=prop["T"],
                        P=prop["P"],
                        H=prop["H"],
                        S=prop["S"],
                        source_tag="A_TP",
                        grid_index=(p_idx, t_idx),
                    )
                )

        # Step 2: generate B nodes by PS projection
        s_low_list = [solver.solve("TP", self.fluid, t_min, px)["S"] for px in p_list]
        s_high_list = [solver.solve("TP", self.fluid, t_max, px)["S"] for px in p_list]
        b_nodes: list[StateNode] = []
        for base in a_nodes:
            for p_idx, px in enumerate(p_list):
                if not (s_low_list[p_idx] - tolerance <= base.S <= s_high_list[p_idx] + tolerance):
                    continue
                if abs(base.P - px) <= tolerance:
                    continue
                prop = solver.solve("PS", self.fluid, px, base.S)
                b_nodes.append(
                    StateNode(
                        node_id=f"N_B_{base.node_id}_{p_idx}",
                        fluid=self.fluid,
                        T=prop["T"],
                        P=prop["P"],
                        H=prop["H"],
                        S=prop["S"],
                        source_tag="B_PS",
                        grid_index=None,
                    )
                )

        # Step 3: deduplicate with tolerance on (P,S)
        node_raw = {}
        seen = {}
        for node in a_nodes + b_nodes:
            p_key = self._normalized_key(node.P, tolerance)
            s_key = self._normalized_key(node.S, tolerance)
            key = (p_key, s_key)
            if key in seen:
                continue
            seen[key] = node.node_id
            node_raw[node.node_id] = node

        self.nodes_raw = node_raw
        self.nodes_working = {nid: n.clone() for nid, n in self.nodes_raw.items()}
        dedup_removed = len(a_nodes) + len(b_nodes) - len(self.nodes_raw)

        node_list = list(self.nodes_working.values())
        group_in_p = self._group_nodes(node_list, key_attr="P", tolerance=tolerance)
        group_in_s = self._group_nodes(node_list, key_attr="S", tolerance=tolerance)
        self.node_groups_in_p = group_in_p
        self.node_groups_in_s = group_in_s

        # Step 4: build global edge registry
        self.edges_raw = {}
        self.edges_working = {}
        self._append_group_edges(group_in_p, sort_attr="S", family="P")
        self._append_group_edges(group_in_s, sort_attr="P", family="S")

        # Step 5: extract subcycles
        self.subcycles, invalid_subcycle_count = self._extract_subcycles(group_in_p, group_in_s)
        self.topology_diagnostics = self._snapshot_topology_diagnostics(
            n_nodes_a=len(a_nodes),
            n_nodes_b=len(b_nodes),
            n_nodes_dedup_removed=dedup_removed,
            n_invalid_subcycles=invalid_subcycle_count,
            n_pressure_groups=len(group_in_p),
            n_entropy_groups=len(group_in_s),
        )
        return

    def calc_cycle_balance(self) -> dict[str, float]:
        total_power = 0.0
        total_heat_in = 0.0
        total_heat_out = 0.0
        for sub in self.subcycles.values():
            p = sub.calc_power()
            q = sub.calc_heat()
            total_power += p["net_power"]
            total_heat_in += q["heat_in_top"]
            total_heat_out += q["heat_out_bottom"]
        return {
            "cycle_power": total_power,
            "cycle_heat_in": total_heat_in,
            "cycle_heat_out": total_heat_out,
        }

    def _snapshot_topology_diagnostics(
        self,
        n_nodes_a: int = 0,
        n_nodes_b: int = 0,
        n_nodes_dedup_removed: int = 0,
        n_invalid_subcycles: int = 0,
        n_pressure_groups: int = 0,
        n_entropy_groups: int = 0,
    ) -> dict:
        return {
            "n_nodes_a": n_nodes_a,
            "n_nodes_b": n_nodes_b,
            "n_nodes_raw": len(self.nodes_raw),
            "n_nodes_working": len(self.nodes_working),
            "n_nodes_dedup_removed": n_nodes_dedup_removed,
            "n_edges_raw": len(self.edges_raw),
            "n_edges_working": len(self.edges_working),
            "n_pressure_groups": n_pressure_groups,
            "n_entropy_groups": n_entropy_groups,
            "n_subcycles": len(self.subcycles),
            "n_invalid_subcycles": n_invalid_subcycles,
        }

    @staticmethod
    def _normalized_key(value: float, tolerance: float) -> float:
        return round(round(value / tolerance) * tolerance, 12)

    def _group_nodes(self, node_list: list[StateNode], key_attr: str, tolerance: float) -> dict[float, list[StateNode]]:
        grouped: dict[float, list[StateNode]] = {}
        for node in node_list:
            raw_val = getattr(node, key_attr)
            key = self._normalized_key(raw_val, tolerance)
            grouped.setdefault(key, []).append(node)
        return grouped

    def _append_group_edges(self, grouped_nodes: dict[float, list[StateNode]], sort_attr: str, family: str) -> None:
        for key, group in grouped_nodes.items():
            if len(group) < 2:
                continue
            sorted_nodes = sorted(group, key=lambda n: getattr(n, sort_attr))
            for i in range(len(sorted_nodes) - 1):
                up = sorted_nodes[i]
                down = sorted_nodes[i + 1]
                edge_id = f"E_{family}_{str(key).replace('.', '_')}_{i}"
                edge = ProcessEdge(
                    edge_id=edge_id,
                    edge_type="unknown",
                    upstream=up,
                    downstream=down,
                    role="unknown",
                    eta=1.0,
                    m_dot=0.0,
                )
                self.edges_raw[edge_id] = edge
                self.edges_working[edge_id] = edge

    def _extract_subcycles(
        self,
        group_in_p: dict[float, list[StateNode]],
        group_in_s: dict[float, list[StateNode]],
    ) -> tuple[dict[str, SubCycle], int]:
        p_keys = sorted(group_in_p.keys())
        s_keys = sorted(group_in_s.keys())
        node_map = {}
        for p_key, group in group_in_p.items():
            for n in group:
                s_key = self._normalized_key(n.S, self.boundary.get("tolerance", DEFAULT_TOLERANCE))
                node_map[(p_key, s_key)] = n

        subcycles: dict[str, SubCycle] = {}
        invalid_subcycle_count = 0
        for p_idx in range(len(p_keys) - 1):
            p_left = p_keys[p_idx]
            p_right = p_keys[p_idx + 1]
            for s_idx in range(len(s_keys) - 1):
                s_low = s_keys[s_idx]
                s_high = s_keys[s_idx + 1]
                n_lb = node_map.get((p_left, s_low))
                n_lt = node_map.get((p_left, s_high))
                n_rt = node_map.get((p_right, s_high))
                n_rb = node_map.get((p_right, s_low))
                if not all([n_lb, n_lt, n_rt, n_rb]):
                    invalid_subcycle_count += 1
                    continue

                sub_id = f"SC_{p_idx}_{s_idx}"
                nodes_raw = {
                    "left_bottom": n_lb,
                    "left_top": n_lt,
                    "right_top": n_rt,
                    "right_bottom": n_rb,
                }
                nodes_working = {k: v.clone() for k, v in nodes_raw.items()}

                edges_raw = self._make_subcycle_edges(sub_id, nodes_working)
                edges_working = {k: e for k, e in edges_raw.items()}
                subcycles[sub_id] = SubCycle(
                    subcycle_id=sub_id,
                    nodes_raw=nodes_raw,
                    nodes_working=nodes_working,
                    edges_raw=edges_raw,
                    edges_working=edges_working,
                    m_dot=DEFAULT_SUBCYCLE_M_DOT,
                    m_dot_min=DEFAULT_SUBCYCLE_M_DOT_MIN,
                    m_dot_max=DEFAULT_SUBCYCLE_M_DOT_MAX,
                    metadata={"p_band": (p_left, p_right), "s_band": (s_low, s_high)},
                )
        return subcycles, invalid_subcycle_count

    @staticmethod
    def _make_subcycle_edges(sub_id: str, nodes_working: dict[str, StateNode]) -> dict[str, ProcessEdge]:
        edge_map = {
            "left": ("compress", "left_bottom", "left_top"),
            "right": ("expand", "right_bottom", "right_top"),
            "top": ("heat_in", "left_top", "right_top"),
            "bottom": ("heat_out", "left_bottom", "right_bottom"),
        }
        result = {}
        for role, (edge_type, up_key, down_key) in edge_map.items():
            result[role] = ProcessEdge(
                edge_id=f"{sub_id}_{role}",
                edge_type=edge_type,
                upstream=nodes_working[up_key],
                downstream=nodes_working[down_key],
                role=role,
                eta=1.0,
                m_dot=DEFAULT_SUBCYCLE_M_DOT,
            )
        return result

    def export_topology_dict(self) -> dict:
        """
        兼容旧版createModel输出风格，便于迁移阶段做结果对比。
        """
        node_list = list(self.nodes_working.values())
        edge_list_in_p = [
            edge for edge in self.edges_working.values()
            if edge.edge_id.startswith("E_P_")
        ]
        edge_list_in_s = [
            edge for edge in self.edges_working.values()
            if edge.edge_id.startswith("E_S_")
        ]
        sub_cycle_list = []
        sub_cycle_in_p_band = {}
        sub_cycle_in_s_band = {}
        for sub in self.subcycles.values():
            sub_item = {
                "id": sub.subcycle_id,
                "nodes": sub.nodes_working,
                "p_band": sub.metadata.get("p_band"),
                "s_band": sub.metadata.get("s_band"),
            }
            sub_cycle_list.append(sub_item)
            p_band = sub_item["p_band"]
            s_band = sub_item["s_band"]
            sub_cycle_in_p_band.setdefault(p_band, []).append(sub_item)
            sub_cycle_in_s_band.setdefault(s_band, []).append(sub_item)

        return {
            "nodeList": node_list,
            "nodeListInP": self.node_groups_in_p,
            "nodeListInS": self.node_groups_in_s,
            "edgeListInP": edge_list_in_p,
            "edgeListInS": edge_list_in_s,
            "subCycleDict": {
                "subCycleList": sub_cycle_list,
                "subCycleInPBand": sub_cycle_in_p_band,
                "subCycleInSBand": sub_cycle_in_s_band,
            },
            "topologyDiagnostics": self.topology_diagnostics,
        }

