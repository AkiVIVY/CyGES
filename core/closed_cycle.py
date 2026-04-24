import math
from dataclasses import dataclass, field

from config import (
    DEFAULT_EDGE_EFF,
    DEFAULT_SUBCYCLE_M_DOT,
    DEFAULT_SUBCYCLE_M_DOT_MAX,
    DEFAULT_SUBCYCLE_M_DOT_MIN,
    DEFAULT_TOLERANCE,
)

from .process_edge import ProcessEdge
from .state_node import StateNode
from .subcycle import SubCycle
from solvers import CoolPropSolver


@dataclass
class ClosedCycle:
    cycle_id: str   # 循环编号
    boundary: dict   # 边界条件
    levels: dict   # 分位点
    fluid: str   # 工质
    subcycles: dict[str, SubCycle] = field(default_factory=dict)    # 子循环
    nodes_raw: dict[str, StateNode] = field(default_factory=dict)   # 原始节点
    nodes_working: dict[str, StateNode] = field(default_factory=dict)   # 工作节点
    edges_raw: dict[str, ProcessEdge] = field(default_factory=dict)   # 原始边
    edges_working: dict[str, ProcessEdge] = field(default_factory=dict)   # 工作边
    node_groups_in_p: dict[float, list[StateNode]] = field(default_factory=dict)   # 压力分组
    node_groups_in_s: dict[float, list[StateNode]] = field(default_factory=dict)   # 熵分组
    topology_diagnostics: dict = field(default_factory=dict)   # 拓扑诊断
    _solver: CoolPropSolver | None = field(default=None, init=False, repr=False)   # 循环内复用的求解器

    def generate_topology(self) -> None:
        """
        根据边界与分位点生成闭式循环拓扑（节点、边、子循环）。
        """
        t_range = self.boundary.get("TRange")
        p_range = self.boundary.get("PRange")
        t_level_input = self.levels.get("TLevel")
        p_level_input = self.levels.get("PLevel")
        tolerance = self.boundary.get("tolerance", DEFAULT_TOLERANCE)

        if not isinstance(t_range, (list, tuple)) or len(t_range) != 2:
            raise ValueError(f"{self.cycle_id}: boundary.TRange must be [Tmin, Tmax]")
        if not isinstance(p_range, (list, tuple)) or len(p_range) != 2:
            raise ValueError(f"{self.cycle_id}: boundary.PRange must be [Pmin, Pmax]")
        if tolerance <= 0:
            raise ValueError(f"{self.cycle_id}: tolerance must be > 0")
        t_level = self._complete_levels(t_level_input, level_name="TLevel")
        p_level = self._complete_levels(p_level_input, level_name="PLevel")

        t_min, t_max = t_range
        p_min, p_max = p_range
        if t_min >= t_max or p_min >= p_max:
            raise ValueError(f"{self.cycle_id}: invalid boundary ranges")

        if self._solver is None:
            # 每个闭式循环只创建一次求解器，后续步骤复用。
            self._solver = CoolPropSolver([self.fluid])
        solver = self._solver

        p_list = [p_min + (p_max - p_min) * x for x in p_level]
        t_list = [t_min + (t_max - t_min) * x for x in t_level]

        # 步骤1：由TP参数对生成一级节点
        a_nodes: list[StateNode] = []
        for p_idx, px in enumerate(p_list):
            for t_idx, tx in enumerate(t_list):
                prop = solver.solve("TP", self.fluid, tx, px)
                a_nodes.append(
                    StateNode(
                        node_id=f"A_{p_idx}_{t_idx}",
                        fluid=self.fluid,
                        T=prop["T"],
                        P=prop["P"],
                        H=prop["H"],
                        S=prop["S"],
                        source_tag="A_TP",
                    )
                )

        # 步骤2：通过PS投影生成二级节点
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
                        node_id=f"{base.node_id}_B_{p_idx}",
                        fluid=self.fluid,
                        T=prop["T"],
                        P=prop["P"],
                        H=prop["H"],
                        S=prop["S"],
                        source_tag="B_PS",
                    )
                )

        # 步骤3：按(P,S)容差去重
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

        self.nodes_raw = node_raw   # 合并后的原始节点
        self.nodes_working = {nid: n.clone() for nid, n in self.nodes_raw.items()}   # 合并后的工作节点
        dedup_removed = len(a_nodes) + len(b_nodes) - len(self.nodes_raw)   # 去重节点数量

        node_list = list(self.nodes_working.values())
        group_in_p = self._group_nodes(node_list, key_attr="P", tolerance=tolerance)   # 压力分组，同时进行了容差处理
        group_in_s = self._group_nodes(node_list, key_attr="S", tolerance=tolerance)   # 熵分组，同时进行了容差处理
        self.node_groups_in_p = group_in_p
        self.node_groups_in_s = group_in_s

        # 步骤4：生成全局边集合
        self.edges_raw = {}
        self.edges_working = {}
        self._append_group_edges(group_in_p, sort_attr="S", family="P")
        self._append_group_edges(group_in_s, sort_attr="P", family="S")

        # 步骤5：提取子循环
        self.subcycles, invalid_subcycle_count = self._extract_subcycles(group_in_p)
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
        normalized = round(value / tolerance) * tolerance
        decimal_digits = max(0, int(math.ceil(-math.log10(tolerance))) + 2)
        return round(normalized, decimal_digits)

    @staticmethod
    def _complete_levels(levels: list | tuple | None, level_name: str) -> list[float]:
        """
        分位点输入仅要求“至少一个”，计算时自动补上0和1。
        """
        if not isinstance(levels, (list, tuple)) or len(levels) < 1:
            raise ValueError(f"{level_name} must have at least 1 value")
        result = [0.0, 1.0]
        for value in levels:
            value_f = float(value)
            if not 0 <= value_f <= 1:
                raise ValueError(f"{level_name} must be within [0, 1]")
            result.append(value_f)
        return sorted(set(result))

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
                    eff=DEFAULT_EDGE_EFF,
                    m_dot=0.0,
                )
                self.edges_raw[edge_id] = edge
                self.edges_working[edge_id] = edge

    def _extract_subcycles(
        self,
        group_in_p: dict[float, list[StateNode]],
    ) -> tuple[dict[str, SubCycle], int]:
        p_keys = sorted(group_in_p.keys())
        s_keys = sorted({
            self._normalized_key(n.S, self.boundary.get("tolerance", DEFAULT_TOLERANCE))
            for group in group_in_p.values()
            for n in group
        })
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
                subcycles[sub_id] = SubCycle(
                    subcycle_id=sub_id,
                    nodes_raw=nodes_raw,
                    nodes_working=nodes_working,
                    m_dot=DEFAULT_SUBCYCLE_M_DOT,
                    m_dot_min=DEFAULT_SUBCYCLE_M_DOT_MIN,
                    m_dot_max=DEFAULT_SUBCYCLE_M_DOT_MAX,
                    metadata={"p_band": (p_left, p_right), "s_band": (s_low, s_high)},
                )
        return subcycles, invalid_subcycle_count

