from dataclasses import dataclass, field

from .closed_cycle import ClosedCycle
from .external_stream import ExternalStream


@dataclass
class SystemModel:
    hot_streams: list[ExternalStream] = field(default_factory=list)
    cold_streams: list[ExternalStream] = field(default_factory=list)
    closed_cycles: list[ClosedCycle] = field(default_factory=list)
    metadata: dict = field(default_factory=dict)

    @classmethod
    def from_spec(cls, spec: dict) -> "SystemModel":
        hot_streams = [
            ExternalStream(
                stream_id=item["stream_id"],
                stream_type="hot",
                fluid=item["fluid"],
                m_dot=item["m_dot"],
                segments=item.get("segments", []),
            )
            for item in spec.get("hot_streams", [])
        ]
        cold_streams = [
            ExternalStream(
                stream_id=item["stream_id"],
                stream_type="cold",
                fluid=item["fluid"],
                m_dot=item["m_dot"],
                segments=item.get("segments", []),
            )
            for item in spec.get("cold_streams", [])
        ]
        closed_cycles = [
            ClosedCycle(
                cycle_id=item["cycle_id"],
                boundary=item["boundary"],
                levels=item["levels"],
                fluid=item["fluid"],
            )
            for item in spec.get("closed_cycles", [])
        ]
        return cls(
            hot_streams=hot_streams,
            cold_streams=cold_streams,
            closed_cycles=closed_cycles,
            metadata=spec.get("metadata", {}),
        )

    def build(self) -> None:
        for cycle in self.closed_cycles:
            cycle.generate_topology()

    def solve(self) -> dict[str, float]:
        cycle_reports = []
        for cycle in self.closed_cycles:
            report = cycle.calc_cycle_balance()
            report["cycle_id"] = cycle.cycle_id
            report["topology"] = cycle.topology_diagnostics
            cycle_reports.append(report)
        return {
            "cycle_reports": cycle_reports,
            "total_cycle_power": sum(item["cycle_power"] for item in cycle_reports),
            "total_cycle_heat_in": sum(item["cycle_heat_in"] for item in cycle_reports),
            "total_cycle_heat_out": sum(item["cycle_heat_out"] for item in cycle_reports),
        }

