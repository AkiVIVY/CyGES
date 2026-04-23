from dataclasses import dataclass


@dataclass
class StateNode:
    node_id: str
    fluid: str
    T: float
    P: float
    H: float
    S: float
    source_tag: str = "unknown"
    grid_index: tuple[int, int] | None = None

    def clone(self, node_id: str | None = None) -> "StateNode":
        return StateNode(
            node_id=node_id or self.node_id,
            fluid=self.fluid,
            T=self.T,
            P=self.P,
            H=self.H,
            S=self.S,
            source_tag=self.source_tag,
            grid_index=self.grid_index,
        )

