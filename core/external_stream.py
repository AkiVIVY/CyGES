from dataclasses import dataclass, field
from typing import Literal


StreamType = Literal["hot", "cold"]


@dataclass
class ExternalStream:
    stream_id: str
    stream_type: StreamType
    fluid: str
    m_dot: float
    segments: list[dict] = field(default_factory=list)

    def to_tq_segments(self) -> list[dict]:
        return self.segments

