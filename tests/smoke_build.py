import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon

sys.path.append(str(Path(__file__).resolve().parents[1]))

from core import SystemModel
from inputs import validate_system_spec


def sample_spec() -> dict:
    return {
        "hot_streams": [
            {"stream_id": "hot-1", "fluid": "He", "m_dot": 2.0, "segments": []},
        ],
        "cold_streams": [
            {"stream_id": "cold-1", "fluid": "Water", "m_dot": 3.0, "segments": []},
        ],
        "closed_cycles": [
            {
                "cycle_id": "cycle-1",
                "fluid": "He",
                "boundary": {"TRange": [200, 1000], "PRange": [10000, 20000]},
                "levels": {"TLevel": [0, 0.5, 1], "PLevel": [0, 0.5, 1]},
            }
        ],
        "metadata": {"case": "smoke"},
    }


def main() -> None:
    spec = sample_spec()
    validate_system_spec(spec)
    model = SystemModel.from_spec(spec)
    model.build()
    report = model.solve()
    print("build-ok")
    print(f"closed_cycles: {len(model.closed_cycles)}")
    print(f"total_cycle_power: {report['total_cycle_power']}")
    plot_closed_cycle_ts(model.closed_cycles[0], Path(__file__).resolve().parents[1] / "smoke_cycle_ts.png")


def plot_closed_cycle_ts(cycle, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))

    # 1) 绘制节点
    nodes = list(cycle.nodes_working.values())
    ax.scatter([n.S for n in nodes], [n.T for n in nodes], s=24, c="black", zorder=3)

    # 2) 绘制全局边（按分组生成的边）
    for edge in cycle.edges_working.values():
        ax.plot(
            [edge.upstream.S, edge.downstream.S],
            [edge.upstream.T, edge.downstream.T],
            color="tab:gray",
            linewidth=0.8,
            alpha=0.45,
            zorder=2,
        )

    # 3) 绘制子循环（绿色填充）以及子循环四条语义边
    for sub in cycle.subcycles.values():
        n = sub.nodes_working
        polygon_points = [
            (n["left_bottom"].S, n["left_bottom"].T),
            (n["left_top"].S, n["left_top"].T),
            (n["right_top"].S, n["right_top"].T),
            (n["right_bottom"].S, n["right_bottom"].T),
        ]
        ax.add_patch(
            Polygon(
                polygon_points,
                closed=True,
                facecolor="tab:green",
                edgecolor="none",
                alpha=0.12,
                zorder=1,
            )
        )
        center_s = sum(p[0] for p in polygon_points) / 4.0
        center_t = sum(p[1] for p in polygon_points) / 4.0
        ax.text(center_s, center_t, sub.subcycle_id, fontsize=6, color="darkgreen")

        for edge in sub.edges_working.values():
            ax.plot(
                [edge.upstream.S, edge.downstream.S],
                [edge.upstream.T, edge.downstream.T],
                color="tab:red" if edge.role in {"left", "right"} else "tab:blue",
                linewidth=1.4,
                alpha=0.85,
                zorder=4,
            )

    ax.set_title(f"Smoke Closed Cycle Topology (T-S): {cycle.cycle_id}")
    ax.set_xlabel("Entropy S [kJ/(kg.K)]")
    ax.set_ylabel("Temperature T [K]")
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="w", markerfacecolor="black", markersize=6, label="nodes"),
            Line2D([0], [0], color="tab:gray", lw=1.2, alpha=0.6, label="global edges"),
            Line2D([0], [0], color="tab:red", lw=1.6, label="left/right edges"),
            Line2D([0], [0], color="tab:blue", lw=1.6, label="top/bottom edges"),
            Line2D([0], [0], color="tab:green", lw=6, alpha=0.2, label="subcycles"),
        ],
        loc="best",
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=170)
    plt.close(fig)
    print(f"plot saved: {output_path.name}")


if __name__ == "__main__":
    main()

