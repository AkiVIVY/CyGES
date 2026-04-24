import sys
import random
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
                "levels": {"TLevel": [0.5], "PLevel": [0.5]},
            }
        ],
        "metadata": {"case": "smoke"},
    }


def main() -> None:
    spec = sample_spec()
    validate_system_spec(spec)
    model = SystemModel.from_spec(spec)
    model.build()
    inject_random_subcycle_flows(model.closed_cycles[0], seed=20260423)
    model.closed_cycles[0].generate_topology()
    report = model.solve()
    print("build-ok")
    print(f"closed_cycles: {len(model.closed_cycles)}")
    print(f"total_cycle_power: {report['total_cycle_power']}")
    print(f"random flows: {model.closed_cycles[0].boundary.get('subcycle_flows', {})}")
    plot_closed_cycle_ts(model.closed_cycles[0], Path(__file__).resolve().parents[1] / "smoke_cycle_ts.png")


def inject_random_subcycle_flows(cycle, seed: int = 0) -> None:
    """
    为当前闭式循环随机注入子循环流量（正负混合）。
    """
    rng = random.Random(seed)
    flow_map = {}
    for sub_id in cycle.subcycles.keys():
        val = 0.0
        while abs(val) < 0.2:
            val = rng.uniform(-2.0, 2.0)
        flow_map[sub_id] = round(val, 3)
    cycle.boundary["subcycle_flows"] = flow_map


def plot_closed_cycle_ts(cycle, output_path: Path, show_local_arrows: bool = False) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))

    # 1) 绘制节点
    nodes = list(cycle.nodes_working.values())
    ax.scatter([n.S for n in nodes], [n.T for n in nodes], s=24, c="black", zorder=3)

    # 2) 绘制全局边（子循环流量矢量叠加后的净方向与净流量）
    for edge in cycle.edge_flow_map.values():
        if abs(edge["net_flow"]) < 1e-9:
            # 净流量为0的边不画方向箭头，避免显示误导
            continue
        start_s, start_t = edge["start_node"].S, edge["start_node"].T
        end_s, end_t = edge["end_node"].S, edge["end_node"].T
        ax.annotate(
            "",
            xy=(end_s, end_t),
            xytext=(start_s, start_t),
            arrowprops=dict(arrowstyle="->", color="tab:gray", lw=1.0, alpha=0.6),
            zorder=2,
        )
        mid_s = (start_s + end_s) / 2.0
        mid_t = (start_t + end_t) / 2.0
        ax.text(mid_s, mid_t, f"{edge['net_flow']:.2f}", fontsize=6, color="dimgray")

    # 3) 绘制子循环区域
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

        if show_local_arrows:
            # 可选：子循环局部箭头（默认关闭，避免与全局净箭头叠加导致双向视觉）
            for contrib in sub.oriented_edges_for_flow():
                start = contrib["start_node"]
                end = contrib["end_node"]
                color = "tab:red" if contrib["role"] in {"left", "right"} else "tab:blue"
                ax.annotate(
                    "",
                    xy=(end.S, end.T),
                    xytext=(start.S, start.T),
                    arrowprops=dict(arrowstyle="->", color=color, lw=1.0, alpha=0.45),
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


