from createModel import workFlow
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


def build_workflow_with_helium():
    """
    依次尝试用户指定工质名，返回可用的workFlow实例和实际工质名。
    """
    candidate_names = ["He", "helium", "Helium"]
    last_error = None
    for name in candidate_names:
        try:
            wf = workFlow([name], input=None)
            return wf, name
        except Exception as exc:
            last_error = exc
    raise RuntimeError(
        "无法初始化氦气工质，请检查CoolProp工质名。"
        f"尝试过: {candidate_names}, 最后错误: {last_error}"
    )


def main():
    wf, substance = build_workflow_with_helium()

    # createModel当前压力单位为kPa，这里把用户给出的MPa换算为kPa
    p_range_kpa = [10_000, 20_000]  # 10-20 MPa
    t_range_k = [200, 1000]         # 200-1000 K
    levels = [0, 1 / 3, 2 / 3, 1]

    result = wf.generateNodesAndEdges(
        TRange=t_range_k,
        PRange=p_range_kpa,
        TLevel=levels,
        PLevel=levels,
        substance=substance,
        tolerance=1e-6,
    )

    print(f"substance used: {substance}")
    print(f"node count: {len(result['nodeList'])}")
    print(f"isobaric edge count: {len(result['edgeListInP'])}")
    print(f"isentropic edge count: {len(result['edgeListInS'])}")
    print(f"subCycle count: {len(result['subCycleDict']['subCycleList'])}")

    if result["nodeList"]:
        n0 = result["nodeList"][0]
        print(
            "sample node -> "
            f"T={n0.T:.3f} K, P={n0.P:.3f} kPa, H={n0.H:.3f} kJ/kg, S={n0.S:.6f} kJ/(kg.K)"
        )

    if result["edgeListInP"]:
        e0 = result["edgeListInP"][0]
        print(
            "sample isobaric edge -> "
            f"dT={e0['deltaT']:.3f} K, dH={e0['deltaH']:.3f} kJ/kg"
        )

    if result["edgeListInS"]:
        e1 = result["edgeListInS"][0]
        print(
            "sample isentropic edge -> "
            f"dT={e1['deltaT']:.3f} K, dH={e1['deltaH']:.3f} kJ/kg"
        )

    # 在T-S图上可视化点和边
    fig, ax = plt.subplots(figsize=(10, 6))

    # 节点散点
    s_values = [n.S for n in result["nodeList"]]
    t_values = [n.T for n in result["nodeList"]]
    ax.scatter(s_values, t_values, s=36, c="black", label="nodes", zorder=3)

    # 等压边（红色）
    for edge in result["edgeListInP"]:
        n1 = edge["from_node"]
        n2 = edge["to_node"]
        ax.plot([n1.S, n2.S], [n1.T, n2.T], color="tab:red", linewidth=1.4, alpha=0.85)

    # 等熵边（蓝色）
    for edge in result["edgeListInS"]:
        n1 = edge["from_node"]
        n2 = edge["to_node"]
        ax.plot([n1.S, n2.S], [n1.T, n2.T], color="tab:blue", linewidth=1.2, alpha=0.75)

    # subCycle（子循环）半透明填充并标注编号
    for sub_cycle in result["subCycleDict"]["subCycleList"]:
        nodes = sub_cycle["nodes"]
        polygon_points = [
            (nodes["left_bottom"].S, nodes["left_bottom"].T),
            (nodes["left_top"].S, nodes["left_top"].T),
            (nodes["right_top"].S, nodes["right_top"].T),
            (nodes["right_bottom"].S, nodes["right_bottom"].T),
        ]
        patch = Polygon(
            polygon_points,
            closed=True,
            facecolor="tab:green",
            edgecolor="none",
            alpha=0.12,
            zorder=1,
        )
        ax.add_patch(patch)

        center_s = sum([p[0] for p in polygon_points]) / 4
        center_t = sum([p[1] for p in polygon_points]) / 4
        ax.text(center_s, center_t, str(sub_cycle["id"]), fontsize=7, color="darkgreen", alpha=0.9)

    ax.set_title(f"CyGES Nodes and Edges on T-S Plane ({substance})")
    ax.set_xlabel("Entropy S [kJ/(kg.K)]")
    ax.set_ylabel("Temperature T [K]")
    ax.grid(True, linestyle="--", alpha=0.4)

    # 手动添加图例项，避免重复
    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="black", markersize=6, label="nodes"),
        Line2D([0], [0], color="tab:red", lw=1.6, label="isobaric edges"),
        Line2D([0], [0], color="tab:blue", lw=1.6, label="isentropic edges"),
        Line2D([0], [0], color="tab:green", lw=6, alpha=0.2, label="subCycle"),
    ]
    ax.legend(handles=legend_handles, loc="best")

    out_path = "nodes_edges_ts_dense4.png"
    plt.tight_layout()
    plt.savefig(out_path, dpi=160)
    print(f"plot saved: {out_path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
