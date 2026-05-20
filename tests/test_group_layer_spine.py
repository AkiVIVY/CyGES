"""组内 layer：主脊分层 + 最长路径并列时最小起点 tie-break。"""

from core.closed_cycle_layer import SimplifiedEdge
from core.non_ideal_closed_cycle_layer import _compute_group_depth_metrics


def _edges(
    pairs: list[tuple[int, int]],
    *,
    kind: str = "heat",
) -> tuple[dict[str, SimplifiedEdge], frozenset[str]]:
    by_key: dict[str, SimplifiedEdge] = {}
    keys: list[str] = []
    for i, (tail, head) in enumerate(pairs, start=1):
        ek = f"SH{i}"
        by_key[ek] = SimplifiedEdge(
            kind=kind,
            tail=tail,
            head=head,
            constituent_edges=(ek,),
            merged_nodes=(),
            mass_flow=1.0,
        )
        keys.append(ek)
    return by_key, frozenset(keys)


def test_layer_merge_to_d_e_is_zero():
    """A→B→C, A→D, E→D：主脊 A-B-C，E 汇入 D=1 → E=0。"""
    by_key, keys = _edges([(1, 2), (2, 3), (1, 4), (5, 4)])
    m = _compute_group_depth_metrics(by_key, keys)
    assert m.layer == {1: 0, 2: 1, 3: 2, 4: 1, 5: 0}


def test_layer_merge_to_c_e_is_one():
    """A→B→C, A→D, E→C：主脊 A-B-C，E 汇入 C=2 → E=1。"""
    by_key, keys = _edges([(1, 2), (2, 3), (1, 4), (5, 3)])
    m = _compute_group_depth_metrics(by_key, keys)
    assert m.layer == {1: 0, 2: 1, 3: 2, 4: 1, 5: 1}


def test_longest_path_tie_break_min_start():
    """两条等长最长路径 A→C 与 B→C：应选起点 A（index 更小）为主脊。"""
    by_key, keys = _edges([(1, 3), (2, 3)])
    m = _compute_group_depth_metrics(by_key, keys)
    assert m.layer[1] == 0
    assert m.layer[2] == 0
    assert m.layer[3] == 1
