r"""CMA-ES HX匹配全参数优化: 2P+1T拓扑下联合优化参数+基函数权重, CoolProp物性。

  搜索空间: t_min/t_max/t_q/p_q + H2流量 + 4×4 高斯基权重,
  目标 hx_unmatched(未匹配功率比), CMA-ES 30代×3重启(σ早停),
  输出 4 张图: TS/PS 状态空间 + HX T-Q概览 + 每配对单元T-Q。
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import pytest

from core import (
    ClosedCycleLayer,
    ClosedCycleTPInput,
    CycleConfig,
    ExternalSourceInput,
    ProcessCategory,
    PropertyRegistry,
    SystemInput,
)
from core.heat_exchanger import match_heat_exchanger_groups
from optimize import Optimizer

TESTS_DIR = Path(__file__).resolve().parent


def _make_system_input(h2_mf: float = 4.3) -> SystemInput:
    hot = ExternalSourceInput(
        fluid="Air", mass_flow=100.0,
        T_in=1250.0, P_in=200.0, T_out=500.0, P_out=180.0,
    )
    cold = ExternalSourceInput(
        fluid="Hydrogen", mass_flow=h2_mf,
        T_in=20.0, P_in=5000.0, T_out=900.0, P_out=4500.0,
    )
    cycle_input = ClosedCycleTPInput(
        fluid="He", t_min=40.0, t_max=1000.0, p_min=2000.0, p_max=10000.0,
        t_quantiles=(0.5,), p_quantiles=(0.33, 0.67),
        subcycle_mass_flow_initial=20.0,
    )
    return SystemInput(
        heat_sources=(hot,), cold_sources=(cold,),
        cycles=(CycleConfig(
            input=cycle_input, use_non_ideal=True,
            delta_T_min=20.0, heat_method=None,
        ),),
        delta_T_min=20.0, heat_method="system_pinch",
    )


def test_cma_full_opt() -> None:
    sys_inp = _make_system_input()
    props = PropertyRegistry()
    hx_dT = 10.0

    run_dir = TESTS_DIR / "run_cma_full"
    run_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 60)
    print("CMA-ES 全参数优化 (2P+1T, H2 3-6, 4x4 basis, CoolProp)")
    print("=" * 60)

    opt = Optimizer(
        base_input=sys_inp, props=props, objective="hx_unmatched",
        mf_step_fraction=0.01, quantile_step=0.01,
        basis_encoding=True, basis_s=4, basis_p=4,
        mf_bounds=(0.0, 50.0), quantile_merge_ratio=0.0,
        n_workers=6, skip_pinch=True, hx_dT_min=hx_dT, hx_max_group_size=3,
        h2_mf_bounds=(3.0, 6.0),
    )

    dim = len(opt.bounds)
    print(f"  维度: {dim} (t_max+t_min+{opt._n_t_q}t_q+{opt._n_p_q}p_q+{opt._mf_dim}basis+1 H2)")
    print(f"  CMA-ES maxiter=30 restarts=3 seed=42, {opt._n_workers} workers")
    print(f"  dT_min={opt._hx_dT_min}K, max_group_size=3, CoolProp\n")

    best_val_log: list[float] = []

    def _on_gen(gen: int, restart: int, best_x: list[float],
                best_val: float, n_evals: int) -> None:
        best_val_log.append(best_val)
        idx = 0
        t_mx, t_mn = best_x[idx], best_x[idx + 1]; idx += 2
        tq = ", ".join(f"{best_x[idx+i]:.3f}" for i in range(opt._n_t_q))
        idx += opt._n_t_q
        pq = ", ".join(f"{best_x[idx+i]:.3f}" for i in range(opt._n_p_q))
        idx += opt._n_p_q
        wm = sum(best_x[idx:idx + opt._mf_dim]) / max(1, opt._mf_dim)
        h2 = best_x[idx + opt._mf_dim] if len(best_x) > idx + opt._mf_dim else 0
        sys.stdout.write(
            f"  gen {gen:3d} | r{restart} | obj={best_val:.5f} | "
            f"t_max={t_mx:.0f} t_min={t_mn:.0f} | "
            f"t_q=[{tq}] p_q=[{pq}] | w_mean={wm:.1f} H2={h2:.2f} | evals={n_evals}\n"
        )
        sys.stdout.flush()

    print("  开始优化...\n")
    t0 = time.time()
    opt_result = opt.run(
        method="cma", maxiter=30, seed=42, early_stop=5, sigma0=0.3,
        restarts=3, callback=_on_gen,
    )
    elapsed = time.time() - t0

    raw = opt_result.system_result

    # 计时
    hots_tmp, colds_tmp = [], []
    for rec in raw.heat_source_records:
        hots_tmp.append(rec)
    for rec in raw.cold_source_records:
        colds_tmp.append(rec)
    for report in raw.cycle_reports:
        for _, rec in report.by_edge:
            if rec.kind == "heat" and rec.power_rate is not None:
                if rec.category == ProcessCategory.HEAT_REJECTION:
                    hots_tmp.append(rec)
                else:
                    colds_tmp.append(rec)

    _t0 = time.time()
    for _ in range(10):
        match_heat_exchanger_groups(hots_tmp, colds_tmp, dT_min=hx_dT, max_group_size=3)
    t_hx_ms = (time.time() - _t0) / 10 * 1000
    print(f"\n  计时: {len(hots_tmp)+len(colds_tmp)}条, 10次HX={t_hx_ms*10:.0f}ms → 单次{t_hx_ms:.2f}ms")
    print(f"  优化总耗时: {elapsed:.1f}s, {opt_result.n_evaluations}次 → {elapsed/opt_result.n_evaluations*1000:.0f}ms/eval")

    # 最优解
    print("\n" + "=" * 55)
    print("最优解")
    print("=" * 55)
    idx = 0
    t_max_o = opt_result.x_opt[idx]; idx += 1
    t_min_o = opt_result.x_opt[idx]; idx += 1
    tq_o = [f"{opt_result.x_opt[idx+i]:.4f}" for i in range(opt._n_t_q)]
    idx += opt._n_t_q
    pq_o = [f"{opt_result.x_opt[idx+i]:.4f}" for i in range(opt._n_p_q)]
    idx += opt._n_p_q
    h2_o = opt_result.x_opt[idx + opt._mf_dim] if len(opt_result.x_opt) > idx + opt._mf_dim else 4.3
    print(f"  t_max / t_min = {t_max_o:.1f} / {t_min_o:.1f} K")
    print(f"  t_q = [{', '.join(tq_o)}]")
    print(f"  p_q = [{', '.join(pq_o)}]")
    print(f"  H2 = {h2_o:.2f} kg/s")
    print(f"  obj = {opt_result.objective:.5f}")
    ct = raw.cycle_reports[0].cycle_totals
    print(f"  net_mech = {ct.net_mechanical_power:.1f} kW")
    print(f"  n_evals = {opt_result.n_evaluations}")

    # HX
    hots, colds = [], []
    for rec in raw.heat_source_records:
        hots.append(rec)
    for rec in raw.cold_source_records:
        colds.append(rec)
    for report in raw.cycle_reports:
        for _, rec in report.by_edge:
            if rec.kind == "heat" and rec.power_rate is not None:
                if rec.category == ProcessCategory.HEAT_REJECTION:
                    hots.append(rec)
                else:
                    colds.append(rec)
    Q_h = sum(abs(float(r.power_rate)) for r in hots if r.power_rate)
    Q_c = sum(abs(float(r.power_rate)) for r in colds if r.power_rate)
    print(f"\n  放热={Q_h:.0f}kW  吸热={Q_c:.0f}kW  不平衡={abs(Q_h-Q_c):.0f}kW")
    print(f"  hot({len(hots)}): {[r.edge_key for r in hots]}")
    print(f"  cold({len(colds)}): {[r.edge_key for r in colds]}")

    hx_result = match_heat_exchanger_groups(hots, colds, dT_min=hx_dT, max_group_size=3)
    print(f"\n  HX(dT_min={hx_dT}K): {hx_result.num_units}组")
    print(f"  matched={hx_result.total_matched:.0f}kW  unmatched={hx_result.total_unmatched:.0f}kW  "
          f"ratio={hx_result.total_unmatched/(Q_h+Q_c):.4f}")
    N_ua = len(hx_result.unassigned_hots) + len(hx_result.unassigned_colds)
    print(f"  原始匹配率 = {hx_result.total_unmatched/(Q_h+Q_c):.4f}  (N+1) = {hx_result.total_unmatched/(Q_h+Q_c)*(N_ua+1):.4f}")
    for ui, unit in enumerate(hx_result.units):
        h_l = ", ".join(r.edge_key for r in unit.hot_records)
        c_l = ", ".join(r.edge_key for r in unit.cold_records)
        print(f"    U{ui+1}: hot=[{h_l}] cold=[{c_l}] matched={unit.matched_heat:.0f}kW resid={unit.residual:.1f}kW pinch={unit.internal_pinch:.0f}K")
    if hx_result.unassigned_hots:
        print(f"    未分配热: {[r.edge_key for r in hx_result.unassigned_hots]}")
    if hx_result.unassigned_colds:
        print(f"    未分配冷: {[r.edge_key for r in hx_result.unassigned_colds]}")

    # 收敛曲线
    if best_val_log:
        matplotlib = pytest.importorskip("matplotlib")
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(range(1, len(best_val_log)+1), best_val_log, "o-", ms=2, lw=1, color="tab:blue")
        ax.set_xlabel("Generation"); ax.set_ylabel("hx_unmatched")
        ax.set_title("CMA-ES 30×3, 2P+1T, 4×4 basis, H2 3-6, CoolProp")
        ax.grid(True, alpha=0.25); ax.set_ylim(bottom=0)
        path = run_dir / "convergence.png"
        fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
        print(f"\n  收敛曲线: {path}")

    assert opt_result.n_evaluations > 0


if __name__ == "__main__":
    test_cma_full_opt()
