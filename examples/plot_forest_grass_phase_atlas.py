"""Japanese phase-diagram atlas; existing forest--grass model, no new ecology.

All rates use model time. N plots an added reflected diffusion; R uses the
already-proved constructed path, not a calibrated environmental response.
"""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from control_carbon.forest_grass_tipping import (
    Parameters, logistic, phi, h, drift, drift_prime, interior_roots, fold_points,
    path_parameters, smooth_step, critical_duration, matched_solutions,
    ramp_solution, potential, mean_first_passage,
)

GREEN, ORANGE, PURPLE, INK = "#147d64", "#cb851b", "#8451ad", "#344250"
PALE_GREEN, PALE_ORANGE = "#d8eee5", "#fae9cc"
BASE = Parameters()


def phase_class(p):
    """Count stable interior equilibria AND the g=1 boundary equilibrium."""
    roots = interior_roots(p)
    stable = roots[drift_prime(roots, p) < 0]
    boundary_stable = h(1, p) > 0
    if len(stable) == 2:
        assert not boundary_stable
        return 2
    if len(stable) == 1:
        return 3 if boundary_stable else 1
    assert boundary_stable
    return 0


def style():
    font = Path("/usr/share/fonts/opentype/ipaexfont-gothic/ipaexg.ttf")
    if font.exists():
        font_manager.fontManager.addfont(str(font))
        family = font_manager.FontProperties(fname=str(font)).get_name()
    else:
        family = font_manager.findfont("IPAGothic", fallback_to_default=False)
        family = font_manager.FontProperties(fname=family).get_name()
    plt.rcParams.update({"font.family": family, "font.size": 11,
                         "axes.titlesize": 13, "axes.labelsize": 11,
                         "axes.unicode_minus": False, "pdf.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "legend.framealpha": .94, "figure.facecolor": "white"})


def page(title):
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.6))
    fig.subplots_adjust(left=.08, right=.965, bottom=.19, top=.82, wspace=.3)
    fig.suptitle(title, fontsize=19, y=.97, fontweight="bold")
    return fig, axes


def footer(fig, explanation):
    fig.text(.08, .075, explanation, fontsize=10, va="center", linespacing=1.6)
    fig.text(.08, .018,
             "出典：Kumar K・Dutta (2026), 式(2.3), (4.1)。文献の被覆率モデルであり、VISIT本体ではない。",
             fontsize=8, color="#555555")


def basin_panel(ax, x, branches):
    ax.fill_between(x, 0, branches[:, 1], color=PALE_GREEN,
                    label="低草本側の吸引域（環境固定時）")
    ax.fill_between(x, branches[:, 1], 1, color=PALE_ORANGE,
                    label="高草本側の吸引域（環境固定時）")
    for j, color, label in [(0, GREEN, "低草本・安定QSE"),
                             (1, INK, "不安定QSE／吸引域境界"),
                             (2, ORANGE, "高草本・安定QSE")]:
        ax.plot(x, branches[:, j], "--" if j == 1 else "-", c=color, lw=2, label=label)
    ax.set_ylim(0, 1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=ROOT / "docs/figures/forest_grass_phase_atlas")
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    style()
    data = {}
    pdf = PdfPages(out / "forest_grass_phase_atlas.pdf")
    def save(fig, name):
        pdf.savefig(fig)
        fig.savefig(out / (name + ".png"), dpi=165)
        fig.savefig(out / (name + ".pdf"))
        plt.close(fig)

    # 1. True two-parameter equilibrium phase diagram, including the boundary.
    fig, (ax, bx) = page("1｜双安定性：同じ環境でも、二つの群集が安定に存在する")
    alphas = np.linspace(.59, 1.5, 200)
    deaths = np.linspace(.28, 1.08, 181)
    classes = np.array([[phase_class(replace(BASE, alpha=float(a), phi1=float(q)))
                         for a in alphas] for q in deaths])
    cmap = ListedColormap(["#e5e5e5", "#dbe9f4", "#d7c2e8", "#f3dab9"])
    ax.pcolormesh(alphas, deaths, classes, shading="nearest", cmap=cmap,
                  norm=BoundaryNorm(np.arange(-.5, 4.5), 4), rasterized=True)
    # Analytic fold locus: delta=phi0/(g L'-L), alpha=delta L'.
    g = np.linspace(.01, .99, 6000)
    l = logistic(g)
    lp = l * (1 - l) / BASE.width
    den = g * lp - l
    valid = den > 0
    aa = np.full_like(g, np.nan)
    qq = aa.copy()
    aa[valid] = BASE.phi0 * lp[valid] / den[valid]
    qq[valid] = BASE.phi0 + BASE.phi0 / den[valid]
    ax.plot(aa, qq, c=INK, lw=1.2, label="鞍点–節点分岐（解析式）")
    ax.plot([float(phi(1, replace(BASE, phi1=float(q)))) for q in deaths],
            deaths, "--", c="#888888", lw=1, label="森林消失境界の安定性変化")
    lam = np.linspace(0, 1, 301)
    path = [path_parameters(float(v)) for v in lam]
    ax.plot([p.alpha for p in path], [p.phi1 for p in path], c=PURPLE, lw=3,
            label="Rの経路：双安定領域内")
    ax.scatter([path[0].alpha, path[-1].alpha], [path[0].phi1, path[-1].phi1],
               c=PURPLE, s=30, zorder=6)
    ax.annotate("開始", (path[0].alpha, path[0].phi1), xytext=(10, -16),
                textcoords="offset points", fontsize=9, color=PURPLE)
    ax.annotate("終了", (path[-1].alpha, path[-1].phi1), xytext=(8, 7),
                textcoords="offset points", fontsize=9, color=PURPLE)
    ax.scatter([BASE.alpha], [BASE.phi1], marker="*", c=INK, s=120, zorder=6)
    ax.annotate("基準値", (BASE.alpha, BASE.phi1), xytext=(8, 8),
                textcoords="offset points", fontsize=10)
    ax.set(xlim=(alphas[0], alphas[-1]), ylim=(deaths[0], deaths[-1]),
           xlabel="森林加入率 α［1/モデル時間］",
           ylabel="高草本側の死亡率 " + r"$\phi_1$" + "［1/モデル時間］", title="凍結系のパラメータ相図")
    handles = [Line2D([], [], marker="s", ls="", ms=10, color=cmap(i), label=t)
               for i, t in enumerate(["森林なしのみ安定", "内部安定平衡が一つ",
                                       "内部安定平衡が二つ", "内部＋森林なしが安定"])]
    ax.legend(handles=handles, fontsize=8, loc="lower right")
    roots = interior_roots()
    for left, right, color in [(0, roots[1], PALE_GREEN), (roots[1], 1, PALE_ORANGE)]:
        bx.axvspan(left, right, color=color)
    bx.axhline(0, c=INK, lw=1)
    bx.plot(g, drift(g), c=INK, lw=2)
    bx.scatter(roots[[0, 2]], [0, 0], c=[GREEN, ORANGE], s=75, zorder=5)
    bx.scatter([roots[1], 1], [0, 0], facecolors="white", edgecolors=INK, s=65, zorder=5)
    for x, dx in [(.07, .05), (.28, -.06), (.54, .07), (.9, -.07)]:
        bx.annotate("", (x + dx, -.012), (x, -.012),
                    arrowprops=dict(arrowstyle="->", color=INK, lw=1.6))
    for r, label, color in zip(roots, ["低草本\n安定", "不安定", "高草本\n安定"], [GREEN, INK, ORANGE]):
        bx.annotate(label + f"\ng={r:.3f}", (r, 0), xytext=(0, 24),
                    textcoords="offset points", ha="center", fontsize=10, color=color)
    bx.set(xlim=(0, 1.02), ylim=(-.07, .2), xlabel="草本被覆率 g",
           ylabel="草本被覆率の変化 dg/dt［1/モデル時間］", title="基準値での位相線：矢印は状態の移動方向")
    footer(fig, "紫の領域は「低草本＋森林多」と「高草本＋森林少」がともに安定。どちらも森林は正。\n"
                "R経路周辺は非常に狭い。分岐を越えないことは解析で保証され、枝の分離は図4で確認できる。")
    data["parameter_phase"] = dict(alpha=alphas.tolist(), phi1=deaths.tolist(),
                                  stable_class=classes.tolist())
    save(fig, "01_bistability")

    # 2. B-tipping and the contrast case of alpha-only forcing without folds.
    fig, (ax, bx) = page("2｜B-tipping と否定結果：平衡が消える場合・消えない場合")
    aa = np.linspace(.93, 1.43, 501)
    folds = fold_points()
    lo, hi = folds[0][1], folds[1][1]
    rr = [interior_roots(replace(BASE, alpha=float(a))) for a in aa]
    boundary = np.array([r[1] if len(r) == 3 else (0 if a < lo else 1)
                         for a, r in zip(aa, rr)])
    ax.fill_between(aa, 0, boundary, color=PALE_GREEN)
    ax.fill_between(aa, boundary, 1, color=PALE_ORANGE)
    # Separate branches to avoid joining the wrong roots across a fold.
    for j, color, label in [(0, GREEN, "低草本・安定"), (1, INK, "不安定"), (2, ORANGE, "高草本・安定")]:
        vals = []
        for a, r in zip(aa, rr):
            vals.append(r[j] if len(r) == 3 else
                        (r[0] if (j == 0 and a > hi) or (j == 2 and a < lo) else np.nan))
        ax.plot(aa, vals, "--" if j == 1 else "-", c=color, lw=2, label=label)
    for gg, a in folds:
        ax.axvline(a, ls=":", c=INK, lw=1)
        ax.scatter(a, gg, c=INK, s=30, zorder=5)
    ax.annotate("加入率を下げる", (.985, .15), (1.22, .15),
                arrowprops=dict(arrowstyle="->", color=GREEN, lw=2), color=GREEN)
    ax.annotate("低位平衡が消失", (.94, .90), (.94, .30),
                arrowprops=dict(arrowstyle="->", color=ORANGE, lw=2), fontsize=9, rotation=90)
    ax.set(xlim=(.93, 1.43), ylim=(0, 1), xlabel="森林加入率 α［1/モデル時間］",
           ylabel="草本被覆率 g", title=f"分岐点 α={lo:.4f}, {hi:.4f} を越える")
    ax.legend(loc="upper right", fontsize=9)
    aa2 = np.linspace(1.02, 1.32, 301)
    branches = np.array([interior_roots(replace(BASE, alpha=float(a))) for a in aa2])
    basin_panel(bx, aa2, branches)
    barrier = branches[0, 1]
    bx.axhline(barrier, c=PURPLE, ls="-.", lw=2, label="全時刻に共通の不変境界")
    u = np.linspace(0, 1, 401)
    alpha_u = 1.32 - .30 * smooth_step(u)
    start = branches[-1, 0]
    for tau, color, name in [(.1, "#3184c3", "非常に速い変化"), (200., "#192d50", "遅い変化")]:
        sol = solve_ivp(lambda v, y: tau * drift(y, replace(
            BASE, alpha=float(1.32 - .30 * smooth_step(v)))),
            (0, 1), [start], rtol=2e-11, atol=2e-13, dense_output=True)
        assert sol.success
        values = sol.sol(u)[0]
        assert np.max(values) < barrier
        bx.plot(alpha_u, values, color=color, lw=2, label=name)
    bx.set(xlim=(1.32, 1.02), ylim=(.09, .6), xlabel="森林加入率 α（右へ進むほど低下）",
           ylabel="草本被覆率 g", title="双安定区間内：加入率だけならRは起こらない")
    bx.legend(fontsize=8, loc="upper left")
    footer(fig, "左：変化が遅くても、低位の安定平衡そのものが消えればB-tippingが起こる。\n"
                "右：死亡率を固定した加入率だけの変化では、紫の境界を上向きに越えられない（低草本から開始）。")
    save(fig, "02_bifurcation_and_no_rate_tipping")

    # 3. Stochastic escape: never present an MFPT as a transition probability.
    fig, (ax, bx) = page("3｜N-tipping：平衡は残っていても、雑音によって障壁を越える")
    v = potential(g) - potential(roots[0])
    ax.plot(g, v, c=INK, lw=2)
    ax.axvspan(0, roots[1], color=PALE_GREEN)
    ax.axvspan(roots[1], 1, color=PALE_ORANGE)
    vr = potential(roots) - potential(roots[0])
    ax.scatter(roots, vr, c=[GREEN, INK, ORANGE], s=65, zorder=4)
    ax.annotate("", (.25, vr[1]), (.25, 0),
                arrowprops=dict(arrowstyle="<->", color=PURPLE, lw=2))
    ax.text(.27, vr[1] / 2, f"障壁 ΔV\n{vr[1]:.5f}", color=PURPLE, fontsize=10)
    ax.annotate("雑音による脱出", (.57, .005), (.21, .011),
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=-.25",
                                color=PURPLE, lw=2), color=PURPLE)
    target = (roots[1] + roots[2]) / 2
    ax.axvline(target, ls=":", c=PURPLE, lw=1)
    ax.text(target + .015, .01, "初到達の目標", fontsize=9, color=PURPLE)
    ax.set(xlim=(0, 1), xlabel="草本被覆率 g", ylabel="V(g) − V(低草本)［1/モデル時間］",
           title="ポテンシャル地形（V′=−f）：谷が安定状態")
    sigmas = np.linspace(.03, .10, 36)
    times = np.array([mean_first_passage(roots[0], target, float(s)) for s in sigmas])
    assert np.all(np.diff(times) < 0) and np.all(np.isfinite(times))
    bx.semilogy(sigmas, times, c=PURPLE, lw=2.5)
    for s in [.03, .04, .05]:
        tt = mean_first_passage(roots[0], target, s)
        bx.scatter(s, tt, color=PURPLE)
        bx.annotate(f"{tt:.2g}", (s, tt), xytext=(7, 6),
                    textcoords="offset points", fontsize=9)
    bx.set(xlabel="雑音の大きさ σ［1/√(モデル時間)］",
           ylabel="低草本からの平均初到達時間［モデル時間］",
           title="雑音が弱いほど、平均的な待ち時間は長い")
    bx.grid(axis="y", alpha=.2)
    footer(fig, "反射境界付きの加法雑音を追加した理論例。右図は平均待ち時間であり、遷移確率ではない。\n"
                "凍結平衡や障壁の位置は変えない。雑音を継続すれば逆方向の遷移も起こり得る。")
    data["noise"] = dict(sigma=sigmas.tolist(), mean_first_passage=times.tolist(), target=float(target))
    save(fig, "03_noise_escape")

    # 4. R-tipping: frozen basins vs actual nonautonomous paths and rate threshold.
    fig, (ax, bx) = page("4｜R-tipping：安定QSEが存続しても、速さだけで到達先が変わる")
    tau_star = critical_duration(rtol=2e-12, atol=2e-14, max_step=.004)
    fw, bw = matched_solutions(tau_star, rtol=2e-12, atol=2e-14)
    assert abs(fw.y[0, -1] - bw.y[0, -1]) < 1e-9
    branches = np.array([interior_roots(p) for p in path])
    basin_panel(ax, lam, branches)
    trajectories = []
    for factor, color, label in [(2., "#2676b7", "遅い変化：τ=2τ*"),
                                 (.5, "#bb3947", "速い変化：τ=0.5τ*")]:
        sol = ramp_solution(factor * tau_star)
        states = sol.sol(u)[0]
        assert np.all((states >= 0) & (states <= 1))
        assert (states[-1] > .32) == (factor < 1)
        ax.plot(smooth_step(u), states, c=color, lw=2.7, label=label)
        trajectories.append(dict(duration=factor * tau_star, lam=smooth_step(u).tolist(),
                                 grass=states.tolist()))
    ax.set(xlim=(0, 1), ylim=(.25, .8), xlabel="環境指標 λ（0→1、全ケースで同じ経路）",
           ylabel="草本被覆率 g", title="凍結吸引域と、実際の軌道の比較")
    ax.legend(fontsize=8, loc="upper left")
    ratios = np.r_[np.geomspace(.25, .98, 18), np.geomspace(1.02, 4, 18)]
    endpoints = np.array([ramp_solution(tau_star / r).y[0, -1] for r in ratios])
    assert np.all(endpoints[ratios < 1] < .32) and np.all(endpoints[ratios > 1] > .32)
    bx.axvspan(.25, 1, color=PALE_GREEN)
    bx.axvspan(1, 4, color=PALE_ORANGE)
    bx.axvline(1, c=INK, ls="--", label="臨界速度 r*")
    bx.plot([.25, 1], [.3, .3], c=GREEN, lw=4, label="最終到達先（解析的判定）")
    bx.plot([1, 4], [branches[-1, 2]] * 2, c=ORANGE, lw=4)
    bx.scatter([1, 1, 1], [.3, .32, branches[-1, 2]], facecolors="white",
               edgecolors=INK, zorder=5)
    bx.scatter(ratios, endpoints, marker="x", s=23, c=PURPLE, label="ランプ終了時（数値）")
    bx.set_xscale("log")
    bx.set_xticks([.25, .5, 1, 2, 4], labels=["0.25", "0.5", "1", "2", "4"])
    bx.set(xlim=(.25, 4), ylim=(.25, .8), xlabel="相対速度 r/r*（右ほど速い、r=1/τ）",
           ylabel="草本被覆率 g", title="速度による最終状態の相図")
    bx.text(.37, .41, "低草本へ\n森林 F=0.700", color=GREEN)
    bx.text(1.28, .59, "高草本へ\n森林 F≈0.258", color=ORANGE)
    bx.legend(fontsize=8, loc="upper left")
    footer(fig, "環境は終了後に固定。τ*≈1741.652 はモデル時間であり、年ではない。r=r* は不安定平衡に到達。\n"
                "加入率と死亡率が連動する、存在証明用の経路。色面は環境固定時の吸引域で、時変系の吸引域そのものではない。")
    data["rate"] = dict(critical_duration=tau_star, relative_rate=ratios.tolist(),
                        ramp_end_grass=endpoints.tolist(), trajectories=trajectories,
                        matched_residual=float(fw.y[0, -1] - bw.y[0, -1]))
    save(fig, "04_rate_tipping")

    # 5. View in forest units, directly addressing the original QSE question.
    fig, (ax, bx) = page("5｜QSEの増加と実状態の増加は同じではない")
    ax.plot(lam, 1 - branches[:, 0], c=GREEN, lw=3, label="追跡対象の森林QSE")
    ax.plot(lam, 1 - branches[:, 1], "--", c=INK, lw=1.5, label="吸引域境界（凍結系）")
    ax.plot(lam, 1 - branches[:, 2], c=ORANGE, lw=2, label="別の安定森林QSE")
    for record, color, label in zip(trajectories, ["#2676b7", "#bb3947"], ["遅い変化", "速い変化"]):
        ax.plot(record["lam"], 1 - np.array(record["grass"]), color=color, lw=2.5, label=label)
    ax.set(xlabel="環境指標 λ", ylabel="森林被覆率 F=1−g", ylim=(.2, .78),
           title="追跡したいQSEは 0.6 → 0.7 に増加する")
    ax.legend(loc="lower left", fontsize=9)
    values = [.6, .7, 1 - branches[-1, 2]]
    bars = bx.bar(["開始時\n平衡化済み", "遅い変化後\n最終到達先", "速い変化後\n最終到達先"],
                  values, color=[INK, "#2676b7", "#bb3947"], width=.6)
    for bar, value in zip(bars, values):
        bx.text(bar.get_x() + bar.get_width()/2, value + .025, f"{value:.3f}",
                ha="center", fontsize=15, fontweight="bold")
    bx.axhline(.7, c=GREEN, ls=":", label="終了環境の追跡対象QSE")
    bx.set(ylim=(0, .9), ylabel="森林被覆率 F", title="同じ開始状態・同じ終了環境でも到達先が異なる")
    bx.legend(loc="upper right", fontsize=9)
    footer(fig, "実状態が失うのは「追跡対象のQSE」であり、速い場合も別のQSEへ収束する。\n"
                "これは森林被覆率の結果。枯死物・土壌の炭素収支を含まないため、炭素放出量を示す図ではない。")
    save(fig, "05_forest_qse")
    pdf.close()
    (out / "plot_data.json").write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")
    sources = [Path(__file__), ROOT / "src/control_carbon/forest_grass_tipping.py",
               ROOT / "docs/forest_grass_nbr_tipping.tex"]
    snapshot = out / "source_snapshot"
    snapshot.mkdir()
    for source in sources:
        shutil.copy2(source, snapshot / source.name)
    (out / "manifest.json").write_text(json.dumps({
        "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in sources},
        "assumptions": "same published cover model; N reflected noise; R constructed path; NOT VISIT",
        "phase_classes": {"0": "boundary only", "1": "one interior",
                          "2": "two interior stable", "3": "one interior plus boundary"},
        "rate_endpoint_checks": "36 rates bracket the analytic threshold; no equality classified",
        "numerical_guards": "no state clipping; bounds and basin checks; matching residual < 1e-9",
        "grid_note": "phase colors sampled numerically; narrow regions near the cusp are under-resolved; black fold curves are analytic",
    }, ensure_ascii=False, indent=2) + "\n")
    print(f"Saved 5 figures, atlas PDF, data, and provenance in {out}")
    print(f"critical_duration={tau_star:.12f}; phase classes={np.unique(classes).tolist()}")


if __name__ == "__main__":
    main()
