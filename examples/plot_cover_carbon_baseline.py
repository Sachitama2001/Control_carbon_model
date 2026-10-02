"""Plot the cover--carbon baseline phase portrait and its prescribed functions."""

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from control_carbon.cover_carbon_baseline import (  # noqa: E402
    CarbonParameters,
    carbon_budget,
    carbon_capacity,
    carbon_rhs,
)
from control_carbon.forest_grass_tipping import (  # noqa: E402
    Parameters as CoverParameters,
    drift,
    drift_prime,
    h,
    interior_roots,
    logistic,
    phi,
)


GREEN = "#147d64"
ORANGE = "#cb851b"
PURPLE = "#8451ad"
INK = "#344250"
PALE_GREEN = "#d8eee5"
PALE_ORANGE = "#fae9cc"
BASE = CoverParameters()
CARBON = CarbonParameters(
    production=[1.1, 0.8],
    allocation=[[0.5, 0.3, 0.2], [0.4, 0.35, 0.25]],
    turnover=[[0.25, 0.12, 0.08], [0.3, 0.15, 0.1]],
    transfer_emission=[[0.2, 0.1, 0.0], [0.0, 0.0, 0.0]],
    decomposition=[0.4, 0.5],
    humification=[0.3, 0.25],
    humus_loss=[0.02, 0.03],
)
POOL_COLORS = ["#3977a8", "#8c6d31", "#25826e", "#707070", "#9a65a5"]
POOL_LABELS = ["葉", "茎・支持組織", "根", "枯死物", "腐植"]


def configure_style():
    candidates = [
        Path("/usr/share/fonts/opentype/ipaexfont-gothic/ipaexg.ttf"),
        Path("/usr/share/fonts/opentype/ipafont-gothic/ipag.ttf"),
    ]
    font_path = next((path for path in candidates if path.exists()), None)
    if font_path is not None:
        font_manager.fontManager.addfont(str(font_path))
        family = font_manager.FontProperties(fname=str(font_path)).get_name()
    else:
        family = font_manager.FontProperties(
            fname=font_manager.findfont("IPAGothic", fallback_to_default=False)
        ).get_name()
    plt.rcParams.update({
        "font.family": family,
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "axes.unicode_minus": False,
        "pdf.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.framealpha": 0.94,
        "figure.facecolor": "white",
    })


def stable_class(parameters):
    """Classify stable interior branches plus the g=1 boundary equilibrium."""
    roots = interior_roots(parameters)
    stable = roots[drift_prime(roots, parameters) < 0]
    boundary_stable = h(1.0, parameters) > 0
    if len(stable) == 2:
        if boundary_stable:
            raise AssertionError("Unexpected third stable equilibrium at g=1.")
        return 2
    if len(stable) == 1:
        return 3 if boundary_stable else 1
    if not boundary_stable:
        raise AssertionError("Parameter cell has no stable equilibrium.")
    return 0


def fold_curve(samples=1200):
    """Return analytic fold locus while varying alpha and phi1."""
    g = np.linspace(0.01, 0.99, samples)
    z = logistic(g, BASE)
    z_prime = z * (1.0 - z) / BASE.width
    denominator = g * z_prime - z
    valid = denominator > 0.0
    alpha = np.full_like(g, np.nan)
    phi1 = np.full_like(g, np.nan)
    alpha[valid] = BASE.phi0 * z_prime[valid] / denominator[valid]
    phi1[valid] = BASE.phi0 + BASE.phi0 / denominator[valid]
    valid &= (alpha > 0) & (phi1 > BASE.phi0)
    return alpha[valid], phi1[valid]


def page_footer(fig, text):
    fig.text(0.075, 0.035, text, fontsize=8.5, color="#4f5962", va="bottom")


def draw_phase_figure():
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(13.8, 6.8))
    fig.subplots_adjust(left=0.08, right=0.97, bottom=0.19, top=0.84, wspace=0.28)
    fig.suptitle(
        "被覆率の相図：同じ固定環境で二つの安定炭素平衡が成立する",
        fontsize=17,
        fontweight="bold",
        y=0.97,
    )

    alpha_axis = np.linspace(0.55, 1.55, 121)
    phi1_axis = np.linspace(0.30, 1.40, 121)
    classes = np.array([
        [
            stable_class(replace(BASE, alpha=float(alpha), phi1=float(phi1)))
            for alpha in alpha_axis
        ]
        for phi1 in phi1_axis
    ])
    cmap = ListedColormap(["#e5e5e5", "#dbe9f4", "#d7c2e8", "#f3dab9"])
    ax.pcolormesh(
        alpha_axis,
        phi1_axis,
        classes,
        shading="nearest",
        cmap=cmap,
        norm=BoundaryNorm(np.arange(-0.5, 4.5), cmap.N),
        rasterized=True,
    )
    fold_alpha, fold_phi1 = fold_curve()
    ax.plot(fold_alpha, fold_phi1, color=INK, lw=1.5, label="内部平衡の折返し（解析式）")
    boundary_phi1 = BASE.phi0 + (alpha_axis - BASE.phi0) / logistic(1.0, BASE)
    visible = (boundary_phi1 >= phi1_axis[0]) & (boundary_phi1 <= phi1_axis[-1])
    ax.plot(
        alpha_axis[visible],
        boundary_phi1[visible],
        "--",
        color="#777777",
        lw=1.1,
        label="g=1境界平衡の安定性境界",
    )
    ax.scatter([BASE.alpha], [BASE.phi1], marker="*", c=INK, s=150, zorder=5)
    ax.annotate(
        "基準値",
        (BASE.alpha, BASE.phi1),
        xytext=(8, 8),
        textcoords="offset points",
        fontsize=9,
    )
    legend_labels = [
        "境界平衡のみ安定",
        "内部安定平衡が一つ",
        "内部安定平衡が二つ",
        "内部平衡と境界平衡が安定",
    ]
    handles = [
        Line2D([], [], marker="s", ls="", ms=9, color=cmap(index), label=label)
        for index, label in enumerate(legend_labels)
    ]
    handles.extend(ax.get_legend_handles_labels()[0])
    ax.legend(handles=handles, fontsize=7.3, loc="upper left")
    ax.set(
        xlim=(alpha_axis[0], alpha_axis[-1]),
        ylim=(phi1_axis[0], phi1_axis[-1]),
        xlabel="草原から森林への加入係数 α［1/モデル時間］",
        ylabel="シグモイド上限 " + r"$\phi_1$" + "［1/モデル時間］",
        title="固定した " + r"$\phi_0$" + "・θ・w のパラメータ相図",
    )

    roots = interior_roots(BASE)
    g = np.linspace(0.0, 1.0, 1000)
    bx.axvspan(0.0, roots[1], color=PALE_GREEN, zorder=0)
    bx.axvspan(roots[1], 1.0, color=PALE_ORANGE, zorder=0)
    bx.axhline(0.0, color=INK, lw=1)
    bx.plot(g, drift(g, BASE), color=INK, lw=2.2)
    bx.scatter(roots[[0, 2]], [0.0, 0.0], c=[GREEN, ORANGE], s=70, zorder=5)
    bx.scatter([roots[1], 1.0], [0.0, 0.0], facecolors="white",
               edgecolors=INK, s=62, zorder=5)
    for x, direction in [
        (0.07, 1),
        (0.28, -1),
        (0.54, 1),
        (0.90, -1),
    ]:
        bx.annotate(
            "",
            xy=(x + direction * 0.055, -0.014),
            xytext=(x, -0.014),
            arrowprops={"arrowstyle": "->", "color": INK, "lw": 1.4},
        )
    labels = [
        ("低草原\n安定", roots[0], GREEN),
        ("不安定\n境界", roots[1], INK),
        ("高草原\n安定", roots[2], ORANGE),
    ]
    for label, root, color in labels:
        bx.annotate(
            f"{label}\ng={root:.3f}",
            (root, 0.0),
            xytext=(0, 21),
            textcoords="offset points",
            ha="center",
            fontsize=8.5,
            color=color,
        )
    bx.set(
        xlim=(0.0, 1.02),
        ylim=(-0.075, 0.20),
        xlabel="草原被覆率 g",
        ylabel="dg/dt［被覆率／モデル時間］",
        title="基準パラメータでの位相線",
    )
    bx.text(
        0.5,
        -0.068,
        "緑・橙の領域：その初期被覆率から安定平衡へ",
        ha="center",
        fontsize=8,
        color=INK,
    )
    page_footer(
        fig,
        "パラメータ相図は内部安定平衡数で分類し、g=1の境界安定性も別に表示。"
        "気象変数はまだ対応づけていない。"
        "全11状態系では被覆率が炭素から独立なので、各安定被覆平衡に安定な炭素平衡が一つずつ対応する。",
    )
    data = {
        "alpha": alpha_axis.tolist(),
        "phi1": phi1_axis.tolist(),
        "stable_class": classes.tolist(),
        "cover_equilibria": roots.tolist(),
        "drift": drift(g, BASE).tolist(),
        "fold_curve": {"alpha": fold_alpha.tolist(), "phi1": fold_phi1.tolist()},
    }
    return fig, data


def draw_function_figure():
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.0))
    fig.subplots_adjust(left=0.085, right=0.97, bottom=0.16, top=0.87,
                        hspace=0.42, wspace=0.28)
    fig.suptitle(
        "基準モデルで採用した関数の概形",
        fontsize=17,
        fontweight="bold",
        y=0.97,
    )
    g = np.linspace(0.0, 1.0, 600)
    roots = interior_roots(BASE)

    ax = axes[0, 0]
    ax.plot(g, phi(g, BASE), color=ORANGE, lw=2.4, label="森林→草原 φ(g)")
    ax.plot(g, BASE.alpha * g, color=GREEN, lw=2.0, label="比較項 αg")
    ax.plot(g, h(g, BASE), color=INK, lw=1.8, ls="--", label="φ(g)−αg")
    ax.axhline(0.0, color="#777777", lw=0.8)
    ax.set(
        xlabel="草原被覆率 g",
        ylabel="転換率／面積率［1/モデル時間］",
        title="ロジスティック応答と平衡条件",
    )
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.axvspan(0.0, roots[1], color=PALE_GREEN, zorder=0)
    ax.axvspan(roots[1], 1.0, color=PALE_ORANGE, zorder=0)
    ax.axhline(0.0, color="#777777", lw=0.8)
    ax.plot(g, drift(g, BASE), color=INK, lw=2.2)
    ax.scatter(roots[[0, 2]], [0.0, 0.0], color=[GREEN, ORANGE], s=55, zorder=5)
    ax.scatter([roots[1]], [0.0], facecolors="white", edgecolors=INK, s=55, zorder=5)
    ax.set(
        xlabel="草原被覆率 g",
        ylabel="dg/dt",
        title="面積制約を含む被覆率の変化",
    )

    ax = axes[1, 0]
    forest_npp = (1.0 - g) * CARBON.production[0]
    grass_npp = g * CARBON.production[1]
    ax.plot(
        g, forest_npp, color=GREEN, lw=2.2,
        label="森林 " + r"$\mu_T=(1-g)P_T$",
    )
    ax.plot(
        g, grass_npp, color=ORANGE, lw=2.2,
        label="草原 " + r"$\mu_G=gP_G$",
    )
    ax.set(
        xlabel="草原被覆率 g",
        ylabel=r"面積当たりNPP［kg C m$^{-2}$／モデル時間］",
        title="群別生産入力（試験係数）",
    )
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    rates = np.vstack((phi(g, BASE), BASE.alpha * (1.0 - g)))
    ax.plot(
        g, rates[0], color=ORANGE, lw=2.2,
        label="森林→草原 " + r"$r_T=\phi(g)$",
    )
    ax.plot(
        g, rates[1], color=GREEN, lw=2.2,
        label="草原→森林 " + r"$r_G=\alpha(1-g)$",
    )
    ax.set(
        xlabel="草原被覆率 g",
        ylabel="個体炭素プールへの面積転換損失［1/モデル時間］",
        title="炭素プールに入る転換損失率",
    )
    ax.legend(fontsize=8)
    page_footer(
        fig,
        "P_T=1.1、P_G=0.8 は図示用の未較正値。"
        "転換率は面積モデルの実効率であり、個体死亡率や燃焼率ではない。",
    )
    data = {
        "grass_cover": g.tolist(),
        "phi": phi(g, BASE).tolist(),
        "alpha_times_cover": (BASE.alpha * g).tolist(),
        "cover_balance": h(g, BASE).tolist(),
        "cover_drift": drift(g, BASE).tolist(),
        "forest_npp": forest_npp.tolist(),
        "grass_npp": grass_npp.tolist(),
        "forest_to_grass_carbon_loss_rate": rates[0].tolist(),
        "grass_to_forest_carbon_loss_rate": rates[1].tolist(),
    }
    return fig, data


def draw_capacity_figure():
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 6.8))
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.20, top=0.82, wspace=0.30)
    fig.suptitle(
        "固定被覆下の炭素capacity：各被覆平衡で炭素平衡へ接続",
        fontsize=17,
        fontweight="bold",
        y=0.96,
    )
    g = np.linspace(0.0, 1.0, 301)
    capacities = np.array([carbon_capacity(value, CARBON, BASE) for value in g])
    roots = interior_roots(BASE)
    stable_roots = roots[drift_prime(roots, BASE) < 0]
    groups = capacities.reshape(len(g), 2, 5)
    group_totals = groups.sum(axis=2)
    total = group_totals.sum(axis=1)

    for group, ax, label in [
        (0, axes[0], "森林由来プール"),
        (1, axes[1], "草原由来プール"),
    ]:
        for pool in range(5):
            ax.plot(
                g,
                groups[:, group, pool],
                color=POOL_COLORS[pool],
                lw=2.0,
                label=POOL_LABELS[pool],
            )
        for root, color in zip(stable_roots, [GREEN, ORANGE]):
            ax.axvline(root, color=color, ls=":", lw=1.3)
        ax.set(
            xlabel="草原被覆率 g",
            ylabel=r"capacity［kg C m$^{-2}$］",
            title=label,
            xlim=(0.0, 1.0),
        )
        ax.legend(fontsize=7.5)

    axes[2].plot(g, group_totals[:, 0], color=GREEN, lw=2.2, label="森林由来合計")
    axes[2].plot(g, group_totals[:, 1], color=ORANGE, lw=2.2, label="草原由来合計")
    axes[2].plot(g, total, color=INK, lw=2.4, label="総capacity")
    for root, color in zip(stable_roots, [GREEN, ORANGE]):
        value = np.interp(root, g, total)
        axes[2].scatter([root], [value], color=color, s=52, zorder=5)
    axes[2].set(
        xlabel="草原被覆率 g",
        ylabel=r"capacity［kg C m$^{-2}$］",
        title="由来別・総capacity",
        xlim=(0.0, 1.0),
    )
    axes[2].legend(fontsize=7.5)
    page_footer(
        fig,
        "係数は未較正のテスト用設定。点線は二つの安定被覆平衡。"
        "その位置のcapacityが全系の炭素平衡となる。"
        "総炭素量の枝間の大小は一般定理ではなく、係数依存の図示例に限る。",
    )
    for root in roots:
        equilibrium = carbon_capacity(float(root), CARBON, BASE)
        residual = np.max(np.abs(carbon_rhs(float(root), equilibrium, CARBON, BASE)))
        if residual > 1e-12:
            raise AssertionError(f"Carbon equilibrium residual is too large: {residual}")
        if np.min(equilibrium) < -1e-12:
            raise AssertionError("Carbon capacity left the nonnegative domain.")
    budget = carbon_budget(0.43, capacities[129], CARBON, BASE)
    if abs(budget["residual"]) > 1e-11:
        raise AssertionError("Carbon budget does not close.")

    data = {
        "grass_cover": g.tolist(),
        "capacity_by_group_and_pool": capacities.tolist(),
        "group_totals": group_totals.tolist(),
        "total_capacity": total.tolist(),
        "interior_cover_equilibria": roots.tolist(),
        "stable_cover_equilibria": stable_roots.tolist(),
        "carbon_budget_check": budget,
    }
    return fig, data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "docs/figures/cover_carbon_baseline",
    )
    args = parser.parse_args()
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    configure_style()

    phase_fig, phase_data = draw_phase_figure()
    function_fig, function_data = draw_function_figure()
    capacity_fig, capacity_data = draw_capacity_figure()
    pages = [
        ("01_cover_phase", phase_fig),
        ("02_function_shapes", function_fig),
        ("03_carbon_capacity", capacity_fig),
    ]
    with PdfPages(output / "cover_carbon_baseline_atlas.pdf") as pdf:
        for name, figure in pages:
            pdf.savefig(figure)
            figure.savefig(output / f"{name}.png", dpi=180)
            figure.savefig(output / f"{name}.pdf")
            plt.close(figure)

    parameter_map = {
        "alpha": {"minimum": 0.55, "maximum": 1.55, "count": 121},
        "phi1": {"minimum": 0.30, "maximum": 1.40, "count": 121},
        "fixed": {
            "phi0": BASE.phi0,
            "theta": BASE.theta,
            "width": BASE.width,
        },
        "stable_class": {
            "0": "interior安定平衡なし（境界平衡のみ）",
            "1": "内部安定平衡が一つ",
            "2": "内部安定平衡が二つ",
            "3": "内部平衡とg=1境界平衡が安定",
        },
    }
    figure_data = {
        "cover_parameter_phase": phase_data,
        "adopted_functions": function_data,
        "carbon_capacity": capacity_data,
    }
    (output / "plot_data.json").write_text(
        json.dumps(figure_data, ensure_ascii=False, indent=2) + "\n"
    )

    sources = [
        Path(__file__),
        ROOT / "src/control_carbon/cover_carbon_baseline.py",
        ROOT / "src/control_carbon/forest_grass_tipping.py",
        ROOT / "docs/cover_carbon_baseline.md",
        ROOT / "docs/forest_grass_cover_climate_research_plan.md",
    ]
    snapshot = output / "source_snapshot"
    snapshot.mkdir(exist_ok=True)
    for source in sources:
        shutil.copy2(source, snapshot / source.name)
    (output / "manifest.json").write_text(
        json.dumps(
            {
                "source_sha256": {
                    str(source.relative_to(ROOT)): hashlib.sha256(
                        source.read_bytes()
                    ).hexdigest()
                    for source in sources
                },
                "cover_parameters": {
                    "alpha": BASE.alpha,
                    "phi0": BASE.phi0,
                    "phi1": BASE.phi1,
                    "theta": BASE.theta,
                    "width": BASE.width,
                },
                "carbon_parameters_are_calibrated": False,
                "carbon_test_parameters": {
                    "production": CARBON.production.tolist(),
                    "allocation": CARBON.allocation.tolist(),
                    "turnover": CARBON.turnover.tolist(),
                    "transfer_emission": CARBON.transfer_emission.tolist(),
                    "decomposition": CARBON.decomposition.tolist(),
                    "humification": CARBON.humification.tolist(),
                    "humus_loss": CARBON.humus_loss.tolist(),
                },
                "parameter_map": parameter_map,
                "scope": "analytical cover-carbon baseline; no climate response or tipping result",
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n"
    )
    print(f"Saved 3 figures, atlas PDF, data, and provenance in {output}")
    print(f"Interior cover equilibria: {interior_roots(BASE).tolist()}")


if __name__ == "__main__":
    main()
