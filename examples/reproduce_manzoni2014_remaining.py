"""Fig. 3 / Appendix A mechanisms and explicitly assumed Fig. 4 sensitivity."""

import argparse
from dataclasses import asdict, replace
from datetime import datetime
import hashlib
import json
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, SymLogNorm
from matplotlib.lines import Line2D
import numpy as np
import scipy
from scipy.optimize import brentq

from reproduce_manzoni2014 import PAPER, ROOT, require
from control_carbon.manzoni_hydraulics import (
    HydraulicParameters, PAPER_DOI, hydraulic_fold, steady_fluxes,
    transpiring_equilibria,
)
from control_carbon.manzoni_extensions import (
    SigmoidShape, SoilRootParameters, deteriorated_traits,
    sigmoid_equilibria, sigmoid_equilibrium_jacobian_scaled,
    sigmoid_log_flux_ratio, sigmoid_steady_fluxes, soil_equilibria,
    soil_equilibrium_state, soil_fold, soil_hydraulic_jacobian,
    soil_hydraulic_rhs,
)


def save_figure(fig, output, name):
    fig.tight_layout()
    fig.savefig(output / f"{name}.png", dpi=160)
    fig.savefig(output / f"{name}.pdf")
    plt.close(fig)


def nullable(array):
    return [[float(value) if np.isfinite(value) else None for value in row] for row in array]


def figure3(output, p):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    roots = np.linspace(-2.999, 0., 240)
    scenarios = []
    for ratio, style in [(39., "--"), (78., "-"), (156., "--")]:
        traits = replace(p, psi50_x=-1.5, gs_max=ratio*p.gx_max)
        folds = np.array([hydraulic_fold(root, traits) for root in roots])
        axes[0].plot(
            roots, folds[:, 0], style,
            label=f"gs/gx={ratio:g} MPa" + (" (paper median)" if ratio == 78 else " (sensitivity)"),
        )
        if ratio == 78:
            axes[0].fill_between(roots, folds[:, 0], .05, color=".8")
        scenarios.append({"conductance_ratio_MPa": ratio, "critical_D": folds[:, 0].tolist()})
    axes[0].set(xlim=(-3., 0.), ylim=(0., .05), xlabel="Root potential (MPa)",
                ylabel="Critical D (mol/mol)", title="Fig. 3a: median + labelled sensitivity")
    axes[0].legend(fontsize=8)

    s_values = np.linspace(-4., -.01, 300)
    boundaries = []
    max_fold_error = 0.
    for d in np.arange(.01, .061, .01):
        t = 78*d
        boundary = []
        for s in s_values:
            b, r = -2*s, .5
            if b <= r:
                boundary.append(-8.)
            elif t >= b-r:
                boundary.append(s)
            else:
                def difference(a):
                    traits = replace(p, psi50_x=-a/2, psi50_s=s, gs_max=78*p.gx_max)
                    return hydraulic_fold(-r, traits)[0]-d
                a = brentq(difference, r+1e-10, b-1e-10)
                traits = replace(p, psi50_x=-a/2, psi50_s=s, gs_max=78*p.gx_max)
                fold_d, leaf = hydraulic_fold(-r, traits)
                max_fold_error = max(max_fold_error, abs(fold_d-d))
                require(len(transpiring_equilibria(-r, .999*d, traits)) == 2, "Fig. 3 pre-fold count failed")
                require(not transpiring_equilibria(-r, 1.001*d, traits), "Fig. 3 post-fold count failed")
                boundary.append(-a/2)
        axes[1].plot(s_values, boundary, label=f"D={d:.2f}")
        if abs(d-.01) < 1e-10:
            axes[1].fill_between(s_values, boundary, 0., color=".85")
        boundaries.append({"D": float(d), "psi50_x_boundary": boundary})
    axes[1].axvline(-.5, color=".5", ls=":", label="root=-0.5 MPa (not closure)")
    axes[1].set(xlim=(-4., 0.), ylim=(-8., 0.), xlabel="Stomatal half-loss potential (MPa)",
                ylabel="Xylem half-loss potential (MPa)", title="Fig. 3b: active-equilibrium-loss boundaries")
    axes[1].legend(fontsize=8)
    save_figure(fig, output, "figure3_thresholds")
    require(max_fold_error < 1e-10, "Fig. 3 boundary residual failed")
    return {
        "root_grid_MPa": roots.tolist(), "ratio_scenarios": scenarios,
        "psi50_s_grid_MPa": s_values.tolist(), "trait_boundaries": boundaries,
        "max_fold_D_error": max_fold_error,
        "limits": "No empirical percentile curves, histogram or observed trait points reproduced.",
    }


def appendix(output, p, shape_values):
    max_log_residual = 0.
    max_potential = 0.
    def equilibria(root, d, traits, shape):
        nonlocal max_log_residual, max_potential
        leaves = sigmoid_equilibria(root, d, traits, shape)
        labels = []
        for leaf in leaves:
            residual = abs(sigmoid_log_flux_ratio(leaf, root, d, traits, shape))
            max_log_residual = max(max_log_residual, residual)
            max_potential = max(max_potential, abs(leaf))
            eig = np.linalg.eigvals(sigmoid_equilibrium_jacobian_scaled(leaf, root, d, traits, shape))
            require(np.max(abs(eig.imag)) < 1e-9, "unexpected complex equilibrium eigenvalues")
            labels.append("stable" if np.max(eig.real) < 0 else "saddle")
        require(labels == ["stable" if i % 2 == 0 else "saddle" for i in range(len(leaves))],
                "Appendix equilibrium stability does not alternate")
        return leaves, labels

    linear_stable = transpiring_equilibria(0., .02, p)[0]
    linear_dry_boundary = 2*p.psi50_x
    errors = np.full((2, len(shape_values), len(shape_values)), np.nan)
    counts_a1 = np.zeros(errors.shape[1:], dtype=int)
    for i, ass in enumerate(shape_values):
        for j, ax in enumerate(shape_values):
            leaves, labels = equilibria(0., .02, p, SigmoidShape(ax, ass))
            counts_a1[i, j] = len(leaves)
            if leaves:
                errors[0, i, j] = linear_stable-leaves[0]
            if len(leaves) >= 2:
                errors[1, i, j] = linear_dry_boundary-leaves[1]

    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    examples = []
    leaf_grid = np.linspace(-8., 0., 600)
    for j, (ax, ass) in enumerate([(4., 12.), (12., 12.), (4., 4.), (12., 4.)]):
        shape = SigmoidShape(ax, ass)
        for traits_name, color in [("linear", "black"), ("sigmoid", "darkorange")]:
            fluxes = np.array([
                steady_fluxes(leaf, 0., .02, p) if traits_name == "linear"
                else sigmoid_steady_fluxes(leaf, 0., .02, p, shape)
                for leaf in leaf_grid
            ])
            axes[0, j].plot(leaf_grid, 1000*fluxes[:, 0], "--", color=color)
            axes[0, j].plot(leaf_grid, 1000*fluxes[:, 1], "-", color=color)
        leaves, labels = equilibria(0., .02, p, shape)
        for leaf, label in zip(leaves, labels):
            axes[0, j].plot(
                leaf, 1000*sigmoid_steady_fluxes(leaf, 0., .02, p, shape)[1],
                "o", color="darkorange",
                markerfacecolor="darkorange" if label == "stable" else "white",
            )
        axes[0, j].set(title=f"{'ABCD'[j]}: ax={ax:g}, as={ass:g}", xlim=(-8., 0.),
                       ylim=(0., 16.), xlabel="Leaf potential (MPa)", ylabel="Flux (mm/day)")
        examples.append({"shape": asdict(shape), "leaves_MPa": list(leaves), "stability": labels})
    for ax in axes[1]:
        ax.remove()
    ax1 = fig.add_subplot(2, 2, 3)
    ax2 = fig.add_subplot(2, 2, 4)
    mesh = ax1.pcolormesh(shape_values, shape_values, np.ma.masked_invalid(errors[0]), shading="nearest")
    fig.colorbar(mesh, ax=ax1, label="Linear stable - sigmoid stable (MPa)")
    mesh = ax2.pcolormesh(
        shape_values, shape_values, np.ma.masked_invalid(errors[1]),
        shading="nearest", norm=SymLogNorm(linthresh=1.),
    )
    fig.colorbar(mesh, ax=ax2, label="Boundary - saddle (MPa; symlog colors)")
    for ax in [ax1, ax2]:
        ax.plot(shape_values, shape_values-1, "k:", label="as=ax-1 (equal tail powers)")
        ax.set(xlim=(1., 15.), ylim=(1., 15.), xlabel="alpha_x", ylabel="alpha_s")
        ax.legend(fontsize=8)
    save_figure(fig, output, "appendix_A1")

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    maps = []
    for row, d in enumerate([.01, .04]):
        for col, traits in enumerate([p, replace(p, psi50_x=-1.5)]):
            counts = np.zeros((len(shape_values), len(shape_values)), dtype=int)
            for i, ass in enumerate(shape_values):
                for j, ax in enumerate(shape_values):
                    counts[i, j] = len(equilibria(-.5, d, traits, SigmoidShape(ax, ass))[0])
            require(np.all(counts <= 3), "more than three Appendix equilibria found")
            cmap = plt.get_cmap("viridis", 4)
            image = axes[row, col].pcolormesh(
                shape_values, shape_values, counts, shading="nearest",
                cmap=cmap, norm=BoundaryNorm(np.arange(-.5, 4., 1.), cmap.N),
            )
            axes[row, col].plot(shape_values, shape_values-1, "k:")
            name = "safe half-loss order" if col == 0 else "unsafe half-loss order"
            axes[row, col].set(title=f"{name}, D={d:.2f}", xlim=(1., 15.), ylim=(1., 15.),
                               xlabel="alpha_x", ylabel="alpha_s")
            fig.colorbar(image, ax=axes[row, col], ticks=[0, 1, 2, 3], label="Finite equilibrium count")
            maps.append({
                "D": d, "psi50_x": traits.psi50_x, "psi50_s": traits.psi50_s,
                "counts": counts.tolist(),
                "histogram": {str(n): int(np.sum(counts == n)) for n in range(4)},
            })
    save_figure(fig, output, "appendix_A2_counts")
    require(max_log_residual < 1e-9, f"Appendix log flux residual too large: {max_log_residual}")
    require(set(n for m in maps for row in m["counts"] for n in row) == {0, 1, 2, 3},
            "Appendix campaign did not exhibit all 0/1/2/3 counts")
    return {
        "shape_grid": shape_values.tolist(), "A1_examples": examples,
        "A1_linear_stable_MPa": linear_stable,
        "A1_linear_nonconducting_boundary_MPa": linear_dry_boundary,
        "A1_counts": counts_a1.tolist(),
        "A1_stable_error_MPa": nullable(errors[0]),
        "A1_boundary_to_saddle_difference_MPa": nullable(errors[1]),
        "A2_maps": maps, "max_abs_log_supply_demand_ratio": max_log_residual,
        "max_finite_equilibrium_magnitude_MPa": max_potential,
        "limits": "A2 low/high D and unsafe traits are explicit scenarios, not caption-specified numerical values; observed points omitted.",
    }


def figure4(output, p):
    damaged = deteriorated_traits(p, .6)
    soil_grid = np.linspace(-3.59, -.1, 60)
    d_grid = np.linspace(0., .03, 65)
    fig = plt.figure(figsize=(13, 5))
    surface_ax = fig.add_subplot(1, 2, 1, projection="3d")
    threshold_ax = fig.add_subplot(1, 2, 2)
    max_residual = 0.
    max_fold_eigenvalue = 0.
    records = []
    surfaces = []
    for gsr_sat in [.03, .1, .3]:
        soil = SoilRootParameters(gsr_sat, -.1)
        folds = np.array([soil_fold(value, damaged, soil) for value in soil_grid])
        threshold_ax.plot(soil_grid, folds[:, 0], label=f"gsr_sat={gsr_sat:g} m/day/MPa")
        for psi_soil, (d, leaf) in zip(soil_grid, folds):
            state = soil_equilibrium_state(leaf, psi_soil, damaged, soil)
            residual = float(np.max(abs(soil_hydraulic_rhs(state, psi_soil, d, damaged, soil))))
            max_residual = max(max_residual, residual)
            eig = np.linalg.eigvals(soil_hydraulic_jacobian(state, psi_soil, d, damaged, soil))
            max_fold_eigenvalue = max(max_fold_eigenvalue, float(min(abs(eig))))
            require(min(eig.real) < 0, "soil fold lacks transverse stable eigenvalue")
            require(len(soil_equilibria(psi_soil, .999*d, damaged, soil)) == 2, "soil pre-fold count failed")
            require(not soil_equilibria(psi_soil, 1.001*d, damaged, soil), "soil post-fold count failed")
        records.append({"soil_parameters": asdict(soil), "critical_D": folds[:, 0].tolist(),
                        "fold_leaf_MPa": folds[:, 1].tolist()})
        if gsr_sat == .1:
            for traits, color, label in [
                (p, "gray", "undamaged, safe order"),
                (damaged, "orange", "damage=0.6, wet stable"),
            ]:
                wet = np.full((len(soil_grid), len(d_grid)), np.nan)
                saddle = wet.copy()
                for i, value in enumerate(soil_grid):
                    for j, d in enumerate(d_grid):
                        leaves = soil_equilibria(value, d, traits, soil)
                        for k, leaf in enumerate(leaves):
                            state = soil_equilibrium_state(leaf, value, traits, soil)
                            residual = float(np.max(abs(soil_hydraulic_rhs(state, value, d, traits, soil))))
                            max_residual = max(max_residual, residual)
                            eig = np.linalg.eigvals(soil_hydraulic_jacobian(state, value, d, traits, soil))
                            require(np.max(eig.real) < 0 if k == 0 else np.prod(eig.real) < 0,
                                    "soil branch stability failed")
                            if k == 0:
                                wet[i, j] = leaf
                            else:
                                saddle[i, j] = leaf
                D, S = np.meshgrid(d_grid, soil_grid)
                surface_ax.plot_surface(D, S, wet, color=color, alpha=.55)
                if traits == damaged:
                    surface_ax.plot_wireframe(D, S, saddle, color="black", linewidth=.4)
                surfaces.append({"label": label, "wet_MPa": nullable(wet), "saddle_MPa": nullable(saddle)})
            surface_ax.plot(folds[:, 0], soil_grid, folds[:, 1], color="red", lw=2)
    surface_ax.set(xlabel="D (mol/mol)", ylabel="Soil potential (MPa)",
                   zlabel="Leaf equilibrium (MPa)", title="Fig. 4 mechanism: assumed soil coefficients")
    surface_ax.legend(handles=[
        Line2D([], [], color="gray", label="undamaged stable"),
        Line2D([], [], color="orange", label="damaged stable"),
        Line2D([], [], color="black", label="damaged saddle"),
        Line2D([], [], color="red", label="damaged fold"),
    ], fontsize=7, loc="upper left")
    threshold_ax.set(xlabel="Soil potential (MPa)", ylabel="Critical D (mol/mol)",
                     title="Fixed damage factor 0.6: soil sensitivity")
    threshold_ax.legend(fontsize=8)
    save_figure(fig, output, "figure4_soil_sensitivity")
    require(max_residual < 1e-10, "soil equilibrium residual exceeds threshold")
    require(max_fold_eigenvalue < 1e-8, "soil fold eigenvalue check failed")
    return {
        "undamaged_traits": asdict(p), "damaged_traits": asdict(damaged),
        "soil_grid_MPa": soil_grid.tolist(), "D_grid": d_grid.tolist(),
        "scenarios": records, "surfaces": surfaces,
        "max_rhs_residual_MPa_per_day": max_residual,
        "max_abs_fold_small_eigenvalue_per_day": max_fold_eigenvalue,
        "limits": "User-approved assumed soil coefficients; not numerical reproduction of original Fig. 4.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "artifacts/manzoni2014" /
                        datetime.now().strftime("remaining_%Y%m%d_%H%M%S"))
    parser.add_argument("--shape-step", type=float, default=.5)
    args = parser.parse_args()
    if not np.isfinite(args.shape_step) or args.shape_step <= 0:
        parser.error("--shape-step must be finite and positive")
    intervals = 14/args.shape_step
    if intervals < 1 or intervals > 140 or abs(intervals-round(intervals)) > 1e-9:
        parser.error("--shape-step must divide 14 and produce 1..140 intervals")
    args.output.mkdir(parents=True, exist_ok=False)
    shapes = np.linspace(1., 15., round(intervals)+1)
    p = HydraulicParameters(.0086, .47, -3., -2., .001, .001)
    results = {
        "paper_doi": PAPER_DOI, "baseline_parameters": asdict(p),
        "figure3": figure3(args.output, p),
        "appendix": appendix(args.output, p, shapes),
        "figure4": figure4(args.output, p),
    }
    sources = [
        PAPER, Path(__file__).resolve(), ROOT / "examples/reproduce_manzoni2014.py",
        ROOT / "src/control_carbon/manzoni_hydraulics.py",
        ROOT / "src/control_carbon/manzoni_extensions.py",
        ROOT / "tests/test_manzoni_extensions.py",
    ]
    results["manifest"] = {
        "python": platform.python_version(), "numpy": np.__version__,
        "scipy": scipy.__version__, "matplotlib": matplotlib.__version__,
        "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                   for path in sources},
    }
    (args.output / "results.json").write_text(
        json.dumps(results, indent=2, allow_nan=False)+"\n", encoding="utf-8",
    )
    print(json.dumps({
        "output": str(args.output),
        "figure3_max_fold_D_error": results["figure3"]["max_fold_D_error"],
        "appendix_max_log_flux_residual": results["appendix"]["max_abs_log_supply_demand_ratio"],
        "appendix_max_finite_equilibrium_magnitude_MPa": results["appendix"]["max_finite_equilibrium_magnitude_MPa"],
        "appendix_A2_histograms": [m["histogram"] for m in results["appendix"]["A2_maps"]],
        "figure4_max_rhs_residual": results["figure4"]["max_rhs_residual_MPa_per_day"],
        "figure4_max_abs_fold_small_eigenvalue": results["figure4"]["max_abs_fold_small_eigenvalue_per_day"],
    }, indent=2))


if __name__ == "__main__":
    main()
