"""Reproduce Fig. 2 mechanisms; no carbon coupling or rate-tipping claim.

Run from the repository: python examples/reproduce_manzoni2014.py --output PATH
PATH must not already exist. See docs/manzoni2014_reproduction.md.
"""

import argparse
from dataclasses import asdict, replace
from datetime import datetime
import hashlib
import json
from pathlib import Path
import platform

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import solve_ivp

from control_carbon.manzoni_hydraulics import (
    HydraulicParameters, PAPER_DOI, conductances, hydraulic_fold,
    hydraulic_jacobian, hydraulic_rhs, steady_fluxes, transpiring_equilibria,
)


ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "docs/Water Resources Research - 2014 - Manzoni - A dynamical system perspective on plant hydraulic failure.pdf"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "artifacts/manzoni2014" / datetime.now().strftime("%Y%m%d_%H%M%S"),
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    p = HydraulicParameters(.0086, .47, -1.5, -2., .001, .001)
    safe = replace(p, psi50_x=-2., psi50_s=-1.5)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    ax = axes[0]
    leaves = np.linspace(-3.999, 0., 1200)
    for traits, d, color, label in [
        (safe, .01, "blue", "safe, D=0.01"),
        (p, .01, "black", "unsafe, D=0.01"),
        (p, .03, "turquoise", "unsafe, D=0.03"),
    ]:
        pairs = [
            (leaf + gs*d/(2*gx), leaf)
            for leaf in leaves
            for gx, gs in [conductances(leaf, traits)] if gx > 1e-7
        ]
        pairs = np.asarray(pairs)
        ax.plot(pairs[:, 0], pairs[:, 1], color=color, label=label)
        for leaf in transpiring_equilibria(0., d, traits):
            state = np.array([leaf/2, leaf])
            eig = np.linalg.eigvals(hydraulic_jacobian(state, 0., d, traits))
            ax.plot(
                *state, "o", color=color,
                markerfacecolor=color if np.max(eig.real) < 0 else "white",
            )
    ax.plot([-2., 0.], [-4., 0.], ":", color="gray", label="conducting xylem nullcline")
    xs, ls = np.linspace(-2., -.001, 30), np.linspace(-3.95, -.001, 35)
    grid = np.array([
        [hydraulic_rhs([x, leaf], 0., .01, p) for x in xs] for leaf in ls
    ])
    ax.streamplot(xs, ls, grid[:, :, 0], grid[:, :, 1], color=".6", density=.7)
    ax.set(xlim=(-2., 0.), ylim=(-4.1, 0.), xlabel="Xylem potential (MPa)",
           ylabel="Leaf potential (MPa)", title="A: conducting nullclines and full ODE")
    ax.legend(fontsize=8)
    ax.axhspan(-4.1, -4., color="lightgray")
    ax.text(-1.9, -4.07, "nonconducting equilibrium continuum", fontsize=7)

    report = {
        "paper_doi": PAPER_DOI,
        "scope": "Independent piecewise-linear Fig. 2 reproduction, root boundary prescribed",
        "unsafe_parameters": asdict(p), "safe_parameters": asdict(safe),
        "numerical_assumptions": {
            "wet_root_MPa": 0.,
            "wet_root_note": "Zero-potential limit; paper uses saturation without a numerical value",
            "capacitances": "Cx=Cl=0.001 m/MPa: numerical fixtures, not paper-calibrated",
        },
        "folds": [], "equilibria": [], "trajectories": [],
    }
    maximum_residual = 0.
    for root, color in [(0., "black"), (-1., "gray")]:
        supply = [steady_fluxes(leaf, root, .01, p)[0]*1000 for leaf in leaves]
        axes[1].plot(leaves, supply, color=color, lw=2, label=f"root={root:g} MPa")
        fold_d, fold_leaf = hydraulic_fold(root, p)
        fold_state = np.array([(root+fold_leaf)/2, fold_leaf])
        fold_eig = np.linalg.eigvals(hydraulic_jacobian(fold_state, root, fold_d, p))
        require(min(abs(fold_eig)) < 1e-9 and min(fold_eig.real) < 0, "invalid full-system fold")
        require(len(transpiring_equilibria(root, .999*fold_d, p)) == 2, "missing pre-fold branches")
        require(not transpiring_equilibria(root, 1.001*fold_d, p), "post-fold active roots remain")
        report["folds"].append({
            "psi_root_MPa": root, "critical_D_mol_per_mol": fold_d,
            "psi_leaf_MPa": fold_leaf, "jacobian_eigenvalues_per_day": fold_eig.real.tolist(),
        })
        for d, style in [(.01, "-"), (.04, ":"), (fold_d, "--")]:
            demand = [steady_fluxes(leaf, root, d, p)[1]*1000 for leaf in leaves]
            axes[1].plot(leaves, demand, style, color=color, lw=.8)
        ds = np.unique(np.append(np.linspace(0., fold_d, 300), [.01, fold_d]))
        branch = np.full((2, len(ds)), np.nan)
        for j, d in enumerate(ds):
            eqs = transpiring_equilibria(root, d, p)
            for k, leaf in enumerate(eqs):
                branch[k, j] = leaf
            if len(eqs) == 1 and abs(d-fold_d) < 1e-14:
                branch[:, j] = eqs[0]
        axes[2].plot(ds, branch[0], color=color, label=f"root={root:g} MPa")
        axes[2].plot(ds, branch[1], "--", color=color)
        axes[2].plot(fold_d, fold_leaf, "o", color=color)
        for d in [.01, .04, fold_d]:
            for leaf in transpiring_equilibria(root, d, p):
                state = np.array([(root+leaf)/2, leaf])
                residual = float(max(abs(hydraulic_rhs(state, root, d, p))))
                maximum_residual = max(maximum_residual, residual)
                eig = np.linalg.eigvals(hydraulic_jacobian(state, root, d, p))
                label = "fold" if min(abs(eig)) < 1e-9 else (
                    "stable" if np.max(eig.real) < 0 else "saddle"
                )
                report["equilibria"].append({
                    "psi_root_MPa": root, "D_mol_per_mol": d,
                    "state_MPa": state.tolist(), "stability": label,
                    "jacobian_eigenvalues_per_day": eig.real.tolist(),
                    "max_rhs_residual_MPa_per_day": residual,
                })
    axes[1].set(xlim=(-4., 0.), ylim=(0., 10.), xlabel="Leaf potential (MPa)",
                ylabel="Water flux (mm/day)", title="B: supply (thick) and demand (thin)")
    axes[1].legend()
    axes[2].axhline(-4., color="black", lw=1)
    axes[2].text(.002, -4.4, "dry leaf limit; not an isolated 2D sink", fontsize=8)
    axes[2].set(xlim=(0., .05), ylim=(-4.6, 0.), xlabel="D (mol/mol)",
                ylabel="Equilibrium leaf potential (MPa)", title="C: solid stable / dashed saddle")
    axes[2].legend()
    fig.tight_layout()
    fig.savefig(args.output / "figure2_reproduction.png", dpi=180)
    fig.savefig(args.output / "figure2_reproduction.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    max_refinement = 0.
    for d, initial in [(.01, [-.2, -.4]), (.01, [-1.7, -3.2]), (.03, [-.2, -.4])]:
        solutions = []
        for rtol, atol in [(1e-10, 1e-12), (1e-12, 1e-14)]:
            sol = solve_ivp(
                lambda t, y: hydraulic_rhs(y, 0., d, p),
                (0., 100.), initial, rtol=rtol, atol=atol, max_step=.1,
                dense_output=True,
            )
            require(sol.success, sol.message)
            require(np.all(sol.y <= 0), "trajectory left nonpositive-potential domain")
            solutions.append(sol)
        ts = np.linspace(0., 100., 1001)
        error = float(np.max(np.abs(solutions[0].sol(ts)-solutions[1].sol(ts))))
        max_refinement = max(max_refinement, error)
        require(error < 2e-7, f"trajectory refinement failed: {error}")
        sol = solutions[-1]
        target = transpiring_equilibria(0., d, p)[0] if d == .01 and initial[0] > -1 else -4.
        require(abs(sol.y[1, -1]-target) < 1e-7, "leaf did not reach expected limiting state")
        report["trajectories"].append({
            "D_mol_per_mol": d, "initial_MPa": initial,
            "final_MPa": sol.y[:, -1].tolist(),
            "max_sampled_refinement_MPa": error, "expected_leaf_limit_MPa": target,
        })
        for i, ax in enumerate(axes):
            ax.plot(ts, sol.sol(ts)[i], label=f"D={d}, initial={initial}")
    for i, ax in enumerate(axes):
        ax.set(xlabel="Time (day; assumed capacitances)",
               ylabel=("Xylem" if i == 0 else "Leaf") + " potential (MPa)")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(args.output / "trajectories.png", dpi=180)
    plt.close(fig)
    require(maximum_residual < 1e-10, "equilibrium residual exceeds threshold")
    require(hydraulic_fold(0., safe) is None, "safe ordering incorrectly has an active fold")
    report["checks"] = {
        "max_equilibrium_rhs_residual_MPa_per_day": maximum_residual,
        "max_sampled_trajectory_refinement_MPa": max_refinement,
        "safe_ordering_has_no_active_fold": True,
    }
    sources = [
        PAPER, Path(__file__).resolve(),
        ROOT / "src/control_carbon/manzoni_hydraulics.py",
        ROOT / "tests/test_manzoni_hydraulics.py",
    ]
    report["manifest"] = {
        "python": platform.python_version(), "numpy": np.__version__,
        "scipy": scipy.__version__, "matplotlib": matplotlib.__version__,
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sources
        },
    }
    (args.output / "results.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8",
    )
    print(json.dumps({"output": str(args.output), "folds": report["folds"], "checks": report["checks"]}, indent=2))


if __name__ == "__main__":
    main()
