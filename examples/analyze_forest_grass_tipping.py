"""Regenerate checks/figures for the Japanese forest--grass analytic report."""
import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import scipy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from control_carbon.forest_grass_tipping import (
    Parameters, drift, drift_prime, interior_roots, fold_points, path_parameters,
    smooth_step, analytic_bounds, critical_duration, matched_solutions,
    ramp_solution, potential, mean_first_passage,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=ROOT / "artifacts/forest_grass_tipping/analytic_20260930")
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    base = Parameters()
    roots = interior_roots(base)
    tau = critical_duration(rtol=2e-12, atol=2e-14, max_step=.004)
    fw, bw = matched_solutions(tau, sensitivity=True, rtol=2e-12, atol=2e-14)
    target = (roots[1] + roots[2]) / 2
    report = dict(
        classification="Published cover model + model-agnostic analytic path, NOT VISIT",
        units="cover fractions; model time, not calibrated years or temperature",
        source="Kumar K & Dutta (2026), doi:10.1098/rspa.2025.0803, (2.3), (4.1)",
        baseline=asdict(base), baseline_roots=roots.tolist(),
        baseline_jacobian=drift_prime(roots).tolist(),
        folds=[dict(grass=g, alpha=alpha) for g, alpha in fold_points()],
        constructed_path=[dict(lam=x, parameters=asdict(path_parameters(x)),
                               roots=interior_roots(path_parameters(x)).tolist())
                          for x in [0., .5, 1.]],
        proof_bounds=analytic_bounds(), critical_duration=tau,
        critical_average_rate=1 / tau,
        matching_residual=float(fw.y[0, -1] - bw.y[0, -1]),
        matching_derivative=float(fw.y[1, -1] - bw.y[1, -1]),
        potential_barriers=(potential(roots[1]) - potential(roots[[0, 2]])).tolist(),
        mfpt=dict(target=target, reflecting_wall=0.,
                  values={str(s): mean_first_passage(roots[0], target, s)
                          for s in [.05, .04, .03]}),
        initialization="Exact past stable equilibrium g=.4; no transient spinup bias",
        state_clipping=False, noise="Reflected additive noise is an explicit added assumption",
        versions=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__),
    )
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    g = np.linspace(0, 1, 601)
    axes[0, 0].plot(g, drift(g), label="f(g)")
    axes[0, 0].axhline(0, c="grey", lw=.8)
    axes[0, 0].scatter(roots, [0, 0, 0], c=["tab:green", "tab:red", "tab:green"])
    axes[0, 0].set(xlabel="Grass cover g", ylabel="dg/dt", title="Two interior stable communities")
    for alpha in np.linspace(.92, 1.44, 301):
        rr = interior_roots(replace(base, alpha=alpha))
        colors = ["tab:red" if drift_prime(x, replace(base, alpha=alpha)) > 0
                  else "tab:green" for x in rr]
        axes[0, 1].scatter(np.full_like(rr, alpha), rr, c=colors, s=3)
    axes[0, 1].set(xlabel="Forest recruitment alpha", ylabel="Equilibrium grass cover",
                   title="B: folds under alpha-only forcing")
    v = potential(g)
    axes[1, 0].plot(g, v - potential(roots[0]))
    axes[1, 0].scatter(roots, potential(roots) - potential(roots[0]),
                       c=["tab:green", "tab:red", "tab:green"])
    axes[1, 0].set(xlabel="Grass cover g", ylabel="V(g) - V(low)",
                   title="N: noise escapes a potential well")
    u = np.linspace(0, 1, 401)
    branches = np.array([interior_roots(path_parameters(float(smooth_step(x)))) for x in u])
    for j, label in enumerate(["Low stable QSE", "Unstable boundary", "High stable QSE"]):
        axes[1, 1].plot(u, branches[:, j], "--" if j == 1 else ":",
                        c="black" if j == 1 else "grey", label=label)
    trajectories = []
    for factor, color in [(.5, "tab:orange"), (2., "tab:blue")]:
        sol = ramp_solution(factor * tau)
        states = sol.sol(u)[0]
        axes[1, 1].plot(u, states, c=color, label=f"duration = {factor} critical")
        trajectories.append(dict(duration=factor * tau, normalized_time=u.tolist(),
                                 grass=states.tolist(), endpoint=float(states[-1]),
                                 eventual_equilibrium=float(branches[-1, 2 if factor < 1 else 0])))
    axes[1, 1].set(xlabel="Normalized ramp time t/duration", ylabel="Grass cover g",
                   title="R: identical parameter path, different speeds")
    axes[1, 1].legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(out / "forest_grass_nbr.pdf")
    fig.savefig(out / "forest_grass_nbr.png", dpi=170)
    plt.close(fig)
    (out / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    (out / "trajectories.json").write_text(json.dumps(trajectories, indent=2) + "\n")
    sources = [Path(__file__), ROOT / "src/control_carbon/forest_grass_tipping.py",
               ROOT / "tests/test_forest_grass_tipping.py",
               ROOT / "docs/forest_grass_nbr_tipping.tex",
               ROOT / "docs/forest_grass_nbr_tipping.md"]
    snapshot = out / "source_snapshot"
    snapshot.mkdir()
    for path in sources:
        shutil.copy2(path, snapshot / path.name)
    visit = ROOT.parent / "VISIT-matrix"
    native = [visit / "visit_local" / name for name in [
        "radiation.c", "photosynthesis.c", "location_proc.c", "plant_proc.c",
        "daily_scheme.c", "experiment.c", "n_flux.c", "disturbance.c",
        "definition.h", "structure.h"]]
    pdfs = [ROOT / "docs" / name for name in ["rspa.2025.0803.pdf", "s12080-009-0067-z.pdf"]]
    hashes = lambda paths: {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (out / "manifest.json").write_text(json.dumps(dict(
        source_sha256=hashes(sources), provided_pdf_sha256=hashes(pdfs),
        inspected_native_sha256=hashes(native),
        visit_commit=subprocess.check_output(["git", "-C", str(visit), "rev-parse", "HEAD"],
                                              text=True).strip(),
        claims="Proofs in TeX; numerical checks are not computer-assisted rigorous interval proofs",
    ), indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
