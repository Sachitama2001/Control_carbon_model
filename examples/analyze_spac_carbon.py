"""Stage A analytic checks and VPD/temperature controls; synthetic, not VISIT."""

import argparse
import csv
from dataclasses import asdict
from datetime import datetime
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.optimize import brentq

from control_carbon.manzoni_hydraulics import HydraulicParameters, hydraulic_fold, transpiring_equilibria
from control_carbon.manzoni_extensions import SigmoidShape, sigmoid_equilibria
from control_carbon.spac_carbon import (
    YEAR_DAYS, CarbonParameters, TemperatureResponse, carbon_rhs, carbon_budget,
    frozen_coefficient_capacity, reproduction_number, positive_carbon_equilibrium,
    hydraulic_factor, combined_rhs, combined_jacobian, simulate,
)

ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT/"artifacts/spac_carbon"/datetime.now().strftime("%Y%m%d_%H%M%S"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    H = HydraulicParameters(.0086, .47, -1.5, -2., .001, .001)
    C = CarbonParameters(1.2/YEAR_DAYS, (.3, .4, .3), tuple(k/YEAR_DAYS for k in [.8, .05, .2, .04]),
                         (.7, .9, .8), 10., .5, .25)
    response = TemperatureResponse(25., 25., 15., 2.)
    fold_d, fold_leaf = hydraulic_fold(0., H)
    fold_r0 = reproduction_number(hydraulic_factor(fold_leaf, H), C)
    extinction_d = (brentq(
        lambda d: reproduction_number(hydraulic_factor(transpiring_equilibria(0., d, H)[0], H), C)-1,
        .01, fold_d,
    ) if fold_r0 <= 1 else None)
    report = {"hydraulic": asdict(H), "carbon": asdict(C), "temperature": asdict(response),
              "classification": "model-agnostic derived analysis; all carbon and temperature coefficients synthetic",
              "hydraulic_fold_D": fold_d, "canopy_extinction_D": extinction_d,
              "R0_at_hydraulic_fold": fold_r0,
              "carbon_threshold_note": "No R0=1 crossing on hydrated branch before fold" if extinction_d is None else "R0=1 precedes hydraulic fold",
              "equilibria": [], "trajectories": []}
    branches = []
    for d in np.linspace(0., fold_d, 180):
        for k, leaf in enumerate(transpiring_equilibria(0., d, H)):
            h = hydraulic_factor(leaf, H)
            positive = positive_carbon_equilibrium(h, C)
            branches.append((d, leaf, reproduction_number(h, C),
                             float(np.sum(positive)) if positive is not None else 0., k))
    branches = np.array(branches)
    shape = SigmoidShape(6., 6.)
    for leaf in sigmoid_equilibria(-.5, .01, H, shape):
        h = hydraulic_factor(leaf, H, shape)
        positive = positive_carbon_equilibrium(h, C)
        stocks = np.zeros(4) if positive is None else positive
        y = np.r_[(-.5+leaf)/2, leaf, stocks]
        eig = np.linalg.eigvals(combined_jacobian(y, -.5, .01, H, C, shape=shape))
        residual = float(np.max(abs(combined_rhs(y, -.5, .01, H, C, shape=shape))))
        require(residual < 1e-10, "six-state equilibrium residual failed")
        report["equilibria"].append({
            "hydraulic_curve": "sigmoid alpha_x=alpha_s=6",
            "state": y.tolist(), "R0": reproduction_number(h, C),
            "eigenvalues_per_day": eig.real.tolist(),
            "stable": bool(np.max(eig.real) < 0), "rhs_residual": residual,
        })

    leaf0 = transpiring_equilibria(0., .01, H)[0]
    h0 = hydraulic_factor(leaf0, H)
    carbon0 = positive_carbon_equilibrium(h0, C)
    initial = np.r_[leaf0/2, leaf0, carbon0]
    temp_critical = brentq(
        lambda t: reproduction_number(h0, C, *response.factors(t))-1, 25., 45.,
    )
    report["temperature_extinction_C"] = temp_critical
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for k, style in [(0, "-"), (1, "--")]:
        rows = branches[branches[:, 4] == k]
        axes[0, 0].plot(rows[:, 0], rows[:, 1], style, color="black")
        axes[0, 1].plot(rows[:, 0], rows[:, 3], style, color="green")
    for ax in axes[0]:
        ax.axvline(fold_d, color="red", ls=":", label="hydraulic fold")
        if extinction_d is not None:
            ax.axvline(extinction_d, color="blue", ls=":", label="R0=1")
        ax.set(xlabel="D (mol/mol)")
        ax.legend(fontsize=8)
    axes[0, 0].set(ylabel="Leaf equilibrium (MPa)")
    axes[0, 1].set(ylabel="Conditional carbon equilibrium (kg C/m2)")

    max_refinement = max_budget = 0.
    for duration in [10., 100.]:
        def forcing(t):
            v = min(max(t/duration, 0.), 1.)
            d = .01+.02*(3*v*v-2*v*v*v)
            return 0., d, 1., 1.
        times = np.unique(np.r_[np.linspace(0., duration, 301),
                                duration+np.linspace(0., 200*YEAR_DAYS, 1001)])
        solutions = [simulate(initial, times, forcing, H, C, rtol=tol) for tol in [1e-8, 1e-10]]
        error = float(np.max(abs(solutions[0]["state"]-solutions[1]["state"])))
        max_refinement = max(max_refinement, error)
        require(error < 1e-5, "VPD refinement failed")
        states = solutions[-1]["state"]
        caps = []
        roots = []
        for t, state in zip(times, states):
            h = hydraulic_factor(state[1], H)
            budget = carbon_budget(state[2:], h, C)
            max_budget = max(max_budget, abs(budget["residual"]))
            caps.append(float(np.sum(frozen_coefficient_capacity(state[2:], h, C))))
            eq = positive_carbon_equilibrium(h, C)
            roots.append(float(np.sum(eq)) if eq is not None else 0.)
        require(abs(states[-1, 1]+4) < 1e-7, "post-fold leaf limit failed")
        require(np.sum(states[-1, 2:]) < .02, "post-failure carbon hold insufficient")
        axes[1, 0].plot(times/YEAR_DAYS, np.sum(states[:, 2:], axis=1), label=f"stocks, ramp={duration:g} d")
        axes[1, 0].plot(times/YEAR_DAYS, caps, "--", label=f"input-frozen capacity, {duration:g} d")
        report["trajectories"].append({
            "name": f"VPD_fold_crossing_{duration:g}days",
            "classification": "B-crossing + carbon threshold, NOT R-tipping",
            "time_days": times.tolist(), "state": states.tolist(), "capacity_total": caps,
            "conditional_attracting_carbon_total": roots, "max_refinement": error,
        })

    for duration in [1., 10.]:
        def forcing(t):
            if t <= duration*YEAR_DAYS:
                temp = 25+20*np.sin(np.pi*t/(duration*YEAR_DAYS))**2
            else:
                temp = 25.
            return 0., .01, *response.factors(temp)
        times = np.unique(np.r_[np.linspace(0., duration*YEAR_DAYS, 301),
                                duration*YEAR_DAYS+np.linspace(0., 500*YEAR_DAYS, 1001)])
        sol = simulate(initial, times, forcing, H, C, rtol=1e-10)
        require(np.max(abs(sol["state"][:, :2]-initial[:2])) < 1e-7, "temperature changed fixed hydraulic subsystem")
        require(np.max(abs(sol["state"][-1, 2:]-carbon0)) < 1e-5, "safe-final-temperature recovery failed")
        axes[1, 1].plot(times/YEAR_DAYS, sol["state"][:, 2], label=f"25-45-25 C pulse, {duration:g} yr")
        report["trajectories"].append({
            "name": f"temperature_pulse_{duration:g}years", "classification": "transient decline and recovery; no R-extinction",
            "time_days": times.tolist(), "state": sol["state"].tolist(),
        })
    def persistent_heat(t):
        v = min(max(t/(10*YEAR_DAYS), 0.), 1.)
        return 0., .01, *response.factors(25+20*(3*v*v-2*v*v*v))
    heat_times = np.unique(np.r_[np.linspace(0., 10*YEAR_DAYS, 301),
                                10*YEAR_DAYS+np.linspace(0., 200*YEAR_DAYS, 1001)])
    heat = simulate(initial, heat_times, persistent_heat, H, C, rtol=1e-10)
    require(np.max(abs(heat["state"][:, :2]-initial[:2])) < 1e-7, "persistent heat changed hydraulics")
    require(np.sum(heat["state"][-1, 2:]) < .02, "persistent high-temperature extinction control failed")
    report["trajectories"].append({
        "name": "persistent_temperature_25_to_45", "classification": "R0=1 boundary crossed; NOT R-tipping",
        "time_days": heat_times.tolist(), "state": heat["state"].tolist(),
    })
    axes[1, 0].set(xlabel="Time (year)", ylabel="Total carbon (kg C/m2)")
    axes[1, 1].set(xlabel="Time (year)", ylabel="Leaf carbon (kg C/m2)", yscale="log")
    for ax in axes[1]:
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(args.output/"stage_a.png", dpi=180)
    fig.savefig(args.output/"stage_a.pdf")
    plt.close(fig)
    temperatures = np.linspace(25., 45., 200)
    r0_values = [reproduction_number(h0, C, *response.factors(t)) for t in temperatures]
    totals = []
    for t in temperatures:
        positive = positive_carbon_equilibrium(h0, C, *response.factors(t))
        totals.append(float(np.sum(positive)) if positive is not None else 0.)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(temperatures, r0_values)
    axes[0].axhline(1., color="gray", ls=":")
    axes[0].set(xlabel="Temperature (C; synthetic response)", ylabel="R0")
    axes[1].plot(temperatures, totals)
    axes[1].set(xlabel="Temperature (C; synthetic response)", ylabel="Attracting carbon equilibrium (kg C/m2)")
    for ax in axes:
        ax.axvline(temp_critical, color="red", ls=":", label="R0=1, not empirical threshold")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(args.output/"temperature_threshold.png", dpi=180)
    fig.savefig(args.output/"temperature_threshold.pdf")
    plt.close(fig)
    report["temperature_branches"] = {"temperature_C": temperatures.tolist(), "R0": r0_values, "carbon_total": totals}
    report["branches"] = branches.tolist()
    report["checks"] = {"max_VPD_refinement": max_refinement, "max_carbon_budget_residual": max_budget}
    sources = [
        Path(__file__).resolve(), ROOT/"src/control_carbon/spac_carbon.py",
        ROOT/"src/control_carbon/manzoni_hydraulics.py",
        ROOT/"src/control_carbon/manzoni_extensions.py", ROOT/"tests/test_spac_carbon.py",
        ROOT/"docs/AGENT_HANDOFF_CWexperiment.md", ROOT/"docs/SPAC_carbon_research_proposal_2026.docx",
    ]
    report["manifest"] = {
        "numpy": np.__version__, "scipy": scipy.__version__,
        "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
    }
    (args.output/"results.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    with (args.output/"equilibria.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["psi_x_MPa", "psi_leaf_MPa", "leaf_kgC_m2", "wood_kgC_m2", "root_kgC_m2", "soil_kgC_m2",
                         "R0", "stable", "max_rhs_residual_mixed_units", "eigenvalues_day-1"])
        for eq in report["equilibria"]:
            writer.writerow([*eq["state"], eq["R0"], eq["stable"], eq["rhs_residual"], json.dumps(eq["eigenvalues_per_day"])])
    with (args.output/"parameters.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["name", "value", "unit", "provenance"])
        for name, value, unit, source in [
            ("GPPmax", C.assimilation_max, "kg C m-2 day-1", "synthetic; proposal value reused as gross input, not its effective NPP"),
            ("allocation", C.allocation, "fraction", "synthetic"),
            ("losses", C.losses, "day-1", "synthetic; proposal annual magnitudes divided by 365.25"),
            ("litter_fractions", C.litter_fractions, "fraction", "synthetic"),
            ("SLA", C.specific_leaf_area, "m2 leaf kg-dry-1", "synthetic proposal"),
            ("fC", C.carbon_dry_fraction, "kg C kg-dry-1", "synthetic proposal"),
            ("kappa", C.extinction, "dimensionless", "synthetic proposal"),
            ("T_reference", response.reference, "C", "synthetic"),
            ("T_optimum", response.optimum, "C", "synthetic"),
            ("T_width", response.width, "C", "synthetic"),
            ("Q10", response.loss_q10, "dimensionless", "synthetic"),
        ]:
            writer.writerow([name, value, unit, source])
    print(json.dumps({key: report[key] for key in [
        "hydraulic_fold_D", "canopy_extinction_D", "temperature_extinction_C", "equilibria", "checks",
    ]}, indent=2))


if __name__ == "__main__":
    main()
