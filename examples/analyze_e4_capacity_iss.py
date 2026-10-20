"""Bounded Phase 1--3 checks for rational E4 capacity tracking, not a ramp sweep."""

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from control_carbon.spac_carbon import YEAR_DAYS, carbon_matrix
from control_carbon.spac_carbon_capacity import reference_scales
from control_carbon.spac_carbon_iss import (
    assimilation_gradient, capacity_diagnostics, carbon_feedback_certificate,
    carbon_iss_certificate, componentwise_iss_bound, equilibrium_tracking_diagnostics,
)
from control_carbon.spac_carbon_rational import RationalSPAC
from prepare_manuscript_phase_portraits import model_family


ROOT = Path(__file__).resolve().parents[1]


def serializable(value):
    if isinstance(value, np.ndarray):
        return serializable(value.tolist())
    if isinstance(value, complex):
        return [value.real, value.imag]
    if isinstance(value, np.generic):
        return serializable(value.item())
    if isinstance(value, dict):
        return {k: serializable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serializable(v) for v in value]
    return value


def save(fig, output, stem):
    for extension in ["png", "pdf", "svg"]:
        fig.savefig(output/f"{stem}.{extension}", dpi=170)
    plt.close(fig)


def manufactured_input_check(model, period_years):
    """Independently integrate C and e with a prescribed sinusoidal A(t)."""
    carbon = model.base.carbon
    cert = carbon_iss_certificate(carbon)
    matrix = cert.matrix
    allocation = np.r_[carbon.allocation, 0.]
    direction = np.linalg.solve(-matrix, allocation)
    mean, amplitude = .2/YEAR_DAYS, .04/YEAR_DAYS
    omega = 2*np.pi/(period_years*YEAR_DAYS)
    day = np.linspace(0., 100*YEAR_DAYS, 1001)
    initial_carbon = direction*mean
    def forcing(time):
        return mean+amplitude*np.sin(omega*time)
    def input_rate(time):
        return direction*amplitude*omega*np.cos(omega*time)
    def solve(tolerance):
        run = solve_ivp(
            lambda t, z: np.r_[matrix@z[:4]+allocation*forcing(t),
                               matrix@z[4:]-input_rate(t)],
            (0., day[-1]), np.r_[initial_carbon, np.zeros(4)], t_eval=day,
            method="Radau", rtol=tolerance, atol=tolerance*.01,
            max_step=period_years*YEAR_DAYS/32)
        if not run.success:
            raise RuntimeError(f"Manufactured-input integration failed: {run.message}")
        return run.y.T
    states = solve(1e-10)
    refined = solve(1e-12)
    capacity = (mean+amplitude*np.sin(omega*day[:, None]))*direction
    capacity_rate = np.cos(omega*day[:, None])*direction*amplitude*omega
    error = states[:, :4]-capacity
    error_independent = states[:, 4:]
    rate_sup = direction*amplitude*omega
    bound = np.array([
        componentwise_iss_bound(carbon, t, np.zeros(4), rate_sup) for t in day
    ])
    l1_bound = -np.expm1(-cert.l1_decay_per_day*day)*sum(rate_sup)/cert.l1_decay_per_day
    l2_bound = (-np.expm1(-cert.euclidean_decay_per_day*day)
                *cert.euclidean_prefactor*np.linalg.norm(rate_sup)
                /cert.euclidean_decay_per_day)
    discrepancy = float(np.max(np.abs(error-error_independent)))
    refinement = float(np.max(np.abs(states-refined)))
    excess = float(np.max(np.abs(error)-bound))
    if discrepancy > 1e-9 or refinement > 1e-8 or excess > 1e-9:
        raise RuntimeError("Manufactured ISS check exceeded declared error tolerances")
    if np.any(np.linalg.norm(error, axis=1) > l2_bound+1e-9):
        raise RuntimeError("Euclidean ISS bound failed")
    if np.any(np.sum(np.abs(error), axis=1) > l1_bound+1e-9):
        raise RuntimeError("l1 ISS bound failed")
    return {
        "day": day, "carbon": states[:, :4], "capacity": capacity,
        "capacity_rate": capacity_rate, "error": error,
        "independently_integrated_error": error_independent,
        "componentwise_bound": bound, "l1_bound": l1_bound, "l2_bound": l2_bound,
        "period_years": period_years,
        "independent_error_max_difference": discrepancy,
        "refinement_max_difference": refinement, "bound_max_excess": excess,
        "scope": "Prescribed assimilation input; NOT a nonlinear E4 trajectory",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    model = RationalSPAC(model_family()["E4"])
    cert = carbon_iss_certificate(model.base.carbon)
    direction = np.linalg.solve(-cert.matrix, np.r_[model.base.carbon.allocation, 0.])
    records = {
        "scope": "Phase 1--3 intermediate report; unchanged E4; no new drought sweep",
        "baseline_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "branch": subprocess.check_output(
            ["git", "branch", "--show-current"], cwd=ROOT, text=True).strip(),
        "parameters": asdict(model.base),
        "units": {"time": "day", "VPD": "MPa", "potentials": "MPa",
                  "carbon": "kg C m^-2", "capacity_rate": "kg C m^-2 day^-1"},
        "certificate": asdict(cert),
        "capacity_direction_days": direction,
        "carbon_feedback_certificates": {
            "dry_region": carbon_feedback_certificate(model, .5625),
            "ratio_below_0_95": carbon_feedback_certificate(model, .95),
        },
        "derivative_checks": [], "equilibria": [], "manufactured_inputs": [],
    }
    for psi, leaf in [(-1.6, .07), (-2.1, .016), (-3., .03), (-3.5, .02), (-4., .02),
                      (-4.2, .01)]:
        state = np.r_[-.9, psi, leaf, 1.4, .27, 3.7]
        diag = capacity_diagnostics(model, state, .00244)
        rate = diag["rhs"]
        step = 1e-6/max(1., np.max(np.abs(rate)))
        numeric = (model.capacity(state+step*rate)-model.capacity(state-step*rate))/(2*step)
        difference = float(np.max(np.abs(numeric-diag["capacity_rate"])))
        if difference > 1e-8 or diag["closed_error_identity_residual"] > 1e-12:
            raise RuntimeError("Capacity derivative or error identity check failed")
        row = {"state": state, **diag, "directional_difference_step_days": step,
               "capacity_rate_finite_difference_max_error": difference}
        if psi not in [-3., -4.]:
            independent_gradient = np.eye(6)[2:]-np.linalg.solve(
                cert.matrix, model.jacobian(state, .00244)[2:])
            gradient_error = float(np.max(np.abs(
                independent_gradient-direction[:, None]*assimilation_gradient(model, state))))
            row["jacobian_capacity_gradient_max_error"] = gradient_error
            if gradient_error > 1e-11:
                raise RuntimeError("Independent capacity gradient disagrees")
        records["derivative_checks"].append(row)

    vpds = np.linspace(.00244, .003, 81)
    for vpd in vpds:
        roots = model.equilibria(vpd)
        if len(roots) != 1:
            raise RuntimeError("Baseline E4 positive root count changed")
        row = equilibrium_tracking_diagnostics(model, roots[0], vpd)
        row["vpd_mpa"] = vpd
        tangent_difference = float(np.max(np.abs(
            row["tangent_per_mpa"]-row["curve_tangent_per_mpa"])))
        if tangent_difference > 1e-7:
            raise RuntimeError("Implicit and eliminated equilibrium tangents disagree")
        row["tangent_independent_max_error"] = tangent_difference
        residual = row["capacity_coordinate_transform"]@row["jacobian"]-(
            row["transformed_jacobian"]@row["capacity_coordinate_transform"])
        row["similarity_residual"] = float(np.max(np.abs(residual)))
        if row["similarity_residual"] > 1e-10:
            raise RuntimeError("Error-coordinate similarity transform failed")
        records["equilibria"].append(row)

    fig, axes = plt.subplots(2, 1, figsize=(8, 7), constrained_layout=True)
    for key, label, color in [
        ("full_eigenvalues_per_day", "Full carbon-water system", "#D62728"),
        ("carbon_feedback_eigenvalues_per_day", "Carbon error; leaf water rate is input", "#0072B2"),
        ("hydraulic_frozen_carbon_eigenvalues_per_day", "Water subsystem; carbon held fixed", ".5"),
    ]:
        leading = [max(row[key].real) for row in records["equilibria"]]
        axes[0].plot(vpds*1000, leading, label=label, color=color)
    axes[0].axhline(0., color=".4", lw=.8)
    axes[0].set(ylabel="Leading real eigenvalue (day$^{-1}$)",
                 yscale="symlog", ylim=(-.5, .04))
    axes[0].set_yscale("symlog", linthresh=1e-4)
    axes[0].legend(fontsize=8)
    axes[1].plot(vpds*1000, [row["leaf_feedback_ratio"] for row in records["equilibria"]],
                 color="#0072B2")
    axes[1].axhline(1., color="#D62728", ls="--")
    axes[1].set(xlabel="VPD (kPa)", ylabel=r"$v_{K,L}\,\partial A/\partial C_L$ (1)")
    for ax in axes:
        ax.grid(alpha=.2)
    save(fig, args.output, "e4_carbon_iss_vs_full_feedback")

    leaves, psi = np.linspace(0., .25, 251), np.linspace(-4.2, 0., 253)
    leaf_grid, psi_grid = np.meshgrid(leaves, psi)
    ratio = (direction[0]*model.base.carbon.assimilation_max
             *model.base.carbon.canopy_coefficient*np.maximum(0., 1+psi_grid/4)
             /(1+model.base.carbon.canopy_coefficient*leaf_grid)**2)
    fig, ax = plt.subplots(figsize=(8, 5.3), constrained_layout=True)
    image = ax.contourf(leaf_grid, psi_grid, ratio, levels=np.linspace(0., 2.25, 19),
                        cmap="viridis")
    ax.contour(leaf_grid, psi_grid, ratio, levels=[1.], colors=["#D62728"], linewidths=2)
    ax.axhline(-3., color="white", ls="--")
    for index, label, color in [(0, "2.44 kPa positive equilibrium", "#0072B2"),
                                (-1, "3.00 kPa positive equilibrium", "#D62728")]:
        eq = records["equilibria"][index]["state"]
        ax.scatter(eq[2], eq[1], marker="*", s=100, color=color, edgecolors="white",
                   label=label)
    ax.set(xlabel=r"$C_L$ (kg C m$^{-2}$)", ylabel=r"$\psi_l$ (MPa)")
    ax.legend(fontsize=8)
    fig.colorbar(image, ax=ax, label=r"$v_{K,L}\,\partial A/\partial C_L$ (1)")
    save(fig, args.output, "e4_carbon_feedback_region")

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for col, period in enumerate([2., 20.]):
        data = manufactured_input_check(model, period)
        np.savez_compressed(args.output/f"manufactured_period_{period:g}yr.npz",
                            **{k: v for k, v in data.items() if isinstance(v, np.ndarray)})
        records["manufactured_inputs"].append({
            k: v for k, v in data.items() if not isinstance(v, np.ndarray)})
        years = data["day"]/YEAR_DAYS
        axes[0, col].plot(years, data["carbon"][:, 0], label=r"$C_L$")
        axes[0, col].plot(years, data["capacity"][:, 0], ls="--", label=r"$K_L$")
        axes[0, col].set(ylabel="Leaf carbon (kg C m$^{-2}$)",
                         title=f"Prescribed input; period = {period:g} years",
                         xlim=(0., 20.))
        axes[0, col].set_xlabel("Time (years; first 20 shown)")
        axes[1, col].plot(years, np.linalg.norm(data["error"], axis=1),
                          label=r"Actual $\|e\|_2$")
        axes[1, col].plot(years, np.linalg.norm(data["componentwise_bound"], axis=1),
                          ls="--", label="Componentwise ISS envelope (2-norm)")
        axes[1, col].plot(years, data["l2_bound"], ls=":", label="Lyapunov ISS bound")
        axes[1, col].set(xlabel="Time (years)", ylabel="Carbon error (kg C m$^{-2}$)")
        axes[1, col].set_yscale("symlog", linthresh=1e-4)
        for row in range(2):
            axes[row, col].grid(alpha=.2)
            axes[row, col].legend(fontsize=8)
    save(fig, args.output, "exogenous_capacity_tracking_validation")

    sources = [
        Path(__file__), ROOT/"src/control_carbon/spac_carbon_iss.py",
        ROOT/"src/control_carbon/spac_carbon_rational.py",
        ROOT/"src/control_carbon/spac_carbon.py",
        ROOT/"src/control_carbon/spac_carbon_tracking.py",
        ROOT/"src/control_carbon/spac_carbon_equilibria.py",
        ROOT/"src/control_carbon/spac_carbon_capacity.py",
        ROOT/"src/control_carbon/spac_carbon_mpa.py",
        ROOT/"src/control_carbon/manzoni_hydraulics.py",
        ROOT/"src/control_carbon/spac_carbon_two_way.py",
        ROOT/"examples/prepare_manuscript_phase_portraits.py",
        ROOT/"configs/spac_carbon_stage_b_synthetic.json",
        ROOT/"configs/spac_carbon_pressure_review.json",
    ]
    records["source_sha256"] = {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    records["created_at_utc"] = datetime.now(timezone.utc).isoformat()
    records["versions"] = {"numpy": np.__version__,
                           "scipy": __import__("scipy").__version__,
                           "matplotlib": matplotlib.__version__}
    (args.output/"results.json").write_text(
        json.dumps(serializable(records), indent=2, allow_nan=False)+"\n")
    provenance = {
        "source_sha256": records["source_sha256"],
        "output_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(args.output.iterdir()) if p.is_file()},
    }
    (args.output/"provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")
    print(json.dumps({
        "output": str(args.output), "equilibrium_samples": len(records["equilibria"]),
        "derivative_checks": len(records["derivative_checks"]),
        "manufactured_inputs": records["manufactured_inputs"],
        "l1_decay_per_year": cert.l1_decay_per_day*YEAR_DAYS,
        "lyapunov_c": cert.euclidean_prefactor,
        "lyapunov_alpha_per_year": cert.euclidean_decay_per_day*YEAR_DAYS,
    }, indent=2))


if __name__ == "__main__":
    main()
