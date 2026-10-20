"""Model A audit and bounded VPD-only six-state local ISS experiments."""

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

from control_carbon.spac_carbon import YEAR_DAYS
from control_carbon.spac_carbon_capacity import reference_scales
from control_carbon.spac_carbon_iss import (
    capacity_diagnostics, equilibrium_tracking_diagnostics,
)
from control_carbon.spac_carbon_rational import RationalSPAC
from control_carbon.spac_carbon_vpd_iss import (
    certify_local_vpd_iss, initial_energy_upper, local_iss_bounds,
    quartic_instability_obstruction,
)
from analyze_e4_capacity_iss import save, serializable
from prepare_manuscript_phase_portraits import model_family


ROOT = Path(__file__).resolve().parents[1]


def experiment(model, certificate, kind, amplitude, certified):
    end, ramp, period = 200*YEAR_DAYS, 30., 10*YEAR_DAYS
    day = np.unique(np.r_[np.linspace(0, 60, 301), np.linspace(0, end, 1601)])
    def vpd(time):
        if kind == "step":
            return certificate.vpd_mpa+amplitude
        if kind == "ramp":
            s = np.clip(time/ramp, 0., 1.)
            return certificate.vpd_mpa+amplitude*s*s*(3-2*s)
        if kind == "periodic":
            return certificate.vpd_mpa+amplitude*np.sin(2*np.pi*time/period)
        raise ValueError(f"Unknown VPD schedule: {kind}")
    def vpd_rate(time):
        if kind == "step":
            return 0.
        if kind == "ramp":
            s = time/ramp
            return amplitude*6*s*(1-s)/ramp if 0 < s < 1 else 0.
        return amplitude*2*np.pi/period*np.cos(2*np.pi*time/period)
    def solve(rtol, atol, max_step):
        run = solve_ivp(lambda t, y: model.rhs(y, vpd(t)), (0., end),
                        certificate.equilibrium, method="Radau", t_eval=day,
                        rtol=rtol, atol=atol, max_step=max_step)
        if not run.success or len(run.t) != len(day):
            raise RuntimeError(f"Six-state {kind} integration failed: {run.message}")
        return run.y.T
    state = solve(2e-12, 1e-14, YEAR_DAYS/4)
    refined = solve(2e-13, 1e-15, YEAR_DAYS/8)
    if (np.any(~np.isfinite(state)) or np.any(state[:, :2] > 0)
            or np.any(state[:, 2:] <= 0)):
        raise RuntimeError("Accepted full-system trajectory left the physical domain")
    scales = certificate.scales
    forcing = np.array([vpd(t) for t in day])
    forcing_rate = np.array([vpd_rate(t) for t in day])
    x = (state-certificate.equilibrium)/scales
    energy = np.einsum("ni,ij,nj->n", x, certificate.metric, x)
    diagnostics = [capacity_diagnostics(model, y, v) for y, v in zip(state, forcing)]
    capacity = np.array([d["capacity"] for d in diagnostics])
    carbon_error = state[:, 2:]-capacity
    cache = {}
    equilibria, equilibrium_rates = [], []
    for v, rate in zip(forcing, forcing_rate):
        if v not in cache:
            equilibrium = model.equilibria(v)[0]
            tangent = equilibrium_tracking_diagnostics(model, equilibrium, v)["tangent_per_mpa"]
            cache[v] = (equilibrium, tangent)
        equilibrium, tangent = cache[v]
        equilibria.append(equilibrium)
        equilibrium_rates.append(tangent*rate)
    equilibria = np.array(equilibria)
    equilibrium_rates = np.array(equilibrium_rates)
    fixed_norm = np.linalg.norm(x, axis=1)
    moving_norm = np.linalg.norm((state-equilibria)/scales, axis=1)
    error_norm = np.linalg.norm(carbon_error/scales[2:], axis=1)
    refinement = float(np.max(np.abs(state-refined)/scales))
    signal = float(np.max(fixed_norm))
    allowed_refinement = max(5e-13, 1e-3*signal)
    if refinement > allowed_refinement:
        raise RuntimeError(f"Tolerance refinement {refinement} exceeds {allowed_refinement}")
    residual = max(d["carbon_identity_residual"] for d in diagnostics)
    if residual > 1e-14:
        raise RuntimeError("Capacity identity failed")
    arrays = {
        "day": day, "state": state, "refined_state": refined,
        "vpd_mpa": forcing, "vpd_rate_mpa_per_day": forcing_rate,
        "capacity": capacity, "capacity_rate": np.array([d["capacity_rate"] for d in diagnostics]),
        "capacity_error": carbon_error, "potential": -carbon_error,
        "fixed_error_scaled": x, "moving_error_scaled": (state-equilibria)/scales,
        "moving_equilibrium": equilibria, "equilibrium_rate": equilibrium_rates,
        "rhs": np.array([d["rhs"] for d in diagnostics]),
        "fixed_norm_scaled": fixed_norm, "moving_norm_scaled": moving_norm,
        "capacity_error_norm_scaled": error_norm, "energy": energy,
        "fluxes": np.array([model.fluxes(y, v) for y, v in zip(state, forcing)]),
        "conductance": np.array([
            [model.base.hydraulic.gx_max*max(0., 1+y[1]/3),
             model.base.hydraulic.gs_max*max(0., 1+y[1]/4)] for y in state]),
    }
    summary = {
        "kind": kind, "certified_input": certified, "amplitude_mpa": amplitude,
        "ramp_days": ramp if kind == "ramp" else None,
        "period_days": period if kind == "periodic" else None,
        "duration_days": end, "samples": len(day),
        "refinement_max_scaled": refinement,
        "allowed_refinement_max_scaled": allowed_refinement,
        "max_fixed_norm_scaled": signal, "max_moving_norm_scaled": float(max(moving_norm)),
        "max_capacity_error_norm_scaled": float(max(error_norm)),
        "min_leaf_water_potential_mpa": float(min(state[:, 1])),
        "max_carbon_identity_residual": residual,
        "step_note": "Initial equilibrium jump; no classical tracking-rate claim at t=0"
                     if kind == "step" else None,
    }
    if certified:
        initial_energy = initial_energy_upper(certificate, certificate.equilibrium)
        bound = local_iss_bounds(
            certificate, day, initial_energy, certificate.max_vpd_perturbation_mpa)
        arrays.update({f"bound_{key}": value for key, value in bound.items()})
        for actual, upper in [(energy, bound["energy"]), (fixed_norm, bound["state_norm_scaled"]),
                              (error_norm, bound["capacity_error_norm_scaled"])]:
            if np.max(actual-upper) > 1e-12:
                raise RuntimeError("Nonlinear trajectory exceeded its certified bound")
        summary["initial_energy_upper"] = initial_energy
        summary["bound_input_sup_mpa"] = certificate.max_vpd_perturbation_mpa
        summary["max_energy_fraction"] = float(max(energy)/certificate.energy_limit)
        summary["max_state_bound_scaled"] = float(max(bound["state_norm_scaled"]))
    return arrays, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    output = parser.parse_args().output
    output.mkdir(parents=True, exist_ok=False)
    model = RationalSPAC(model_family()["E4"])
    certificate = certify_local_vpd_iss(
        model, .00244, reference_scales(model.base).state, .0015)
    records = {
        "scope": "plan1020 Phase 1-2, Model A; no soil or regrowth extension",
        "certificate": asdict(certificate),
        "unstable_3kpa": quartic_instability_obstruction(model, .003),
        "parameters": asdict(model.base),
        "units": {"time": "day", "VPD": "MPa", "potentials": "MPa", "carbon": "kg C m^-2",
                  "norms": "dimensionless, scaled by reference_scales"},
        "experiments": [],
        "baseline_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "versions": {"numpy": np.__version__, "scipy": __import__("scipy").__version__,
                     "matplotlib": matplotlib.__version__},
        "integrator": {"method": "Radau", "rtol": 2e-12, "atol": 1e-14,
                       "max_step_days": YEAR_DAYS/4, "refined_rtol": 2e-13,
                       "refined_atol": 1e-15, "refined_max_step_days": YEAR_DAYS/8},
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    runs = {}
    cases = [("certified_step", "step", .8*certificate.max_vpd_perturbation_mpa, True),
             ("certified_ramp", "ramp", .8*certificate.max_vpd_perturbation_mpa, True),
             ("certified_periodic", "periodic", .8*certificate.max_vpd_perturbation_mpa, True),
             ("exploratory_ramp", "ramp", 5e-6, False)]
    for name, kind, amplitude, certified in cases:
        arrays, summary = experiment(model, certificate, kind, amplitude, certified)
        runs[name] = arrays
        np.savez_compressed(output/f"{name}.npz", **arrays)
        records["experiments"].append({"name": name, **summary})

    fig, axes = plt.subplots(3, 3, figsize=(13, 10), constrained_layout=True)
    for col, (name, _, _, _) in enumerate(cases[:3]):
        run = runs[name]
        axes[0, col].plot(run["day"], (run["vpd_mpa"]-certificate.vpd_mpa)*1e6)
        axes[0, col].set(xlim=(0, 60), ylim=(-6e-5, 6e-5),
                         xlabel="Time (day)", ylabel="VPD anomaly (Pa)")
        axes[0, col].text(.04, .9, name.replace("certified_", ""), transform=axes[0, col].transAxes)
        for key, label in [("fixed_norm_scaled", "Fixed equilibrium"),
                           ("moving_norm_scaled", "Moving equilibrium"),
                           ("capacity_error_norm_scaled", "Capacity error")]:
            axes[1, col].plot(run["day"]/YEAR_DAYS, run[key], label=label)
        axes[1, col].set(xlabel="Time (year)", ylabel="Scaled error (-)")
        axes[1, col].legend(fontsize=7)
        axes[2, col].semilogy(run["day"]/YEAR_DAYS, run["bound_state_norm_scaled"],
                            label="Certified state bound")
        axes[2, col].semilogy(run["day"]/YEAR_DAYS, run["fixed_norm_scaled"],
                            label="Actual fixed error")
        axes[2, col].set(xlabel="Time (year)", ylabel="Scaled norm (-)")
        axes[2, col].legend(fontsize=7)
    save(fig, output, "local_vpd_iss_and_conservatism")

    fig, axes = plt.subplots(2, 3, figsize=(12, 6), constrained_layout=True)
    for col, (name, _, _, _) in enumerate(cases[:3]):
        run = runs[name]
        anomaly = (run["vpd_mpa"]-certificate.vpd_mpa)*1e6
        axes[0, col].plot(run["day"]/YEAR_DAYS, anomaly)
        axes[1, col].plot(run["day"], anomaly)
        axes[0, col].set(xlim=(0, 200), ylim=(-6e-5, 6e-5),
                         xlabel="Time (year)", ylabel="VPD anomaly (Pa)")
        axes[1, col].set(xlim=(0, 60), ylim=(-6e-5, 6e-5),
                         xlabel="Time (day)", ylabel="VPD anomaly (Pa)")
        axes[0, col].text(.04, .9, name.replace("certified_", ""),
                          transform=axes[0, col].transAxes)
    save(fig, output, "vpd_forcing_schedules")

    run = runs["exploratory_ramp"]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    axes[0, 0].plot(run["day"], run["vpd_mpa"]*1000)
    axes[0, 0].set(xlim=(0, 60), xlabel="Time (day)", ylabel="VPD (kPa)")
    for index, label in [(0, "Xylem"), (1, "Leaf")]:
        axes[0, 1].plot(run["day"]/YEAR_DAYS, run["state"][:, index], label=label)
    axes[0, 1].set(xlabel="Time (year)", ylabel="Water potential (MPa)")
    axes[0, 1].legend()
    for values, label in [(run["state"][:, 2], "Actual leaf carbon"),
                          (run["capacity"][:, 0], "Instantaneous capacity"),
                          (run["moving_equilibrium"][:, 2], "Full equilibrium")]:
        axes[1, 0].plot(run["day"]/YEAR_DAYS, values, label=label)
    axes[1, 0].set(xlabel="Time (year)", ylabel="Leaf carbon (kg C m$^{-2}$)")
    axes[1, 0].legend(fontsize=7)
    for key, label in [("fixed_norm_scaled", "Fixed equilibrium"),
                       ("moving_norm_scaled", "Moving equilibrium"),
                       ("capacity_error_norm_scaled", "Capacity error")]:
        axes[1, 1].plot(run["day"]/YEAR_DAYS, run[key], label=label)
    axes[1, 1].set(xlabel="Time (year)", ylabel="Scaled error (-)")
    axes[1, 1].legend(fontsize=7)
    save(fig, output, "outside_certificate_exploratory_response")

    sources = [
        Path(__file__), ROOT/"examples/analyze_e4_capacity_iss.py",
        ROOT/"src/control_carbon/spac_carbon_vpd_iss.py",
        ROOT/"src/control_carbon/spac_carbon_stability.py",
        ROOT/"src/control_carbon/spac_carbon_iss.py",
        ROOT/"src/control_carbon/spac_carbon_rational.py",
        ROOT/"src/control_carbon/spac_carbon.py",
        ROOT/"src/control_carbon/spac_carbon_capacity.py",
        ROOT/"src/control_carbon/spac_carbon_equilibria.py",
        ROOT/"src/control_carbon/spac_carbon_mpa.py",
        ROOT/"src/control_carbon/manzoni_hydraulics.py",
        ROOT/"src/control_carbon/spac_carbon_two_way.py",
        ROOT/"examples/prepare_manuscript_phase_portraits.py",
        ROOT/"configs/spac_carbon_stage_b_synthetic.json",
        ROOT/"configs/spac_carbon_pressure_review.json",
        ROOT/"tests/test_spac_carbon_vpd_iss.py",
        ROOT/"docs/plan1020.md",
    ]
    records["source_sha256"] = {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    (output/"results.json").write_text(json.dumps(serializable(records), indent=2, allow_nan=False)+"\n")
    provenance = {"source_sha256": records["source_sha256"], "output_sha256": {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(output.iterdir()) if p.is_file()}}
    (output/"provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")
    print(json.dumps(serializable({
        "output": str(output), "delta_vpd_kpa": certificate.max_vpd_perturbation_mpa*1000,
        "experiments": records["experiments"],
        "instability_delta3": records["unstable_3kpa"]["third_hurwitz_determinant"],
    }), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
