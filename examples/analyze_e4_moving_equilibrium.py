"""VPD-rate tracking bounds and prescribed-soil/root-node comparisons."""

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

from analyze_e4_capacity_iss import save, serializable
from prepare_manuscript_phase_portraits import model_family
from control_carbon.manzoni_extensions import SoilRootParameters, soil_root_conductance
from control_carbon.spac_carbon import YEAR_DAYS, carbon_matrix
from control_carbon.spac_carbon_capacity import reference_scales
from control_carbon.spac_carbon_iss import assimilation_gradient
from control_carbon.spac_carbon_moving_iss import (
    certify_moving_equilibrium, moving_initial_energy_upper, moving_ramp_bounds,
    ramp_speed_upper,
)
from control_carbon.spac_carbon_rational import RationalSPAC
from control_carbon.spac_carbon_soil_boundary import SoilBoundaryRationalSPAC
from control_carbon.spac_carbon_tracking import SmoothVPDRamp, forcing_derivative


ROOT = Path(__file__).resolve().parents[1]


def run_ramp(model, ramp, *, tiny=False):
    plant = model.plant if isinstance(model, SoilBoundaryRationalSPAC) else model
    scales = reference_scales(plant.base).state
    initial_roots = model.equilibria(ramp.start_mpa)
    if len(initial_roots) != 1:
        raise RuntimeError("These experiments require a unique initial positive equilibrium")
    initial = initial_roots[0]
    end = ramp.duration_days+200*YEAR_DAYS
    day = np.unique(np.r_[
        np.linspace(0, min(30., ramp.duration_days), 151),
        np.linspace(0, ramp.duration_days, 241),
        ramp.duration_days+np.r_[0., np.geomspace(.001, 200*YEAR_DAYS, 601)]])
    def solve(refine):
        rtol, atol = ((2e-13, 1e-15) if refine else (2e-12, 1e-14)) if tiny else (
            (1e-11, 1e-13) if refine else (1e-9, 1e-11))
        max_step = YEAR_DAYS*(.125 if refine else .25) if tiny else YEAR_DAYS/(2 if refine else 1)
        ramp_run = solve_ivp(
            lambda t, y: model.rhs(y, ramp.value_rate(t)[0]),
            (0, ramp.duration_days), initial, method="Radau", dense_output=True,
            rtol=rtol, atol=atol, max_step=min(max_step, ramp.duration_days/(60 if refine else 30)))
        if not ramp_run.success:
            raise RuntimeError(f"VPD ramp integration failed: {ramp_run.message}")
        hold = solve_ivp(
            lambda t, y: model.rhs(y, ramp.end_mpa),
            (ramp.duration_days, end), ramp_run.y[:, -1], method="Radau", dense_output=True,
            rtol=rtol, atol=atol, max_step=max_step)
        if not hold.success:
            raise RuntimeError(f"Post-ramp integration failed: {hold.message}")
        state = np.empty((len(day), 6))
        during = day <= ramp.duration_days
        state[during], state[~during] = ramp_run.sol(day[during]).T, hold.sol(day[~during]).T
        return state
    state, refined = solve(False), solve(True)
    if (np.any(~np.isfinite(state)) or np.any(state[:, :2] > 0) or np.any(state[:, 2:] <= 0)):
        raise RuntimeError("Accepted trajectory left the physical domain")
    forcing = np.array([ramp.value_rate(t) for t in day])
    cache = {}
    target, tangent, spectral = [], [], []
    for v, _ in forcing:
        if v not in cache:
            roots = model.equilibria(float(v))
            if len(roots) != 1:
                raise RuntimeError("Moving reference lost its unique positive branch")
            equilibrium = roots[0]
            jac = model.jacobian(equilibrium, v)
            derivative = np.linalg.solve(jac, -forcing_derivative(equilibrium, v, plant.base))
            cache[v] = (equilibrium, derivative, max(np.linalg.eigvals(jac).real))
        equilibrium, derivative, eigenvalue = cache[v]
        target.append(equilibrium)
        tangent.append(derivative)
        spectral.append(eigenvalue)
    target, tangent = np.array(target), np.array(tangent)
    target_rate = tangent*forcing[:, 1, None]
    rhs = np.array([model.rhs(y, v) for y, (v, _) in zip(state, forcing)])
    capacity = np.array([model.capacity(y) for y in state])
    matrix = carbon_matrix(plant.base.carbon)
    direction = np.linalg.solve(-matrix, np.r_[plant.base.carbon.allocation, 0.])
    capacity_rate = np.array([
        direction*(assimilation_gradient(plant, y)@f) for y, f in zip(state, rhs)])
    moving_error = (state-target)/scales
    error = (state[:, 2:]-capacity)/scales[2:]
    fluxes = np.array([model.fluxes(y, v) for y, (v, _) in zip(state, forcing)])
    h = plant.base.hydraulic
    water_residual = np.max(np.abs(h.cx*rhs[:, 0]+h.cl*rhs[:, 1]-fluxes[:, 0]+fluxes[:, 2]))
    carbon_residual = np.max(np.abs(rhs[:, 2:]-(state[:, 2:]-capacity)@matrix.T))
    if max(water_residual, carbon_residual) > 1e-14:
        raise RuntimeError("Water/carbon balance identity failed")
    difference = float(np.max(np.abs(state-refined)/scales))
    signal = float(np.max(np.linalg.norm((state-initial)/scales, axis=1)))
    tolerance = max(1e-12, 1e-4*signal)
    if difference > tolerance:
        raise RuntimeError(f"Refinement {difference} exceeds declared limit {tolerance}")
    root = (np.array([model.root_potential(y) for y in state])
            if isinstance(model, SoilBoundaryRationalSPAC) else np.zeros(len(day)))
    root_residual = 0.
    if isinstance(model, SoilBoundaryRationalSPAC):
        gsr = soil_root_conductance(model.psi_soil_mpa, model.soil)
        root_residual = float(np.max(np.abs(gsr*(model.psi_soil_mpa-root)-fluxes[:, 0])))
        if root_residual > 1e-14:
            raise RuntimeError("Root-node continuity failed")
    quasi_static = np.array([
        np.linalg.solve(model.jacobian(y, v), rate)/scales
        for y, (v, _), rate in zip(target, forcing, target_rate)])
    arrays = {
        "day": day, "state": state, "refined_state": refined, "vpd_mpa": forcing[:, 0],
        "vpd_rate_mpa_per_day": forcing[:, 1], "equilibrium": target,
        "equilibrium_tangent_per_mpa": tangent, "equilibrium_rate": target_rate,
        "rhs": rhs, "capacity": capacity, "capacity_rate": capacity_rate,
        "potential": capacity-state[:, 2:], "capacity_error_scaled": error,
        "moving_error_scaled": moving_error, "fixed_error_scaled": (state-initial)/scales,
        "moving_error_rate_scaled": (rhs-target_rate)/scales,
        "quasi_static_lag_scaled": quasi_static,
        "equilibrium_speed_scaled_per_day": np.linalg.norm(target_rate/scales, axis=1),
        "state_speed_scaled_per_day": np.linalg.norm(rhs/scales, axis=1),
        "capacity_speed_scaled_per_day": np.linalg.norm(capacity_rate/scales[2:], axis=1),
        "moving_error_norm_scaled": np.linalg.norm(moving_error, axis=1),
        "capacity_error_norm_scaled": np.linalg.norm(error, axis=1),
        "spectral_abscissa_per_day": np.array(spectral), "root_potential_mpa": root,
        "fluxes": fluxes,
    }
    peak = int(np.argmax(arrays["moving_error_norm_scaled"]))
    summary = {
        "ramp": asdict(ramp), "samples": len(day), "end_day": end,
        "max_vpd_rate_mpa_per_day": ramp_speed_upper(ramp),
        "max_equilibrium_speed_scaled_per_day": float(max(arrays["equilibrium_speed_scaled_per_day"])),
        "max_moving_error_norm_scaled": float(max(arrays["moving_error_norm_scaled"])),
        "peak_error_day": float(day[peak]),
        "peak_after_ramp": bool(day[peak] > ramp.duration_days),
        "max_capacity_error_norm_scaled": float(max(arrays["capacity_error_norm_scaled"])),
        "max_spectral_abscissa_per_day": float(max(spectral)),
        "refinement_max_scaled": difference, "refinement_limit_scaled": tolerance,
        "min_leaf_potential_mpa": float(min(state[:, 1])),
        "water_balance_max_residual": float(water_residual),
        "carbon_capacity_identity_max_residual": float(carbon_residual),
        "root_node_max_residual": root_residual,
        "evidence": "C: ordinary nonlinear ODE integration, not interval integration",
    }
    return arrays, summary


def certified_ramp(model, certificate, ramp):
    arrays, summary = run_ramp(model, ramp, tiny=True)
    initial = moving_initial_energy_upper(model, certificate, arrays["state"][0], ramp.start_mpa)
    bounds = moving_ramp_bounds(certificate, ramp, arrays["day"], initial)
    arrays.update({f"bound_{key}": value for key, value in bounds.items()})
    error = arrays["moving_error_scaled"]
    arrays["metric_radius"] = np.sqrt(np.einsum("ni,ij,nj->n", error, certificate.metric, error))
    for actual, bound in [(arrays["metric_radius"], bounds["metric_radius"]),
                          (arrays["moving_error_norm_scaled"], bounds["moving_norm_scaled"]),
                          (arrays["capacity_error_norm_scaled"], bounds["capacity_error_norm_scaled"])]:
        if np.max(actual-bound) > 1e-11:
            raise RuntimeError("Certified moving-ramp envelope failed")
    return arrays, {
        **summary, "initial_energy_upper": initial,
        "rate_condition": ramp_speed_upper(ramp) <= certificate.max_vpd_rate_mpa_per_day,
        "max_metric_radius": float(max(arrays["metric_radius"])),
        "max_bound_metric_radius": float(max(bounds["metric_radius"])),
        "max_tube_fraction": float(max(arrays["metric_radius"])**2/certificate.energy_limit),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    output = parser.parse_args().output
    output.mkdir(parents=True, exist_ok=False)
    plant = RationalSPAC(model_family()["E4"])
    scales = reference_scales(plant.base).state
    certificate = certify_moving_equilibrium(plant, (.0024399, .0024401), scales)
    records = {
        "scope": "Phase 3 plus prescribed-soil/root-node connection; no soil/root storage ODE",
        "certificate": asdict(certificate), "plant_parameters": asdict(plant.base),
        "baseline_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "units": {"time": "day", "VPD": "MPa", "potential": "MPa",
                  "carbon": "kg C m^-2", "scaled_speed": "day^-1"},
        "normal_integrator": {"method": "Radau", "rtol": 1e-9, "atol": 1e-11,
                              "refined_rtol": 1e-11, "refined_atol": 1e-13,
                              "max_step_day": YEAR_DAYS, "refined_max_step_day": YEAR_DAYS/2,
                              "ramp_max_step": "min(max_step, ramp_duration/30); refinement /60"},
        "tiny_integrator": {"rtol": 2e-12, "atol": 1e-14, "refined_rtol": 2e-13,
                            "refined_atol": 1e-15, "max_step_day": YEAR_DAYS/4,
                            "refined_max_step_day": YEAR_DAYS/8},
        "experiments": [], "soil_comparisons": [],
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "versions": {"numpy": np.__version__, "scipy": __import__("scipy").__version__,
                     "matplotlib": matplotlib.__version__},
    }
    cases = [(10., "10 day"), (YEAR_DAYS, "1 year"), (10*YEAR_DAYS, "10 year"),
             (100*YEAR_DAYS, "100 year")]
    runs = []
    for duration, label in cases:
        ramp = SmoothVPDRamp(.00244, .00254, duration)
        arrays, summary = run_ramp(plant, ramp)
        if summary["max_spectral_abscissa_per_day"] >= 0:
            raise RuntimeError("Exploratory branch unexpectedly became unstable")
        name = f"ramp_{duration:g}day"
        np.savez_compressed(output/f"{name}.npz", **arrays)
        records["experiments"].append({"name": name, "label": label, "certified": False, **summary})
        runs.append(arrays)

    certified_runs = []
    records["certified_experiments"] = []
    for duration, label in [(10., "10 day"), (100*YEAR_DAYS, "100 year")]:
        ramp = SmoothVPDRamp(.00244, .00244+5e-11, duration)
        tiny, summary = certified_ramp(plant, certificate, ramp)
        np.savez_compressed(output/f"certified_{duration:g}day.npz", **tiny)
        certified_runs.append((tiny, label))
        records["certified_experiments"].append({"label": label, **summary})

    fig, axes = plt.subplots(3, 4, figsize=(16, 10), constrained_layout=True)
    for col, (run, (_, label)) in enumerate(zip(runs, cases)):
        years = run["day"]/YEAR_DAYS
        for values, name, style in [(run["state"][:, 2], "Actual leaf carbon", "-"),
                                     (run["capacity"][:, 0], "Instantaneous capacity", "--"),
                                     (run["equilibrium"][:, 2], "Full equilibrium", ":")]:
            axes[0, col].plot(years, values, style, label=name)
        axes[0, col].text(.06, .9, label, transform=axes[0, col].transAxes)
        axes[0, col].set(ylabel="Leaf carbon (kg C m$^{-2}$)")
        axes[0, col].set_ylim(.055, .070)
        axes[0, col].legend(fontsize=7)
        axes[1, col].plot(years, run["moving_error_norm_scaled"], label="Moving equilibrium error")
        axes[1, col].plot(years, run["capacity_error_norm_scaled"], label="Capacity error")
        axes[1, col].set(ylabel="Scaled error (-)")
        axes[1, col].set_ylim(0, .115)
        axes[1, col].legend(fontsize=7)
        for key, name in [("state_speed_scaled_per_day", "State"),
                          ("equilibrium_speed_scaled_per_day", "Equilibrium"),
                          ("capacity_speed_scaled_per_day", "Capacity")]:
            values = run[key]
            axes[2, col].semilogy(years, np.where(values > 0, values, np.nan), label=name)
        axes[2, col].set(ylabel="Scaled speed (day$^{-1}$)")
        axes[2, col].set_ylim(1e-9, 1e-1)
        axes[2, col].legend(fontsize=7)
        for row in range(3):
            axes[row, col].set(xlim=(0, 300), xlabel="Time (year)")
    save(fig, output, "moving_equilibrium_capacity_and_speeds")

    fig, axes = plt.subplots(2, 1, figsize=(9, 6), constrained_layout=True)
    for run, (_, label) in zip(runs, cases):
        axes[0].plot(run["day"]/YEAR_DAYS, run["vpd_mpa"]*1000, label=label)
        axes[1].plot(run["day"], run["vpd_mpa"]*1000, label=label)
    axes[0].set(xlim=(0, 120), xlabel="Time (year)")
    axes[1].set(xlim=(0, 30), xlabel="Time (day)")
    for axis in axes:
        axis.set(ylim=(2.435, 2.545), ylabel="VPD (kPa)")
        axis.legend(fontsize=8)
    save(fig, output, "vpd_ramp_schedules")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    peaks, errors, cap_errors = [], [], []
    for run, (duration, label) in zip(runs, cases):
        inside = (run["day"] > 0) & (run["day"] < duration)
        axes[0].loglog(run["equilibrium_speed_scaled_per_day"][inside],
                       run["state_speed_scaled_per_day"][inside], label=label)
        peaks.append(max(run["equilibrium_speed_scaled_per_day"]))
        errors.append(max(run["moving_error_norm_scaled"]))
        cap_errors.append(max(run["capacity_error_norm_scaled"]))
    axes[0].set(xlabel="Equilibrium speed (day$^{-1}$)", ylabel="State speed (day$^{-1}$)")
    axes[0].legend(fontsize=8)
    axes[1].loglog(peaks, errors, "o-", label="Maximum moving error")
    axes[1].loglog(peaks, cap_errors, "s--", label="Maximum capacity error")
    axes[1].set(xlabel="Maximum equilibrium speed (day$^{-1}$)", ylabel="Maximum scaled error (-)")
    axes[1].legend(fontsize=8)
    save(fig, output, "speed_plane_and_peak_errors")

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    axes[0, 0].set(ylabel="VPD anomaly (Pa)")
    for (tiny, label), color in zip(certified_runs, ["tab:blue", "tab:orange"]):
        years = tiny["day"]/YEAR_DAYS
        axes[0, 0].plot(years, (tiny["vpd_mpa"]-.00244)*1e6, color=color, label=label)
        for axis, actual, upper, ylabel in [
                (axes[0, 1], "metric_radius", "bound_metric_radius", "Metric error radius (-)"),
                (axes[1, 0], "moving_error_norm_scaled", "bound_moving_norm_scaled", "Scaled moving error (-)"),
                (axes[1, 1], "capacity_error_norm_scaled", "bound_capacity_error_norm_scaled", "Scaled capacity error (-)")]:
            axis.semilogy(years, np.where(tiny[actual] > 0, tiny[actual], np.nan),
                           color=color, label=f"{label}: actual")
            axis.semilogy(years, tiny[upper], "--", color=color, label=f"{label}: bound")
            axis.set(ylabel=ylabel)
    for axis in axes.ravel():
        axis.set(xlabel="Time (year)", xlim=(0, 300))
        axis.legend(fontsize=7)
    save(fig, output, "certified_moving_tube")

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    soil_parameters = SoilRootParameters(.1, -.1)
    for psi, color in [(-.1, "tab:orange"), (-.3, "tab:green")]:
        coupled = SoilBoundaryRationalSPAC(plant, soil_parameters, psi)
        arrays, summary = run_ramp(coupled, SmoothVPDRamp(.00244, .00254, 10*YEAR_DAYS))
        name = f"soil_{abs(psi):g}mpa"
        np.savez_compressed(output/f"{name}.npz", **arrays)
        records["soil_comparisons"].append({
            "name": name, "soil_potential_mpa": psi, "soil_parameters": asdict(soil_parameters),
            "coefficient_scope": "Existing user-approved sensitivity coefficients, NOT paper calibration",
            "storage_scope": "Prescribed soil potential; algebraic root, no added storage state",
            "initial_equilibrium": arrays["state"][0], **summary})
        years = arrays["day"]/YEAR_DAYS
        label = f"Soil {psi:g} MPa"
        axes[0, 0].plot(years, arrays["root_potential_mpa"], color=color, label=label)
        axes[0, 0].axhline(psi, color=color, ls=":", lw=.8)
        axes[0, 1].plot(years, arrays["state"][:, 1], color=color, label=label)
        axes[1, 0].plot(years, arrays["state"][:, 2], color=color, label=label)
        axes[1, 0].plot(years, arrays["equilibrium"][:, 2], color=color, ls=":")
        axes[1, 1].plot(years, arrays["moving_error_norm_scaled"], color=color, label=label)
    baseline = runs[2]
    for axis, values in [(axes[0, 0], baseline["root_potential_mpa"]),
                         (axes[0, 1], baseline["state"][:, 1]), (axes[1, 0], baseline["state"][:, 2]),
                         (axes[1, 1], baseline["moving_error_norm_scaled"])]:
        axis.plot(baseline["day"]/YEAR_DAYS, values, color="tab:blue", label="Model A, root 0 MPa")
    for axis, ylabel in zip(axes.ravel(), ["Root potential (MPa)", "Leaf potential (MPa)",
                                         "Leaf carbon (kg C m$^{-2}$)", "Scaled moving error (-)"]):
        axis.set(xlabel="Time (year)", ylabel=ylabel, xlim=(0, 210))
        axis.legend(fontsize=7)
    save(fig, output, "prescribed_soil_root_connection")

    sources = [
        Path(__file__), ROOT/"examples/analyze_e4_capacity_iss.py",
        ROOT/"examples/prepare_manuscript_phase_portraits.py",
        ROOT/"examples/reproduce_manzoni2014_remaining.py",
        *[ROOT/f"src/control_carbon/{name}.py" for name in [
            "spac_carbon_moving_iss", "spac_carbon_soil_boundary", "spac_carbon_vpd_iss",
            "spac_carbon_rational", "spac_carbon_stability", "spac_carbon_tracking",
            "spac_carbon_iss", "spac_carbon_capacity", "spac_carbon_equilibria",
            "spac_carbon_mpa", "spac_carbon_two_way", "spac_carbon", "manzoni_hydraulics",
            "manzoni_extensions"]],
        ROOT/"configs/spac_carbon_stage_b_synthetic.json",
        ROOT/"configs/spac_carbon_pressure_review.json",
        ROOT/"tests/test_spac_carbon_moving_iss.py", ROOT/"tests/test_spac_carbon_soil_boundary.py",
        ROOT/"docs/plan1020.md",
    ]
    records["source_sha256"] = {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    (output/"results.json").write_text(json.dumps(serializable(records), indent=2, allow_nan=False)+"\n")
    provenance = {"source_sha256": records["source_sha256"], "output_sha256": {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(output.iterdir()) if p.is_file()}}
    (output/"provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")
    print(json.dumps(serializable({
        "output": str(output), "certified_rate_mpa_per_day": certificate.max_vpd_rate_mpa_per_day,
        "certified_total_variation_mpa": certificate.max_total_vpd_variation_mpa,
        "rate_comparison": [
            {key: row[key] for key in ["label", "max_equilibrium_speed_scaled_per_day",
                                       "max_moving_error_norm_scaled", "max_capacity_error_norm_scaled",
                                       "refinement_max_scaled", "peak_after_ramp"]}
            for row in records["experiments"]],
        "certified_experiments": [
            {key: row[key] for key in ["label", "rate_condition", "max_metric_radius",
                                       "max_bound_metric_radius", "refinement_max_scaled"]}
            for row in records["certified_experiments"]],
        "soil_comparisons": [
            {key: row[key] for key in ["soil_potential_mpa", "root_node_max_residual",
                                       "max_moving_error_norm_scaled"]}
            for row in records["soil_comparisons"]],
    }), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
