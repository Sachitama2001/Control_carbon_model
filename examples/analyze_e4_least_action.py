"""The four bounded first experiments in docs/leastaction.md, without changing E4."""

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
from time import perf_counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad, simpson, solve_ivp
from scipy.linalg import expm, solve_triangular
from scipy.optimize import LinearConstraint, minimize

from analyze_e4_capacity_iss import save, serializable
from prepare_manuscript_phase_portraits import model_family
from control_carbon.carbon_compensation import compensation_gramian
from control_carbon.spac_carbon import YEAR_DAYS
from control_carbon.spac_carbon_potential import (
    certify_drying_step_potential, endpoint_peak_lower, potential_diagnostics, potential_structure,
)
from control_carbon.spac_carbon_rational import RationalSPAC


ROOT = Path(__file__).resolve().parents[1]
PREVIOUS = ROOT/"artifacts/spac_carbon/rational_e4_moving_20261020/phase3_v3"


def discrete_optimum(matrix, inputs, weight, duration, terminal, intervals):
    """Independent piecewise-constant convex QP, solved by SLSQP."""
    start = perf_counter()
    n, controls = inputs.shape
    dt = duration/intervals
    augmented = np.zeros((n+controls, n+controls))
    augmented[:n, :n], augmented[:n, n:] = matrix, inputs
    gain = expm(augmented*dt)[:n, n:]
    cholesky = np.linalg.cholesky(weight)
    reach = np.concatenate([
        solve_triangular(cholesky, (expm(matrix*(duration-(i+1)*dt))@gain).T, lower=True).T
        /np.sqrt(dt) for i in range(intervals)], axis=1)
    left, values, right = np.linalg.svd(reach, full_matrices=False)
    if min(values) <= np.finfo(float).eps*max(reach.shape)*max(values):
        raise RuntimeError("Independent discrete benchmark lost resolved controllability")
    constraint = LinearConstraint(right, -(left.T@terminal)/values, -(left.T@terminal)/values)
    run = minimize(lambda x: .5*x@x, np.zeros(intervals*controls),
                   jac=lambda x: x, constraints=[constraint], method="SLSQP",
                   options={"ftol": 1e-12, "maxiter": 30})
    if not run.success:
        raise RuntimeError(f"Independent compensation QP failed: {run.message}")
    residual = float(np.linalg.norm(reach@run.x+terminal))
    if residual > 1e-8*max(1., np.linalg.norm(terminal)):
        raise RuntimeError("Discrete compensation terminal constraint failed")
    return {"intervals": intervals, "action": float(run.fun), "terminal_residual": residual,
            "optimizer_iterations": run.nit, "formation_and_optimization_seconds": perf_counter()-start}


def linear_benchmark(structure, duration, inputs, weight):
    matrix, direction = structure.matrix, structure.capacity_direction_days
    amplitude = .05/YEAR_DAYS
    def input_rate(t):
        s = t/duration
        return direction*amplitude*6*s*(1-s)/duration
    start = perf_counter()
    natural = solve_ivp(lambda t, p: matrix@p+input_rate(t), (0, duration), np.zeros(4),
                        method="DOP853", dense_output=True, rtol=2e-12, atol=1e-14)
    if not natural.success:
        raise RuntimeError(f"External-capacity integration failed: {natural.message}")
    natural_time = perf_counter()-start
    terminal = natural.y[:, -1]
    start = perf_counter()
    problem = compensation_gramian(matrix, inputs, weight, duration)
    preprocessing = perf_counter()-start
    start = perf_counter()
    minimum = problem.minimum_action(terminal)
    query = perf_counter()-start
    refined = compensation_gramian(matrix, inputs, weight, duration, quadrature_order=128).minimum_action(terminal)
    action_error = abs(refined.action-minimum.action)/minimum.action
    integral, quadrature_error = quad(
        lambda t: .5*minimum.control(t)@weight@minimum.control(t), 0, duration,
        epsabs=1e-13, epsrel=1e-10)
    day = np.linspace(0, duration, 401)
    controls = np.array([minimum.control(t) for t in day])
    compensated = solve_ivp(lambda t, p: matrix@p+input_rate(t)+inputs@minimum.control(t),
                            (0, duration), np.zeros(4), method="DOP853", t_eval=day,
                            rtol=2e-12, atol=1e-13)
    if not compensated.success:
        raise RuntimeError(f"Compensated integration failed: {compensated.message}")
    terminal_residual = float(np.linalg.norm(compensated.y[:, -1]))
    discrete = [discrete_optimum(matrix, inputs, weight, duration, terminal, count) for count in [32, 64, 128]]
    gap = [(row["action"]/minimum.action-1) for row in discrete]
    if (action_error > 1e-8 or abs(integral/minimum.action-1) > 1e-8 or terminal_residual > 1e-8
            or min(gap) < -1e-8 or not gap[-1] < gap[0] or gap[-1] > .005):
        raise RuntimeError("Gramian, independent optimization or direct integration did not agree")
    budget = minimum.action*np.unique(np.r_[0., np.geomspace(1e-10, 1, 301), np.linspace(0, 1, 101)])
    budget_lower = np.array([
        problem.output_budget_interval(terminal, np.ones(4), energy)["unavoidable_absolute_lower"]
        for energy in budget])
    adjoint_outputs = {}
    for label, weights in [
            ("total_potential", np.ones(4)), ("leaf_potential", np.eye(4)[0]),
            ("soil_potential", np.eye(4)[3]), ("total_carbon_rate", -matrix.T@np.ones(4))]:
        integrated, _ = quad(
            lambda t: weights@expm(matrix*(duration-t))@input_rate(t),
            0, duration, epsabs=1e-13, epsrel=1e-10)
        actual = float(weights@terminal)
        if abs(integrated-actual) > 1e-9*max(abs(actual), 1e-8):
            raise RuntimeError("Rank-one adjoint identity disagrees with independent natural ODE")
        adjoint_outputs[label] = {"terminal_output": actual, "adjoint_integral": integrated}
    arrays = {"day": day, "assimilation_kg_per_day": .2/YEAR_DAYS+amplitude*(day/duration)**2*(3-2*day/duration),
              "natural_potential": natural.sol(day).T, "compensated_potential": compensated.y.T,
              "control": controls, "gramian": problem.gramian,
              "energy_budget": budget, "total_output_budget_lower": budget_lower}
    compensated_carbon = arrays["assimilation_kg_per_day"][:, None]*direction-compensated.y.T
    total_output = float(sum(terminal))
    total_variance = float(np.ones(4)@problem.gramian@np.ones(4))
    record = {
        "duration_days": duration, "B": inputs, "R": weight, "action": minimum.action,
        "natural_terminal_potential": terminal, "direct_terminal_residual_l2": terminal_residual,
        "action_direct_quadrature": integral, "quadrature_error_estimate": quadrature_error,
        "quadrature_refinement_relative_error": action_error, "resolved_rank": problem.resolved_rank,
        "gramian_condition": float((problem.singular_values[0]/problem.singular_values[-1])**2),
        "discrete_optimizations": discrete, "discrete_relative_gaps": gap,
        "adjoint_outputs": adjoint_outputs,
        "min_compensated_carbon_per_pool": compensated_carbon.min(axis=0),
        "max_absolute_control_per_input": np.abs(controls).max(axis=0),
        "minimum_control_per_input": controls.min(axis=0),
        "action_to_cancel_total_output_only": total_output**2/(2*total_variance),
        "timing_seconds": {"natural_four_pool_ode": natural_time,
                           "gramian_preprocessing": preprocessing, "terminal_query": query},
        "evidence": "C: external capacity only; signed, unconstrained virtual fluxes",
    }
    return arrays, record


def analyze_ramp(model, path, previous_provenance):
    start = perf_counter()
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != previous_provenance["output_sha256"][path.name]:
        raise RuntimeError(f"Reused trajectory hash mismatch: {path.name}")
    with np.load(path) as source:
        day, state, vpd = source["day"], source["state"], source["vpd_mpa"]
        refined_state = source["refined_state"]
        old_potential, old_capacity = source["potential"], source["capacity"]
    loading_time = perf_counter()-start
    start = perf_counter()
    structure = potential_structure(model)
    diagnostics = [potential_diagnostics(model, y, v) for y, v in zip(state, vpd)]
    potential = np.array([row["potential"] for row in diagnostics])
    capacity = np.array([row["capacity"] for row in diagnostics])
    assimilation = np.array([row["assimilation"] for row in diagnostics])
    assimilation_rate = np.array([row["assimilation_rate"] for row in diagnostics])
    if max(np.max(np.abs(potential-old_potential)), np.max(np.abs(capacity-old_capacity))) > 1e-13:
        raise RuntimeError("Current E4 diagnostic disagrees with reused trajectory")
    lower = np.array([endpoint_peak_lower(structure.matrix, t, k-capacity[0], potential[0])
                      for t, k in zip(day, capacity)])
    general_lower = np.array([endpoint_peak_lower(structure.matrix, t, k-capacity[0])
                              for t, k in zip(day, capacity)])
    norm = np.abs(potential).sum(axis=1)
    observed_peak = np.maximum.accumulate(norm)
    history_lower = np.maximum.accumulate(lower)
    if np.any(lower > observed_peak+1e-10):
        raise RuntimeError("Endpoint lower bound exceeds observed peak; investigate output resolution")
    duration = float(path.stem.removeprefix("ramp_").removesuffix("day"))
    indices = np.flatnonzero(day == duration)
    if len(indices) != 1:
        raise RuntimeError("Reused ramp does not contain its exact end time")
    end_index = int(indices[0])
    windows = []
    refined_potential = np.array([model.capacity(y)-y[2:] for y in refined_state])
    for label, index in [("ramp_end", end_index), ("hold_end", len(day)-1)]:
        windows.append({
            "window": label, "T_days": day[index], "capacity_endpoint_change_l1": np.abs(capacity[index]-capacity[0]).sum(),
            "peak_lower_l1": lower[index], "unknown_initial_peak_lower_l1": general_lower[index],
            "sampled_peak_l1": observed_peak[index], "final_potential_l1": norm[index],
            "fraction_of_sampled_peak": lower[index]/observed_peak[index],
            "history_peak_lower_l1": history_lower[index],
            "history_fraction_of_sampled_peak": history_lower[index]/observed_peak[index],
        })
    evaluation_time = perf_counter()-start
    # This is a posteriori signed source analysis, not a sign certificate.
    threshold = 1e-12
    record = {"source": str(path.relative_to(ROOT)), "source_sha256": digest,
              "ramp_duration_days": duration, "windows": windows,
              "initial_potential_l1": norm[0], "max_potential_identity_residual": max(
                  row["potential_identity_residual"] for row in diagnostics),
              "min_assimilation_rate": min(assimilation_rate), "max_assimilation_rate": max(assimilation_rate),
              "rate_sign_threshold_kg_per_day2": threshold,
              "positive_rate_samples": int(sum(assimilation_rate > threshold)),
              "negative_rate_samples": int(sum(assimilation_rate < -threshold)),
              "positive_rate_samples_during_ramp": int(sum((assimilation_rate > threshold) & (day <= duration))),
              "positive_rate_samples_during_hold": int(sum((assimilation_rate > threshold) & (day > duration))),
              "positive_rate_samples_while_vpd_increasing": int(sum(
                  (assimilation_rate > threshold) & (day > 0) & (day < duration))),
              "maximum_assimilation_rate_day": float(day[np.argmax(assimilation_rate)]),
              "potential_refinement_max_l1": float(np.abs(refined_potential-potential).sum(axis=1).max()),
              "new_six_state_integration_count": 0,
              "timing_seconds": {"cached_trajectory_loading": loading_time, "diagnostics_and_bounds": evaluation_time,
                                 "original_integration_time": None},
              "evidence": "C: endpoint-informed bound; no VPD-only prediction and no continuous-time maximum certificate"}
    arrays = {"day": day, "vpd_mpa": vpd, "state": state, "capacity": capacity, "potential": potential,
              "assimilation": assimilation, "assimilation_rate": assimilation_rate,
              "capacity_rate": np.array([row["capacity_rate"] for row in diagnostics]),
              "potential_rate": np.array([row["potential_rate"] for row in diagnostics]),
              "rhs": np.array([row["rhs"] for row in diagnostics]),
              "endpoint_peak_lower_l1": lower, "unknown_initial_peak_lower_l1": general_lower,
              "history_peak_lower_l1": history_lower,
              "sampled_running_peak_l1": observed_peak, "potential_l1": norm}
    return arrays, record


def step_experiment(model):
    start = perf_counter()
    certificate = certify_drying_step_potential(model, .00244, .00254, .02)
    certificate_time = perf_counter()-start
    initial = model.equilibria(.00244)[0]
    day = np.linspace(0, .02, 201)
    states = []
    start = perf_counter()
    for tolerance, max_step in [(2e-11, .001), (2e-12, .0005)]:
        run = solve_ivp(lambda t, y: model.rhs(y, .00254), (0, .02), initial, t_eval=day,
                        method="Radau", rtol=tolerance, atol=tolerance*.01, max_step=max_step)
        if not run.success:
            raise RuntimeError(f"Step validation integration failed: {run.message}")
        states.append(run.y.T)
    validation_time = perf_counter()-start
    diag = [potential_diagnostics(model, y, .00254) for y in states[1]]
    potential = np.array([row["potential"] for row in diag])
    rate = np.array([row["assimilation_rate"] for row in diag])
    lower, upper = [np.array(certificate[key]) for key in [
        "absolute_potential_component_lower", "absolute_potential_component_upper"]]
    box = np.array(certificate["state_box"], dtype=float)
    if (np.any(states[1] < box[:, 0]) or np.any(states[1] > box[:, 1])
            or np.any(lower > -potential[-1]) or np.any(-potential[-1] > upper)):
        raise RuntimeError("Independent step trajectory failed the certified state/potential enclosure")
    cumulative = float(simpson(np.abs(potential).sum(axis=1), x=day))
    if cumulative < certificate["cumulative_potential_l1_lower"]:
        raise RuntimeError("Step cumulative lower bound failed")
    certificate.update({
        "observed_final_absolute_potential": -potential[-1],
        "sampled_peak_l1": np.abs(potential).sum(axis=1).max(),
        "observed_cumulative_potential_l1": cumulative,
        "fraction_of_sampled_peak": certificate["final_potential_l1_lower"]/np.abs(potential).sum(axis=1).max(),
        "numerical_state_refinement_max": float(np.abs(states[0]-states[1]).max()),
        "timing_seconds": {"interval_preprocessing_and_scenario": certificate_time,
                           "two_validation_six_state_integrations": validation_time},
    })
    arrays = {"day": day, "state": states[1], "potential": potential, "assimilation_rate": rate,
              "component_lower_at_end": lower, "component_upper_at_end": upper}
    return arrays, certificate


def figures(output, linear, ramps, step, certificate):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    case = linear["identity_10year"]
    axes[0, 0].plot(case["day"]/YEAR_DAYS, case["assimilation_kg_per_day"]*YEAR_DAYS)
    axes[0, 0].set(xlabel="Time (year)", ylabel="Prescribed assimilation (kg C m$^{-2}$ year$^{-1}$)")
    for key, label in [("identity_10year", "B = I"), ("allocation_10year", "B = b")]:
        arrays = linear[key]
        norm = np.abs(arrays["compensated_potential"]).sum(axis=1)
        axes[0, 1].semilogy(arrays["day"]/YEAR_DAYS, np.where(norm > 0, norm, np.nan), label=label)
        positive = arrays["energy_budget"] > 0
        axes[1, 1].semilogx(
            arrays["energy_budget"][positive]/arrays["energy_budget"][-1],
            arrays["total_output_budget_lower"][positive]/abs(sum(arrays["natural_potential"][-1])), label=label)
    natural_norm = np.abs(case["natural_potential"]).sum(axis=1)
    axes[0, 1].semilogy(case["day"]/YEAR_DAYS, np.where(natural_norm > 0, natural_norm, np.nan),
                        "k--", label="Natural")
    axes[0, 1].set(xlabel="Time (year)", ylabel="Potential L1 norm (kg C m$^{-2}$)")
    names = list(linear)
    actions = [linear[key]["energy_budget"][-1] for key in names]
    axes[1, 0].bar(range(len(names)), actions)
    axes[1, 0].set(yscale="log", ylabel="Minimum action ((kg C m$^{-2}$)$^2$ day$^{-1}$)")
    axes[1, 0].set_xticks(range(len(names)), ["I, 1 yr", "I, 10 yr", "b, 10 yr", "I, weighted, 1 yr"], rotation=15)
    axes[1, 1].set(xlabel="Energy budget / minimum terminal action (-)",
                    ylabel="Unavoidable total output / natural output (-)", ylim=(0, 1.05), xlim=(1e-8, 1))
    axes[0, 1].legend()
    axes[1, 1].legend()
    save(fig, output, "exogenous_capacity_compensation")

    labels = ["10 day", "1 year", "10 year", "100 year"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    bars = []
    for arrays, label in zip(ramps, labels):
        line, = axes[0].plot(arrays["day"]/YEAR_DAYS, arrays["sampled_running_peak_l1"], label=label)
        axes[0].plot(arrays["day"]/YEAR_DAYS, arrays["history_peak_lower_l1"], "--", color=line.get_color())
        duration = [10, YEAR_DAYS, 10*YEAR_DAYS, 100*YEAR_DAYS][len(bars)]
        i = np.flatnonzero(arrays["day"] == duration)[0]
        bars.append([arrays["endpoint_peak_lower_l1"][i], arrays["sampled_running_peak_l1"][i], arrays["potential_l1"][i]])
    for j, label in enumerate(["Endpoint lower", "Sampled peak", "Final norm"]):
        axes[1].bar(np.arange(4)+(j-1)*.25, np.array(bars)[:, j], width=.25, label=label)
    axes[0].set(xlabel="Time (year)", ylabel="Peak potential L1 norm (kg C m$^{-2}$)", xlim=(0, 300))
    axes[0].legend()
    axes[1].set_xticks(range(4), labels)
    axes[1].set(ylabel="Potential over ramp window (kg C m$^{-2}$)")
    axes[1].legend()
    save(fig, output, "ramp_peak_lower_bounds")

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True, constrained_layout=True)
    for arrays, label in zip(ramps, labels):
        axes[0].plot(arrays["day"]/YEAR_DAYS, arrays["assimilation"]*YEAR_DAYS, label=label)
        axes[1].plot(arrays["day"]/YEAR_DAYS, arrays["assimilation_rate"], label=label)
    axes[0].set(ylabel="Assimilation (kg C m$^{-2}$ year$^{-1}$)")
    axes[0].legend()
    axes[1].axhline(0, color="k", lw=.7)
    axes[1].set(yscale="symlog", ylim=(-1e-5, 1e-7), xlabel="Time (year)", xlim=(0, 300),
                ylabel="Assimilation rate (kg C m$^{-2}$ day$^{-2}$)")
    axes[1].set_yscale("symlog", linthresh=1e-11)
    inset = axes[1].inset_axes([.56, .52, .42, .44])
    for arrays in ramps:
        inset.plot(arrays["day"]/YEAR_DAYS, arrays["assimilation_rate"])
    inset.axhline(0, color="k", lw=.5)
    inset.set_yscale("symlog", linthresh=1e-11)
    inset.set(xlim=(0, 2), ylim=(-1e-6, 1e-7), xlabel="Time (year)")
    inset.tick_params(labelsize=7)
    inset.xaxis.label.set_fontsize(8)
    save(fig, output, "ramp_capacity_direction")

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True)
    minutes = step["day"]*1440
    box, rate_box = np.array(certificate["state_box"], dtype=float), np.array(certificate["assimilation_rate_interval"], dtype=float)
    axes[0, 0].axhspan(*box[1], alpha=.2, label="Verified state box")
    axes[0, 0].plot(minutes, step["state"][:, 1], label="ODE check")
    axes[0, 0].set(xlabel="Time after VPD step (min)", ylabel="Leaf water potential (MPa)")
    axes[0, 0].legend()
    axes[0, 1].axhspan(*rate_box, alpha=.2, label="Verified derivative interval")
    axes[0, 1].plot(minutes, step["assimilation_rate"], label="ODE check")
    axes[0, 1].set(xlabel="Time after VPD step (min)", ylabel="Assimilation rate (kg C m$^{-2}$ day$^{-2}$)")
    axes[0, 1].legend()
    for i, label in enumerate(["Leaf", "Wood", "Root", "Soil"]):
        line, = axes[1, 0].plot(minutes, -step["potential"][:, i], label=label)
        axes[1, 0].plot(minutes[-1], step["component_lower_at_end"][i], "o", color=line.get_color())
    axes[1, 0].set(xlabel="Time after VPD step (min)", ylabel="Negative potential magnitude (kg C m$^{-2}$)")
    axes[1, 0].legend()
    for j, (values, label) in enumerate([(step["component_lower_at_end"], "Certified lower"),
                                        (-step["potential"][-1], "ODE check")]):
        axes[1, 1].bar(np.arange(4)+(j-.5)*.35, values, width=.35, label=label)
    axes[1, 1].set_xticks(range(4), ["Leaf", "Wood", "Root", "Soil"])
    axes[1, 1].set(yscale="log", ylabel="Final negative potential magnitude (kg C m$^{-2}$)")
    axes[1, 1].legend()
    save(fig, output, "certified_short_step_lower_bounds")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    start = perf_counter()
    model = RationalSPAC(model_family()["E4"])
    structure = potential_structure(model)
    structure_time = perf_counter()-start
    previous_provenance = json.loads((PREVIOUS/"provenance.json").read_text())
    for name, digest in previous_provenance["source_sha256"].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
            raise RuntimeError(f"Reused experiment source changed: {name}")
    records = {"created_utc": datetime.now(timezone.utc).isoformat(),
               "git_head": subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
               "model": asdict(model), "potential_structure": asdict(structure),
               "model_and_structure_extraction_seconds": structure_time,
               "norm": "physical carbon L1, kg C m^-2; M induced L1, day^-1",
               "linear_benchmarks": [], "reused_ramps": []}
    linear = {}
    for name, duration, inputs, weight in [
            ("identity_1year", YEAR_DAYS, np.eye(4), np.eye(4)),
            ("identity_10year", 10*YEAR_DAYS, np.eye(4), np.eye(4)),
            ("allocation_10year", 10*YEAR_DAYS, structure.allocation[:, None], np.eye(1)),
            ("weighted_1year", YEAR_DAYS, np.eye(4), np.diag([1., 4., 1., .25]))]:
        arrays, row = linear_benchmark(structure, duration, inputs, weight)
        linear[name] = arrays
        records["linear_benchmarks"].append({"name": name, **row})
        np.savez_compressed(args.output/f"{name}.npz", **arrays)
    ramps = []
    for duration in [10., YEAR_DAYS, 10*YEAR_DAYS, 100*YEAR_DAYS]:
        path = PREVIOUS/f"ramp_{duration:g}day.npz"
        arrays, row = analyze_ramp(model, path, previous_provenance)
        ramps.append(arrays)
        records["reused_ramps"].append(row)
        np.savez_compressed(args.output/path.name, **arrays)
    step, certificate = step_experiment(model)
    records["short_step_certificate"] = certificate
    np.savez_compressed(args.output/"short_step.npz", **step)
    figures(args.output, linear, ramps, step, certificate)
    sources = [Path(__file__), ROOT/"docs/leastaction.md", ROOT/"tests/test_spac_carbon_potential.py",
               ROOT/"examples/analyze_e4_capacity_iss.py", ROOT/"examples/prepare_manuscript_phase_portraits.py",
               *[ROOT/f"src/control_carbon/{name}.py" for name in [
                   "carbon_compensation", "spac_carbon_potential", "spac_carbon_iss", "spac_carbon_rational",
                   "spac_carbon_vpd_iss", "spac_carbon_stability", "spac_carbon", "spac_carbon_mpa",
                   "spac_carbon_two_way", "spac_carbon_equilibria", "manzoni_hydraulics"]],
               ROOT/"configs/spac_carbon_stage_b_synthetic.json", ROOT/"configs/spac_carbon_pressure_review.json",
               PREVIOUS/"provenance.json"]
    records["source_sha256"] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    (args.output/"results.json").write_text(json.dumps(serializable(records), indent=2, allow_nan=False)+"\n")
    provenance = {"source_sha256": records["source_sha256"],
                  "output_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                    for p in sorted(args.output.iterdir()) if p.is_file()}}
    (args.output/"provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")
    print(json.dumps(serializable({
        "output": str(args.output), "action_checks": records["linear_benchmarks"],
        "ramp_bounds": records["reused_ramps"], "step_certificate": certificate}), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
