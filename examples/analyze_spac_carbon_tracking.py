"""Stage C2/C3 bounded rate tracking and local saddle-cycle sections."""

import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.linalg import expm

from analyze_spac_carbon_capacity import encode, record, require, write_csv
from control_carbon.capacity_feedback import (
    CubicCapacityFeedback, periodic_attractor, periodic_capacity, simulate_capacity_ramp,
)
from control_carbon.manzoni_hydraulics import HydraulicParameters, transpiring_equilibria
from control_carbon.spac_carbon import YEAR_DAYS, CarbonParameters, hydraulic_factor, positive_carbon_equilibrium
from control_carbon.spac_carbon_two_way import StructuralFeedback
from control_carbon.spac_carbon_mpa import PressureSPAC
from control_carbon.spac_carbon_capacity import diagnose_capacity, reference_scales
from control_carbon.spac_carbon_tracking import (
    FrozenBranch, SmoothVPDRamp, frozen_lyapunov_metric, linearized_tracking, tracking_diagnostics,
)
from control_carbon.spac_carbon_e4 import hopf_point, periodic_orbit
from control_carbon.spac_carbon_section import local_cycle_section, section_return

ROOT = Path(__file__).resolve().parents[1]


def sampled_outcome(solution, rows, hold_days):
    if solution["status"] == "vegetation_domain_exit":
        return "vegetation_domain_exit_not_alternate_attractor"
    tail = [row for row in rows if row["day"] >= solution["time_days"][-1]-.1*hold_days]
    require(len(tail) >= 30, "Equilibrium-neighborhood classification needs a resolved final hold window")
    if (max(row["equilibrium_error_norm_scaled"] for row in tail) < 1e-4
            and max(row["full_residual_norm_scaled_per_day"] for row in tail) < 1e-7):
        return "sampled_equilibrium_neighborhood_not_global_basin_proof"
    return "finite_time_unresolved"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT/"configs/spac_carbon_pressure_review.json")
    parser.add_argument("--output", type=Path, default=ROOT/"artifacts/spac_carbon"/datetime.now().strftime("stage_c23_tracking_%Y%m%d_%H%M%S"))
    args = parser.parse_args()
    settings = json.loads(args.config.read_text())
    require(settings["units"]["vpd"] == "MPa", "Tracking runner requires MPa settings")
    require(settings["routine_vpd_mpa_bounds"] == [0., .003], "Campaign main range must be 0--.003 MPa")
    config_path = args.config.parent/settings["model_config"]
    cfg = json.loads(config_path.read_text())
    h = HydraulicParameters(**cfg["hydraulic"])
    k = cfg["carbon"]
    c = CarbonParameters(k["assimilation_max_per_year"]/YEAR_DAYS, tuple(k["allocation"]),
                         tuple(v/YEAR_DAYS for v in k["losses_per_year"]), tuple(k["litter_fractions"]),
                         k["specific_leaf_area"], k["carbon_dry_fraction"], k["extinction"])
    reference_psi = transpiring_equilibria(cfg["root_potential"], cfg["reference_vpd"], h)[0]
    stocks = positive_carbon_equilibrium(hydraulic_factor(reference_psi, h), c)
    require(stocks is not None, "Canonical structural references require a positive root")
    model = PressureSPAC(h, c, StructuralFeedback(stocks[0], stocks[1], cfg["leaf_floor"], cfg["minimum_leaf"],
                                                1., 1., 1., cfg["wood_exponent"]),
                         cfg["root_potential"], settings["atmospheric_pressure_mpa"])
    scales = reference_scales(model)
    scale = np.array(scales.state)
    branch = FrozenBranch(model, (-2.05, -.3))
    start = settings["reference_vpd_mpa"]
    endpoints = [.0024, .00266]
    durations = [10., 10*YEAR_DAYS, 100*YEAR_DAYS]
    hold = 300*YEAR_DAYS
    initial = branch.state(start)
    args.output.mkdir(parents=True, exist_ok=False)
    all_rows, summaries, metrics, derivative_checks = [], [], [], []

    for v in np.linspace(start, max(endpoints), 41):
        y, tangent, _ = branch.tangent(v, scales)
        step = 1e-8
        fd = (branch.state(v+step)-branch.state(v-step))/(2*step)
        error = float(np.max(abs((fd-tangent)/scale)))
        relative = error/max(1., np.max(abs(tangent/scale)))
        require(relative < 1e-6, "Full equilibrium tangent finite-difference refinement failed")
        d = tracking_diagnostics(y, v, 0., branch, scales)
        require(d["spectral_abscissa_per_day"] < -1e-10, "Rate campaign has an unstable sampled target")
        metric = frozen_lyapunov_metric(v, branch, scales)
        refined_metric = frozen_lyapunov_metric(v, branch, scales, step_mpa=5e-8)
        q_error = float(np.linalg.norm(metric["metric_derivative_per_mpa"]
                                     -refined_metric["metric_derivative_per_mpa"], 2)
                        /np.linalg.norm(refined_metric["metric_derivative_per_mpa"], 2))
        require(q_error < 1e-4, "Lyapunov metric derivative refinement failed")
        gains = [float(np.linalg.norm(expm(d["scaled_jacobian"]*t), 2))
                 for t in np.r_[0., np.geomspace(.001, 1e5, 35)]]
        metrics.append({
            "vpd_mpa": v, "spectral_abscissa_per_day": d["spectral_abscissa_per_day"],
            "numerical_abscissa_per_day": d["numerical_abscissa_per_day"],
            "slowest_decay_years": d["slowest_decay_days"]/YEAR_DAYS,
            "eigenvector_condition_2": d["eigenvector_condition_2"],
            "sampled_frozen_transient_gain": max(gains),
            "minimum_metric_eigenvalue_days": metric["minimum_metric_eigenvalue_days"],
            "maximum_metric_eigenvalue_days": metric["maximum_metric_eigenvalue_days"],
            "lyapunov_residual": metric["lyapunov_residual"],
            "metric_derivative_relative_refinement_error": q_error,
        })
        derivative_checks.append({"vpd_mpa": v, "tangent_relative_fd_error": relative})

    for endpoint in endpoints:
        for duration in durations:
            ramp = SmoothVPDRamp(start, endpoint, duration)
            label = f"E4_v{endpoint:g}_ramp{duration/YEAR_DAYS:g}year"
            times = np.unique(np.r_[np.linspace(0., duration, 101),
                                    duration+np.geomspace(.01, hold, 151),
                                    duration+np.linspace(.9*hold, hold, 31)])
            forcing = lambda t: ramp.value_rate(t)[0]
            solution = model.simulate(initial, times, forcing, rtol=1e-9)
            refined = model.simulate(initial, times, forcing, rtol=1e-10)
            require(solution["status"] == refined["status"], f"{label}: status changed with tolerance")
            common = times[times <= min(solution["time_days"][-1], refined["time_days"][-1])]
            interpolated = [
                np.column_stack([np.interp(common, run["time_days"], run["state"][:, j]) for j in range(6)])
                for run in [solution, refined]
            ]
            state_error = float(np.max(abs(interpolated[0]-interpolated[1])/scale))
            event_error = abs(solution["time_days"][-1]-refined["time_days"][-1])
            require(state_error < 2e-6 and event_error < .1, f"{label}: trajectory/event refinement failed")
            linear_times = refined["time_days"]
            linear = linearized_tracking(linear_times, ramp, branch, scales, rtol=1e-9)
            linear_refined = linearized_tracking(linear_times, ramp, branch, scales, rtol=1e-10)
            linear_error = float(np.max(abs(linear["error_scaled"]-linear_refined["error_scaled"])))
            require(linear_error < 1e-7, f"{label}: linearized tracking refinement failed")
            rows = []
            for i, (t, y) in enumerate(zip(refined["time_days"], refined["state"])):
                v, rate = ramp.value_rate(t)
                td = tracking_diagnostics(y, v, rate, branch, scales)
                data = diagnose_capacity(y, v, model, scales, equilibrium=td["equilibrium"])
                require(data["identity_residual_max"] < 1e-14
                        and data["potential_identity_residual_max"] < 1e-11
                        and abs(data["water_budget_residual"]) < 1e-14
                        and abs(data["carbon_budget_residual"]) < 1e-14,
                        f"{label}: capacity identity or conservation budget failed")
                row = record(label, t, y, v, data, td["equilibrium"])
                row.update({
                    "vpd_rate_mpa_per_day": rate,
                    "branch_speed_scaled_per_day": td["branch_speed_norm_scaled_per_day"],
                    "spectral_abscissa_per_day": td["spectral_abscissa_per_day"],
                    "numerical_abscissa_per_day": td["numerical_abscissa_per_day"],
                    "slowest_decay_years": td["slowest_decay_days"]/YEAR_DAYS,
                    "nonlinear_remainder_scaled_per_day": td["nonlinear_remainder_norm_scaled_per_day"],
                    "linear_error_norm_scaled": float(np.linalg.norm(linear_refined["error_scaled"][i])),
                    "nonlinear_vs_linear_error_norm_scaled": float(np.linalg.norm(
                        td["error_scaled"]-linear_refined["error_scaled"][i])),
                })
                for j in range(6):
                    row[f"linear_error_scaled_{j}"] = linear_refined["error_scaled"][i, j]
                rows.append(row)
            # Q(v) drift is evaluated along this forcing schedule at stated sample times.
            drift_margins = []
            for t in np.linspace(0., duration, 41):
                v, rate = ramp.value_rate(t)
                q = frozen_lyapunov_metric(v, branch, scales, step_mpa=5e-8)
                drift_margins.append(float(min(np.linalg.eigvalsh(
                    np.eye(6)-q["metric_derivative_per_mpa"]*rate))))
            summary = {
                "case": label, "start_mpa": start, "end_mpa": endpoint,
                "duration_days": duration, "hold_days": hold,
                "status": refined["status"], "outcome": sampled_outcome(refined, rows, hold),
                "end_day": float(refined["time_days"][-1]),
                "max_error_scaled": max(row["equilibrium_error_norm_scaled"] for row in rows),
                "final_error_scaled": rows[-1]["equilibrium_error_norm_scaled"],
                "max_P_norm_scaled": max(row["potential_norm_scaled"] for row in rows),
                "max_linearization_error_scaled": max(row["nonlinear_vs_linear_error_norm_scaled"] for row in rows),
                "min_sampled_metric_drift_margin": min(drift_margins),
                "state_refinement_error_scaled": state_error, "event_refinement_error_days": event_error,
                "linear_refinement_error_scaled": linear_error,
                "final_window_samples": sum(row["day"] >= refined["time_days"][-1]-.1*hold for row in rows),
            }
            all_rows.extend(rows)
            summaries.append(summary)
            print(f"{label}: {summary['outcome']}; max error={summary['max_error_scaled']:.6g}, "
                  f"final={summary['final_error_scaled']:.6g}", flush=True)

    small_rows, small_summaries = [], []
    for delta in [1e-6, 5e-7]:
        ramp = SmoothVPDRamp(.0024, .0024+delta, 10.)
        small_times = np.unique(np.r_[np.linspace(0., 10., 41), 10.+np.geomspace(.01, 100*YEAR_DAYS, 121)])
        small_initial = branch.state(.0024)
        nonlinear = model.simulate(small_initial, small_times, lambda t: ramp.value_rate(t)[0], rtol=1e-11)
        linear = linearized_tracking(small_times, ramp, branch, scales, rtol=1e-11)
        require(nonlinear["status"] == "completed", "Small-amplitude linearity verification exited domain")
        differences = []
        for i, (t, y) in enumerate(zip(nonlinear["time_days"], nonlinear["state"])):
            v, rate = ramp.value_rate(t)
            td = tracking_diagnostics(y, v, rate, branch, scales)
            error = float(np.linalg.norm(td["error_scaled"]-linear["error_scaled"][i]))
            differences.append(error)
            small_rows.append({
                "delta_vpd_mpa": delta, "day": t, "vpd_mpa": v,
                "nonlinear_error_scaled": td["error_norm_scaled"],
                "linear_error_scaled": float(np.linalg.norm(linear["error_scaled"][i])),
                "nonlinear_vs_linear_error_scaled": error,
            })
        small_summaries.append({"delta_vpd_mpa": delta, "max_nonlinear_vs_linear_error_scaled": max(differences)})
    quadratic_ratio = small_summaries[0]["max_nonlinear_vs_linear_error_scaled"]/small_summaries[1]["max_nonlinear_vs_linear_error_scaled"]
    require(3.8 < quadratic_ratio < 4.2, "Small-ramp full trajectory linearity does not refine quadratically")

    hopf = hopf_point((-2.01, -1.99), model)
    cycle = periodic_orbit(hopf, .002, model, tol=1e-9)
    section = local_cycle_section(cycle, model, scales, equilibrium=branch.state(cycle["vpd_mpa"]))
    section_refined = local_cycle_section(cycle, model, scales, rtol=1e-12, steps_per_period=240,
                                          equilibrium=branch.state(cycle["vpd_mpa"]))
    multiplier_error = float(max(abs(section["multipliers"]-section_refined["multipliers"])))
    require(multiplier_error < 1e-7, "Section transverse multipliers failed tolerance/step refinement")
    section = section_refined
    vectors, covectors = section["right_vectors"], section["left_covectors"]
    returns, fd_errors = [], []
    for delta in [5e-7, 2e-7]:
        minus = section_return(section, -delta*vectors[:, 0], model, rtol=1e-12, steps_per_period=240)
        plus = section_return(section, delta*vectors[:, 0], model, rtol=1e-12, steps_per_period=240)
        fd = (plus["coordinates"]-minus["coordinates"])/(2*delta)
        fd_errors.append(float(max(abs(fd-section["return_derivative"]@vectors[:, 0]))))
    require(fd_errors[1] < 1e-5 and fd_errors[1] < .4*fd_errors[0],
            "Nonlinear section finite-difference refinement failed")
    zero_return = section_return(section, np.zeros(3), model, rtol=1e-12, steps_per_period=240)
    for stable in [-2e-8, 0., 2e-8]:
        for unstable in [-2e-7, 0., 2e-7]:
            start_coordinates = unstable*vectors[:, 0]+stable*vectors[:, 1]
            returned = section_return(section, start_coordinates, model, rtol=1e-12, steps_per_period=240)
            baseline = section_return(section, start_coordinates, model, rtol=1e-11, steps_per_period=160)
            difference = float(max(abs(returned["coordinates"]-baseline["coordinates"])))
            require(difference < 1e-8, "Nonlinear section return tolerance refinement failed")
            centered = returned["coordinates"]-zero_return["coordinates"]
            modes = covectors@centered
            if unstable:
                require(modes[0]*unstable > 0, "Signed local sides not resolved at the chosen perturbation")
            prediction = section["return_derivative"]@start_coordinates
            returns.append({
                "stable_initial_scaled": stable, "unstable_initial_scaled": unstable,
                "stable_return_scaled": float(modes[1]), "unstable_return_scaled": float(modes[0]),
                "return_days": returned["time_days"], "section_residual": returned["section_residual"],
                "nonlinear_vs_tangent_prediction_norm_scaled": float(np.linalg.norm(centered-prediction)),
                "return_refinement_error_scaled": difference,
                "interpretation": "one_period_local_side_not_recovery_exit_label",
            })
    print(f"E4 section: transverse multipliers={section['multipliers']}; "
          f"radial e-folding={section['radial_e_folding_days']/YEAR_DAYS:.6g} years", flush=True)

    scalar = CubicCapacityFeedback()
    scalar_runs = [simulate_capacity_ramp(duration, scalar) for duration in [20., 120.]]
    require([run["outcome"] for run in scalar_runs] == ["low", "high"],
            "Existing positive-carbon rate-transition benchmark failed")
    benchmark_rows = []
    for duration, run in zip([20., 120.], scalar_runs):
        for t, stock, driver, capacity in zip(run["time"], run["stock"], run["driver"], run["capacity"]):
            low, saddle, high = scalar.equilibria(driver)
            benchmark_rows.append({"duration_years": duration, "year": t, "stock": stock,
                                   "K": capacity, "low_equilibrium": low,
                                   "boundary_equilibrium": saddle, "high_equilibrium": high})
    periodic_times = np.linspace(0., 3., 401)
    linear_benchmark = [{
        "year": t, "instantaneous_capacity": float(periodic_capacity(t)),
        "exact_attracting_non_autonomous_trajectory": float(periodic_attractor(t)),
    } for t in periodic_times]
    legacy_path = {
        "start_vpd_mpa": .00101325, "end_vpd_mpa": .00303975,
        "main_range_exceeded": True, "description": "Archived E2/E3 10-day exit versus 10-year recovery, not this E4 path.",
        "source": "docs/spac_carbon_stage_b_results.md", "rerun_in_this_campaign": False,
    }
    files = [
        ("tracking.csv", all_rows), ("rates.csv", summaries), ("branch_metrics.csv", metrics),
        ("section_returns.csv", returns), ("scalar_benchmark.csv", benchmark_rows),
        ("linear_benchmark.csv", linear_benchmark), ("small_ramps.csv", small_rows),
    ]
    for name, rows in files:
        write_csv(args.output/name, rows)
    fig, axes = plt.subplots(2, 3, figsize=(16, 8))
    for summary in summaries:
        rows = [row for row in all_rows if row["case"] == summary["case"]]
        t = np.array([row["day"] for row in rows])/YEAR_DAYS
        ax = axes[0, endpoints.index(summary["end_mpa"])]
        errors = np.array([row["equilibrium_error_norm_scaled"] for row in rows])
        selected = t > 0
        duration_label = "10 day" if summary["duration_days"] == 10 else f"{summary['duration_days']/YEAR_DAYS:g} yr"
        ax.semilogy(t[selected], errors[selected], label=f"ramp {duration_label}")
        if summary["status"] == "vegetation_domain_exit":
            ax.plot(t[-1], errors[-1], "rx", ms=9, label="L_min domain exit")
        ax.set(xlabel="Time (year)", ylabel="Fixed-scale full target error",
               title=f"E4: same path to {summary['end_mpa']:g} MPa", ylim=(1e-5, 10.))
    slow_rows = [row for row in all_rows if row["case"] == summaries[2]["case"]]
    for key, label in [("L", "actual L"), ("K_L", "K_L"), ("equilibrium_L", "moving equilibrium L*")]:
        axes[0, 2].plot(np.array([row["day"] for row in slow_rows])/YEAR_DAYS,
                        [row[key] for row in slow_rows], label=label)
    axes[0, 2].set(xlabel="Time (year)", ylabel="kg C/m2", title="100-year ramp: map != target")
    axes[1, 0].plot([row["vpd_mpa"] for row in metrics], [row["sampled_frozen_transient_gain"] for row in metrics],
                    label="sampled ||exp(A t)||_2")
    axes[1, 0].set(xlabel="VPD (MPa)", ylabel="Fixed-scale gain", title="Non-normal amplification, frozen diagnostic")
    for stable in [-2e-8, 0., 2e-8]:
        rows = [row for row in returns if row["stable_initial_scaled"] == stable]
        axes[1, 1].plot([row["unstable_initial_scaled"] for row in rows],
                        [row["unstable_return_scaled"]-row["unstable_initial_scaled"] for row in rows],
                        "o-", label=f"stable shift {stable:g}")
    axes[1, 1].axhline(0., color="grey", lw=.8)
    axes[1, 1].set(xlabel="Initial unstable coordinate (scaled)", ylabel="Unstable increment per return (scaled)",
                   title="Signed local growth: NOT global basin labels")
    for duration, run in zip([20., 120.], scalar_runs):
        axes[1, 2].plot(run["time"], run["stock"], label=f"scalar ramp {duration:g} yr -> {run['outcome']}")
    axes[1, 2].set(xlabel="Time (year)", ylabel="Arbitrary carbon units", title="True branch transition: mathematical benchmark")
    for ax in axes.ravel():
        ax.grid(alpha=.25)
        ax.legend(fontsize=8)
    fig.suptitle("Stage C2/C3: moving capacity targets and LOCAL saddle-cycle geometry | synthetic, no mortality")
    fig.tight_layout()
    for extension in ["png", "pdf"]:
        fig.savefig(args.output/f"tracking_sections.{extension}", dpi=180)
    plt.close(fig)

    sources = [
        args.config, config_path, ROOT/"examples/analyze_spac_carbon_tracking.py",
        ROOT/"examples/analyze_spac_carbon_capacity.py", ROOT/"tests/test_spac_carbon_tracking.py",
        *[ROOT/f"src/control_carbon/{name}.py" for name in [
            "spac_carbon_tracking", "spac_carbon_section", "spac_carbon_capacity",
            "spac_carbon_e4", "spac_carbon_bifurcation", "spac_carbon_two_way",
            "spac_carbon_mpa", "spac_carbon", "manzoni_hydraulics", "capacity_feedback",
        ]],
        ROOT/"docs/SPAC_nonlinear_carbon_storage_capacity_research_plan.md",
    ]
    git_head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
                              capture_output=True, text=True).stdout.strip()
    result = {
        "classification": "analytic_tracking_equations_and_precision_checked_finite_experiments_not_global_basin_or_R_tipping",
        "baseline_git_SHA": "757f359", "run_git_HEAD": git_head,
        "settings": settings, "baseline_config": cfg, "state_scales": scales,
        "frozen_branch_leaf_bounds_mpa": branch.leaf_bounds_mpa,
        "ramp_shape": "C1_cubic_smoothstep_then_fixed_endpoint", "rates": summaries,
        "small_ramp_linearity_checks": small_summaries, "small_ramp_quadratic_ratio": quadratic_ratio,
        "derivative_checks": derivative_checks, "branch_metrics": metrics,
        "metric_scope": "sampled_frozen_Q_and_Qv_drift_not_uniform_interval_certificate_or_nonlinear_bound",
        "cycle": cycle, "section": section, "section_multiplier_refinement_error": multiplier_error,
        "section_FD_errors": fd_errors, "section_FD_steps_scaled": [5e-7, 2e-7],
        "section_returns": returns,
        "scalar_benchmark_outcomes": [run["outcome"] for run in scalar_runs],
        "linear_periodic_benchmark": "exact_non_autonomous_attractor_not_an_autonomous_stable_cycle",
        "legacy_rate_path": legacy_path,
        "solver_settings": {
            "nonlinear": {"method": "Radau_log_carbon", "rtol": [1e-9, 1e-10], "atol_over_rtol": .01},
            "linearized": {"method": "Radau", "rtol": [1e-9, 1e-10], "atol_over_rtol": .01},
            "cycle": {"collocation_tol": 1e-9, "amplitude": .002},
            "section": {"method": "Radau_variational_plus_positive_oriented_event",
                        "rtol": [1e-10, 1e-12], "steps_per_period": [160, 240]},
            "metric": {"equation": "A.T Q+Q A=-I", "A_central_step_mpa": [1e-7, 5e-8]},
        },
        "acceptance": {
            "nonlinear_refinement_scaled": 2e-6, "exit_time_refinement_days": .1,
            "linear_refinement_scaled": 1e-7, "section_multiplier_refinement": 1e-7,
            "section_FD_scaled": 1e-5, "section_return_refinement_scaled": 1e-8,
            "sampled_target_spectral_abscissa_per_day": -1e-10,
            "small_ramp_quadratic_ratio_bounds": [3.8, 4.2],
            "sampled_equilibrium_neighborhood_error_scaled": 1e-4,
            "sampled_equilibrium_neighborhood_residual_per_day": 1e-7,
            "sampled_neighborhood_window": "last_10_percent_of_300_year_hold",
            "minimum_final_window_samples_for_neighborhood_label": 30,
            "capacity_identity_residual": 1e-14, "potential_identity_residual": 1e-11,
            "water_carbon_budget_residual": 1e-14,
        },
        "versions": {"numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "limits": [
            "Six predefined E4 ramps, not a critical-rate or global basin search.",
            "Local stable tangent plane is not the nonlinear stable manifold.",
            "One-return signed sides are not long-term recovery/exit classifications.",
            "No interval-certified uniform tracking inequality or new physiological regime.",
            "C4 new structural laws and C5 manuscript are deferred.",
        ],
        "source_sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
    }
    (args.output/"results.json").write_text(json.dumps(result, default=encode, indent=2, allow_nan=False)+"\n")
    print(f"Saved {len(all_rows)} tracking samples and {len(returns)} local returns to {args.output}", flush=True)


if __name__ == "__main__":
    main()
