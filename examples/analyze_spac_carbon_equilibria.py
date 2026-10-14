"""Scalar-equilibrium classification and minimal fixed-VPD hydraulic outcomes."""

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
from scipy.optimize import brentq

from analyze_spac_carbon_capacity import encode, record, require, write_csv
from control_carbon.manzoni_hydraulics import HydraulicParameters, conductances, transpiring_equilibria
from control_carbon.spac_carbon import (
    YEAR_DAYS, CarbonParameters, carbon_matrix, hydraulic_factor, positive_carbon_equilibrium,
    assimilation, frozen_coefficient_capacity, carbon_budget,
)
from control_carbon.spac_carbon_two_way import StructuralFeedback, feedback_factors
from control_carbon.spac_carbon_mpa import PressureSPAC
from control_carbon.spac_carbon_capacity import diagnose_capacity, reference_scales
from control_carbon.spac_carbon_equilibria import (
    stock_ratios, positive_domain, leaf_curve, baseline_landmarks, positive_equilibria,
)
from control_carbon.spac_carbon_e4 import hopf_point

ROOT = Path(__file__).resolve().parents[1]
BASELINE = "629314706aaac5f2f78ff37f6e19f01a01d58b42"


def flux_capacity_record(label, day, state, vpd, model, scales, equilibrium):
    """RHS/flux diagnostics remain defined at kinks; full Jacobians do not."""
    y = np.asarray(state)
    scale = np.asarray(scales.state)
    rhs = model.rhs(y, vpd)
    h = hydraulic_factor(y[1], model.hydraulic)
    k = frozen_coefficient_capacity(y[2:], h, model.carbon)
    erx, _, e = model.fluxes(y, vpd)
    gx, gs = conductances(y[1], model.hydraulic)
    return {
        "case": label, "day": day, "vpd_mpa": vpd,
        **dict(zip(["xylem_mpa", "leaf_mpa", "L", "W", "R", "S"], y)),
        **{f"K_{name}": value for name, value in zip(["L", "W", "R", "S"], k)},
        "gx": gx, "gs": gs, "GPP": assimilation(y[2], h, model.carbon),
        "stock_total": float(sum(y[2:])), "capacity_total": float(sum(k)),
        "potential_norm_scaled": float(np.linalg.norm((k-y[2:])/scale[2:])),
        "equilibrium_error_norm_scaled": float(np.linalg.norm((y-equilibrium)/scale)),
        "full_residual_norm_scaled_per_day": float(np.linalg.norm(rhs/scale)),
        "identity_residual_max": float(max(abs(carbon_matrix(model.carbon)@(y[2:]-k)-rhs[2:]))),
        "water_budget_residual": float(model.hydraulic.cx*rhs[0]+model.hydraulic.cl*rhs[1]-erx+e),
        "carbon_budget_residual": carbon_budget(y[2:], h, model.carbon)["residual"],
        "jacobian_status": ("undefined_at_conductance_kink" if y[1] in
                            [2*model.hydraulic.psi50_x, 2*model.hydraulic.psi50_s]
                            else "not_evaluated_in_flux_only_timeseries"),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT/"configs/spac_carbon_pressure_review.json")
    parser.add_argument("--output", type=Path, default=ROOT/"artifacts/spac_carbon"/datetime.now().strftime("scalar_equilibria_%Y%m%d_%H%M%S"))
    args = parser.parse_args()
    settings = json.loads(args.config.read_text())
    require(settings["units"]["vpd"] == "MPa", "Equilibrium runner requires MPa boundary")
    config_path = args.config.parent/settings["model_config"]
    cfg = json.loads(config_path.read_text())
    h = HydraulicParameters(**cfg["hydraulic"])
    p = cfg["carbon"]
    c = CarbonParameters(p["assimilation_max_per_year"]/YEAR_DAYS, tuple(p["allocation"]),
                         tuple(v/YEAR_DAYS for v in p["losses_per_year"]), tuple(p["litter_fractions"]),
                         p["specific_leaf_area"], p["carbon_dry_fraction"], p["extinction"])
    psi = transpiring_equilibria(cfg["root_potential"], cfg["reference_vpd"], h)[0]
    refs = positive_carbon_equilibrium(hydraulic_factor(psi, h), c)
    require(refs is not None, "Canonical reference must have positive stocks")
    models = {
        label: PressureSPAC(
            h, c, StructuralFeedback(refs[0], refs[1], cfg["leaf_floor"], cfg["minimum_leaf"],
                                     *cfg["variants"][label], wood_exponent=cfg["wood_exponent"]),
            cfg["root_potential"], settings["atmospheric_pressure_mpa"],
        ) for label in ["E0", "E1", "E4"]
    }
    scales = reference_scales(models["E0"])
    landmarks = {label: baseline_landmarks(m) for label, m in models.items()}
    e4 = models["E4"]
    def abscissa(leaf):
        y, v, _ = leaf_curve(leaf, e4)
        return float(max(np.linalg.eigvals(e4.jacobian(y, v)).real))
    scan = np.unique(np.r_[np.geomspace(1e-9, positive_domain(e4)["leaf_upper_closed_at_zero_vpd"], 321),
                           np.linspace(.001, positive_domain(e4)["leaf_upper_closed_at_zero_vpd"], 241)])
    values = [abscissa(leaf) for leaf in scan]
    crossings = []
    for left, right, fl, fr in zip(scan[:-1], scan[1:], values[:-1], values[1:]):
        if fl*fr < 0:
            leaf = brentq(abscissa, left, right, xtol=1e-14)
            y, v, _ = leaf_curve(leaf, e4)
            eig = np.linalg.eigvals(e4.jacobian(y, v))
            require(max(abs(eig.imag)) > 1e-8, "E4 sampled stability crossing lost its complex pair")
            steps = [leaf*1e-4, leaf*5e-5]
            slopes = [(abscissa(leaf+s)-abscissa(leaf-s))/(2*s*leaf_curve(leaf, e4)[2]) for s in steps]
            require(abs(slopes[0]-slopes[1]) < 1e-3*max(abs(slopes[-1]), 1.),
                    "E4 stability crossing derivative refinement failed")
            fd_checks = []
            jac = e4.jacobian(y, v)
            for relative_step in [1e-4, 5e-5]:
                columns = []
                for j in range(6):
                    delta = relative_step*max(abs(y[j]), 1e-8)
                    plus, minus = y.copy(), y.copy()
                    plus[j] += delta
                    minus[j] -= delta
                    columns.append((e4.rhs(plus, v)-e4.rhs(minus, v))/(2*delta))
                fd = np.column_stack(columns)
                eigen_error = float(max(np.min(abs(eig[:, None]-np.linalg.eigvals(fd)), axis=1)))
                require(eigen_error < 1e-8, "E4 crossed-pair finite-difference verification failed")
                fd_checks.append({
                    "relative_step": relative_step, "max_matched_eigenvalue_error_per_day": eigen_error,
                    "relative_jacobian_infinity_norm_error": float(np.linalg.norm(fd-jac, np.inf)/np.linalg.norm(jac, np.inf)),
                })
            crossings.append({
                "leaf": leaf, "vpd_mpa": v, "eigenvalues": eig, "slope_per_mpa_refinement": slopes,
                "finite_difference_jacobian_checks": fd_checks,
                "scope": "numerical_linear_stability_crossing_not_new_l1_or_periodic_analysis",
                "in_main_vpd_range": bool(v <= .003),
            })
    args.output.mkdir(parents=True, exist_ok=False)
    curves, equilibria, counts, trajectories, outcomes, exact_closed = [], [], [], [], [], []
    vgrid = np.unique(np.r_[np.linspace(0., .0051, 121), .00244, .00266, .003,
                           [point["zero_leaf_vpd_limit_mpa"] for point in landmarks.values()],
                           [point["fold_vpd_mpa"] for point in landmarks.values() if point["fold_leaf"] is not None],
                           [point["vpd_mpa"] for point in crossings]])
    for label, model in models.items():
        top = positive_domain(model)["leaf_upper_closed_at_zero_vpd"]
        leaf_grid = np.unique(np.r_[np.geomspace(1e-9, top, 161), np.linspace(top/100, top, 141),
                                   [] if landmarks[label]["fold_leaf"] is None else [landmarks[label]["fold_leaf"]],
                                   [point["leaf"] for point in crossings] if label == "E4" else []])
        for leaf in leaf_grid:
            y, v, slope = leaf_curve(leaf, model)
            data = diagnose_capacity(y, v, model, scales, equilibrium=y)
            row = record(label, 0., y, v, data, y)
            row.update({"vpd_slope_per_leaf": slope,
                        "in_main_vpd_range": bool(v <= .003),
                        "D_mole_fraction": model.mole_fraction(v)})
            curves.append(row)
        for v in vgrid:
            result = positive_equilibria(v, model, scales)
            roots = result["roots"]
            counts.append({
                "model": label, "vpd_mpa": v, "positive_conducting_roots": len(roots),
                "stable_roots_numerical_spectrum": sum(root["diagnostics"]["equilibrium_stability"] == "stable" for root in roots),
                "nonhyperbolic_roots_numerical_spectrum": sum(root["diagnostics"]["equilibrium_center_dimensions"] > 0 for root in roots),
                "zero_leaf_limit_hit_excluded": result["zero_leaf_limit_hit"],
                "in_main_vpd_range": bool(v <= .003),
                "formal_zero_nonconducting_sets": "enumerated_separately_not_counted_as_isolated_roots",
            })
            for i, root in enumerate(roots):
                row = record(f"{label}:root{i}", 0., root["state"], v, root["diagnostics"], root["state"])
                row["vpd_slope_per_leaf"] = root["vpd_slope_per_leaf"]
                equilibria.append(row)
    print("E0/E1 single folds and interval-checked E4 monotonicity: scalar roots reconstructed.", flush=True)
    experiment_vpd = .00244
    require(len(positive_equilibria(experiment_vpd, models["E0"], scales)["roots"]) == 2,
            "Representative VPD must have two positive E0 equilibria")
    times = np.unique(np.r_[0., np.geomspace(.0001, 30., 121), np.linspace(31., 100*YEAR_DAYS, 161)])
    for label in ["E0", "E4"]:
        model = models[label]
        roots = positive_equilibria(experiment_vpd, model, scales)["roots"]
        sink = [root["state"] for root in roots if root["diagnostics"]["equilibrium_stability"] == "stable"]
        require(len(sink) == 1, "Representative initial carbon requires one stable positive root")
        for name, psi, stocks in [
            ("wet_initial", -1.5, sink[0][2:]),
            ("dry_initial", -2.5, sink[0][2:]),
            ("both_closed_initial", -4.2, sink[0][2:]),
            ("xylem_closed_low_leaf", -3.5, stock_ratios(model)*(2*model.feedback.minimum_leaf)),
        ]:
            fl, _, fw, *_ = feedback_factors(stocks[0], stocks[1], model.feedback)
            x = (fw*model.root_potential_mpa+fl*psi)/(fw+fl)
            initial = np.r_[x, psi, stocks]
            run = model.simulate(initial, times, lambda t: experiment_vpd, rtol=1e-9)
            refined = model.simulate(initial, times, lambda t: experiment_vpd, rtol=1e-10)
            require(run["status"] == refined["status"], "Representative outcome changed on tolerance refinement")
            common = times[times <= min(run["time_days"][-1], refined["time_days"][-1])]
            states = [np.column_stack([np.interp(common, result["time_days"], result["state"][:, j])
                                       for j in range(6)]) for result in [run, refined]]
            error = float(np.max(abs(states[0]-states[1])/np.array(scales.state)))
            time_error = float(abs(run["time_days"][-1]-refined["time_days"][-1]))
            require(error < 1e-6 and time_error < .1, "Representative IVP failed state/event refinement")
            for t, y in zip(refined["time_days"], refined["state"]):
                row = flux_capacity_record(f"{label}:{name}", t, y, experiment_vpd, model, scales, sink[0])
                require(row["identity_residual_max"] < 1e-14 and abs(row["water_budget_residual"]) < 1e-14
                        and abs(row["carbon_budget_residual"]) < 1e-14, "Representative capacity/budget identity failed")
                trajectories.append(row)
            last = refined["state"][-1]
            gx, gs = conductances(last[1], h)
            outcomes.append({
                "model": label, "case": name, "initial": initial, "status": refined["status"],
                "last_day": float(refined["time_days"][-1]), "final_state": last,
                "final_gx": gx, "final_gs": gs,
                "final_equilibrium_error_scaled": trajectories[-1]["equilibrium_error_norm_scaled"],
                "state_refinement_scaled": error, "event_refinement_days": time_error,
                "kink_samples": sum(row["jacobian_status"] == "undefined_at_conductance_kink"
                                    for row in trajectories if row["case"] == f"{label}:{name}"),
                "interpretation": "finite_fixed_vpd_outcome_domain_exit_not_mortality_or_global_basin",
            })
            if name == "both_closed_initial":
                require(np.max(abs(refined["state"][:, :2]-initial[:2])) < 1e-12,
                        "Closed water potentials failed invariance")
                exact = np.array([expm(carbon_matrix(c)*t)@stocks for t in refined["time_days"]])
                require(np.max(abs(refined["state"][:, 2:]-exact)/np.array(scales.state[2:])) < 1e-8,
                        "Closed carbon IVP failed exact exponential comparison")
                for t in np.linspace(0., 750*YEAR_DAYS, 301):
                    future = expm(carbon_matrix(c)*t)@stocks
                    require(np.all(future > 0), "Analytical carbon comparison reached numerical underflow")
                    exact_closed.append({
                        "model": label, "day": t, **dict(zip(["L", "W", "R", "S"], future)),
                        "stock_total": float(sum(future)), "capacity_total": 0., "GPP": 0.,
                        "in_integrator_leaf_domain": bool(future[0] > model.feedback.minimum_leaf),
                        "interpretation": "analytic_carbon_subsystem_limit_not_post_exit_full_model_simulation",
                    })
            print(f"{label}/{name}: {refined['status']}, last day={refined['time_days'][-1]:.6g}, "
                  f"gx={gx:.6g}, gs={gs:.6g}", flush=True)
    hopf = hopf_point((-2.01, -1.99), models["E4"])
    for name, rows in [
        ("leaf_curves.csv", curves), ("equilibria.csv", equilibria), ("root_counts.csv", counts),
        ("trajectories.csv", trajectories), ("analytic_closed_carbon.csv", exact_closed),
    ]:
        write_csv(args.output/name, rows)
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for j, label in enumerate(models):
        rows = [row for row in curves if row["case"] == label]
        for status, style, color in [("stable", "-", "tab:blue"), ("unstable", "--", "tab:orange"),
                                     ("nonhyperbolic_within_tolerance", ":", "grey")]:
            segments, segment = [], []
            for row in rows:
                if row["equilibrium_stability"] == status:
                    segment.append(row)
                elif segment:
                    segments.append(segment)
                    segment = []
            if segment:
                segments.append(segment)
            for i, segment in enumerate(segments):
                axes[0, j].plot([row["vpd_mpa"] for row in segment], [row["L"] for row in segment],
                                 style, color=color, label=status if i == 0 else None)
        if label == "E4":
            zoom = axes[0, j].inset_axes([.53, .38, .43, .43])
            for status, color in [("stable", "tab:blue"), ("unstable", "tab:orange"),
                                   ("nonhyperbolic_within_tolerance", "grey")]:
                selected = [row for row in rows if row["vpd_mpa"] > .00478 and row["equilibrium_stability"] == status]
                zoom.scatter([1e3*row["vpd_mpa"] for row in selected], [row["L"] for row in selected], s=3, color=color)
            zoom.set_yscale("log")
            zoom.axhline(models[label].feedback.minimum_leaf, color="black", ls=":", lw=.7)
            zoom.set(title="Formal high-VPD zoom", xlabel="VPD (kPa)", ylabel="L", xlim=(4.788, 4.8075))
            zoom.set_xticks([4.79, 4.806])
            zoom.tick_params(labelsize=7)
            zoom.title.set_fontsize(8)
            zoom.xaxis.label.set_fontsize(8)
            zoom.yaxis.label.set_fontsize(8)
        axes[0, j].set(xlabel="VPD (MPa)", ylabel="Leaf carbon (kg C/m2)", title=f"{label}: complete positive conducting branch")
        selected = [row for row in counts if row["model"] == label]
        boundaries = sorted(set([0., max(row["vpd_mpa"] for row in selected),
                                 landmarks[label]["zero_leaf_vpd_limit_mpa"]] +
                                ([] if landmarks[label]["fold_leaf"] is None else [landmarks[label]["fold_vpd_mpa"]])))
        for i, (left, right) in enumerate(zip(boundaries[:-1], boundaries[1:])):
            count = len(positive_equilibria((left+right)/2, models[label], scales)["roots"])
            axes[1, j].plot([left, right], [count, count], color="tab:blue",
                            label="positive roots (open intervals)" if i == 0 else None)
        axes[1, j].scatter(boundaries, [len(positive_equilibria(v, models[label], scales)["roots"]) for v in boundaries],
                           color="tab:blue", s=14, zorder=3)
        axes[1, j].scatter([row["vpd_mpa"] for row in selected], [row["stable_roots_numerical_spectrum"] for row in selected],
                           color="tab:orange", s=6, label="stable roots (sampled spectrum)", zorder=4)
        axes[1, j].set(xlabel="VPD (MPa)", ylabel="Isolated positive root count", ylim=(-.1, 2.2))
        for ax in axes[:, j]:
            ax.axvline(.003, color="grey", ls=":", label="main range upper")
    for ax in axes.ravel():
        ax.grid(alpha=.25)
        ax.legend(fontsize=8, loc="lower left" if ax is axes[0, 2] else "best")
    fig.suptitle("Scalar equilibrium elimination | closed/zero-carbon continua excluded from positive root counts")
    fig.tight_layout()
    for extension in ["png", "pdf"]:
        fig.savefig(args.output/f"equilibrium_classification.{extension}", dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for ax, name in zip(axes.ravel(), ["L", "W", "R", "S"]):
        for label in models:
            rows = [row for row in curves if row["case"] == label]
            ax.plot([row["vpd_mpa"] for row in rows], [row[f"K_{name}"] for row in rows], label=label)
        ax.set(xlabel="VPD (MPa)", ylabel=f"Equilibrium K_{name}=C_{name} (kg C/m2)")
        ax.axvline(.003, color="grey", ls=":")
        ax.legend()
        ax.grid(alpha=.25)
    fig.suptitle("Self-consistent capacity VECTORS on all positive branches (stability in classification figure)")
    fig.tight_layout()
    for extension in ["png", "pdf"]:
        fig.savefig(args.output/f"equilibrium_capacity.{extension}", dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(2, 4, figsize=(19, 8))
    for j, label in enumerate(["E0", "E4"]):
        for name, color in zip(["wet_initial", "dry_initial", "both_closed_initial", "xylem_closed_low_leaf"],
                               ["tab:green", "tab:red", "tab:purple", "tab:orange"]):
            rows = [row for row in trajectories if row["case"] == f"{label}:{name}"]
            t = np.array([row["day"] for row in rows])
            axes[j, 0].plot(t, [row["leaf_mpa"] for row in rows], color=color, label=name)
            axes[j, 1].plot(t, [row["gx"]/h.gx_max for row in rows], color=color, label=f"{name}: gx")
            axes[j, 1].plot(t, [row["gs"]/h.gs_max for row in rows], color=color, ls="--")
            axes[j, 2].plot(t, [row["GPP"] for row in rows], color=color, label=name)
            axes[j, 3].plot(t/YEAR_DAYS, [row["capacity_total"] for row in rows], color=color, label=f"{name}: K")
            axes[j, 3].plot(t/YEAR_DAYS, [row["stock_total"] for row in rows], color=color, ls="--")
        for ax in axes[j, :3]:
            ax.set_xscale("symlog", linthresh=.01)
            ax.set_xlabel("Time (day), early/late scales")
        axes[j, 0].set(ylabel="Leaf potential (MPa)", title=f"{label}: fixed VPD={experiment_vpd} MPa")
        axes[j, 1].set(ylabel="Conductance / own maximum", title="Solid gx; dashed gs")
        axes[j, 2].set(ylabel="GPP (kg C/m2/day)", title="Input is not stock")
        axes[j, 3].set(xlabel="Time (year)", ylabel="kg C/m2", title="Solid total K; dashed total C")
    for ax in axes.ravel():
        ax.grid(alpha=.25)
        ax.legend(fontsize=6)
    fig.suptitle("Eight representative initial states only | stop at L_min; no new basin/ramp/mortality mechanism")
    fig.tight_layout()
    for extension in ["png", "pdf"]:
        fig.savefig(args.output/f"hydraulic_capacity_timeseries.{extension}", dpi=180)
    plt.close(fig)
    sources = [
        args.config, config_path, ROOT/"examples/analyze_spac_carbon_equilibria.py",
        ROOT/"examples/analyze_spac_carbon_capacity.py", ROOT/"tests/test_spac_carbon_equilibria.py",
        *[ROOT/f"src/control_carbon/{name}.py" for name in [
            "spac_carbon_equilibria", "spac_carbon_capacity", "spac_carbon_stability",
            "spac_carbon_mpa", "spac_carbon_two_way", "spac_carbon", "manzoni_hydraulics",
            "spac_carbon_e4", "spac_carbon_bifurcation",
        ]],
        ROOT/"docs/research_plan_hydraulic_failure_capacity_20261014.md",
    ]
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True).stdout.strip()
    result = {
        "baseline_git_SHA": BASELINE, "run_git_HEAD": head,
        "settings": settings, "baseline_config": cfg, "state_scales": scales,
        "landmarks_and_shape_proofs": landmarks, "stock_ratios": stock_ratios(models["E0"]),
        "positive_domain_scope": "global_positive_conducting_branch_under_checked_baseline_shape_hypotheses",
        "stability_scope": "full_J_numerical_spectrum_E0_E1_analytical_reasoning_in_note_E4_not_interval_certified",
        "hopf_E4_vpd_mpa": hopf["vpd_mpa"], "representative_vpd_mpa": experiment_vpd,
        "E4_detected_stability_crossings": crossings,
        "E4_stability_scan_scope": {"leaf_samples": len(scan), "lower_leaf": 1e-9,
                                   "coverage": "numerical_scan_not_interval_stability_certificate"},
        "outcomes": outcomes,
        "solver_settings": {"root_method": "brentq_on_proved_monotone_intervals", "leaf_xtol": 1e-14,
                            "IVP": "Radau_log_carbon", "rtol": [1e-9, 1e-10], "atol_over_rtol": .01,
                            "exact_closed_comparison": "expm(M*t), 0--750year, formal_carbon_subsystem_only",
                            "Decimal_constant_margin_precision": 40},
        "acceptance": {"full_scaled_root_residual_per_day": 1e-8, "IVP_state_refinement_scaled": 1e-6,
                       "exit_time_refinement_days": .1, "water_carbon_budget_residual": 1e-14,
                       "closed_exact_carbon_scaled_error": 1e-8},
        "versions": {"numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "limits": [
            "Formal zero-carbon/nonconducting continua are not isolated hyperbolic sinks.",
            "gx=0 alone does not imply instantaneous gs=0, GPP=0 or K=0.",
            "Analytic closed-carbon long-time limit is not a post-exit full-model simulation.",
            "Values above .003 MPa are formal branch classification, not new routine forcing.",
            "No nonautonomous ramp, global basin, new mortality or same-slice periodic continuation.",
        ],
        "source_sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
    }
    (args.output/"results.json").write_text(json.dumps(result, default=encode, indent=2, allow_nan=False)+"\n")
    print(f"Saved {len(equilibria)} reconstructed roots and {len(outcomes)} representative outcomes to {args.output}", flush=True)


if __name__ == "__main__":
    main()
