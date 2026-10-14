"""Stage C0/C1: capacity identities, counterexamples and E4 diagnostics."""

import argparse
import csv
from dataclasses import asdict, is_dataclass
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
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from control_carbon.capacity_feedback import CubicCapacityFeedback
from control_carbon.manzoni_hydraulics import HydraulicParameters, transpiring_equilibria
from control_carbon.spac_carbon import (
    YEAR_DAYS, CarbonParameters, hydraulic_factor, positive_carbon_equilibrium,
)
from control_carbon.spac_carbon_two_way import StructuralFeedback
from control_carbon.spac_carbon_mpa import PressureSPAC
from control_carbon.spac_carbon_e4 import core_rhs, core_jacobian, hopf_point, periodic_orbit
from control_carbon.spac_carbon_capacity import (
    conditional_carbon_roots, diagnose_capacity, frozen_equilibria, reference_scales,
)

ROOT = Path(__file__).resolve().parents[1]
BASELINE = "2ef9438c080d8183edc1efeff6e3a3b9fad62eb1"
STATE_NAMES = ("xylem_mpa", "leaf_mpa", "L", "W", "R", "S")
CARBON_NAMES = STATE_NAMES[2:]


def encode(value):
    if is_dataclass(value):
        return asdict(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, complex):
        return [value.real, value.imag]
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Unsupported result: {type(value)}")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def write_csv(path, records):
    require(bool(records), f"No diagnostic records for {path}")
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def record(label, day, y, vpd, data, equilibrium=None):
    result = {"case": label, "day": day, "vpd_mpa": vpd}
    result.update(dict(zip(STATE_NAMES, y)))
    for prefix, values in [
        ("K", data["capacity"]), ("P", data["potential"]),
        ("conditional", data["conditional_roots"][-1].state),
        ("equilibrium", equilibrium[2:] if equilibrium is not None else [None]*4),
    ]:
        result.update({f"{prefix}_{name}": value for name, value in zip(CARBON_NAMES, values)})
    for key in [
        "potential_sum", "potential_norm_scaled", "conditional_distance_norm_scaled",
        "full_residual_norm_scaled_per_day", "carbon_residual_norm_scaled_per_day",
        "water_residual_norm_scaled_per_day", "equilibrium_error_norm_scaled",
        "potential_energy_rate_scaled_per_day", "capacity_rate_norm_scaled_per_day",
        "identity_residual_max", "potential_identity_residual_max",
        "carbon_budget_residual", "water_budget_residual", "matrix_condition_2",
        "reproduction_number", "conditional_target_kind",
        "conditional_target_in_positive_leaf_domain", "equilibrium_reference_status",
        "equilibrium_stability", "equilibrium_type", "equilibrium_unstable_dimensions",
        "equilibrium_center_dimensions", "state_above_leaf_exit_boundary",
    ]:
        result[key] = data[key]
    result["local_carbon_max_real_per_day"] = float(max(data["local_carbon_eigenvalues"].real))
    result["local_full_max_real_per_day"] = float(max(data["local_full_eigenvalues"].real))
    for key in ["local_carbon_eigenvalues", "local_full_eigenvalues", "equilibrium_eigenvalues"]:
        result[f"{key}_per_day"] = json.dumps(data[key], default=encode, allow_nan=False)
    return result


def legacy_audit():
    results = {}
    for folder in [
        "stage_a_final_20261009", "stage_b1_final_20261009", "stage_b2_verified_20261012",
        "pressure_review_final_v2_20261013", "stage_b3_e4_final_v2_20261013",
    ]:
        path = ROOT/"artifacts/spac_carbon"/folder/"results.json"
        if not path.exists():
            results[folder] = {"status": "not_present_locally"}
            continue
        saved = json.loads(path.read_text())
        if "source_sha256" in saved:
            hashes = saved["source_sha256"]
        elif "sha256" in saved:
            hashes = saved["sha256"]
        elif "manifest" in saved and "sha256" in saved["manifest"]:
            hashes = saved["manifest"]["sha256"]
        else:
            raise RuntimeError(f"Unrecognized archived source-hash schema: {path}")
        changed = [key for key, digest in hashes.items()
                   if not (ROOT/key).exists() or hashlib.sha256((ROOT/key).read_bytes()).hexdigest() != digest]
        results[folder] = {
            "status": "saved_results_present_not_campaign_rerun",
            "results_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "source_hash_count": len(hashes), "changed_or_missing_sources": changed,
            "hash_note": "checked against current files; recognized original manifest/source_sha256/sha256 schema",
        }
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT/"configs/spac_carbon_pressure_review.json")
    parser.add_argument("--output", type=Path, default=ROOT/"artifacts/spac_carbon"/datetime.now().strftime("stage_c_capacity_%Y%m%d_%H%M%S"))
    args = parser.parse_args()
    settings = json.loads(args.config.read_text())
    require(settings["units"]["vpd"] == "MPa", "Capacity runner requires MPa settings")
    baseline_path = args.config.parent/settings["model_config"]
    cfg = json.loads(baseline_path.read_text())
    h = HydraulicParameters(**cfg["hydraulic"])
    k = cfg["carbon"]
    c = CarbonParameters(
        k["assimilation_max_per_year"]/YEAR_DAYS, tuple(k["allocation"]),
        tuple(value/YEAR_DAYS for value in k["losses_per_year"]), tuple(k["litter_fractions"]),
        k["specific_leaf_area"], k["carbon_dry_fraction"], k["extinction"],
    )
    leaf = transpiring_equilibria(cfg["root_potential"], cfg["reference_vpd"], h)[0]
    stocks = positive_carbon_equilibrium(hydraulic_factor(leaf, h), c)
    require(stocks is not None, "Fixed structural reference must have positive carbon")
    models = {
        label: PressureSPAC(
            h, c, StructuralFeedback(stocks[0], stocks[1], cfg["leaf_floor"], cfg["minimum_leaf"],
                                     *strengths, wood_exponent=cfg["wood_exponent"]),
            cfg["root_potential"], settings["atmospheric_pressure_mpa"],
        ) for label, strengths in cfg["variants"].items()
    }
    scales = reference_scales(models["E0"])
    require(settings["routine_vpd_mpa_bounds"] == [0., .003],
            "This bounded C1 campaign is specified for 0--.003 MPa")
    args.output.mkdir(parents=True, exist_ok=False)
    bounds = cfg["curve_bounds"]
    vgrid = np.unique(np.r_[np.linspace(.00001, .003, 45), .0013, .00244, .0028])
    branches, counters, conditional, checks = [], [], [], []
    branch_inventory = {}
    for label, model in models.items():
        folds = model.stationary_points(bounds)
        inventory = {
            "leaf_bounds_mpa": bounds, "stationary_scan_samples": 300,
            "detected_folds": [model.equilibrium_curve(z)[1] for z in folds],
            "coverage": "sampled_partition_of_conducting_positive_branch_not_global_certificate",
            "zero_vpd": "outside_open_root_to_leaf_drop_parameterization_not_absent_equilibrium",
            "formal_zero_full_roots": "not_defined_by_positive_leaf_model",
        }
        # Doubling the scan checks detected stationary points, not global completeness.
        doubled = model.stationary_points(bounds, samples=600)
        require(len(folds) == len(doubled) and np.allclose(folds, doubled, atol=1e-10, rtol=0),
                f"{label} stationary scan refinement failed")
        local_grid = np.unique(np.r_[vgrid, [model.equilibrium_curve(z)[1] for z in folds]])
        require(np.all((local_grid >= 0) & (local_grid <= .003)),
                f"{label} fold diagnostic exceeds this campaign's main VPD range")
        inventory["sampled_vpd_mpa"] = local_grid
        for v in local_grid:
            roots = model.equilibria(v, bounds, folds)
            for y in roots:
                data = diagnose_capacity(y, v, model, scales, equilibrium=y)
                branches.append(record(label, 0., y, v, data, y))
                checks.append(data)
        inventory["sampled_root_records"] = sum(row["case"] == label for row in branches)
        branch_inventory[label] = inventory
        for psi in [-.5, -1.2, -2., -3., -4.2]:
            for root in conditional_carbon_roots(psi, model):
                conditional.append({
                    "case": label, "leaf_potential_mpa": psi, "kind": root.kind,
                    "carbon_stability": root.stability,
                    "in_positive_leaf_domain": root.in_positive_leaf_domain,
                    **dict(zip(CARBON_NAMES, root.state)),
                    "carbon_max_real_per_day": float(max(root.eigenvalues.real)),
                })
        target = positive_carbon_equilibrium(hydraulic_factor(-1.2, h), c)
        require(target is not None, "Counterexample requires a positive conditional carbon root")
        nonwater_root = np.r_[-.8, -1.2, target]
        cancelled = nonwater_root.copy()
        cancelled[3] += .05
        cancelled[4] -= .05
        input_off = np.r_[-4.1, -4.2, stocks]
        for name, y in [
            ("P_zero_H_nonzero", nonwater_root), ("sum_P_zero_P_nonzero", cancelled),
            ("input_off_positive_stock", input_off),
        ]:
            data = diagnose_capacity(y, .002, model, scales)
            counters.append(record(f"{label}:{name}", 0., y, .002, data))
            checks.append(data)
            if name == "P_zero_H_nonzero":
                require(data["potential_norm_scaled"] < 1e-12 and
                        data["water_residual_norm_scaled_per_day"] > .1,
                        f"{label}: P=0 / H!=0 counterexample failed")
            elif name == "sum_P_zero_P_nonzero":
                require(abs(data["potential_sum"]) < 1e-12 and data["potential_norm_scaled"] > .01,
                        f"{label}: potential cancellation counterexample failed")
            else:
                require(max(abs(data["capacity"])) == 0 and
                        data["carbon_residual_norm_scaled_per_day"] > 0,
                        f"{label}: input-off capacity is not instantaneous equilibrium")
    print("E0-E4 branch and counterexample diagnostics verified.", flush=True)

    model = models["E4"]
    reference = frozen_equilibria(.0024, bounds, model, scales)["roots"]
    require(len(reference) == 1 and reference[0]["stability"] == "stable",
            "Transient requires a uniquely selected stable root in the stated domain")
    eq = reference[0]["state"]
    initial = eq.copy()
    initial[1] -= .02
    initial[2] *= .95
    initial[3] *= .99
    times = np.unique(np.r_[0., np.geomspace(.01, 5., 31), np.linspace(10., 20*YEAR_DAYS, 241)])
    transient = model.simulate(initial, times, lambda t: .0024, rtol=1e-9)
    refined = model.simulate(initial, times, lambda t: .0024, rtol=1e-10)
    require(transient["status"] == refined["status"] == "completed", "Bounded transient left the vegetation domain")
    transient_error = float(np.max(abs(transient["state"]-refined["state"])/np.array(scales.state)))
    require(transient_error < 1e-7, "Transient tolerance refinement failed")
    transient_rows = []
    for t, y in zip(refined["time_days"], refined["state"]):
        data = diagnose_capacity(y, .0024, model, scales, equilibrium=eq)
        transient_rows.append(record("E4:finite_transient_not_basin_label", t, y, .0024, data, eq))
        checks.append(data)

    hopf = hopf_point((-2.01, -1.99), model)
    cycle = periodic_orbit(hopf, .001, model, tol=1e-9)
    cycle_eq = model.equilibrium_curve(brentq(
        lambda z: model.equilibrium_curve(z)[1]-cycle["vpd_mpa"], -2.12, -1.9, xtol=1e-13,
    ))[0]
    periodic_rows = []
    for phase, y in zip(cycle["phase"], cycle["state"]):
        data = diagnose_capacity(y, cycle["vpd_mpa"], model, scales, equilibrium=cycle_eq)
        periodic_rows.append(record("E4:repelling_cycle_not_attractor",
                                    phase*cycle["period_days"], y, cycle["vpd_mpa"], data, cycle_eq))
        checks.append(data)
    flow_refinement = []
    for steps, tol in [(80, 1e-10), (160, 1e-11)]:
        direct = solve_ivp(
            lambda t, y: core_rhs(y, cycle["vpd_mpa"], model),
            (0., cycle["period_days"]), cycle["state"][0, :4], method="Radau",
            jac=lambda t, y: core_jacobian(y, cycle["vpd_mpa"], model),
            rtol=tol, atol=tol*.01, max_step=cycle["period_days"]/steps,
            t_eval=cycle["phase"]*cycle["period_days"],
        )
        require(direct.success, f"E4 independent cycle IVP failed: {direct.message}")
        error = float(np.max(abs(direct.y.T-cycle["state"][:, :4])/np.array(scales.state[:4])))
        require(error < 1e-7, "E4 cycle timestep/tolerance check failed")
        flow_refinement.append({"max_step_days": cycle["period_days"]/steps, "rtol": tol,
                                "max_scaled_core_state_error": error})
    energy_rates = [row["potential_energy_rate_scaled_per_day"] for row in periodic_rows]
    require(min(energy_rates) < -1e-14 and max(energy_rates) > 1e-14,
            "Nonmonotone potential not resolved on this periodic orbit")
    multipliers = cycle["full_floquet"]
    neutral = int(np.argmin(abs(multipliers-1)))
    transverse = np.delete(multipliers, neutral)
    require(sum(abs(transverse) > 1+1e-6) == 1 and sum(abs(transverse) < 1-1e-6) == 4,
            "Saddle-cycle Floquet dimensions unresolved")
    radial = float(max(abs(transverse)))
    saddle_diagnostics = {
        "stable_transverse_dimensions": 4, "unstable_transverse_dimensions": 1,
        "stable_manifold_dimension_in_full_system": 5,
        "stable_manifold_codimension": 1,
        "radial_e_folding_years": cycle["period_days"]/np.log(radial)/YEAR_DAYS,
        "claim": "local_saddle_cycle_boundary_candidate_not_verified_global_recovery_exit_boundary",
    }
    scalar = CubicCapacityFeedback()
    scalar_result = {
        "stock": 2.8, "driver": 0., "time_unit": "year_arbitrary_carbon_unit_not_SPAC",
        "capacity_slope": scalar.capacity_derivative(2.8, 0.),
        "potential_energy_rate": scalar.potential_energy_derivative(2.8, 0.),
    }
    require(scalar_result["potential_energy_rate"] > 0, "Scalar potential-growth counterexample failed")
    maxima = {
        key: float(max(abs(data[key]) for data in checks)) for key in [
            "identity_residual_max", "potential_identity_residual_max",
            "water_budget_residual", "carbon_budget_residual",
        ]
    }
    require(maxima["identity_residual_max"] < 1e-14 and maxima["potential_identity_residual_max"] < 1e-11
            and maxima["water_budget_residual"] < 1e-14 and maxima["carbon_budget_residual"] < 1e-14,
            "Capacity identity or conservation budget failed")
    print(f"E4 cycle P energy rate range: [{min(energy_rates):.6g}, {max(energy_rates):.6g}]/day; "
          f"radial e-folding {saddle_diagnostics['radial_e_folding_years']:.3f} years.", flush=True)

    for name, rows in [("branches.csv", branches), ("counterexamples.csv", counters),
                       ("conditional_roots.csv", conditional), ("transient.csv", transient_rows),
                       ("periodic_capacity.csv", periodic_rows)]:
        write_csv(args.output/name, rows)
    fig, axes = plt.subplots(2, 3, figsize=(16, 8))
    years = np.array([row["day"] for row in transient_rows])/YEAR_DAYS
    for key, label in [("L", "actual L"), ("K_L", "instantaneous K_L"),
                       ("conditional_L", "conditional C*_L"), ("equilibrium_L", "coupled equilibrium L*")]:
        axes[0, 0].plot(years, [row[key] for row in transient_rows], label=label)
    axes[0, 0].set(xlabel="Time (year)", ylabel="kg C/m2", title="E4 frozen transient: four different objects")
    for key, label in [("potential_norm_scaled", "||P||"), ("conditional_distance_norm_scaled", "conditional error"),
                       ("equilibrium_error_norm_scaled", "full equilibrium error")]:
        axes[0, 1].semilogy(years, [row[key] for row in transient_rows], label=label)
    axes[0, 1].set(xlabel="Time (year)", ylabel="Fixed-scale norm", title="No common notion of distance")
    for key, label in [("full_residual_norm_scaled_per_day", "full F"), ("carbon_residual_norm_scaled_per_day", "carbon RHS")]:
        axes[0, 2].semilogy(years, [row[key] for row in transient_rows], label=label)
    axes[0, 2].set(xlabel="Time (year)", ylabel="Scaled rate (day^-1)", title="Residual is not equilibrium distance")
    e4_branch = [row for row in branches if row["case"] == "E4"]
    for key, label in [("local_full_max_real_per_day", "full J"), ("local_carbon_max_real_per_day", "conditional J_CC")]:
        axes[1, 0].plot([row["vpd_mpa"] for row in e4_branch], [row[key] for row in e4_branch], label=label)
    axes[1, 0].axhline(0, color="black", lw=.8)
    axes[1, 0].axvline(hopf["vpd_mpa"], color="grey", ls=":", label="Hopf")
    axes[1, 0].set_yscale("symlog", linthresh=1e-6)
    axes[1, 0].set(xlabel="VPD (MPa)", ylabel="Largest Re(lambda) (day^-1)", title="Stable carbon block != stable full system")
    phase = cycle["phase"]
    for key, label in [("potential_norm_scaled", "||P||"), ("equilibrium_error_norm_scaled", "full equilibrium error")]:
        axes[1, 1].plot(phase, [row[key] for row in periodic_rows], label=label)
    axes[1, 1].set(xlabel="Cycle phase", ylabel="Fixed-scale norm", title="Repelling cycle: NOT an attracting capacity")
    axes[1, 2].plot(phase, energy_rates, label="d(||P||^2/2)/dt")
    axes[1, 2].axhline(0, color="black", lw=.8)
    axes[1, 2].set(xlabel="Cycle phase", ylabel="Scaled potential energy rate (day^-1)", title="Potential both increases and decreases")
    for ax in axes.ravel():
        ax.grid(alpha=.25)
        ax.legend(fontsize=8)
    fig.suptitle("Stage C1: capacity map / conditional root / coupled root / trajectory | synthetic E0-E4")
    fig.tight_layout()
    for extension in ["png", "pdf"]:
        fig.savefig(args.output/f"capacity_diagnostics.{extension}", dpi=180)
    plt.close(fig)

    sources = [
        args.config, baseline_path, ROOT/"examples/analyze_spac_carbon_capacity.py",
        ROOT/"tests/test_spac_carbon_capacity.py",
        *[ROOT/f"src/control_carbon/{name}.py" for name in [
            "spac_carbon_capacity", "spac_carbon", "spac_carbon_two_way", "spac_carbon_mpa",
            "spac_carbon_e4", "spac_carbon_bifurcation", "spac_carbon_stability",
            "manzoni_hydraulics", "manzoni_extensions", "visitc_hydrology", "capacity_feedback",
        ]],
        ROOT/"docs/SPAC_nonlinear_carbon_storage_capacity_research_plan.md",
        ROOT/"docs/spac_carbon_capacity_definition_ledger.md",
        ROOT/"docs/capacity_qse_attractor_exploration.tex",
    ]
    git_head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
                              capture_output=True, text=True).stdout.strip()
    result = {
        "classification": "analytic_identities_and_precision_checked_numerics_not_global_basin_or_R_tipping",
        "baseline_git_SHA": BASELINE, "run_git_HEAD": git_head,
        "settings": settings, "baseline_config": cfg, "state_scales": scales,
        "scale_rule": "fixed_1_MPa_and_original_structural_carbon_references_not_moving_targets",
        "versions": {"numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "branch_inventory": branch_inventory, "branch_vpd_samples_mpa": vgrid,
        "spectral_tolerance_per_day": 1e-10, "norm_convention": "Euclidean_after_fixed_diagonal_scaling",
        "checks": maxima, "legacy_artifact_audit": legacy_audit(),
        "transient": {"vpd_mpa": .0024, "duration_years": 20., "rtol": [1e-9, 1e-10],
                      "scaled_refinement_error": transient_error, "status": refined["status"],
                      "initial": initial, "equilibrium": eq,
                      "final_equilibrium_error_scaled": transient_rows[-1]["equilibrium_error_norm_scaled"],
                      "interpretation": "finite_positive_transient_not_asymptotic_basin_classification"},
        "cycle": cycle, "cycle_solver": {"amplitude": .001, "collocation_tolerance": 1e-9},
        "cycle_IVP_refinement": flow_refinement,
        "solver_settings": {
            "conditional_positive_root": {"method": "brentq_divided_by_leaf_to_exclude_zero", "xtol": 1e-14, "rtol": 1e-14},
            "stationary_and_full_branch_roots": {"method": "brentq", "xtol": 1e-12,
                                                 "stationary_scan_samples": [300, 600]},
            "transient": {"method": "Radau_log_carbon", "rtol": [1e-9, 1e-10], "atol_over_rtol": .01,
                          "stop": "minimum_leaf_domain_exit_not_mortality"},
            "cycle": {"method": "periodic_solve_bvp_plus_passive_R_S_lift", "tol": 1e-9,
                      "max_nodes": 5000, "initial_nodes": 161},
            "IVP_refinement": "Radau_with_analytic_core_J_and_both_tolerance_and_max_step_refined",
        },
        "acceptance_limits": {
            "identity_residual_max": 1e-14, "potential_identity_residual_max": 1e-11,
            "water_budget_residual": 1e-14, "carbon_budget_residual": 1e-14,
            "scaled_trajectory_refinement_error": 1e-7,
            "stationary_point_refinement_error_mpa": 1e-10,
            "potential_energy_rate_sign_resolution_per_day": 1e-14,
            "Floquet_transverse_resolution": 1e-6,
        },
        "periodic_potential_energy_rate_range": [min(energy_rates), max(energy_rates)],
        "periodic_potential_norm_range": [min(row["potential_norm_scaled"] for row in periodic_rows),
                                          max(row["potential_norm_scaled"] for row in periodic_rows)],
        "saddle_cycle_diagnostics": saddle_diagnostics, "scalar_counterexample": scalar_result,
        "complex_encoding": "[real, imaginary]",
        "limits": [
            "No new physiological mechanism or change to Stage A/B equations.",
            "Formal zero carbon is a conditional root, not a defined full-model mortality equilibrium.",
            "Root enumeration is confined to the declared conducting positive branch and sampled stationary partition.",
            "Stable manifold dimension is inferred from numerical Floquet modes, not a computed global manifold.",
            "No forcing-speed comparison, alternate attractor, true R-tipping, or global recovery/exit boundary.",
        ],
        "source_sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
    }
    (args.output/"results.json").write_text(json.dumps(result, default=encode, indent=2, allow_nan=False)+"\n")
    print(f"Saved {len(branches)} branch records, {len(counters)} counterexamples, "
          f"{len(periodic_rows)} periodic capacity samples to {args.output}", flush=True)


if __name__ == "__main__":
    main()
