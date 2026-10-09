"""Reproducible synthetic Stage B1 branches, stability and domain exits."""

import argparse
import csv
from dataclasses import asdict, replace
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
from scipy.optimize import brentq, least_squares

from control_carbon.manzoni_hydraulics import HydraulicParameters, conductances, transpiring_equilibria
from control_carbon.spac_carbon import (
    YEAR_DAYS, CarbonParameters, carbon_budget, frozen_coefficient_capacity,
    hydraulic_factor, positive_carbon_equilibrium, reproduction_number,
)
from control_carbon.spac_carbon_two_way import (
    StructuralFeedback, feedback_factors, two_way_rhs, two_way_jacobian,
    equilibrium_curve, stationary_points, equilibria_at_vpd,
    continue_curve, simulate_two_way, structural_fluxes,
)

ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def stability_boundaries(bounds, root_potential, hydraulic, carbon, feedback, samples=300):
    def alpha(leaf):
        state, d, _ = equilibrium_curve(leaf, root_potential, hydraulic, carbon, feedback)
        return float(np.max(np.linalg.eigvals(two_way_jacobian(
            state, root_potential, d, hydraulic, carbon, feedback,
        )).real))
    grid = np.linspace(*bounds, samples)
    alphas = [alpha(z) for z in grid]
    crossings = []
    for left, right, fl, fr in zip(grid[:-1], grid[1:], alphas[:-1], alphas[1:]):
        if fl*fr < 0:
            leaf = brentq(alpha, left, right, xtol=1e-12)
            state, d, _ = equilibrium_curve(leaf, root_potential, hydraulic, carbon, feedback)
            eig = np.linalg.eigvals(two_way_jacobian(state, root_potential, d, hydraulic, carbon, feedback))
            dominant = eig[np.argmax(eig.real)]
            crossings.append({
                "leaf_MPa": float(leaf), "D": d,
                "type": "oscillatory_stability_boundary_candidate" if abs(dominant.imag) > 1e-6 else "real_zero_candidate",
                "dominant_imaginary_per_day": float(abs(dominant.imag)),
            })
    return crossings


def fold_checks(leaf, root_potential, hydraulic, carbon, feedback):
    y, d, _ = equilibrium_curve(leaf, root_potential, hydraulic, carbon, feedback)
    jac = two_way_jacobian(y, root_potential, d, hydraulic, carbon, feedback)
    schur = jac[:2, :2]-jac[:2, 2:]@np.linalg.solve(jac[2:, 2:], jac[2:, :2])
    u, singular, vh = np.linalg.svd(jac)
    left, right = u[:, -1], vh[-1]
    gs = conductances(leaf, hydraulic)[1]
    fe = feedback_factors(y[2], y[3], feedback)[1]
    fd = np.r_[0., -gs*fe/hydraulic.cl, np.zeros(4)]
    quadratic = []
    for step in [1e-3, 5e-4]:
        value = (two_way_rhs(y+step*right, root_potential, d, hydraulic, carbon, feedback)
                 -2*two_way_rhs(y, root_potential, d, hydraulic, carbon, feedback)
                 +two_way_rhs(y-step*right, root_potential, d, hydraulic, carbon, feedback))/step**2
        quadratic.append(float(left@value))
    eig = np.linalg.eigvals(jac)
    require(np.count_nonzero(abs(eig) < 1e-9) == 1, "fold zero eigenvalue is not simple")
    require(abs(left@fd) > 1e-8 and abs(quadratic[-1]) > 1e-8, "fold nondegeneracy failed")
    require(abs(quadratic[0]-quadratic[1]) < 1e-3*abs(quadratic[-1]), "fold Hessian refinement failed")
    return {
        "D": d, "leaf_MPa": leaf, "state": y.tolist(),
        "R0": reproduction_number(hydraulic_factor(leaf, hydraulic), carbon),
        "eigenvalues_per_day": [[float(z.real), float(z.imag)] for z in eig],
        "schur_smallest_abs_eigenvalue": float(min(abs(np.linalg.eigvals(schur)))),
        "smallest_singular_value": float(singular[-1]),
        "transversality": float(left@fd),
        "quadratic_directional_coefficient_refinement": quadratic,
        "classification": "numerically validated simple saddle-node; not interval-certified",
    }


def multistart_roots(vpd, bounds, root_potential, hydraulic, carbon, feedback):
    """Independent bounded full-six-state solves; retain failed attempts too."""
    scale = np.r_[1., 1., YEAR_DAYS/np.array([
        feedback.leaf_reference, feedback.wood_reference,
        feedback.leaf_reference*4, feedback.wood_reference*3,
    ])]
    lower = np.r_[bounds[0], bounds[0], feedback.minimum_leaf*(1+1e-6), [1e-12]*3]
    upper = np.r_[root_potential, bounds[1], [np.inf]*4]
    attempts, roots = [], []
    for seed in np.linspace(bounds[0]+.005, bounds[1]-.005, 7):
        initial = equilibrium_curve(seed, root_potential, hydraulic, carbon, feedback)[0]
        solution = least_squares(
            lambda y: two_way_rhs(y, root_potential, vpd, hydraulic, carbon, feedback)*scale,
            initial, bounds=(lower, upper),
            jac=lambda y: two_way_jacobian(y, root_potential, vpd, hydraulic, carbon, feedback)*scale[:, None],
            ftol=1e-12, xtol=1e-12, gtol=1e-12, max_nfev=500,
        )
        y = solution.x
        conditional = positive_carbon_equilibrium(hydraulic_factor(y[1], hydraulic), carbon)
        residual = float(max(abs(two_way_rhs(y, root_potential, vpd, hydraulic, carbon, feedback))))
        carbon_error = (float(max(abs(y[2:]-conditional)/np.maximum(conditional, 1e-12)))
                        if conditional is not None else None)
        accepted = bool(solution.success and residual < 1e-10 and carbon_error is not None and carbon_error < 1e-6)
        attempts.append({
            "D": vpd, "seed_leaf": float(seed), "accepted_root": accepted,
            "solver_message": solution.message, "nfev": solution.nfev,
            "full_rhs_residual": residual, "conditional_carbon_relative_error": carbon_error,
            "final_state": y.tolist(),
        })
        if accepted and not any(abs(y[1]-other[1]) < 1e-7 for other in roots):
            roots.append(y)
    return sorted(roots, key=lambda y: -y[1]), attempts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT/"configs/spac_carbon_stage_b_synthetic.json")
    parser.add_argument("--output", type=Path, default=ROOT/"artifacts/spac_carbon"/datetime.now().strftime("stage_b1_%Y%m%d_%H%M%S"))
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    args.output.mkdir(parents=True, exist_ok=False)
    h = HydraulicParameters(**cfg["hydraulic"])
    cdata = cfg["carbon"]
    c = CarbonParameters(
        cdata["assimilation_max_per_year"]/YEAR_DAYS, tuple(cdata["allocation"]),
        tuple(k/YEAR_DAYS for k in cdata["losses_per_year"]), tuple(cdata["litter_fractions"]),
        cdata["specific_leaf_area"], cdata["carbon_dry_fraction"], cdata["extinction"],
    )
    root_potential = cfg["root_potential"]
    leaf0 = transpiring_equilibria(root_potential, cfg["reference_vpd"], h)[0]
    stocks0 = positive_carbon_equilibrium(hydraulic_factor(leaf0, h), c)
    require(stocks0 is not None, "Stage A reference must be positive")
    initial = np.r_[(root_potential+leaf0)/2, leaf0, stocks0]
    base = StructuralFeedback(
        stocks0[0], stocks0[1], cfg["leaf_floor"], cfg["minimum_leaf"],
        wood_exponent=cfg["wood_exponent"],
    )
    variants = {
        name: replace(base, leaf_supply=values[0], leaf_demand=values[1], wood_supply=values[2])
        for name, values in cfg["variants"].items()
    }
    bounds = tuple(cfg["curve_bounds"])
    report = {
        "classification": "model-agnostic synthetic theory; not VISIT or physiological mortality",
        "git_HEAD": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "config": cfg, "hydraulic": asdict(h), "carbon_per_day": asdict(c),
        "fixed_reference": {"state": initial.tolist(), "leaf": base.leaf_reference, "wood": base.wood_reference},
        "vegetation_domain": "L > L_min; event exits are NOT demonstrated alternative attractors",
        "branches": {}, "trajectories": [], "basin": [], "sensitivities": [],
        "checks": {},
    }
    branch_rows, spectra, root_rows, continuation_rows = [], [], [], []
    max_residual = max_arc_refinement = 0.
    stationary = {}
    for name, feedback in variants.items():
        points = stationary_points(bounds, root_potential, h, c, feedback)
        refined = stationary_points(bounds, root_potential, h, c, feedback, samples=600)
        require(len(points) == len(refined), "stationary point grid refinement changed count")
        require(all(abs(a-b) < 1e-9 for a, b in zip(points, refined)), "stationary point refinement failed")
        stationary[name] = points
        boundaries = stability_boundaries(bounds, root_potential, h, c, feedback)
        boundary_refined = stability_boundaries(bounds, root_potential, h, c, feedback, samples=600)
        require(len(boundaries) == len(boundary_refined), "stability boundary grid count changed")
        require(all(abs(a["D"]-b["D"]) < 1e-9 for a, b in zip(boundaries, boundary_refined)), "stability boundary refinement failed")
        continuation = []
        for step in [.02, .01]:
            records, status = continue_curve(-.01, step, 4000, (-2.12, bounds[1]), root_potential, h, c, feedback)
            require(status == "domain_boundary", f"{name} continuation exhausted prescribed steps")
            continuation.append(records)
        coarse, fine = continuation
        for step, records in zip([.02, .01], continuation):
            for y, d in records:
                continuation_rows.append({"variant": name, "arclength_step": step,
                                          "D_mol_mol": d, "leaf_MPa": y[1],
                                          "xylem_MPa": y[0], "leaf_carbon": y[2],
                                          "total_carbon": float(sum(y[2:]))})
        fine_leaf = np.array([y[1] for y, d in fine])
        fine_d = np.array([d for y, d in fine])
        errors = [abs(d-np.interp(y[1], fine_leaf[::-1], fine_d[::-1]))
                  for y, d in coarse if y[1] >= fine_leaf[-1]]
        arc_error = float(max(errors))
        max_arc_refinement = max(max_arc_refinement, arc_error)
        require(arc_error < 2e-5, "pseudo-arclength interpolation refinement failed")
        # Plot at uniform leaf intervals too, without depending on D monotonicity.
        for leaf in np.linspace(*bounds, 600):
            y, d, slope = equilibrium_curve(leaf, root_potential, h, c, feedback)
            residual = float(max(abs(two_way_rhs(y, root_potential, d, h, c, feedback))))
            max_residual = max(max_residual, residual)
            eig = np.linalg.eigvals(two_way_jacobian(y, root_potential, d, h, c, feedback))
            branch_rows.append({
                "variant": name, "D_mol_mol": d, "leaf_MPa": leaf, "xylem_MPa": y[0],
                "leaf_carbon": y[2], "total_carbon": float(sum(y[2:])),
                "R0": reproduction_number(hydraulic_factor(leaf, h), c),
                "D_slope_per_MPa": slope, "max_real_eigenvalue_per_day": float(max(eig.real)),
            })
            spectra.append({"variant": name, "D_mol_mol": d, "leaf_MPa": leaf,
                            **{f"eig{j}_{part}": float(getattr(z, part)) for j, z in enumerate(eig) for part in ["real", "imag"]}})
        multistart = []
        for d in [.01, .015, .02, .03, .05]:
            # Separate monotone-interval brackets are distinct root searches.
            roots = equilibria_at_vpd(d, bounds, points, root_potential, h, c, feedback)
            independent, attempts = multistart_roots(d, bounds, root_potential, h, c, feedback)
            multistart.extend(attempts)
            require(len(independent) == len(roots), f"{name}, D={d}: multistart root count disagrees")
            require(all(abs(a[1]-b[1]) < 1e-7 for a, b in zip(independent, roots)),
                    "multistart full-system roots disagree with reduced curve")
            for y in roots:
                root_rows.append({"variant": name, "D_mol_mol": d, "leaf_MPa": y[1],
                                  "total_carbon": float(sum(y[2:])),
                                  "max_real_eigenvalue_per_day": float(max(np.linalg.eigvals(two_way_jacobian(
                                      y, root_potential, d, h, c, feedback,
                                  )).real))})
        report["branches"][name] = {
            "feedback": asdict(feedback), "folds": [fold_checks(z, root_potential, h, c, feedback) for z in points],
            "stability_boundaries": boundaries, "continuation_points": [len(coarse), len(fine)],
            "continuation_max_D_refinement_error": arc_error,
            "continuation_bounds_MPa": [-2.12, bounds[1]],
            "root_search_bounds_MPa": list(bounds),
            "full_six_state_multistart_attempts": multistart,
        }
    require(max_residual < 1e-10, "full equilibrium residual failed")
    write_csv(args.output/"branches.csv", branch_rows)
    write_csv(args.output/"spectrum.csv", spectra)
    write_csv(args.output/"equilibria.csv", root_rows)
    write_csv(args.output/"continuation.csv", continuation_rows)

    trajectory_rows = []
    max_budget = max_water = max_trajectory_refinement = 0.
    trajectories = {}
    for name in ["E0", "E1", "E2", "E3"]:
        feedback = variants[name]
        for duration in [10., 10*YEAR_DAYS]:
            def forcing(t):
                s = min(max(t/duration, 0.), 1.)
                return .01+.02*(3*s*s-2*s*s*s)
            times = np.unique(np.r_[np.linspace(0., duration, 201),
                                    duration+np.linspace(0., 500*YEAR_DAYS, 801)])
            solutions = [simulate_two_way(initial, times, forcing, root_potential, h, c, feedback, rtol=tol)
                         for tol in [1e-8, 1e-10]]
            coarse, fine = solutions
            require(coarse["status"] == fine["status"], "tolerance changed vegetation exit classification")
            common = min(len(coarse["time_days"]), len(fine["time_days"]))-1
            error = float(np.max(abs(coarse["state"][:common]-fine["state"][:common])))
            end_error = abs(coarse["time_days"][-1]-fine["time_days"][-1])
            require(error < 2e-5 and end_error < .01, "trajectory/event refinement failed")
            max_trajectory_refinement = max(max_trajectory_refinement, error)
            final = fine["state"][-1]
            final_roots = equilibria_at_vpd(.03, bounds, stationary[name], root_potential, h, c, feedback)
            attracting_roots = [y for y in final_roots if max(np.linalg.eigvals(
                two_way_jacobian(y, root_potential, .03, h, c, feedback),
            ).real) < 0]
            relative_error = (float(min(np.max(abs(final-y)/np.maximum(abs(y), 1e-6)) for y in attracting_roots))
                              if fine["status"] == "completed" and attracting_roots else None)
            if fine["status"] == "completed":
                require(relative_error is not None and relative_error < 1e-5, "500-year final hold did not verify attracting equilibrium")
            frozen_path = []
            for d in np.linspace(.01, .03, 201):
                roots = equilibria_at_vpd(d, bounds, stationary[name], root_potential, h, c, feedback)
                frozen_path.append(sum(max(np.linalg.eigvals(two_way_jacobian(
                    y, root_potential, d, h, c, feedback,
                )).real) < 0 for y in roots))
            has_fold = any(.01 < fold["D"] < .03 for fold in report["branches"][name]["folds"])
            classification = ("frozen_fold_crossing_control" if has_fold else
                              "rate_comparison_on_numerically_persistent_attracting_branch"
                              if min(frozen_path) >= 1 else "frozen_stability_not_verified")
            item = {
                "variant": name, "ramp_days": duration, "path": "D=.01 -> .03 smoothstep, final hold 500 years",
                "status": fine["status"], "classification": classification,
                "min_attracting_equilibrium_count_on_201_path_nodes": int(min(frozen_path)),
                "last_day": float(fine["time_days"][-1]), "final_state": final.tolist(),
                "final_equilibrium_relative_error": relative_error,
                "refinement_max_state_error": error, "refinement_end_day_error": end_error,
                "claim": "vegetation-domain exit versus positive-equilibrium recovery; NOT a verified R-tipping claim",
            }
            report["trajectories"].append(item)
            for t, y in zip(fine["time_days"], fine["state"]):
                stress = hydraulic_factor(y[1], h)
                capacity = frozen_coefficient_capacity(y[2:], stress, c)
                budget = carbon_budget(y[2:], stress, c)
                erx, exl, e = structural_fluxes(y, root_potential, forcing(t), h, feedback)
                rhs = two_way_rhs(y, root_potential, forcing(t), h, c, feedback)
                max_budget = max(max_budget, abs(budget["residual"]))
                max_water = max(max_water, abs(h.cx*rhs[0]+h.cl*rhs[1]-erx+e))
                trajectory_rows.append({
                    "variant": name, "ramp_days": duration, "day": t, "D_mol_mol": forcing(t),
                    "xylem_MPa": y[0], "leaf_MPa": y[1], "leaf_carbon": y[2],
                    "total_carbon": float(sum(y[2:])), "capacity_map_total": float(sum(capacity)),
                    "GPP": budget["gpp"], "NEP": budget["nep"],
                })
            trajectories[(name, duration)] = fine
    require(max_budget < 1e-14 and max_water < 1e-14, "trajectory mass balance failed")
    write_csv(args.output/"trajectories.csv", trajectory_rows)
    report["minimum_leaf_sensitivity"] = []
    for minimum in [1e-5, 1e-6, 1e-7]:
        f = replace(variants["E3"], minimum_leaf=minimum)
        def fast_forcing(t):
            s = min(max(t/10., 0.), 1.)
            return .01+.02*(3*s*s-2*s*s*s)
        result = simulate_two_way(initial, [0., 500*YEAR_DAYS], fast_forcing, root_potential, h, c, f)
        report["minimum_leaf_sensitivity"].append({
            "minimum_leaf": minimum, "status": result["status"],
            "last_day": float(result["time_days"][-1]),
            "final_total": float(sum(result["state"][-1, 2:])),
        })

    basin_rows = []
    stable_e1 = equilibria_at_vpd(.01, bounds, stationary["E1"], root_potential, h, c, variants["E1"])[0]
    for leaf in np.linspace(-2.4, -.6, 6):
        for scale in np.linspace(.2, 1.2, 6):
            y = initial.copy()
            y[1], y[2] = leaf, base.leaf_reference*scale
            fl = feedback_factors(y[2], y[3], variants["E1"])[0]
            y[0] = (root_potential+fl*leaf)/(1+fl)
            result = simulate_two_way(y, [0., 500*YEAR_DAYS], lambda t: .01, root_potential, h, c, variants["E1"])
            label = "vegetation_domain_exit"
            distance = None
            if result["status"] == "completed":
                distance = float(max(abs(result["state"][-1]-stable_e1)/np.maximum(abs(stable_e1), 1e-6)))
                require(distance < 1e-5, "basin hold did not reach verified wet equilibrium")
                label = "verified_wet_equilibrium"
            basin_rows.append({"initial_leaf_MPa": leaf, "initial_leaf_scale": scale,
                               "outcome": label, "last_day": float(result["time_days"][-1]),
                               "final_relative_distance": distance})
    report["basin"] = basin_rows
    write_csv(args.output/"basin.csv", basin_rows)

    sensitivity_rows = []
    cases = []
    for root_value in [0., -.2, -.5]:
        for strength in [0., .25, .5, .75, 1.]:
            for name in ["E1", "E2", "E3"]:
                f = variants[name]
                cases.append((f"{name}_strength_{strength:g}", root_value,
                              replace(f, leaf_supply=f.leaf_supply*strength, leaf_demand=f.leaf_demand*strength)))
    for floor in [.001, .01, .05]:
        cases.append((f"E3_floor_{floor:g}", root_potential, replace(variants["E3"], leaf_floor=floor)))
    for factor in [.75, 1., 1.25]:
        cases.append((f"E3_fixed_reference_x{factor:g}", root_potential,
                      replace(variants["E3"], leaf_reference=factor*base.leaf_reference,
                              wood_reference=factor*base.wood_reference)))
    for exponent in [.5, 1., 2.]:
        cases.append((f"E4_exponent_{exponent:g}", root_potential, replace(variants["E4"], wood_exponent=exponent)))
    for label, root_value, f in cases:
        local_bounds = (bounds[0], root_value-.001)
        points = stationary_points(local_bounds, root_value, h, c, f)
        boundaries = stability_boundaries(local_bounds, root_value, h, c, f)
        row = {
            "case": label, "root_MPa": root_value,
            "leaf_supply_strength": f.leaf_supply, "leaf_demand_strength": f.leaf_demand,
            "wood_supply_strength": f.wood_supply, "floor": f.leaf_floor,
            "leaf_ref": f.leaf_reference, "wood_ref": f.wood_reference, "wood_exponent": f.wood_exponent,
            "fold_D": json.dumps([equilibrium_curve(z, root_value, h, c, f)[1] for z in points]),
            "stability_boundaries": json.dumps(boundaries),
        }
        sensitivity_rows.append(row)
    report["sensitivities"] = sensitivity_rows
    write_csv(args.output/"sensitivities.csv", sensitivity_rows)

    # E2 invariant-box sufficient condition: no hydraulic/leaf-domain exit
    # for any D(t) in the stated interval, independent of forcing rate.
    leaf_barrier = -1.8
    upper_leaf = positive_carbon_equilibrium(hydraulic_factor(root_potential, h), c)[0]
    lower_leaf = positive_carbon_equilibrium(hydraulic_factor(leaf_barrier, h), c)[0]
    gx, gs = conductances(leaf_barrier, h)
    safe_d = gx*(root_potential-leaf_barrier)/(gs*(upper_leaf/base.leaf_reference))
    require(lower_leaf > base.minimum_leaf and initial[2] >= lower_leaf and initial[2] <= upper_leaf, "invariant carbon box invalid")
    require(initial[1] >= leaf_barrier and initial[0] >= (root_potential+leaf_barrier)/2, "initial state outside hydraulic box")
    report["E2_rate_independent_non_exit_sufficient_condition"] = {
        "psi_leaf_lower": leaf_barrier, "psi_x_lower": (root_potential+leaf_barrier)/2,
        "leaf_carbon_interval": [float(lower_leaf), float(upper_leaf)],
        "maximum_D": float(safe_d),
        "statement": "Any nonnegative D(t)<=maximum_D preserves this box; not a global convergence/R-tipping theorem",
    }
    report["checks"] = {
        "max_full_equilibrium_rhs_residual": max_residual,
        "max_pseudo_arclength_D_refinement_error": max_arc_refinement,
        "max_trajectory_state_refinement_error": max_trajectory_refinement,
        "max_carbon_budget_residual": max_budget, "max_water_budget_residual": max_water,
    }
    paths = [
        Path(__file__), args.config, ROOT/"src/control_carbon/spac_carbon_two_way.py",
        ROOT/"tests/test_spac_carbon_two_way.py", ROOT/"src/control_carbon/spac_carbon.py",
        ROOT/"src/control_carbon/manzoni_hydraulics.py", ROOT/"src/control_carbon/manzoni_extensions.py",
    ]
    report["source_sha256"] = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                             for path in paths}
    report["versions"] = {"numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__}

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    leaves = np.linspace(.001, .5, 150)
    factors = np.array([feedback_factors(z, base.wood_reference*z/base.leaf_reference, variants["E4"])[:3] for z in leaves])
    for j, label in enumerate(["FL: E1/E3/E4", "FE: E2/E3/E4", "FW: E4 (equilibrium W/L)"]):
        axes[0, 0].plot(leaves, factors[:, j], label=label)
    axes[0, 0].axvline(base.leaf_reference, color="gray", ls=":")
    axes[0, 0].set(xlabel="Leaf carbon (kg C/m2 ground)", ylabel="Dimensionless F")
    sla_factor = c.specific_leaf_area/c.carbon_dry_fraction
    secondary = axes[0, 0].secondary_xaxis("top", functions=(lambda z: sla_factor*z, lambda z: z/sla_factor))
    secondary.set_xlabel("LAI (m2 leaf/m2 ground)")
    for name in variants:
        rows = [r for r in branch_rows if r["variant"] == name and r["D_mol_mol"] <= .065]
        ds = [r["D_mol_mol"] for r in rows]
        stable = np.array([r["max_real_eigenvalue_per_day"] < 0 for r in rows])
        for ax, key in [(axes[0, 1], "leaf_MPa"), (axes[0, 2], "total_carbon")]:
            values = np.array([r[key] for r in rows])
            line, = ax.plot(ds, np.where(stable, values, np.nan), label=name)
            ax.plot(ds, np.where(stable, np.nan, values), ls="--", color=line.get_color())
        axes[1, 0].plot(ds, [YEAR_DAYS*r["max_real_eigenvalue_per_day"] for r in rows], label=name)
    axes[0, 1].set(ylabel="Equilibrium leaf potential (MPa)")
    axes[0, 2].set(ylabel="Equilibrium total carbon (kg C/m2)")
    axes[1, 0].set(ylabel="Full Jacobian max Re(lambda) (1/year)", yscale="symlog")
    axes[1, 0].set_yscale("symlog", linthresh=.03)
    axes[1, 0].axhline(0, color="gray", ls=":")
    for ax in [axes[0, 1], axes[0, 2], axes[1, 0]]:
        ax.set_xlabel("D (mol/mol)")
    for row in basin_rows:
        color = "green" if row["outcome"] == "verified_wet_equilibrium" else "orange"
        axes[1, 1].scatter(row["initial_leaf_MPa"], row["initial_leaf_scale"], color=color, marker="s", s=150)
    axes[1, 1].set(xlabel="Initial leaf potential (MPa)", ylabel="Initial L/L_ref",
                   title="E1, D=.01: green=wet; orange=domain exit")
    for duration in [10., 10*YEAR_DAYS]:
        rows = [r for r in trajectory_rows if r["variant"] == "E3" and r["ramp_days"] == duration and r["day"]/YEAR_DAYS <= 30]
        years = [r["day"]/YEAR_DAYS for r in rows]
        line, = axes[1, 2].plot(years, [r["total_carbon"] for r in rows], label=f"stock, ramp {duration/YEAR_DAYS:.3g} yr")
        axes[1, 2].plot(years, [r["capacity_map_total"] for r in rows], ls="--", color=line.get_color(), label="capacity map")
    axes[1, 2].set(xlabel="Time (year)", ylabel="E3 carbon (kg C/m2)")
    for ax in axes.flat:
        if ax != axes[1, 1]:
            ax.legend(fontsize=7)
        ax.grid(alpha=.2)
    fig.suptitle("Stage B1 | synthetic | fixed capacities/references | solid=stable, dashed=unstable | NOT mortality")
    fig.tight_layout()
    for suffix in ["png", "pdf"]:
        fig.savefig(args.output/f"stage_b1.{suffix}", dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    for ax, name in zip(axes.flat, variants):
        rows = [r for r in spectra if r["variant"] == name and r["D_mol_mol"] <= .065]
        ds = [r["D_mol_mol"] for r in rows]
        for j in range(6):
            ax.plot(ds, [YEAR_DAYS*r[f"eig{j}_real"] for r in rows], lw=.7, label=f"Re eig {j+1}")
        ax.set(title=name, xlabel="D (mol/mol)", ylabel="Re(lambda) (1/year)", yscale="symlog")
        ax.axhline(0, color="black", ls=":")
        ax.grid(alpha=.2)
    axes[1, 2].axis("off")
    axes[1, 2].text(0, .5, "Full 6x6 spectra (CSV also stores Im).\nIndices are not continuous eigenmode labels.\nOscillatory crossings are candidates only.\nAll structural/carbon coefficients synthetic.")
    fig.suptitle("Stage B1 full spectra | no fold inference from hydraulic block alone")
    fig.tight_layout()
    for suffix in ["png", "pdf"]:
        fig.savefig(args.output/f"spectrum.{suffix}", dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, name in zip(axes.flat, ["E0", "E1", "E2", "E3"]):
        for duration in [10., 10*YEAR_DAYS]:
            solution = trajectories[(name, duration)]
            ax.plot(solution["time_days"]/YEAR_DAYS, np.sum(solution["state"][:, 2:], axis=1),
                    label=f"{duration/YEAR_DAYS:.3g} yr ramp: {solution['status']}")
            if solution["status"] == "vegetation_domain_exit":
                ax.scatter(solution["time_days"][-1]/YEAR_DAYS, sum(solution["state"][-1, 2:]), marker="x")
        ax.set(title=name, xlabel="Time (year), planned final hold 500 years", ylabel="Total carbon (kg C/m2)")
        ax.legend(fontsize=7)
        ax.grid(alpha=.2)
    fig.suptitle("Same D=.01 -> .03 path | speed only | exit is not a death/alternate-attractor claim")
    fig.tight_layout()
    for suffix in ["png", "pdf"]:
        fig.savefig(args.output/f"rate_comparison.{suffix}", dpi=180)
    plt.close(fig)
    (args.output/"results.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps({
        "fixed_reference": report["fixed_reference"],
        "branches": {k: {"folds": [f["D"] for f in v["folds"]], "stability": v["stability_boundaries"]}
                     for k, v in report["branches"].items()},
        "trajectories": [{"variant": r["variant"], "ramp_days": r["ramp_days"], "status": r["status"],
                          "last_year": r["last_day"]/YEAR_DAYS, "final_total": sum(r["final_state"][2:])}
                         for r in report["trajectories"]],
        "basin_wet_count": sum(r["outcome"] == "verified_wet_equilibrium" for r in basin_rows),
        "E2_non_exit_condition": report["E2_rate_independent_non_exit_sufficient_condition"],
        "checks": report["checks"],
    }, indent=2))


if __name__ == "__main__":
    main()
