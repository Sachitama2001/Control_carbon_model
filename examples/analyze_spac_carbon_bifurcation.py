"""Stage B2 continuous-path certificates and local subcritical cycles."""

import argparse
import csv
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

from control_carbon.manzoni_hydraulics import HydraulicParameters, transpiring_equilibria
from control_carbon.spac_carbon import (
    YEAR_DAYS, CarbonParameters, hydraulic_factor, positive_carbon_equilibrium, carbon_budget,
)
from control_carbon.spac_carbon_two_way import (
    StructuralFeedback, equilibrium_curve, equilibria_at_vpd, two_way_rhs,
    two_way_jacobian, structural_fluxes,
)
from control_carbon.spac_carbon_bifurcation import (
    core_jacobian, core_tensors, cubic_hurwitz, first_lyapunov, hopf_point, periodic_orbit,
)
from control_carbon.spac_carbon_stability import certify_stable_path

ROOT = Path(__file__).resolve().parents[1]


def encode(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, complex):
        return [value.real, value.imag]
    raise TypeError(f"unserializable result type: {type(value)}")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT/"configs/spac_carbon_stage_b_synthetic.json")
    parser.add_argument("--output", type=Path, default=ROOT/"artifacts/spac_carbon"/datetime.now().strftime("stage_b2_%Y%m%d_%H%M%S"))
    args = parser.parse_args()
    args.config = args.config.resolve()
    cfg = json.loads(args.config.read_text())
    args.output.mkdir(parents=True, exist_ok=False)
    h = HydraulicParameters(**cfg["hydraulic"])
    k = cfg["carbon"]
    c = CarbonParameters(k["assimilation_max_per_year"]/YEAR_DAYS, tuple(k["allocation"]),
                         tuple(v/YEAR_DAYS for v in k["losses_per_year"]), tuple(k["litter_fractions"]),
                         k["specific_leaf_area"], k["carbon_dry_fraction"], k["extinction"])
    root = cfg["root_potential"]
    leaf0 = transpiring_equilibria(root, cfg["reference_vpd"], h)[0]
    stocks = positive_carbon_equilibrium(hydraulic_factor(leaf0, h), c)
    require(stocks is not None, "B1 fixed reference is not positive")
    report = {
        "classification": "synthetic B2 mathematical analysis, unchanged B1 field; not mortality or R-tipping",
        "git_HEAD": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "config": cfg, "reference_stocks": stocks, "models": {},
        "normal_form_convention": "l1=Re(G21)/(2 omega), ||q||2=1, <p,q>=1, original mixed-unit coordinates and day time",
        "complex_encoding": "[real,imaginary]",
    }
    path_rows, cycle_rows, state_rows = [], [], []
    cycles = {}
    max_water = max_carbon = 0.
    for name, supply, leaf_bounds in [("E2", 0., (.19, .30)), ("E3", 1., (.145, .30))]:
        feedback = StructuralFeedback(
            stocks[0], stocks[1], cfg["leaf_floor"], cfg["minimum_leaf"],
            leaf_supply=supply, leaf_demand=1.,
        )
        certificate = certify_stable_path(leaf_bounds, (.01, .03), root, h, c, feedback)
        (args.output/f"certificate_{name}.json").write_text(json.dumps(certificate, indent=2)+"\n")
        hopf = hopf_point((-2.1, -1.9), root, h, c, feedback)
        fd_l1 = []
        for step in [1e-4, 5e-5]:
            y, d = hopf["state"][:3], hopf["D"]
            b = np.stack([(core_jacobian(y+step*e, root, d, h, c, feedback)
                           -core_jacobian(y-step*e, root, d, h, c, feedback))/(2*step)
                          for e in np.eye(3)], axis=2)
            third = np.stack([(core_tensors(y+step*e, root, d, h, c, feedback)[0]
                               -core_tensors(y-step*e, root, d, h, c, feedback)[0])/(2*step)
                              for e in np.eye(3)], axis=3)
            value = first_lyapunov(core_jacobian(y, root, d, h, c, feedback), b, third)["l1"]
            fd_l1.append(value)
        require(max(abs(v/hopf["l1"]-1) for v in fd_l1) < 1e-5, "Lyapunov derivative refinement failed")
        for d in np.linspace(.01, .03, 101):
            y = equilibria_at_vpd(d, tuple(cfg["curve_bounds"]), (), root, h, c, feedback)[0]
            coefficients = cubic_hurwitz(core_jacobian(y[:3], root, d, h, c, feedback))
            full_eigen = np.linalg.eigvals(two_way_jacobian(y, root, d, h, c, feedback))
            path_rows.append({"variant": name, "D": d, "leaf_MPa": y[1], "leaf_carbon": y[2],
                              **dict(zip(["a1", "a2", "a3", "delta"], coefficients)),
                              "full_max_real_eigenvalue": float(max(full_eigen.real))})
        item = {
            "continuous_certificate_file": f"certificate_{name}.json",
            "certificate_cells": certificate["cell_count"],
            "certificate_strict_lower_bounds": certificate["strict_lower_bounds"],
            "certificate_Dprime_L_upper_bound": certificate["Dprime_L_upper_bound"],
            "certificate_endpoint_D": certificate["endpoint_D"],
            "hopf": hopf, "finite_difference_l1_refinement": fd_l1,
            "local_classification": "numerically verified subcritical Hopf and repelling cycles",
            "cycles": [], "refinement": {},
        }
        for amplitude in [.0005, .001, .002]:
            cycle = periodic_orbit(hopf, amplitude, root, h, c, feedback,
                                   leaf_bracket=(-2.1, -1.9), rtol=1e-10, steps_per_period=100)
            require(cycle["D"] < hopf["D"] and cycle["radial_floquet"].real > 1,
                    "computed cycle does not verify the predicted subcritical side/stability")
            state_at_cycle = equilibria_at_vpd(cycle["D"], tuple(cfg["curve_bounds"]), (), root, h, c, feedback)[0]
            require(max(np.linalg.eigvals(two_way_jacobian(state_at_cycle, root, cycle["D"], h, c, feedback)).real) < 0,
                    "cycle is not on the attracting-equilibrium side")
            summary = {key: value for key, value in cycle.items() if key not in ["state", "time_days"]}
            summary["D_shift_ratio_to_normal_form"] = (cycle["D"]-hopf["D"])/cycle["predicted_D_shift"]
            require(abs(summary["D_shift_ratio_to_normal_form"]-1) < .03, "small-cycle branch disagrees with normal form")
            item["cycles"].append(summary)
            cycle_rows.append({"variant": name, "amplitude": amplitude, "D": cycle["D"],
                               "D_H_minus_D": hopf["D"]-cycle["D"], "period_days": cycle["period_days"],
                               "radial_multiplier": float(cycle["radial_floquet"].real),
                               "neutral_error": float(min(abs(cycle["core_floquet"]-1))),
                               "closure": cycle["full_state_closure"],
                               "leaf_potential_half_range": float(np.ptp(cycle["state"][:, 1])/2)})
            for t, y in zip(cycle["time_days"], cycle["state"]):
                f = two_way_rhs(y, root, cycle["D"], h, c, feedback)
                erx, exl, eatm = structural_fluxes(y, root, cycle["D"], h, feedback)
                budget = carbon_budget(y[2:], hydraulic_factor(y[1], h), c)
                max_water = max(max_water, abs(h.cx*f[0]+h.cl*f[1]-erx+eatm))
                max_carbon = max(max_carbon, abs(sum(f[2:])-budget["nep"]))
                state_rows.append({"variant": name, "amplitude": amplitude, "day": t,
                                   **dict(zip(["xylem_MPa", "leaf_MPa", "L", "W", "R", "S"], y))})
            cycles[(name, amplitude)] = cycle
            print(f"{name}: amplitude={amplitude:g}, D={cycle['D']:.12g}, "
                  f"T={cycle['period_days']:.6g} day, radial Floquet={cycle['radial_floquet'].real:.9g}", flush=True)
        refined = periodic_orbit(hopf, .001, root, h, c, feedback,
                                 leaf_bracket=(-2.1, -1.9), rtol=1e-11, steps_per_period=160)
        original = cycles[(name, .001)]
        refinement = {
            "D_error": abs(refined["D"]-original["D"]),
            "period_days_error": abs(refined["period_days"]-original["period_days"]),
            "radial_multiplier_error": float(abs(refined["radial_floquet"]-original["radial_floquet"])),
            "phase_sample_max_state_error": float(np.max(abs(refined["state"]-original["state"]))),
            "refined_shooting_residual": refined["shooting_residual"],
            "refined_neutral_multiplier_error": float(min(abs(refined["core_floquet"]-1))),
        }
        require(refinement["D_error"] < 1e-9 and refinement["period_days_error"] < 5e-4
                and refinement["radial_multiplier_error"] < 1e-6
                and refinement["phase_sample_max_state_error"] < 1e-6, "cycle accuracy refinement failed")
        item["refinement"] = refinement
        report["models"][name] = item
    require(max_water < 1e-14 and max_carbon < 1e-14, "full periodic orbit mass balance failed")
    report["checks"] = {"max_water_budget_residual": max_water, "max_carbon_budget_residual": max_carbon}
    write_csv(args.output/"path_hurwitz.csv", path_rows)
    write_csv(args.output/"periodic_summary.csv", cycle_rows)
    write_csv(args.output/"periodic_states.csv", state_rows)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for ax, name in zip(axes[0], ["E2", "E3"]):
        rows = [r for r in path_rows if r["variant"] == name]
        for key in ["a1", "a2", "a3", "delta"]:
            ax.plot([r["D"] for r in rows], [r[key]/rows[0][key] for r in rows], label=key)
        ax.set(title=f"{name}: continuous stability certificate ({report['models'][name]['certificate_cells']} cells)",
               xlabel="D (mol/mol)", ylabel="Hurwitz margins / values at D=.01", yscale="log")
        ax.legend()
    for name in ["E2", "E3"]:
        rows = [r for r in cycle_rows if r["variant"] == name]
        axes[1, 0].plot([r["D_H_minus_D"] for r in rows], [r["leaf_potential_half_range"] for r in rows],
                        "o--", label=f"{name}: repelling cycles")
    axes[1, 0].set(xlabel="D_H - D (mol/mol)", ylabel="Leaf potential half-range (MPa)", xscale="log", yscale="log")
    axes[1, 0].legend()
    for amplitude in [.0005, .001, .002]:
        y = cycles[("E3", amplitude)]["state"]
        axes[1, 1].plot(y[:, 1]-np.mean(y[:, 1]), 1e6*(y[:, 2]-np.mean(y[:, 2])), label=f"r={amplitude:g}")
    axes[1, 1].set(title="E3: local cycle phase portraits (centered by sampled mean)",
                   xlabel="Leaf potential deviation (MPa)", ylabel="Leaf C deviation (1e-6 kg C/m2)")
    axes[1, 1].legend()
    for ax in axes.flat:
        ax.grid(alpha=.2)
    fig.suptitle("Stage B2 | unchanged synthetic B1 | stable frozen path != no-R | local cycles are NOT attractors")
    fig.tight_layout()
    for suffix in ["png", "pdf"]:
        fig.savefig(args.output/f"stage_b2.{suffix}", dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for j, name in enumerate(["E2", "E3"]):
        cycle = cycles[(name, .002)]
        phase = cycle["time_days"]/cycle["period_days"]
        y = cycle["state"]
        axes[0, j].plot(phase, y[:, 1])
        axes[0, j].set(title=f"{name}, r=.002: T={cycle['period_days']:.3f} day",
                       xlabel="Time / period", ylabel="Leaf potential (MPa)")
        for k, label in enumerate(["L", "W", "R", "S"]):
            axes[1, j].plot(phase, 1e6*(y[:, k+2]-np.mean(y[:, k+2])), label=label)
        axes[1, j].set(xlabel="Time / period", ylabel="C deviation from sampled mean (1e-6 kg C/m2)")
        axes[1, j].legend()
    for ax in axes.flat:
        ax.grid(alpha=.2)
    fig.suptitle("Stage B2 full six-state periodic lift | positive carbon, fixed capacities | repelling cycles")
    fig.tight_layout()
    for suffix in ["png", "pdf"]:
        fig.savefig(args.output/f"periodic_states.{suffix}", dpi=180)
    plt.close(fig)
    sources = [Path(__file__), args.config, ROOT/"src/control_carbon/spac_carbon_bifurcation.py",
               ROOT/"src/control_carbon/spac_carbon_stability.py", ROOT/"tests/test_spac_carbon_bifurcation.py",
               ROOT/"src/control_carbon/spac_carbon_two_way.py", ROOT/"src/control_carbon/spac_carbon.py",
               ROOT/"src/control_carbon/manzoni_hydraulics.py", ROOT/"docs/SPAC_carbon_stage_B2_plan.md"]
    report["source_sha256"] = {str(path.relative_to(ROOT) if path.is_relative_to(ROOT) else path): hashlib.sha256(path.read_bytes()).hexdigest()
                             for path in sources}
    report["versions"] = {"numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__}
    (args.output/"results.json").write_text(json.dumps(report, indent=2, default=encode, allow_nan=False)+"\n")
    print(json.dumps({
        name: {"certificate_cells": item["certificate_cells"], "Hopf_D": item["hopf"]["D"],
               "omega_per_day": item["hopf"]["omega"], "l1": item["hopf"]["l1"],
               "refinement": item["refinement"]} for name, item in report["models"].items()
    }, indent=2, default=encode))


if __name__ == "__main__":
    main()
