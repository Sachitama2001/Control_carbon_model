"""Bounded E4 local Hopf/cycle analysis, MPa forcing and synthetic sensitivities."""

import argparse
import csv
from dataclasses import replace
from datetime import datetime
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from control_carbon.manzoni_hydraulics import HydraulicParameters, transpiring_equilibria
from control_carbon.spac_carbon import YEAR_DAYS, CarbonParameters, positive_carbon_equilibrium, hydraulic_factor, carbon_budget
from control_carbon.spac_carbon_two_way import StructuralFeedback
from control_carbon.spac_carbon_mpa import PressureSPAC
from control_carbon.spac_carbon_e4 import core_jacobian, core_tensors, hopf_point, periodic_orbit
from control_carbon.spac_carbon_bifurcation import first_lyapunov

ROOT = Path(__file__).resolve().parents[1]


def encode(value):
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
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT/"configs/spac_carbon_pressure_review.json")
    parser.add_argument("--output", type=Path, default=ROOT/"artifacts/spac_carbon"/datetime.now().strftime("stage_b3_e4_%Y%m%d_%H%M%S"))
    args = parser.parse_args()
    settings = json.loads(args.config.read_text())
    require(settings["units"]["vpd"] == "MPa", "E4 runner requires MPa-facing settings")
    baseline_path = args.config.parent/settings["model_config"]
    cfg = json.loads(baseline_path.read_text())
    h = HydraulicParameters(**cfg["hydraulic"])
    k = cfg["carbon"]
    c = CarbonParameters(k["assimilation_max_per_year"]/YEAR_DAYS, tuple(k["allocation"]),
                         tuple(v/YEAR_DAYS for v in k["losses_per_year"]), tuple(k["litter_fractions"]),
                         k["specific_leaf_area"], k["carbon_dry_fraction"], k["extinction"])
    leaf = transpiring_equilibria(cfg["root_potential"], cfg["reference_vpd"], h)[0]
    stocks = positive_carbon_equilibrium(hydraulic_factor(leaf, h), c)
    require(stocks is not None, "Canonical fixed reference must have positive carbon")
    f = StructuralFeedback(stocks[0], stocks[1], cfg["leaf_floor"], cfg["minimum_leaf"], 1., 1., 1.,
                           wood_exponent=cfg["wood_exponent"])
    model = PressureSPAC(h, c, f, cfg["root_potential"], settings["atmospheric_pressure_mpa"])
    args.output.mkdir(parents=True, exist_ok=False)
    hopf = hopf_point((-2.01, -1.99), model)
    diagnostic = model.air_diagnostic(hopf["vpd_mpa"], settings["diagnostic_air_temperature_c"])
    require(diagnostic["admissible_at_diagnostic_temperature"], "Baseline E4 candidate outside diagnostic climate range")
    fd = []
    for step in [1e-4, 5e-5]:
        y, v = hopf["state"][:4], hopf["vpd_mpa"]
        b = np.stack([(core_jacobian(y+step*e, v, model)-core_jacobian(y-step*e, v, model))/(2*step)
                      for e in np.eye(4)], axis=2)
        third = np.stack([(core_tensors(y+step*e, v, model)[0]-core_tensors(y-step*e, v, model)[0])/(2*step)
                          for e in np.eye(4)], axis=3)
        fd.append(first_lyapunov(core_jacobian(y, v, model), b, third)["l1"])
    require(max(abs(value/hopf["l1"]-1) for value in fd) < 1e-5, "E4 l1 finite difference refinement failed")
    cycles, summary, states = [], [], []
    max_water = max_carbon = 0.
    for amplitude in [.0005, .001, .002]:
        cycle = periodic_orbit(hopf, amplitude, model, tol=1e-9)
        cycles.append(cycle)
        neutral = np.argmin(abs(cycle["core_floquet"]-1))
        radial = max(abs(np.delete(cycle["core_floquet"], neutral)))
        require(radial > 1 and cycle["vpd_mpa"] < hopf["vpd_mpa"], "Cycle does not match subcritical prediction")
        summary.append({"amplitude": amplitude, "vpd_mpa": cycle["vpd_mpa"], "period_days": cycle["period_days"],
                        "radial_multiplier": float(radial), "neutral_error": float(min(abs(cycle["core_floquet"]-1))),
                        "closure": cycle["full_state_closure"], "flow_error": cycle["independent_flow_max_state_error"]})
        for t, y in zip(cycle["phase"]*cycle["period_days"], cycle["state"]):
            rhs = model.rhs(y, cycle["vpd_mpa"])
            erx, exl, e = model.fluxes(y, cycle["vpd_mpa"])
            max_water = max(max_water, abs(h.cx*rhs[0]+h.cl*rhs[1]-erx+e))
            max_carbon = max(max_carbon, abs(sum(rhs[2:])-carbon_budget(y[2:], hydraulic_factor(y[1], h), c)["nep"]))
            states.append({"amplitude": amplitude, "day": t,
                           **dict(zip(["xylem_mpa", "leaf_mpa", "L", "W", "R", "S"], y))})
        print(f"E4 r={amplitude}: VPD={cycle['vpd_mpa']:.12g} MPa, T={cycle['period_days']:.6g} d, multiplier={radial:.9g}", flush=True)
    refined = periodic_orbit(hopf, .001, model, tol=1e-10)
    old = cycles[1]
    refinement = {"vpd_mpa_error": abs(refined["vpd_mpa"]-old["vpd_mpa"]),
                  "period_days_error": abs(refined["period_days"]-old["period_days"]),
                  "max_phase_aligned_state_error": float(max(abs(refined["state"]-old["state"]).ravel())),
                  "radial_multiplier_error": float(abs(max(abs(refined["core_floquet"]))-max(abs(old["core_floquet"])))),
                  "refined_flow_error": refined["independent_flow_max_state_error"]}
    require(refinement["vpd_mpa_error"] < 1e-10 and refinement["period_days_error"] < 1e-3
            and refinement["max_phase_aligned_state_error"] < 1e-6 and refinement["radial_multiplier_error"] < 1e-6,
            "E4 periodic tolerance refinement failed")
    cases = [("baseline", model)]
    cases += [(f"capacity_x{factor:g}", replace(model, hydraulic=replace(h, cx=h.cx*factor, cl=h.cl*factor)))
              for factor in [.5, 2.]]
    cases += [(f"wood_exponent_{beta:g}", replace(model, feedback=replace(f, wood_exponent=beta)))
              for beta in [.75, 1.25]]
    for factor in [.5, 2.]:
        losses = list(c.losses)
        losses[1] *= factor
        cases.append((f"wood_loss_x{factor:g}", replace(model, carbon=replace(c, losses=tuple(losses)))))
    sensitivities = []
    for label, variant in cases:
        def pair_real(psi):
            state, vpd, _ = variant.equilibrium_curve(psi)
            eigen = np.linalg.eigvals(core_jacobian(state[:4], vpd, variant))
            require(max(eigen.imag) > 1e-8, "Sensitivity bracket lost its complex pair")
            return float(eigen[np.argmax(eigen.imag)].real)
        psi = brentq(pair_real, -2.01, -1.99, xtol=1e-13)
        state, vpd, _ = variant.equilibrium_curve(psi)
        eigen = np.linalg.eigvals(core_jacobian(state[:4], vpd, variant))
        other_unstable = bool(np.any(eigen[abs(eigen.imag) < 1e-8].real >= 0))
        candidate = None if other_unstable else hopf_point((-2.01, -1.99), variant)
        sensitivity_folds = [variant.equilibrium_curve(z)[1] for z in variant.stationary_points(cfg["curve_bounds"])]
        sensitivities.append({
            "case": label, "vpd_mpa": vpd,
            "status": "pair_crossing_on_already_unstable_branch" if other_unstable else "nondegenerate_attracting_side_Hopf",
            "l1": None if candidate is None else candidate["l1"],
            "omega_per_day": float(max(eigen.imag)),
            "slope_per_mpa": None if candidate is None else candidate["slope_refinement_per_mpa"][-1],
            "noncritical_real_eigenvalues": json.dumps(eigen[abs(eigen.imag) < 1e-8].real.tolist()),
            "detected_fold_vpd_mpa": json.dumps(sensitivity_folds),
            "in_routine_range": bool(settings["routine_vpd_mpa_bounds"][0] <= vpd <= settings["routine_vpd_mpa_bounds"][1]),
            **variant.air_diagnostic(vpd, settings["diagnostic_air_temperature_c"]),
        })
    require(max_water < 1e-14 and max_carbon < 1e-14, "E4 cycle budget failed")
    write_csv(args.output/"cycles.csv", summary)
    write_csv(args.output/"periodic_states.csv", states)
    write_csv(args.output/"sensitivities.csv", sensitivities)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for cycle in cycles:
        delta = cycle["state"]-hopf["state"]
        axes[0].plot(delta[:, 1], 1e6*delta[:, 2], label=f"r={cycle['amplitude']:g}")
        axes[1].plot(cycle["phase"], 1e6*delta[:, 3], label=f"r={cycle['amplitude']:g}")
    axes[0].set(xlabel="Leaf potential - Hopf equilibrium (MPa)", ylabel="Leaf C deviation (1e-6 kg C/m2)")
    axes[1].set(xlabel="Time / period", ylabel="Wood C deviation (1e-6 kg C/m2)")
    bars = axes[2].barh([r["case"] for r in sensitivities], [r["vpd_mpa"] for r in sensitivities])
    for bar, row in zip(bars, sensitivities):
        if row["status"] == "pair_crossing_on_already_unstable_branch":
            bar.set_color("darkorange")
            bar.set_hatch("//")
            axes[2].plot(json.loads(row["detected_fold_vpd_mpa"]), [bar.get_y()+bar.get_height()/2],
                         "kD", label="fold; orange bar is already unstable")
    axes[2].axvline(.003, color="black", ls=":", label="routine upper bound")
    axes[2].axvline(diagnostic["saturation_vapor_pressure_mpa"], color="gray", ls="--", label="saturation at 25 C")
    axes[2].set(xlabel="Complex-pair crossing VPD (MPa)")
    for ax in axes:
        ax.legend(fontsize=7)
        ax.grid(alpha=.2)
    fig.suptitle("E4: numerically verified local subcritical Hopf | synthetic coefficients, fixed references | NOT mortality/R-tipping")
    fig.tight_layout()
    for suffix in ["png", "pdf"]:
        fig.savefig(args.output/f"e4_nonlinear.{suffix}", dpi=180)
    plt.close(fig)
    report = {"classification": "model-agnostic synthetic theory; local numerical verification, not interval-certified Hopf",
              "settings": settings, "baseline_config": cfg, "reference_stocks": stocks, "hopf": hopf,
              "air_diagnostic": diagnostic, "finite_difference_l1": fd,
              "cycles": [{key: value for key, value in cycle.items() if key not in ["state", "phase"]} for cycle in cycles],
              "refinement": refinement, "sensitivities": sensitivities,
              "sensitivity_rule": "7 predefined one-at-a-time cases; fixed original references, no recalibration",
              "budgets": {"water": max_water, "carbon": max_carbon}}
    sources = [Path(__file__), args.config, baseline_path, ROOT/"src/control_carbon/spac_carbon_e4.py",
               ROOT/"src/control_carbon/spac_carbon_bifurcation.py", ROOT/"src/control_carbon/spac_carbon_mpa.py",
               ROOT/"src/control_carbon/spac_carbon_two_way.py", ROOT/"src/control_carbon/spac_carbon.py",
               ROOT/"src/control_carbon/manzoni_hydraulics.py", ROOT/"tests/test_spac_carbon_e4.py"]
    report["source_sha256"] = {str(path.resolve().relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                             for path in sources}
    (args.output/"results.json").write_text(json.dumps(report, default=encode, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"hopf_vpd_mpa": hopf["vpd_mpa"], "l1": hopf["l1"],
                      "refinement": refinement, "sensitivities": sensitivities}, indent=2))


if __name__ == "__main__":
    main()
