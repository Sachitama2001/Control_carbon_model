"""MPa-facing unit/scale audit; unchanged synthetic physiology."""

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
from control_carbon.spac_carbon import YEAR_DAYS, CarbonParameters, hydraulic_factor, positive_carbon_equilibrium
from control_carbon.spac_carbon_two_way import StructuralFeedback
from control_carbon.spac_carbon_mpa import PressureSPAC

ROOT = Path(__file__).resolve().parents[1]


def e4_local_audit(model, threshold, bounds):
    """Exact four-state block; R/S remain passive stable filters."""
    samples = []
    center = threshold["vpd_mpa"]
    step = 1e-8
    for offset in [-step, -step/2, 0., step/2, step]:
        vpd = center+offset
        states = model.equilibria(vpd, bounds, model.stationary_points(bounds))
        if len(states) != 1:
            raise RuntimeError("E4 local check requires a unique detected positive equilibrium")
        full = model.jacobian(states[0], vpd)
        if np.any(full[:4, 4:] != 0):
            raise RuntimeError("E4 four-state block is not closed")
        core = full[:4, :4]
        a1, a2, a3, a4 = np.poly(core)[1:]
        factored = np.polymul(np.poly(core), np.polymul([1., model.carbon.losses[2]],
                                                      [1., model.carbon.losses[3]]))
        np.testing.assert_allclose(np.poly(full), factored, rtol=1e-10, atol=1e-12)
        eigen = np.linalg.eigvals(core)
        pair = eigen[np.argmax(eigen.imag)]
        if pair.imag <= 0 or min(a1, a2, a3, a4, a1*a2-a3) <= 0:
            raise RuntimeError("E4 check did not find the expected oscillatory Hurwitz boundary")
        samples.append({
            "vpd_mpa": vpd, "offset_mpa": offset, "quartic_coefficients": [a1, a2, a3, a4],
            "delta2": a1*a2-a3, "delta3": a1*a2*a3-a3*a3-a1*a1*a4,
            "critical_pair_real_per_day": float(pair.real),
            "critical_pair_imaginary_per_day": float(pair.imag),
        })
    if not (samples[0]["delta3"] > 0 > samples[-1]["delta3"]
            and samples[0]["critical_pair_real_per_day"] < 0 < samples[-1]["critical_pair_real_per_day"]):
        raise RuntimeError("E4 local stability change was not resolved")
    slopes = [
        (samples[4]["critical_pair_real_per_day"]-samples[0]["critical_pair_real_per_day"])/(2*step),
        (samples[3]["critical_pair_real_per_day"]-samples[1]["critical_pair_real_per_day"])/step,
    ]
    if abs(slopes[0]-slopes[1]) > 1e-4*abs(slopes[1]) or slopes[1] <= 0:
        raise RuntimeError("E4 critical eigenvalue slope refinement failed")
    return {
        "active_states": ["xylem_mpa", "leaf_mpa", "leaf_carbon", "wood_carbon"],
        "passive_eigenvalues_per_day": [-model.carbon.losses[2], -model.carbon.losses[3]],
        "samples": samples,
        "critical_pair_real_slope_per_vpd_mpa_refinement": slopes,
        "claim": "numerically resolved local oscillatory stability change; exact block factorization",
        "not_claimed": "nondegenerate Hopf, Lyapunov coefficient, periodic orbit, global basin, or mortality",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT/"configs/spac_carbon_pressure_review.json")
    parser.add_argument("--output", type=Path, default=ROOT/"artifacts/spac_carbon"/datetime.now().strftime("pressure_review_%Y%m%d_%H%M%S"))
    args = parser.parse_args()
    settings = json.loads(args.config.read_text())
    if settings["units"]["vpd"] != "MPa" or settings["units"]["water_potential"] != "MPa":
        raise ValueError("pressure-review config requires explicitly declared MPa inputs")
    baseline_path = args.config.parent/settings["model_config"]
    cfg = json.loads(baseline_path.read_text())
    h = HydraulicParameters(**cfg["hydraulic"])
    k = cfg["carbon"]
    c = CarbonParameters(k["assimilation_max_per_year"]/YEAR_DAYS, tuple(k["allocation"]),
                         tuple(v/YEAR_DAYS for v in k["losses_per_year"]), tuple(k["litter_fractions"]),
                         k["specific_leaf_area"], k["carbon_dry_fraction"], k["extinction"])
    p = settings["atmospheric_pressure_mpa"]
    # Bootstrap the fixed stocks with the archived Stage A reference, not a new climate calibration.
    leaf = transpiring_equilibria(cfg["root_potential"], cfg["reference_vpd"], h)[0]
    stocks = positive_carbon_equilibrium(hydraulic_factor(leaf, h), c)
    if stocks is None:
        raise ValueError("archived reference must have positive carbon stocks")
    base = StructuralFeedback(stocks[0], stocks[1], cfg["leaf_floor"], cfg["minimum_leaf"],
                              wood_exponent=cfg["wood_exponent"])
    models = {name: PressureSPAC(
        h, c, replace(base, leaf_supply=values[0], leaf_demand=values[1], wood_supply=values[2]),
        cfg["root_potential"], p,
    ) for name, values in cfg["variants"].items()}
    reference = models["E0"]
    if not np.isclose(reference.mole_fraction(settings["reference_vpd_mpa"]), cfg["reference_vpd"], rtol=1e-14, atol=0):
        raise ValueError("MPa reference must match the archived fixed-stock reference")
    temperature = settings["diagnostic_air_temperature_c"]
    routine = settings["routine_vpd_mpa_bounds"]
    for vpd in routine:
        if not reference.air_diagnostic(vpd, temperature)["admissible_at_diagnostic_temperature"]:
            raise ValueError("routine range is not admissible at the diagnostic temperature/pressure")
    if not 0 <= routine[0] < routine[1]:
        raise ValueError("routine VPD bounds must be ordered")
    formal = settings["formal_vpd_mpa_bounds"]
    if not 0 <= formal[0] < formal[1]:
        raise ValueError("formal VPD bounds must be ordered")
    rows, thresholds, curves = [], [], {}
    bounds = cfg["curve_bounds"]
    for name, model in models.items():
        records = []
        for leaf in np.linspace(*bounds, 801):
            y, vpd, slope = model.equilibrium_curve(leaf)
            eig = np.linalg.eigvals(model.jacobian(y, vpd))
            alpha = float(max(eig.real))
            row = {"variant": name, "vpd_mpa": vpd, "xylem_mpa": y[0], "leaf_mpa": y[1],
                   "leaf_carbon_kgC_m2": y[2], "wood_carbon_kgC_m2": y[3],
                   "root_carbon_kgC_m2": y[4], "soil_carbon_kgC_m2": y[5],
                   "vpd_slope_per_leaf_mpa": slope, "spectral_abscissa_per_day": alpha,
                   **model.air_diagnostic(vpd, temperature)}
            rows.append(row)
            records.append(row)
        curves[name] = records
        candidates = [(z, "fold") for z in model.stationary_points(bounds)]
        if name in ["E2", "E3"]:
            hopf = model.hopf_point((-2.1, -1.9))
            candidates.append((hopf["state"][1], "subcritical_Hopf"))
        if name == "E4":
            def alpha_at(leaf):
                y, vpd, _ = model.equilibrium_curve(leaf)
                return max(np.linalg.eigvals(model.jacobian(y, vpd)).real)
            for left, right in zip(records[:-1], records[1:]):
                if left["spectral_abscissa_per_day"]*right["spectral_abscissa_per_day"] < 0:
                    z = brentq(alpha_at, left["leaf_mpa"], right["leaf_mpa"], xtol=1e-12)
                    candidates.append((z, "oscillatory_candidate_not_Hopf_certified"))
        for z, kind in candidates:
            y, vpd, _ = model.equilibrium_curve(z)
            thresholds.append({"variant": name, "kind": kind, "vpd_mpa": vpd,
                               "leaf_mpa": z, "leaf_carbon_kgC_m2": y[2],
                               "in_routine_range": bool(routine[0] <= vpd <= routine[1]),
                               **model.air_diagnostic(vpd, temperature)})
    certificates = {}
    for name, leaf_bounds in [("E2", (.19, .30)), ("E3", (.145, .30))]:
        certificates[name] = {
            "routine_reference_to_max": models[name].certify_stable_path(
                leaf_bounds, (settings["reference_vpd_mpa"], routine[1])),
            "mapped_legacy_path": models[name].certify_stable_path(
                leaf_bounds, settings["legacy_path_vpd_mpa_bounds"]),
        }
    y = max(reference.equilibria(settings["reference_vpd_mpa"], bounds, reference.stationary_points(bounds)),
            key=lambda state: state[1])
    water_eig = np.linalg.eigvals(reference.jacobian(y, settings["reference_vpd_mpa"])[:2, :2])
    source_paths = [Path(__file__), args.config, baseline_path,
                    *[ROOT/"src/control_carbon"/name for name in [
                        "manzoni_hydraulics.py", "spac_carbon.py", "spac_carbon_two_way.py",
                        "spac_carbon_bifurcation.py", "spac_carbon_stability.py",
                        "spac_carbon_mpa.py", "visitc_hydrology.py"]]]
    report = {"settings": settings, "thresholds": thresholds, "certificates": certificates,
              "e4_local_four_state_check": e4_local_audit(
                  models["E4"], next(row for row in thresholds if row["variant"] == "E4"), bounds),
              "saturation_diagnostic": reference.air_diagnostic(routine[1], temperature),
              "scales": {"gs_max_m_water_per_day_per_mpa": reference.gs_max_per_mpa,
                         "gx_max_m_water_per_day_per_mpa": h.gx_max,
                         "capacitances_m_water_per_mpa": [h.cx, h.cl],
                         "reference_stocks_kgC_m2": stocks.tolist(),
                         "reference_leaf_area_index": c.canopy_coefficient*stocks[0]/k["extinction"],
                         "carbon_loss_times_days": (1/np.array(c.losses)).tolist(),
                         "reference_hydraulic_decay_times_days": (-1/water_eig.real).tolist()},
              "limits": ["synthetic uncalibrated coefficients", "temperature used only for RH diagnosis",
                         "certificates cover stated subintervals, not all VPD starting at zero",
                         "sampled curves do not certify global absence of other equilibria",
                         "unit conversion alone introduces no new dynamics or climate thresholds"],
              "sha256": {str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest() for path in source_paths}}
    args.output.mkdir(parents=True, exist_ok=False)
    for filename, data in [("curves.csv", rows), ("thresholds.csv", thresholds)]:
        with (args.output/filename).open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    (args.output/"results.json").write_text(json.dumps(report, indent=2)+"\n")
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
    saturation = report["saturation_diagnostic"]["saturation_vapor_pressure_mpa"]
    for ax, limits, title in zip(axes, [routine, formal], [
            f"Routine range at {temperature:g} C",
            f"Formal exploration (not climate forcing at {temperature:g} C)"]):
        for name, records in curves.items():
            vpd = np.array([r["vpd_mpa"] for r in records])
            carbon = np.array([r["leaf_carbon_kgC_m2"] for r in records])
            stable = np.array([r["spectral_abscissa_per_day"] < 0 for r in records])
            line, = ax.plot(vpd, np.where(stable, carbon, np.nan), label=name)
            ax.plot(vpd, np.where(stable, np.nan, carbon), color=line.get_color(), ls="--")
        for row in thresholds:
            if limits[0] <= row["vpd_mpa"] <= limits[1]:
                ax.scatter(row["vpd_mpa"], row["leaf_carbon_kgC_m2"],
                           marker="x", color="black", s=40)
                ax.axvline(row["vpd_mpa"], color="black", alpha=.15)
        if limits[1] > saturation:
            ax.axvspan(saturation, limits[1], color="gray", alpha=.15, label=f"RH < 0 at {temperature:g} C")
        ax.set(xlim=limits, ylim=(0, .4), xlabel="VPD (MPa)", ylabel="Leaf carbon (kg C/m2)", title=title)
        ax.ticklabel_format(axis="x", style="plain", useOffset=False)
        ax.legend(fontsize=8)
    fig.suptitle(f"Unchanged synthetic SPAC | p_atm={p} MPa | gas deficit is NOT water potential\n"
                 "Solid: stable equilibrium; dashed: unstable | x: fold / oscillatory boundary")
    fig.savefig(args.output/"pressure_scales.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(json.dumps({k: report[k] for k in ["thresholds", "saturation_diagnostic", "scales", "e4_local_four_state_check"]}, indent=2))


if __name__ == "__main__":
    main()
