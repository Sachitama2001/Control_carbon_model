"""Projected six-state trajectories for the Japanese manuscript revision."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from control_carbon.manzoni_hydraulics import HydraulicParameters, transpiring_equilibria
from control_carbon.spac_carbon import (
    YEAR_DAYS, CarbonParameters, hydraulic_factor, positive_carbon_equilibrium,
)
from control_carbon.spac_carbon_capacity import reference_scales
from control_carbon.spac_carbon_equilibria import (
    baseline_landmarks, baseline_shape_proof, leaf_curve, positive_domain,
    positive_equilibria, stock_ratios,
)
from control_carbon.spac_carbon_mpa import PressureSPAC
from control_carbon.spac_carbon_two_way import StructuralFeedback, feedback_factors


ROOT = Path(__file__).resolve().parents[1]


def model_family():
    cfg = json.loads((ROOT / "configs/spac_carbon_stage_b_synthetic.json").read_text())
    pressure = json.loads((ROOT / "configs/spac_carbon_pressure_review.json").read_text())
    hydraulic = HydraulicParameters(**cfg["hydraulic"])
    p = cfg["carbon"]
    carbon = CarbonParameters(
        p["assimilation_max_per_year"] / YEAR_DAYS, tuple(p["allocation"]),
        tuple(k / YEAR_DAYS for k in p["losses_per_year"]), tuple(p["litter_fractions"]),
        p["specific_leaf_area"], p["carbon_dry_fraction"], p["extinction"],
    )
    psi = transpiring_equilibria(cfg["root_potential"], cfg["reference_vpd"], hydraulic)[0]
    reference = positive_carbon_equilibrium(hydraulic_factor(psi, hydraulic), carbon)
    if reference is None:
        raise RuntimeError("Reference carbon equilibrium must be positive")
    return {
        label: PressureSPAC(
            hydraulic, carbon,
            StructuralFeedback(
                reference[0], reference[1], cfg["leaf_floor"], cfg["minimum_leaf"],
                *strengths, wood_exponent=cfg["wood_exponent"],
            ),
            cfg["root_potential"], pressure["atmospheric_pressure_mpa"],
        )
        for label, strengths in cfg["variants"].items()
    }


def roots_at(vpd, model, scales):
    if model.feedback.leaf_demand == 0 or model.feedback.wood_supply:
        return [r["state"] for r in positive_equilibria(vpd, model, scales)["roots"]]
    top = positive_domain(model)["leaf_upper_closed_at_zero_vpd"]
    if vpd == 0:
        return [leaf_curve(top, model)[0]]
    low = top / 2
    while leaf_curve(low, model)[1] < vpd:
        low /= 2
        if low <= np.finfo(float).tiny:
            raise RuntimeError("Could not bracket the positive E2/E3 equilibrium")
    leaf = brentq(lambda value: leaf_curve(value, model)[1] - vpd, low, top,
                  xtol=np.finfo(float).tiny)
    return [leaf_curve(leaf, model)[0]]


def write_equilibrium_matrices(path, equilibria):
    def number(value):
        if value == 0:
            return "0"
        mantissa, exponent = f"{value:.4e}".split("e")
        return mantissa if int(exponent) == 0 else rf"{mantissa}\times10^{{{int(exponent)}}}"

    blocks = []
    for index, rec in enumerate(equilibria):
        if index and index % 3 == 0:
            blocks.append(r"\clearpage")
        vector = r"\\".join(number(v) for v in rec["state"])
        matrix = r"\\".join("&".join(number(v) for v in row) for row in rec["jacobian"])
        blocks.append(
            rf"\paragraph{{{rec['model']}, $v={rec['vpd_mpa'] * 1000:.8g}$ kPa.}}"
            "\n" + r"\begin{equation*}" + "\n"
            + r"\resizebox{.96\linewidth}{!}{$\displaystyle "
            + r"y_*=\begin{pmatrix}" + vector + r"\end{pmatrix},\qquad "
            + r"J_*=\begin{pmatrix}" + matrix + r"\end{pmatrix}.$}"
            + "\n" + r"\end{equation*}" + "\n"
        )
    path.write_text("\n".join(blocks))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    models = model_family()
    scales = reference_scales(models["E0"])
    proof = baseline_shape_proof(models["E4"])
    upper = float(proof["positive_log_D0_derivative_upper_interval"][1])
    lower2 = 1 / .05
    lower3 = models["E3"].feedback.leaf_floor / (
        .05 * (.05 + models["E3"].feedback.leaf_floor)
    )
    if not upper < min(lower2, lower3):
        raise RuntimeError("E2/E3 strict monotonicity bound failed")
    cases = []
    for label, single, double, none in [
        ("E0", .002, .00244, .0026), ("E1", .0001, .001, .002),
    ]:
        fold = baseline_landmarks(models[label])["fold_vpd_mpa"]
        cases.extend((label, v, count, kind) for v, count, kind in [
            (single, 1, "one positive root"), (double, 2, "two positive roots"),
            (fold, 1, "fold: one double root"), (none, 0, "no positive root"),
        ])
    cases += [
        ("E2", .00244, 1, "stable positive root"),
        ("E2", .0155, 1, "unstable positive root; formal VPD"),
        ("E3", .00244, 1, "stable positive root"),
        ("E3", .005, 1, "unstable positive root; formal VPD"),
        ("E4", .00244, 1, "stable positive root"),
        ("E4", .0028, 1, "unstable positive root"),
        ("E4", .0049, 0, "no positive root; formal VPD"),
    ]
    trajectories, equilibria, panels, refinements = [], [], [], []
    fig_fold, axes_fold = plt.subplots(2, 4, figsize=(15, 7), constrained_layout=True)
    fig_other, axes_other = plt.subplots(3, 3, figsize=(12, 10), constrained_layout=True)
    other_positions = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2)]
    times = np.unique(np.r_[
        0., np.geomspace(1e-7, 30, 300),
        np.linspace(30, 10 * YEAR_DAYS, 1200),
        np.linspace(10 * YEAR_DAYS, 400 * YEAR_DAYS, 5000),
    ])
    max_residual = 0.
    for index, (label, vpd, expected, kind) in enumerate(cases):
        model = models[label]
        roots = roots_at(vpd, model, scales)
        if len(roots) != expected:
            raise RuntimeError(f"{label} {vpd}: root count {len(roots)} != {expected}")
        ax = axes_fold.flat[index] if index < 8 else axes_other[other_positions[index - 8]]
        reference = roots[0] if roots else leaf_curve(.18, model)[0]
        initials = []
        for psi, leaf_multiplier in [(-.8, .9), (-2.5, 1.1)]:
            stocks = stock_ratios(model) * max(reference[2] * leaf_multiplier, .002)
            fl, _, fw, *_ = feedback_factors(stocks[0], stocks[1], model.feedback)
            initials.append(np.r_[fl * psi / (fw + fl), psi, stocks])
        for root in roots:
            for sign in [-1, 1]:
                initial = root.copy()
                initial[1] += sign * .015
                initial[2] *= 1 + sign * .03
                initials.append(initial)
        panel_id = f"{label}_{vpd:.15g}"
        root_records = []
        for root in roots:
            residual = float(np.max(np.abs(model.rhs(root, vpd) / scales.state)))
            max_residual = max(max_residual, residual)
            eigenvalues = np.linalg.eigvals(model.jacobian(root, vpd))
            abscissa = float(max(eigenvalues.real))
            stable = abscissa < -1e-10
            center = abs(abscissa) <= 1e-10
            marker = "D" if center else "o" if stable else "X"
            ax.scatter(root[1], root[2], marker=marker, s=60, color="black",
                       zorder=5, facecolors="white" if stable else "black")
            rec = {"panel": panel_id, "model": label, "vpd_mpa": vpd,
                   "state": root.tolist(), "max_real_eigenvalue_per_day": abscissa,
                   "scaled_residual_per_day": residual,
                   "jacobian": model.jacobian(root, vpd).tolist()}
            equilibria.append(rec)
            root_records.append(rec)
        for number, initial in enumerate(initials):
            run = model.simulate(initial, times, lambda t, v=vpd: v, rtol=1e-9)
            ys = run["state"]
            ax.plot(ys[:, 1], ys[:, 2], lw=1, alpha=.8)
            ax.plot(ys[0, 1], ys[0, 2], "+", color="grey", ms=5)
            if run["status"] == "vegetation_domain_exit":
                ax.scatter(ys[-1, 1], ys[-1, 2], marker="v", s=23, color="tab:red")
            projected = ys[:, [1, 2]]
            distance = np.linalg.norm(np.diff(projected, axis=0) / [4., .4], axis=1)
            length = np.r_[0., np.cumsum(distance)]
            for fraction in [.3, .65]:
                at = min(int(np.searchsorted(length, fraction * length[-1])), len(ys) - 2)
                before = max(0, at - 2)
                if np.linalg.norm(projected[at + 1] - projected[before]) > 1e-9:
                    ax.annotate("", xy=projected[at + 1], xytext=projected[before],
                                arrowprops={"arrowstyle": "->", "color": "grey", "lw": .7})
            trajectories.append({"panel": panel_id, "initial_index": number,
                                 "initial_state": initial.tolist(), "status": run["status"],
                                 "end_year": float(run["time_days"][-1] / YEAR_DAYS)})
            if number == 0:
                fine = model.simulate(initial, times, lambda t, v=vpd: v, rtol=1e-10)
                # Refined exit times can differ slightly; compare at common times.
                common = min(len(ys), len(fine["state"])) - 1
                error = float(np.max(np.abs(
                    (ys[:common] - fine["state"][:common]) / scales.state
                )))
                if fine["status"] != run["status"] or error > 1e-5:
                    raise RuntimeError(f"Trajectory refinement failed: {panel_id}: {error}")
                refinements.append({"panel": panel_id, "max_scaled_error": error,
                                    "status_agrees": True})
            np.savez_compressed(args.output / f"{panel_id}_trajectory_{number}.npz",
                                day=run["time_days"], state=ys)
        ax.set(xlabel=r"$\psi_l$ (MPa)", ylabel=r"$C_L$ (kg C m$^{-2}$)",
               title=f"{label}: {vpd * 1000:.6g} kPa\n{kind}", ylim=(0, None))
        ax.grid(alpha=.25)
        panels.append({"id": panel_id, "positive_root_count": expected, "roots": root_records,
                       "classification": kind})
        print(panel_id, kind, flush=True)
    if max_residual > 1e-9:
        raise RuntimeError(f"Reconstructed equilibria fail full RHS: {max_residual}")
    for pos in [(0, 2), (1, 2)]:
        axes_other[pos].axis("off")
        axes_other[pos].text(
            .05, .8, "Projected full 6D trajectories\n"
            "o stable; X unstable; D fold\n"
            "+ initial; red v: L_min exit\n"
            "400-year hold or L_min exit\n"
            "No claim of global basin boundary\n"
            "Formal VPD: above 25 C saturation",
            transform=axes_other[pos].transAxes, va="top", fontsize=10,
        )
    for fig, name in [(fig_fold, "phase_E0_E1"), (fig_other, "phase_E2_E3_E4")]:
        fig.savefig(args.output / f"{name}.pdf")
        fig.savefig(args.output / f"{name}.png", dpi=220)
        plt.close(fig)
    paths = [Path(__file__), *list((ROOT / "src/control_carbon").glob("spac_carbon*.py")),
             ROOT / "configs/spac_carbon_stage_b_synthetic.json",
             ROOT / "configs/spac_carbon_pressure_review.json"]
    results = {
        "scope": "Fixed VPD full-six-state trajectories projected to psi_leaf,C_L; not a closed 2D vector field",
        "root_count_scope": "Isolated conducting strictly positive carbon equilibria only; zero-carbon continua excluded",
        "E2_E3_shape_bounds": {"inherited_outward_upper_log_D0": upper,
                              "E2_negative_lower": lower2, "E3_negative_lower": lower3},
        "sample_days": times.tolist(),
        "panels": panels, "equilibria": equilibria, "trajectories": trajectories,
        "refinements": refinements, "max_equilibrium_residual_scaled_per_day": max_residual,
        "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in paths},
    }
    (args.output / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    write_equilibrium_matrices(args.output / "equilibrium_jacobians.tex", equilibria)
    with (args.output / "equilibria.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(["panel", "model", "vpd_mpa", "psi_x", "psi_l", "C_L", "C_W",
                         "C_R", "C_S", "max_real_eigenvalue_per_day"])
        for rec in equilibria:
            writer.writerow([rec["panel"], rec["model"], rec["vpd_mpa"], *rec["state"],
                             rec["max_real_eigenvalue_per_day"]])
    print(f"Verified {len(panels)} panels, {len(equilibria)} equilibria, "
          f"{len(trajectories)} trajectories, {len(refinements)} refinements.")


if __name__ == "__main__":
    main()
