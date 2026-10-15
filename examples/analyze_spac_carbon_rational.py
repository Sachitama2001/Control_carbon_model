"""Reproduce algebraic E0--E4 equilibria, stability and fixed-VPD projections."""

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import Polynomial
from scipy.optimize import brentq

from control_carbon.spac_carbon import YEAR_DAYS
from control_carbon.spac_carbon_capacity import reference_scales
from control_carbon.spac_carbon_equilibria import baseline_landmarks, leaf_curve, positive_domain, stock_ratios
from control_carbon.spac_carbon_rational import RationalSPAC
from prepare_manuscript_phase_portraits import model_family, write_equilibrium_matrices


ROOT = Path(__file__).resolve().parents[1]


def spectrum(model, leaf):
    state, vpd, _ = model.curve(leaf)
    return np.linalg.eigvals(model.jacobian(state, vpd))


def expanded_equation(model, vpd):
    """Independent, baseline-unit polynomial from the displayed equations."""
    leaf = Polynomial([0., 1.])
    h, f = model.base.hydraulic, model.base.feedback
    left = 27*h.gs_max*model.base.mole_fraction(vpd)*(1+5*leaf)
    right = h.gx_max*(35+260*leaf-1600*leaf**2)
    if model.label == "E1":
        left *= leaf+f.leaf_reference+2*f.leaf_floor
        right *= 2*(leaf+f.leaf_floor)
    elif model.label == "E2":
        left *= leaf
        right *= f.leaf_reference
    elif model.label == "E3":
        left *= leaf*(leaf+f.leaf_reference+2*f.leaf_floor)
        right *= 2*f.leaf_reference*(leaf+f.leaf_floor)
    elif model.label == "E4":
        left *= ((64/3)*leaf*(f.leaf_reference+f.leaf_floor)
                 +f.wood_reference*(leaf+f.leaf_floor))
        right *= 2*(64/3)*f.leaf_reference*(leaf+f.leaf_floor)
    return left-right


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    bases = model_family()
    models = {label: RationalSPAC(base) for label, base in bases.items()}
    scales = reference_scales(bases["E0"]).state
    records = {"scope": "baseline rational canopy; fixed VPD; no new cycle or mortality claim",
               "parameters": {name: asdict(base) for name, base in bases.items()},
               "shape_bounds": models["E4"].shape_bounds(), "models": {},
               "equilibria": [], "panels": [], "trajectories": [], "refinements": []}
    fig, axes = plt.subplots(2, 3, figsize=(12, 8), constrained_layout=True)
    max_jac_error, max_coefficient_error = 0., 0.
    leaves = np.unique(np.r_[np.geomspace(1e-8, .25, 350), np.linspace(.001, .25, 351)])
    for ax, (label, model) in zip(axes.flat, models.items()):
        mark = model.landmarks()
        eig = np.array([spectrum(model, leaf) for leaf in leaves])
        abscissa = eig.real.max(axis=1)
        vpd = np.array([model.curve(leaf)[1] for leaf in leaves])
        found = []
        for left, right, a, b in zip(leaves[:-1], leaves[1:], abscissa[:-1], abscissa[1:]):
            if a*b >= 0:
                continue
            leaf = brentq(lambda value: max(spectrum(model, value).real), left, right, xtol=1e-15)
            state, v, _ = model.curve(leaf)
            eigenvalues = spectrum(model, leaf)
            neutral = eigenvalues[np.argmax(eigenvalues.real)]
            step = leaf*1e-4
            slope = ((max(spectrum(model, leaf+step).real)-max(spectrum(model, leaf-step).real))
                     /(model.curve(leaf+step)[1]-model.curve(leaf-step)[1])
                     if label not in ["E0", "E1"] else None)
            found.append({"leaf": leaf, "vpd_mpa": v, "state": state.tolist(),
                          "eigenvalues_per_day": [[z.real, z.imag] for z in eigenvalues],
                          "frequency_per_day": abs(neutral.imag),
                          "real_eigenvalue_slope_per_mpa": slope,
                          "classification": "fold" if abs(neutral.imag) < 1e-8
                          else "numerical imaginary-pair crossing; l1 not computed"})
            ax.scatter(v*1000, leaf, marker="D", s=25, color="black")
        for v in [0., .00244, .005]:
            a, b = model.equilibrium_polynomial(v).coef, expanded_equation(model, v).coef
            error = np.max(np.abs(a/a[-1]-b/b[-1]))
            max_coefficient_error = max(max_coefficient_error, float(error))
            if error > 1e-11:
                raise RuntimeError(f"Displayed polynomial disagrees for {label}")
        base = bases[label]
        old_top = positive_domain(base)["leaf_upper_closed_at_zero_vpd"]
        old_leaves = np.unique(np.r_[np.geomspace(1e-8, old_top, 300),
                                    np.linspace(.001, old_top, 300)])
        old_vpd = np.array([leaf_curve(leaf, base)[1] for leaf in old_leaves])
        ax.plot(old_vpd*1000, old_leaves, "--", color=".6", lw=1.2, label="exponential")
        for stable, color, name in [(True, "tab:blue", "rational: stable"),
                                     (False, "tab:red", "rational: unstable")]:
            ax.plot(vpd*1000, np.where((abscissa < 0) == stable, leaves, np.nan),
                    color=color, lw=1.5, label=name)
        ax.axvline(base.air_diagnostic(.00244, 25)["saturation_vapor_pressure_mpa"]*1000,
                   color=".7", ls=":", lw=1)
        ax.set(title=label, xlabel="VPD (kPa)", ylabel=r"$C_L$ (kg C m$^{-2}$)",
               xlim=(0, {"E2": 35, "E3": 8}.get(label, 5.2)), yscale="log", ylim=(1e-6, .5))
        ax.grid(alpha=.2)
        records["models"][label] = {
            "rational_landmarks": mark,
            "exponential_landmarks": baseline_landmarks(base) if label in ["E0", "E1", "E4"] else None,
            "crossings": found, "scan_leaf_interval": [float(leaves[0]), float(leaves[-1])],
            "scan_samples": len(leaves),
            "polynomial_coefficients_ascending_at_2_44kpa": model.equilibrium_polynomial(.00244).coef.tolist(),
        }
        for leaf in [.0001, .01, .05, .15, .24]:
            state, v, _ = model.curve(leaf)
            jac = model.jacobian(state, v)
            fd = np.empty_like(jac)
            for j in range(6):
                step = 1e-6*max(1., abs(state[j])) if j != 2 else min(1e-6, leaf*1e-3)
                delta = np.eye(6)[j]*step
                fd[:, j] = (model.rhs(state+delta, v)-model.rhs(state-delta, v))/(2*step)
            error = np.max(np.abs(fd-jac)/(1+np.abs(jac)))
            max_jac_error = max(max_jac_error, float(error))
            if error > 1e-6:
                raise RuntimeError(f"Jacobian mismatch: {label} {leaf} {error}")
    axes.flat[0].legend(fontsize=7, loc="lower left")
    z = np.linspace(0, 4, 401)
    axes.flat[5].plot(z, 1-np.exp(-z), "--", label="exponential")
    axes.flat[5].plot(z, z/(1+z), label="rational")
    axes.flat[5].plot(z, z, ":", label="linear")
    axes.flat[5].set(xlabel=r"$q C_L$", ylabel="canopy response", ylim=(0, 1.2))
    axes.flat[5].legend()
    fig.savefig(args.output/"equilibrium_comparison.pdf")
    fig.savefig(args.output/"equilibrium_comparison.png", dpi=160)
    plt.close(fig)

    cases = []
    for label, single, double, none in [
        ("E0", .002, .00244, .0026), ("E1", .0001, .0009, .0012),
    ]:
        cases.extend([(label, single, 1), (label, double, 2),
                      (label, models[label].landmarks()["fold_vpd_mpa"], 1), (label, none, 0)])
    cases += [("E2", .00244, 1), ("E2", .030, 1),
              ("E3", .00244, 1), ("E3", .0061, 1),
              ("E4", .00244, 1), ("E4", .003, 1), ("E4", .0049, 0)]
    fig1, axes1 = plt.subplots(2, 4, figsize=(14, 7), constrained_layout=True)
    fig2, axes2 = plt.subplots(3, 3, figsize=(11, 10), constrained_layout=True)
    times = np.unique(np.r_[0., np.geomspace(1e-7, 30, 180),
                           np.linspace(30, 10*YEAR_DAYS, 800),
                           np.linspace(10*YEAR_DAYS, 400*YEAR_DAYS, 2500)])
    max_residual = 0.
    other = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2)]
    for index, (label, vpd, expected) in enumerate(cases):
        model = models[label]
        roots = model.equilibria(vpd)
        if len(roots) != expected:
            raise RuntimeError(f"Unexpected root count: {label} {vpd}")
        ax = axes1.flat[index] if index < 8 else axes2[other[index-8]]
        panel_id = f"{index:02d}_{label}"
        reference = roots[0] if roots else model.curve(.15)[0]
        initials = []
        for psi, factor in [(-.8, .9), (-2.5, 1.1)]:
            stocks = stock_ratios(model.base)*max(reference[2]*factor, .002)
            fl, _, fw = model.factors(np.r_[0., psi, stocks])
            initials.append(np.r_[fl*psi/(fl+fw), psi, stocks])
        for root in roots:
            eigenvalues = np.linalg.eigvals(model.jacobian(root, vpd))
            abscissa = float(max(eigenvalues.real))
            residual = float(max(abs(model.rhs(root, vpd)/scales)))
            max_residual = max(max_residual, residual)
            if residual > 1e-10 or max(abs((model.capacity(root)-root[2:])/scales[2:])) > 1e-10:
                raise RuntimeError("Equilibrium RHS/capacity validation failed")
            ax.scatter(root[1], root[2], color="black", s=45, zorder=5,
                       marker="D" if abs(abscissa) < 1e-10 else "o" if abscissa < 0 else "X",
                       facecolors="white" if abscissa < -1e-10 else "black")
            records["equilibria"].append({
                "panel": panel_id, "model": label, "vpd_mpa": vpd, "state": root.tolist(),
                "jacobian": model.jacobian(root, vpd).tolist(),
                "eigenvalues_per_day": [[z.real, z.imag] for z in eigenvalues],
                "max_real_eigenvalue_per_day": abscissa, "scaled_residual_per_day": residual,
                "capacity": model.capacity(root).tolist(),
            })
            for sign in [-1, 1]:
                initial = root.copy()
                initial[1] += sign*.015
                initial[2] *= 1+sign*.03
                initials.append(initial)
        records["panels"].append({"id": panel_id, "model": label, "vpd_mpa": vpd,
                                  "positive_roots": len(roots), "horizon_years": 400,
                                  "formal_at_25c": not model.base.air_diagnostic(vpd, 25)[
                                      "admissible_at_diagnostic_temperature"]})
        for number, initial in enumerate(initials):
            run = model.simulate(initial, times, vpd)
            ys = run["state"]
            ax.plot(ys[:, 1], ys[:, 2], lw=.9, alpha=.85)
            ax.plot(ys[0, 1], ys[0, 2], "+", color=".5", ms=5)
            if run["status"] == "leaf_domain_exit":
                ax.plot(ys[-1, 1], ys[-1, 2], "v", color="black", ms=4)
            filename = f"{panel_id}_{number}.npz"
            np.savez_compressed(args.output/filename, **run)
            records["trajectories"].append({
                "file": filename, "panel": panel_id, "initial": initial.tolist(),
                "status": run["status"], "last_day": float(run["day"][-1]),
                "final": ys[-1].tolist(),
            })
            if number == 0:
                refined = model.simulate(initial, times, vpd, rtol=2e-10)
                common = run["day"][run["day"] <= refined["day"][-1]]
                original = np.column_stack([np.interp(common, run["day"], ys[:, j]) for j in range(6)])
                tight = np.column_stack([np.interp(common, refined["day"], refined["state"][:, j])
                                         for j in range(6)])
                error = float(np.max(abs((original-tight)/scales)))
                event_error = abs(run["day"][-1]-refined["day"][-1])
                if refined["status"] != run["status"] or error > 2e-5 or event_error > .02:
                    raise RuntimeError(f"Trajectory refinement failed: {panel_id}, {error}, {event_error}")
                records["refinements"].append({
                    "panel": panel_id, "max_scaled_state_error": error,
                    "last_time_difference_days": event_error,
                    "same_status": refined["status"] == run["status"],
                })
        ax.set(title=f"{label}, {vpd*1000:.6g} kPa; {len(roots)} positive roots",
               xlabel=r"$\psi_l$ (MPa)", ylabel=r"$C_L$ (kg C m$^{-2}$)")
        ax.grid(alpha=.2)
    axes2[0, 2].axis("off")
    axes2[1, 2].axis("off")
    for name, fig in [("phase_E0_E1", fig1), ("phase_E2_E3_E4", fig2)]:
        fig.savefig(args.output/f"{name}.pdf")
        fig.savefig(args.output/f"{name}.png", dpi=160)
        plt.close(fig)
    write_equilibrium_matrices(args.output/"equilibrium_jacobians.tex", records["equilibria"])
    records["validation"] = {
        "max_scaled_equilibrium_residual_per_day": max_residual,
        "max_normalized_fd_jacobian_error": max_jac_error,
        "max_normalized_polynomial_coefficient_error": max_coefficient_error,
        "max_trajectory_refinement_error": max(r["max_scaled_state_error"] for r in records["refinements"]),
    }
    sources = [Path(__file__), ROOT/"examples/prepare_manuscript_phase_portraits.py",
               ROOT/"tests/test_spac_carbon_rational.py",
               *sorted((ROOT/"src/control_carbon").glob("*.py")),
               ROOT/"configs/spac_carbon_stage_b_synthetic.json",
               ROOT/"configs/spac_carbon_pressure_review.json"]
    records["source_sha256"] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sources}
    (args.output/"results.json").write_text(json.dumps(records, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"panels": len(records["panels"]), "equilibria": len(records["equilibria"]),
                      "trajectories": len(records["trajectories"]),
                      "validation": records["validation"]}, indent=2))


if __name__ == "__main__":
    main()
