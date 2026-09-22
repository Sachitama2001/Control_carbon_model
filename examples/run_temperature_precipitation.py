"""Reproduce WP1 / first WP2 T/P controls without asserting tipping.

Run: PYTHONPATH=src python examples/run_temperature_precipitation.py
All parameter choices are uncalibrated assumptions; see the output manifest.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import platform
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import cumulative_trapezoid, solve_ivp

from control_carbon.temperature_precipitation import (
    Climate, Controls, Parameters, STATE_NAMES, STATE_UNITS, budget_residuals,
    domain_violation, equilibrium, fluxes, frozen_carbon_capacity, integrate,
    jacobian, parameter_manifest, rhs,
)


THRESHOLDS = {"budget": 1e-10, "equilibrium": 1e-8, "stable_max_real": -1e-8,
              "refinement": 1e-5, "domain": 1e-6}


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def equilibrium_record(eq, climate, p, control, **labels):
    f = fluxes(eq.state, climate, p, control)
    return {**labels, "T_degC": climate.temperature, "P_mm_day": climate.precipitation,
            **dict(zip(STATE_NAMES, eq.state)), "C_total": eq.state[:4].sum(),
            "GPP": f.gpp, "NEP": f.nep, "mortality_flux": f.mortality.sum(),
            "max_real_day_inv": eq.max_real_part,
            "residual": eq.residual, "condition_number": eq.condition_number,
            "eigenvalues": json.dumps([[v.real, v.imag] for v in eq.eigenvalues]),
            "solver_evaluations": eq.evaluations}


def diagnostic_series(times, states, climate, p, control):
    rows = []
    for t, x in zip(times, states):
        f = fluxes(x, climate, p, control)
        cap = frozen_carbon_capacity(x, climate, p, control)
        bc, bw = budget_residuals(x, climate, p, control)
        rows.append({"day": t, "T_degC": climate.temperature,
                     "P_mm_day": climate.precipitation, **dict(zip(STATE_NAMES, x)),
                     "C_total": x[:4].sum(), "capacity": cap.sum(),
                     "potential": cap.sum()-x[:4].sum(), "GPP": f.gpp,
                     "Ra": f.respiration.sum(), "Rh": f.soil_respiration,
                     "NEP": f.nep, "C_export": f.carbon_export,
                     "mortality_L": f.mortality[0], "mortality_S": f.mortality[1],
                     "mortality_R": f.mortality[2], "mortality_total": f.mortality.sum(),
                     "E_leaf": f.transpiration, "E_soil": f.evaporation,
                     "drainage": f.drainage, "runoff": f.runoff,
                     "water_export": f.water_export,
                     "release_L": f.water_release[0], "release_S": f.water_release[1],
                     "release_R": f.water_release[2],
                     "q_OR": f.transport[0], "q_RS": f.transport[1], "q_SL": f.transport[2],
                     "carbon_budget_residual": bc, "water_budget_residual": bw,
                     "heat_hazard": f.heat_hazard, "drought_hazard": f.drought_hazard,
                     "carbon_deficit_rate": max(0., -f.nep)})
    for src, dst in (("heat_hazard", "heat_dose"), ("drought_hazard", "drought_dose"),
                     ("carbon_deficit_rate", "cumulative_carbon_deficit")):
        vals = cumulative_trapezoid([r[src] for r in rows], times, initial=0.)
        for r, val in zip(rows, vals):
            r[dst] = val
    return rows


def plot_experiment(out, stage, series, equilibria, branches):
    colors = {"0": "#777777", "T": "#cc6b22", "P": "#288ac5", "TP": "#943da2"}
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    for key, rows in series.items():
        t = np.array([r["day"] for r in rows])/365.25
        c = np.array([r["C_total"] for r in rows])
        soil = np.array([r["W_O"] for r in rows])
        plant = np.array([r["W_L"]+r["W_S"]+r["W_R"] for r in rows])
        eq = equilibria[key].state
        axes[0, 0].plot(t, c, color=colors[key], label=key)
        axes[0, 0].axhline(eq[:4].sum(), color=colors[key], linestyle=":", alpha=.6)
        axes[0, 1].plot(soil, c, color=colors[key], label=key)
        axes[0, 1].scatter(soil[0], c[0], color=colors[key], marker="o", s=15)
        axes[0, 1].scatter(soil[-1], c[-1], color=colors[key], marker=">", s=30)
        axes[0, 1].scatter(eq[7], eq[:4].sum(), color=colors[key], marker="x", s=40)
        axes[0, 2].plot(plant, c, color=colors[key])
        axes[0, 2].scatter(plant[-1], c[-1], color=colors[key], marker=">", s=30)
        for ax, horizontal in ((axes[0, 1], soil), (axes[0, 2], plant)):
            a, b = max(1, len(c)//5), max(2, len(c)//4)
            ax.annotate("", xy=(horizontal[b], c[b]), xytext=(horizontal[a], c[a]),
                        arrowprops={"arrowstyle": "->", "color": colors[key]})
        axes[1, 0].plot(t, [r["capacity"] for r in rows], color=colors[key], label=key)
        axes[1, 1].plot(t, [r["NEP"]*365.25 for r in rows], color=colors[key])
    t = np.array([r["day"] for r in series["0"]])/365.25
    interaction = np.array([r["C_total"] for r in series["TP"]])-np.array(
        [r["C_total"] for r in series["T"]])-np.array([r["C_total"] for r in series["P"]])+np.array(
        [r["C_total"] for r in series["0"]])
    axes[1, 2].plot(t, interaction, color=colors["TP"])
    axes[1, 2].axhline(0, color="black", lw=.5)
    extent = max(1., float(np.max(np.abs(interaction)))*1.1)
    axes[1, 2].set_ylim(-extent, extent)
    axes[1, 2].set_title(f"max |I_C| = {np.max(np.abs(interaction)):.3g}")
    axes[0, 0].set(xlabel="Time (years)", ylabel="Total C (Mg C/ha)", title="Stocks; dotted: coupled equilibrium")
    axes[0, 1].set(xlabel="Soil water (mm)", ylabel="Total C (Mg C/ha)", title="8D trajectory projected to 2D; x: QSE")
    axes[0, 2].set(xlabel="Plant water (mm)", ylabel="Total C (Mg C/ha)", title="Projected trajectory; >: final sample")
    axes[1, 0].set(xlabel="Time (years)", ylabel="Frozen C capacity (Mg C/ha)")
    axes[1, 1].set(xlabel="Time (years)", ylabel="NEP (Mg C/ha/year)")
    axes[1, 2].set(xlabel="Time (years)", ylabel="I_C = TP - T - P + 0 (Mg C/ha)")
    axes[0, 0].legend()
    fig.suptitle(f"T/P theoretical model: {stage}; assumptions, not site calibration")
    fig.savefig(out/f"{stage}_trajectories.png", dpi=160)
    fig.savefig(out/f"{stage}_trajectories.pdf")
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for j, driver in enumerate(("T", "P")):
        for direction, style in (("forward", "-"), ("reverse", "--")):
            rows = [r for r in branches if r["driver"] == driver and r["direction"] == direction]
            xkey = "T_degC" if driver == "T" else "P_mm_day"
            axes[0, j].plot([r[xkey] for r in rows], [r["C_total"] for r in rows], style, label=direction)
            axes[1, j].plot([r[xkey] for r in rows], [r["max_real_day_inv"] for r in rows], style)
        axes[0, j].set(ylabel="Equilibrium C (Mg C/ha)", xlabel=xkey)
        axes[1, j].set(ylabel="max Re(lambda) (1/day)", xlabel=xkey)
        axes[1, j].axhline(0, color="black", lw=.6)
        axes[0, j].legend()
    fig.suptitle(f"{stage}: local parameter continuation; not exhaustive root search")
    fig.savefig(out/f"{stage}_single_driver_branches.png", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("artifacts/temperature_precipitation/wp1_wp2"))
    parser.add_argument("--years", type=float, default=20.)
    parser.add_argument("--precipitation-end", type=float, default=4., help="mm/day; use 2 for water-limited comparison")
    args = parser.parse_args()
    if args.years <= 0:
        parser.error("years must be positive")
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    p = Parameters()
    # Constant forcing is a special case of monthly piecewise-constant forcing.
    # The control step is imposed at day 0; this is NOT a finite-rate experiment.
    climates = {"0": Climate(27, 6), "T": Climate(31, 6),
                "P": Climate(27, args.precipitation_end), "TP": Climate(31, args.precipitation_end)}
    times = np.unique(np.r_[np.arange(0, args.years*365.25, 5.), args.years*365.25])
    summary, all_eq, all_branches, interaction_rows, clamp_rows = {}, [], [], [], []
    solver = {"method": "Radau", "max_step": 10., "rtol": 1e-8, "atol": 1e-10}
    refined_solver = {"method": "BDF", "max_step": 5., "rtol": 2e-9, "atol": 2e-11}
    for stage in ("none", "heat", "water", "additive"):
        control = Controls(mortality=stage)
        base = equilibrium(climates["0"], p, control)
        branches = []
        max_reversal = 0.
        for driver, values in (("T", np.linspace(27, 33, 13)), ("P", np.linspace(6, 2, 13))):
            guess = base.state
            forward = []
            for direction, seq in (("forward", values), ("reverse", values[::-1])):
                for i, value in enumerate(seq):
                    climate = Climate(value, 6) if driver == "T" else Climate(27, value)
                    eq = equilibrium(climate, p, control, guess)
                    if eq.state[:4].min() < 1e-6:
                        raise RuntimeError("continuation left the positive branch; requires WP3 branch analysis")
                    guess = eq.state
                    branches.append(equilibrium_record(eq, climate, p, control,
                                                      stage=stage, driver=driver, direction=direction))
                    if direction == "forward":
                        forward.append(eq.state)
                    else:
                        max_reversal = max(max_reversal, float(np.max(np.abs(eq.state-forward[-i-1])/np.maximum(1., abs(eq.state)))))
        all_branches.extend(branches)
        series, equilibria, refinements = {}, {}, {}
        max_budget = np.zeros(2)
        max_domain = 0.
        for key, climate in climates.items():
            eq = equilibrium(climate, p, control, base.state)
            equilibria[key] = eq
            all_eq.append(equilibrium_record(eq, climate, p, control, stage=stage, experiment=key))
            x = integrate(base.state, times, climate, p, control, **solver)
            series[key] = diagnostic_series(times, x, climate, p, control)
            write_csv(out/f"{stage}_{key}_trajectory.csv", series[key])
            refined = integrate(base.state, times, climate, p, control, **refined_solver)
            refinements[key] = float(np.max(abs(x-refined)/np.maximum(1., abs(refined))))
            max_budget = np.maximum(max_budget, np.max(np.abs([[r["carbon_budget_residual"], r["water_budget_residual"]] for r in series[key]]), axis=0))
            max_domain = max(max_domain, max(domain_violation(z, p) for z in x))
            np.savez_compressed(out/f"{stage}_{key}.npz", time_day=times, states=x, refined_states=refined,
                                equilibrium=eq.state, eigenvalues=eq.eigenvalues,
                                jacobian=jacobian(eq.state, climate, p, control),
                                state_names=STATE_NAMES, state_units=STATE_UNITS)
        for j, t in enumerate(times):
            interaction_rows.append({"stage": stage, "day": t, **{
                f"I_{name}": series["TP"][j][name]-series["T"][j][name]-series["P"][j][name]+series["0"][j][name]
                for name in ("C_total", "capacity", "NEP", "mortality_total")}})
        initial_hazard = fluxes(base.state, climates["0"], p, control).drought_hazard
        clamps = {"heat_clamped": replace(control, heat_temperature=27),
                  "evap_temperature_only": replace(control, heat_temperature=27,
                      photo_temperature=27, respiration_temperature=27),
                  "drought_clamped": replace(control, drought_rate=initial_hazard)}
        clamp_series = {}
        for name, clamped in clamps.items():
            climate = climates["P"] if name == "drought_clamped" else climates["T"]
            states = integrate(base.state, times, climate, p, clamped, **solver)
            rows = diagnostic_series(times, states, climate, p, clamped)
            write_csv(out/f"{stage}_{name}.csv", rows)
            clamp_series[name] = rows
            clamp_rows.append({"stage": stage, "clamp": name,
                               "C_initial": base.state[:4].sum(), "C_end": states[-1, :4].sum(),
                               "max_domain_violation": max(domain_violation(x, p) for x in states)})
        # Fixed absolute water intervention: carbon-only ODE, not water-conserving.
        fixed = solve_ivp(lambda t, c: rhs(t, np.r_[c, base.state[4:]], climates["T"], p, control)[:4],
                          (times[0], times[-1]), base.state[:4], t_eval=times, **solver)
        if not fixed.success:
            raise RuntimeError(fixed.message)
        fixed_x = np.c_[fixed.y.T, np.tile(base.state[4:], (len(times), 1))]
        fixed_rows = diagnostic_series(times, fixed_x, climates["T"], p, control)
        for row, x in zip(fixed_rows, fixed_x):
            intervention = -rhs(0, x, climates["T"], p, control)[4:]
            row.update(dict(zip(("water_intervention_L", "water_intervention_S", "water_intervention_R", "water_intervention_O"), intervention)))
            row["max_domain_violation"] = domain_violation(x, p)
            # Here dW/dt=0; include the externally imposed four water fluxes.
            external = (row["P_mm_day"]-row["E_leaf"]-row["E_soil"]
                        -row["drainage"]-row["runoff"]-row["water_export"])
            row["water_budget_residual"] = -(external+intervention.sum())
        write_csv(out/f"{stage}_fixed_water_intervention.csv", fixed_rows)
        clamp_series["fixed_water_intervention"] = fixed_rows
        clamp_rows.append({"stage": stage, "clamp": "fixed_water_intervention",
                           "C_initial": base.state[:4].sum(), "C_end": fixed.y[:, -1].sum(),
                           "max_domain_violation": max(domain_violation(x, p) for x in fixed_x)})
        fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
        for name, rows in {"T_total": series["T"], "P_total": series["P"], **clamp_series}.items():
            c = [r["C_total"] for r in rows]
            axes[0].plot(times/365.25, c, label=name)
            axes[1].plot([r["W_O"] for r in rows], c, label=name)
        axes[0].set(xlabel="Time (years)", ylabel="Total C (Mg C/ha)")
        axes[1].set(xlabel="Soil W (mm)", ylabel="Total C (Mg C/ha)", title="Projected trajectories")
        axes[0].legend(fontsize=7)
        fig.suptitle(f"{stage}: pathway controls; fixed water is an intervention")
        fig.savefig(out/f"{stage}_clamps.png", dpi=160)
        plt.close(fig)
        eq_interaction = (equilibria["TP"].state[:4].sum()-equilibria["T"].state[:4].sum()
                          -equilibria["P"].state[:4].sum()+base.state[:4].sum())
        eq_spectral_interaction = (equilibria["TP"].max_real_part-equilibria["T"].max_real_part
                                  -equilibria["P"].max_real_part+base.max_real_part)
        summary[stage] = {"baseline_state": base.state.tolist(), "baseline_residual": base.residual,
                          "slowest_baseline_years": -1/base.max_real_part/365.25,
                          "baseline_max_real": base.max_real_part,
                          "C_equilibria": {k: e.state[:4].sum() for k, e in equilibria.items()},
                          "equilibrium_I_C": eq_interaction,
                          "equilibrium_I_max_real": eq_spectral_interaction,
                          "transient_final_I_C": [r for r in interaction_rows if r["stage"] == stage][-1]["I_C_total"],
                          "max_budget_residual": max_budget.tolist(), "max_domain_violation": max_domain,
                          "refinement": refinements, "forward_reverse_max_scaled_difference": max_reversal,
                          "max_axis_eigenvalue_real": max(r["max_real_day_inv"] for r in branches),
                          "continuation_points_including_reverse": len(branches),
                          "claim": "sampled positive branches and quantitative interaction only; uniqueness and R-tipping untested"}
        if (max(max_budget) > THRESHOLDS["budget"] or max(refinements.values()) > THRESHOLDS["refinement"]
                or max_domain > THRESHOLDS["domain"] or base.max_real_part >= THRESHOLDS["stable_max_real"]):
            raise RuntimeError(f"acceptance thresholds failed: {stage}")
        plot_experiment(out, stage, series, equilibria, branches)
        print(stage, json.dumps(summary[stage], ensure_ascii=False), flush=True)
    write_csv(out/"equilibria.csv", all_eq)
    write_csv(out/"single_driver_continuation.csv", all_branches)
    write_csv(out/"interactions.csv", interaction_rows)
    write_csv(out/"clamp_summary.csv", clamp_rows)
    write_csv(out/"parameters.csv", [{"name": k, **v} for k, v in parameter_manifest(p).items()])
    root = Path(__file__).resolve().parents[1]
    paths = [root/"src/control_carbon/temperature_precipitation.py", Path(__file__).resolve(),
             root/"docs/temperature_precipitation_tipping_research_plan.md",
             root/"docs/temperature_precipitation_implementation.md"]
    manifest = {"model_level": "model-agnostic derived analysis; source-informed low-dimensional theoretical model",
                "parameters": parameter_manifest(p), "solver": solver, "refined_solver": refined_solver,
                "equilibrium_solver": {"method": "bounded least_squares on scaled states",
                                       "ftol": 1e-13, "xtol": 1e-13, "gtol": None,
                                       "raw_residual_acceptance": THRESHOLDS["equilibrium"]},
                "thresholds": THRESHOLDS, "forcing": {k: asdict(v) for k, v in climates.items()},
                "forcing_protocol": "step at day 0, constant for 20 years by default; no finite-rate test",
                "horizon_years": args.years, "diagnostic_step_days": 5,
                "dose_quadrature": "trapezoidal on saved 5-day samples; diagnostic, not exact exposure integral",
                "state_names": STATE_NAMES, "state_units": STATE_UNITS,
                "units": {"T": "degC", "P": "mm day-1", "time": "day",
                          "carbon_flux": "Mg C ha-1 day-1", "water_flux": "mm day-1"},
                "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
                "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
                "git_status": subprocess.check_output(["git", "status", "--short"], cwd=root, text=True),
                "source_sha256": {str(f.relative_to(root)): hashlib.sha256(f.read_bytes()).hexdigest() for f in paths},
                "summary": summary,
                "remaining": ["WP2 complete basin/boundary search and no-tipping scope", "WP3 two-parameter continuation and bifurcations",
                              "WP4 matched rate/path families and long holds", "WP5 ecological ranges and alternative closures",
                              "seasonal periodic reference", "state expansion only after evidence"]}
    (out/"manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
