"""Explicit spinup, nonnegative guards, configurable VISIT-informed closures.

PYTHONPATH=src python examples/run_tp_v2.py --config configs/tp_v2_spinup.json
Output directory must be new; prior v1 and v2 results are never overwritten.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, replace
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import cumulative_trapezoid

from control_carbon import temperature_precipitation as v1
from control_carbon.tp_processes import ProcessOptions, process_provenance
from control_carbon.tp_experiments_v2 import Model, Solver, SpinupSettings, integrate, spinup


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")


def table(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def configure(path):
    cfg = json.loads(path.read_text())
    allowed = {"schema_version", "baseline", "target", "initial_state", "years", "sample_days",
               "solver", "spinup", "parameters", "models"}
    if set(cfg)-allowed or cfg.get("schema_version") != 2:
        raise ValueError("unknown configuration keys or schema version")
    p = v1.Parameters(**cfg.get("parameters", {}))
    baseline, target = v1.Climate(**cfg["baseline"]), v1.Climate(**cfg["target"])
    solver, setting = Solver(**cfg.get("solver", {})), SpinupSettings(**cfg.get("spinup", {}))
    initial = np.array(cfg["initial_state"], dtype=float)
    if initial.shape != (8,) or not np.all(np.isfinite(initial)) or initial.min() < 0:
        raise ValueError("initial_state must have eight finite nonnegative pools")
    for name in ("years", "sample_days"):
        if not np.isfinite(cfg[name]) or cfg[name] <= 0:
            raise ValueError(f"positive {name} required")
    models, names = [], set()
    for spec in cfg["models"]:
        if set(spec)-{"name", "mortality", "processes", "parameters", "refine_spinup"}:
            raise ValueError("unknown model option")
        name = spec["name"]
        if not re.fullmatch(r"[a-zA-Z0-9_-]+", name) or name in names:
            raise ValueError("model names must be unique safe directory names")
        names.add(name)
        parameters = replace(p, **spec.get("parameters", {}))
        if v1.domain_violation(initial, parameters) > solver.bound_tolerance:
            raise ValueError(f"initial water exceeds {name} capacity")
        model = Model(parameters, v1.Controls(mortality=spec["mortality"]), ProcessOptions(**spec.get("processes", {})))
        models.append((spec, model))
    if not models:
        raise ValueError("at least one model required")
    return cfg, initial, baseline, target, solver, setting, models


def records(model, times, states, climate):
    rows = []
    for t, x in zip(times, states):
        f = model.fluxes(x, climate)
        bc, bw = model.budgets(x, climate)
        try:
            capacity = float(np.linalg.solve(-f.carbon_matrix, f.carbon_input).sum())
        except np.linalg.LinAlgError:
            capacity = None
        rows.append({"day": float(t), "T_degC": climate.temperature, "P_mm_day": climate.precipitation,
            **dict(zip(v1.STATE_NAMES, x)), "C_total": float(x[:4].sum()), "GPP": f.gpp,
            "Ra": float(f.respiration.sum()), "Rh": f.soil_respiration, "NEP": float(f.nep),
            "C_export": f.carbon_export, "capacity": capacity,
            "potential": capacity-x[:4].sum() if capacity is not None else None,
            "mortality_L": f.mortality[0], "mortality_S": f.mortality[1], "mortality_R": f.mortality[2],
            "E_leaf": f.transpiration, "E_soil": f.evaporation, "drainage": f.drainage, "runoff": f.runoff,
            "W_export": f.water_export, "release_L": f.water_release[0], "release_S": f.water_release[1],
            "release_R": f.water_release[2], "q_OR": f.transport[0], "q_RS": f.transport[1], "q_SL": f.transport[2],
            "heat_hazard": f.heat_hazard, "drought_hazard": f.drought_hazard,
            "carbon_deficit_rate": float(max(0., -f.nep)), "C_budget": bc, "W_budget": bw})
    for src, dest in (("carbon_deficit_rate", "carbon_deficit"), ("heat_hazard", "heat_dose"), ("drought_hazard", "drought_dose")):
        for row, value in zip(rows, cumulative_trapezoid([r[src] for r in rows], times, initial=0)):
            row[dest] = float(value)
    return rows


def plots(folder, spin, results, baseline, eq):
    fig, axes = plt.subplots(2, 4, figsize=(15, 7), constrained_layout=True)
    t = (spin["times"]-spin["times"][-1])/365.25
    for k, ax in enumerate(axes.flat):
        ax.plot(t, spin["states"][:, k], color="#777777")
        ax.axhline(eq.state[k], color="black", linestyle=":", label="direct root")
        ax.set(xlabel="Years before experiment", ylabel=f"{v1.STATE_NAMES[k]} ({v1.STATE_UNITS[k]})", ylim=(0, None))
        ax.axvline(0, color="black", lw=.5)
    fig.suptitle(f"Explicit constant-climate spinup: {baseline.temperature} degC, {baseline.precipitation} mm/day; {spin['years']} years")
    fig.savefig(folder/"spinup_all_pools.png", dpi=145); plt.close(fig)
    fig, axes = plt.subplots(3, 4, figsize=(15, 10), constrained_layout=True)
    colors = {"0": "gray", "T": "#c76b22", "P": "#258bcb", "TP": "#8f39a1"}
    for key, rows in results.items():
        te = np.array([r["day"] for r in rows])/365.25
        for k, ax in enumerate(axes[:2].flat):
            keep = t >= -10
            ax.plot(t[keep], spin["states"][keep, k], color="gray", alpha=.4)
            ax.plot(te, [r[v1.STATE_NAMES[k]] for r in rows], label=key, color=colors[key])
            ax.axvline(0, color="black", lw=.5)
            ax.set(xlabel="Years (forcing changes at 0)", ylabel=f"{v1.STATE_NAMES[k]} ({v1.STATE_UNITS[k]})", ylim=(0, None))
        axes[2, 0].plot(np.r_[-10., -1e-8, te], np.r_[baseline.temperature, baseline.temperature, [r["T_degC"] for r in rows]], color=colors[key])
        axes[2, 1].plot(np.r_[-10., -1e-8, te], np.r_[baseline.precipitation, baseline.precipitation, [r["P_mm_day"] for r in rows]], color=colors[key])
        axes[2, 2].plot([r["W_O"] for r in rows], [r["C_total"] for r in rows], color=colors[key])
        axes[2, 2].scatter(rows[0]["W_O"], rows[0]["C_total"], s=12, color=colors[key])
        a, b = len(rows)//4, len(rows)//3
        axes[2, 2].annotate("", xy=(rows[b]["W_O"], rows[b]["C_total"]), xytext=(rows[a]["W_O"], rows[a]["C_total"]), arrowprops={"arrowstyle": "->", "color": colors[key]})
    contrast = np.array([r["C_total"] for r in results["TP"]])-np.array([r["C_total"] for r in results["T"]])-np.array([r["C_total"] for r in results["P"]])+np.array([r["C_total"] for r in results["0"]])
    axes[2, 3].plot(te, contrast)
    axes[2, 3].axhline(0, color="black", lw=.5)
    axes[2, 0].set(xlabel="Years", ylabel="Temperature (degC)")
    axes[2, 1].set(xlabel="Years", ylabel="Precipitation (mm/day)")
    axes[2, 2].set(xlabel="Soil water (mm)", ylabel="Total C (Mg C/ha)", title="8D trajectory projected to 2D")
    axes[2, 3].set(xlabel="Years", ylabel="I_C: signed contrast (Mg C/ha)", title="This panel is not a carbon pool")
    axes[0, 0].legend(ncol=4)
    fig.suptitle(folder.name+": explicit spinup followed by shared-initial-state factorial runs")
    fig.savefig(folder/"states_forcing_and_contrast.png", dpi=145); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("configs/tp_v2_spinup.json"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cfg, seed, baseline, target, solver, settings, models = configure(args.config)
    root = Path(__file__).resolve().parents[1]
    out = args.output or root/"artifacts/temperature_precipitation/v2"/datetime.now(timezone.utc).strftime("run_%Y%m%dT%H%M%S")
    if (root/"artifacts/temperature_precipitation/v2").resolve() not in out.resolve().parents:
        raise ValueError("v2 output must be inside artifacts/temperature_precipitation/v2")
    out.mkdir(parents=True, exist_ok=False)
    dump(out/"configuration.json", cfg)
    climates = {"0": baseline, "T": v1.Climate(target.temperature, baseline.precipitation),
                "P": v1.Climate(baseline.temperature, target.precipitation), "TP": target}
    refined_solver = replace(solver, method="BDF" if solver.method != "BDF" else "Radau",
                             max_step=solver.max_step/2, rtol=solver.rtol/5, atol=solver.atol/5)
    times = np.unique(np.r_[np.arange(0., cfg["years"]*365.25, cfg["sample_days"]), cfg["years"]*365.25])
    summary = []
    for spec, model in models:
        folder = out/spec["name"]; folder.mkdir()
        print(f"{spec['name']}: integrating spinup from user seed", flush=True)
        spun = spinup(model, seed, baseline, settings, solver)
        np.savez_compressed(folder/"spinup.npz", time_day=spun["times"]-spun["times"][-1], states=spun["states"])
        table(folder/"spinup_convergence.csv", spun["checks"])
        table(folder/"spinup.csv", records(model, spun["times"]-spun["times"][-1], spun["states"], baseline))
        dump(folder/"spinup_audit.json", {"converged": spun["converged"], "years": spun["years"], "audits": spun["audits"]})
        if not spun["converged"]:
            summary.append({"model": spec["name"], "status": "spinup_failed_no_experiment", "years": spun["years"]})
            print(spec["name"], "spinup failed; no forcing experiment", flush=True)
            continue
        initial = spun["states"][-1].copy()
        eq = model.equilibrium(baseline, initial)
        distance = float(np.max(abs(initial-eq.state)/np.maximum(1., abs(eq.state))))
        if distance > 1e-5 or eq.max_real_part >= -1e-8:
            raise RuntimeError("spinup did not match a verified stable root")
        spin_ref = None
        if spec.get("refine_spinup", False):
            other = spinup(model, seed, baseline, settings, refined_solver)
            if not other["converged"]:
                raise RuntimeError("refined spinup failed")
            spin_ref = float(np.max(abs(initial-other["states"][-1])/np.maximum(1., abs(initial))))
            if spin_ref > 1e-5:
                raise RuntimeError("spinup refinement failed")
            dump(folder/"spinup_refinement.json", {"years": other["years"], "scaled_error": spin_ref,
                                                  "solver": asdict(refined_solver), "audit": other["audits"]})
        all_rows, experiment_info = {}, {}
        for key, climate in climates.items():
            result = integrate(model, initial, times, climate, solver)
            check = integrate(model, initial, times, climate, refined_solver)
            err = float(np.max(abs(result.states-check.states)/np.maximum(1., abs(check.states))))
            if err > 1e-5:
                raise RuntimeError("experiment refinement failed")
            eq_end = model.equilibrium(climate, result.states[-1])
            rows = records(model, times, result.states, climate)
            budget = max(abs(r[n]) for r in rows for n in ("C_budget", "W_budget"))
            if budget > 1e-10:
                raise RuntimeError("mass budget failed")
            all_rows[key] = rows
            table(folder/f"{key}_trajectory.csv", rows)
            np.savez_compressed(folder/f"{key}.npz", time_day=times, states=result.states,
                refined_states=check.states, equilibrium=eq_end.state, eigenvalues=eq_end.eigenvalues,
                jacobian=model.jacobian(eq_end.state, climate), state_names=v1.STATE_NAMES)
            experiment_info[key] = {"audit": result.audit, "refined_audit": check.audit,
                "refinement": err, "budget": budget, "initial": result.states[0].tolist(),
                "equilibrium": eq_end.state.tolist(), "equilibrium_residual": eq_end.residual,
                "max_real": eq_end.max_real_part, "condition_number": eq_end.condition_number}
        interaction = [{"day": float(t), "I_C": float(all_rows["TP"][i]["C_total"]-all_rows["T"][i]["C_total"]-all_rows["P"][i]["C_total"]+all_rows["0"][i]["C_total"])} for i, t in enumerate(times)]
        table(folder/"interaction.csv", interaction)
        legacy_distance = None
        if model.processes == ProcessOptions() and model.parameters == v1.Parameters():
            legacy = root/"artifacts/temperature_precipitation/wp1_wp2_dry"
            if all((legacy/f"{model.control.mortality}_{k}.npz").exists() for k in climates):
                diffs = []
                for k in climates:
                    old = np.load(legacy/f"{model.control.mortality}_{k}.npz")
                    if np.array_equal(old["time_day"], times) and asdict(baseline) == {"temperature": 27, "precipitation": 6} and asdict(target) == {"temperature": 31, "precipitation": 2}:
                        current = np.array([[r[n] for n in v1.STATE_NAMES] for r in all_rows[k]])
                        diffs.append(float(np.max(abs(current-old["states"])/np.maximum(1., abs(old["states"])))))
                if diffs:
                    legacy_distance = max(diffs)
        info = {"name": spec["name"], "model_level": "source-informed theoretical ODE; not native VISIT",
            "parameters": v1.parameter_manifest(model.parameters), "control": asdict(model.control),
            "processes": process_provenance(model.processes), "spinup_years": spun["years"],
            "spinup_end": initial.tolist(), "root_check": eq.state.tolist(), "root_distance": distance,
            "spinup_refinement": spin_ref, "experiments": experiment_info,
            "legacy_max_scaled_difference": legacy_distance,
            "I_C_final": interaction[-1]["I_C"], "explicit_clipping_mass": 0.}
        dump(folder/"manifest.json", info)
        plots(folder, spun, all_rows, baseline, eq)
        summary.append({"model": spec["name"], "status": "passed", "years": spun["years"],
            "root_distance": distance, "legacy_difference": legacy_distance,
            "I_C_final": interaction[-1]["I_C"],
            "minimum_pool": min(min(a["minima"].values()) for a in spun["audits"]),
            "max_refinement": max(i["refinement"] for i in experiment_info.values())})
        print(json.dumps(summary[-1]), flush=True)
    dump(out/"summary.json", summary)
    source_paths = [root/"src/control_carbon/tp_experiments_v2.py", root/"src/control_carbon/tp_processes.py",
        root/"src/control_carbon/temperature_precipitation.py", Path(__file__).resolve(),
        root/"docs/temperature_precipitation_tipping_research_plan.md", args.config.resolve()]
    dump(out/"manifest.json", {"schema_version": 2, "configuration": cfg, "solver": asdict(solver),
        "refined_solver": asdict(refined_solver), "spinup": asdict(settings), "summary": summary,
        "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        "git_status": subprocess.check_output(["git", "status", "--short"], cwd=root, text=True),
        "source_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths},
        "numpy_version": np.__version__, "scipy_version": scipy.__version__,
        "dose_quadrature": "trapezoidal over output samples; spinup annual and experiment sample_days",
        "units": {"state": dict(zip(v1.STATE_NAMES, v1.STATE_UNITS)), "time": "day", "T": "degC", "P": "mm/day"}})
    print("Output:", out, flush=True)


if __name__ == "__main__":
    main()
