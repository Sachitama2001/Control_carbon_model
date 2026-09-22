"""WP2 check: distinguish a 20-year interaction from its frozen equilibrium.

Consumes the saved dry-comparison manifest; does not introduce a rate family.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from control_carbon.temperature_precipitation import Climate, Controls, Parameters, integrate
from run_temperature_precipitation import diagnostic_series, write_csv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("artifacts/temperature_precipitation/wp1_wp2_dry/manifest.json"))
    parser.add_argument("--output", type=Path, default=Path("artifacts/temperature_precipitation/long_hold"))
    args = parser.parse_args()
    parent = json.loads(args.manifest.read_text())
    p = Parameters(**{k: v["value"] for k, v in parent["parameters"].items()})
    control = Controls(mortality="water")
    initial = np.array(parent["summary"]["water"]["baseline_state"])
    times = np.unique(np.r_[np.arange(0, 200*365.25, 10.), 20*365.25, 200*365.25])
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    series, errors, residuals = {}, {}, {}
    for key, spec in parent["forcing"].items():
        climate = Climate(**spec)
        states = integrate(initial, times, climate, p, control, **parent["solver"])
        series[key] = diagnostic_series(times, states, climate, p, control)
        write_csv(out/f"water_{key}_200years.csv", series[key])
        np.savez_compressed(out/f"water_{key}_200years.npz", time_day=times, states=states)
        eq = np.load(args.manifest.parent/f"water_{key}.npz")["equilibrium"]
        errors[key] = float(np.max(np.abs(states[-1]-eq)/np.maximum(1., np.abs(eq))))
        residuals[key] = max(abs(r[n]) for r in series[key]
                             for n in ("carbon_budget_residual", "water_budget_residual"))
    totals = {k: np.array([r["C_total"] for r in rows]) for k, rows in series.items()}
    interaction = totals["TP"]-totals["T"]-totals["P"]+totals["0"]
    crossing = np.flatnonzero((interaction[:-1] < 0) & (interaction[1:] >= 0))
    # Ignore numerical sign noise while the physical contrast is initially zero.
    crossing = [i for i in crossing if np.min(interaction[:i+1]) < -1e-3]
    brackets = [[times[i]/365.25, times[i+1]/365.25] for i in crossing]
    summary = {"stage": "water", "years": 200, "protocol": "constant post-step climate; not a rate experiment",
               "initial_state": initial.tolist(), "parent_manifest": str(args.manifest),
               "parent_manifest_sha256": hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
               "solver": parent["solver"], "diagnostic_step_days": 10,
               "dose_quadrature": "trapezoidal on saved 10-day samples",
               "final_scaled_distances_to_equilibrium": errors, "max_budget_residuals": residuals,
               "I_C_20_years": float(interaction[np.where(times == 20*365.25)[0][0]]),
               "I_C_200_years": float(interaction[-1]),
               "I_C_equilibrium": parent["summary"]["water"]["equilibrium_I_C"],
               "negative_to_positive_crossing_year_brackets": brackets}
    # Cross-solver refinement specifically for the long combined trajectory.
    ref = integrate(initial, times, Climate(**parent["forcing"]["TP"]), p, control, **parent["refined_solver"])
    main_states = np.array([[r[n] for n in ("C_L", "C_S", "C_R", "C_O", "W_L", "W_S", "W_R", "W_O")] for r in series["TP"]])
    summary["TP_refinement_scaled_error"] = float(np.max(abs(main_states-ref)/np.maximum(1, abs(ref))))
    if max(errors.values()) > 1e-3 or summary["TP_refinement_scaled_error"] > 1e-5:
        raise RuntimeError(f"long-hold checks failed: {summary}")
    root = Path(__file__).resolve().parents[1]
    paths = [Path(__file__).resolve(), root/"examples/run_temperature_precipitation.py",
             root/"src/control_carbon/temperature_precipitation.py"]
    summary["source_sha256"] = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    summary["git_commit"] = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    (out/"manifest.json").write_text(json.dumps(summary, indent=2)+"\n")
    write_csv(out/"interaction.csv", [{"day": t, "I_C": y} for t, y in zip(times, interaction)])
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), constrained_layout=True)
    for key, total in totals.items():
        axes[0].plot(times/365.25, total, label=key)
        axes[2].plot([r["W_O"] for r in series[key]], total, label=key)
    axes[0].set(xlabel="Time (years)", ylabel="Total C (Mg C/ha)")
    axes[0].legend()
    axes[1].plot(times/365.25, interaction)
    axes[1].axhline(summary["I_C_equilibrium"], linestyle=":", color="black", label="Frozen equilibrium")
    axes[1].axhline(0, color="black", linewidth=.5)
    axes[1].axvline(20, color="gray", linestyle="--", label="20-year checkpoint")
    axes[1].set(xlabel="Time (years)", ylabel="I_C (Mg C/ha)")
    axes[1].legend(fontsize=8)
    axes[2].set(xlabel="Soil water (mm)", ylabel="Total C (Mg C/ha)", title="8D trajectory projected to 2D")
    fig.suptitle("Water mortality: long transient interaction changes sign; no rate comparison")
    fig.savefig(out/"water_interaction_long_hold.png", dpi=160)
    fig.savefig(out/"water_interaction_long_hold.pdf")
    plt.close(fig)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
