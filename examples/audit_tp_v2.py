"""Audit saved v2 results and draw stocks separately from signed contrasts."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from control_carbon.temperature_precipitation import STATE_NAMES, STATE_UNITS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, default=Path("artifacts/temperature_precipitation/v2/spinup_comparison"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    top = json.loads((args.run/"manifest.json").read_text())
    for path, expected in top["source_sha256"].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == expected, path
    legacy = root/"artifacts/temperature_precipitation/wp1_wp2/manifest.json"
    original = json.loads(legacy.read_text())
    for path, expected in original["source_sha256"].items():
        file = root/path
        if path == "docs/temperature_precipitation_tipping_research_plan.md":
            file = root/"docs/temperature_precipitation_tipping_research_plan_v1.md"
        assert hashlib.sha256(file.read_bytes()).hexdigest() == expected, path
    rows, comparisons = [], []
    for result in top["summary"]:
        if result["status"] != "passed":
            continue
        folder = args.run/result["model"]
        info = json.loads((folder/"manifest.json").read_text())
        spin = np.load(folder/"spinup.npz")
        spun = json.loads((folder/"spinup_audit.json").read_text())
        groups = {"spinup": spun["audits"]}
        groups.update({k: [i["audit"], i["refined_audit"]] for k, i in info["experiments"].items()})
        for phase, audits in groups.items():
            minima = {n: min(a["minima"][n] for a in audits) for n in STATE_NAMES}
            assert min(minima.values()) >= 0
            assert all(a["mass_added_by_clipping"] == 0 for a in audits)
            rows.append({"model": result["model"], "phase": phase, **minima,
                         "inspected_points": sum(a["inspected_points"] for a in audits),
                         "capacity_excess": max(a["max_capacity_excess"] for a in audits)})
        for key, info_exp in info["experiments"].items():
            z = np.load(folder/f"{key}.npz")
            assert np.array_equal(z["states"][0], spin["states"][-1])
            assert np.all(np.isfinite(z["states"])) and z["states"].min() >= 0
            assert info_exp["refinement"] < 1e-5 and info_exp["budget"] < 1e-10
        comparisons.append({"model": result["model"], "spinup_years": info["spinup_years"],
            "C0": sum(info["spinup_end"][:4]), "I_C_20years": info["I_C_final"],
            "root_distance": info["root_distance"], "legacy_difference": info["legacy_max_scaled_difference"]})
        # A zero origin PLUS headroom makes the nonnegative baseline visible.
        fig, axes = plt.subplots(2, 4, figsize=(15, 7), constrained_layout=True)
        colors = {"0": "#777777", "T": "#c46b22", "P": "#208bcc", "TP": "#923da3"}
        years = spin["time_day"]/365.25
        keep = years >= -10
        for k, ax in enumerate(axes.flat):
            ymax = float(spin["states"][keep, k].max())
            ax.plot(years[keep], spin["states"][keep, k], color="#555555", linestyle=":", label="spinup")
            for key in info["experiments"]:
                z = np.load(folder/f"{key}.npz")
                ax.plot(z["time_day"]/365.25, z["states"][:, k], color=colors[key], label=key)
                ymax = max(ymax, float(z["states"][:, k].max()))
            ax.set(xlabel="Years; climate changes at 0", ylabel=f"{STATE_NAMES[k]} ({STATE_UNITS[k]})", ylim=(0, max(1e-5, ymax*1.12)))
            ax.axvline(0, color="black", linewidth=.5)
        axes[0, 0].legend(ncol=3, fontsize=8)
        fig.suptitle(result["model"]+": all eight pools are nonnegative (I_C is a separate diagnostic)")
        fig.savefig(folder/"nonnegative_pools.png", dpi=150); plt.close(fig)
    for filename, data in (("nonnegative_minima.csv", rows), ("model_comparison.csv", comparisons)):
        with (args.run/filename).open("w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0])); w.writeheader(); w.writerows(data)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for result in top["summary"]:
        folder = args.run/result["model"]
        if result["status"] != "passed":
            continue
        info = json.loads((folder/"manifest.json").read_text())
        if info["control"]["mortality"] != "water":
            continue
        z = np.load(folder/"TP.npz")
        axes[0].plot(z["time_day"]/365.25, z["states"][:, :4].sum(axis=1), label=result["model"])
        data = list(csv.DictReader((folder/"interaction.csv").open(encoding="utf-8-sig")))
        axes[1].plot([float(r["day"])/365.25 for r in data], [float(r["I_C"]) for r in data], label=result["model"])
    axes[0].set(xlabel="Years since climate change", ylabel="Total C (Mg C/ha)", ylim=(0, None))
    axes[1].set(xlabel="Years since climate change", ylabel="I_C: signed contrast (Mg C/ha)")
    axes[1].axhline(0, color="black", linewidth=.5)
    axes[0].legend(fontsize=7)
    fig.suptitle("Function comparison after each model's own spinup; all parameters uncalibrated")
    fig.savefig(args.run/"function_comparison.png", dpi=160)
    fig.savefig(args.run/"function_comparison.pdf"); plt.close(fig)
    native = root.parent/"VISIT-matrix/visit_local"
    audit = {"v2_source_hashes_match": True, "v1_source_hashes_match_using_archived_plan": True,
             "passed_models": len(comparisons), "nonnegative": True, "clipping_mass_added": 0,
             "inspected_points": sum(r["inspected_points"] for r in rows),
             "minima": {n: min(r[n] for r in rows) for n in STATE_NAMES},
             "max_capacity_excess": max(r["capacity_excess"] for r in rows),
             "native_source_sha256": {n: hashlib.sha256((native/n).read_bytes()).hexdigest()
                                      for n in ("photosynthesis.c", "hydro_balance.c", "definition.h")},
             "plot_audit_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (args.run/"audit.json").write_text(json.dumps(audit, indent=2)+"\n")
    print(json.dumps(audit, indent=2))
    print(json.dumps(comparisons, indent=2))


if __name__ == "__main__":
    main()
