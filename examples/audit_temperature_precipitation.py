"""Read saved WP1/WP2 results, including failed fixed-water interventions."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np


def main():
    root = Path(__file__).resolve().parents[1]
    for folder in ("wp1_wp2", "wp1_wp2_dry"):
        out = root/"artifacts/temperature_precipitation"/folder
        manifest = json.loads((out/"manifest.json").read_text())
        for filename, digest in manifest["source_sha256"].items():
            if hashlib.sha256((root/filename).read_bytes()).hexdigest() != digest:
                raise RuntimeError(f"source changed since experiment: {filename}")
        audit = {"source_hashes_match": True, "fixed_absolute_water_interventions": {}}
        for stage in manifest["summary"]:
            filename = out/f"{stage}_fixed_water_intervention.csv"
            rows = list(csv.DictReader(filename.open(encoding="utf-8-sig")))
            invalid = [r for r in rows if float(r["max_domain_violation"]) > 1e-6]
            audit["fixed_absolute_water_interventions"][stage] = {
                "admissible_for_full_duration": not invalid,
                "first_invalid_saved_day": float(invalid[0]["day"]) if invalid else None,
                "max_excess_storage_mm": max(float(r["max_domain_violation"]) for r in rows),
                "interpretation": "reject ecological attribution after capacity violation; keep as failed intervention diagnostic",
            }
            for key in manifest["forcing"]:
                z = np.load(out/f"{stage}_{key}.npz")
                if not np.all(np.isfinite(z["states"])) or z["states"].shape[1] != 8:
                    raise RuntimeError("invalid state archive")
                if np.max(z["eigenvalues"].real) >= -1e-8:
                    raise RuntimeError("unstable sampled endpoint")
                if not np.allclose(np.sort_complex(np.linalg.eigvals(z["jacobian"])),
                                   np.sort_complex(z["eigenvalues"]), atol=1e-8):
                    raise RuntimeError("saved Jacobian and eigenvalues disagree")
        (out/"audit.json").write_text(json.dumps(audit, indent=2)+"\n")
        print(folder, json.dumps(audit), flush=True)


if __name__ == "__main__":
    main()
