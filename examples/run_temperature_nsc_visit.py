"""First VISIT-informed target/GPP/respiration comparison for the 8-state model."""
import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from control_carbon.nsc_analysis import equilibrium, spinup
from control_carbon.temperature_nsc import NAMES, Parameters
from control_carbon.temperature_nsc_visit import (
    FixedEnvironment, SOURCE_COMMIT, SOURCE_REPOSITORY, evaluate_variant,
    leaf_area, load_tky_tree_parameters, optimum_leaf_carbon,
)


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def seed(biomass, chi):
    x = np.zeros(8)
    x[:3] = biomass * np.array([0.03, 0.80, 0.17])
    x[4:7] = biomass * chi / (1-chi) * np.array([0.10, 0.75, 0.15])
    x[[3, 7]] = [100, 1]
    return x


def finite_linear_info(x, temperature, p, evaluator):
    j = np.empty((8, 8))
    for i in range(8):
        h = min(1e-5*max(1.0, x[i]), 0.1*x[i])
        d = np.eye(8)[i]*h
        j[:, i] = (evaluator(x+d, temperature, p)[0]
                   - evaluator(x-d, temperature, p)[0])/(2*h)
    eig = np.linalg.eigvals(j)
    residual = evaluator(x, temperature, p)[0]
    return dict(state=x.tolist(), max_real=float(eig.real.max()),
                eigen_real=eig.real.tolist(), eigen_imag=eig.imag.tolist(),
                residual=float(np.max(np.abs(residual))))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=str(ROOT / "configs/temperature_nsc_visit_v1.json"))
    parser.add_argument("--visit-config", default="/mnt/d/ct/VISIT-matrix/visit_local/INPUT/Config_TKY.xlsx")
    parser.add_argument("--output")
    args = parser.parse_args()
    config_path = Path(args.config)
    cfg = json.loads(config_path.read_text())
    out = Path(args.output or ROOT / "artifacts/temperature_nsc_visit_allocation" / cfg["experiment_id"])
    out.mkdir(parents=True, exist_ok=False)

    visit_config = Path(args.visit_config)
    visit = load_tky_tree_parameters(visit_config)
    env = FixedEnvironment(**cfg["fixed_environment"])
    base = replace(Parameters(), nsc_mortality=False, mu_max=0.0)
    sources = [Path(__file__), config_path, ROOT / "src/control_carbon/temperature_nsc_visit.py",
               ROOT / "src/control_carbon/nsc_analysis.py", visit_config]
    snapshot = out / "source_snapshot"
    snapshot.mkdir()
    for source in sources[:-1]:
        shutil.copy2(source, snapshot / source.name)

    save(out / "manifest.json", dict(
        status="exploratory first comparison; no R-tipping classification",
        phase="WP-B and first part of WP-C (V0 through V3)",
        config=cfg, base_parameters=asdict(base), visit_tree_parameters=asdict(visit),
        fixed_environment=asdict(env), state_names=NAMES,
        stock_unit="Mg C ha-1", time_unit="day",
        visit_source_repository=SOURCE_REPOSITORY, visit_source_commit=SOURCE_COMMIT,
        visit_config=str(visit_config),
        fixed_environment_provenance="explicit experiment assumptions, not Config_TKY.xlsx",
        commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        dirty=subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True),
        hashes={str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
        python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__))

    targets = []
    for temperature in cfg["temperatures_c"]:
        c_leaf, lai, status = optimum_leaf_carbon(temperature, visit, env)
        targets.append(dict(temperature_c=temperature, leaf_carbon=c_leaf, lai=lai,
                            status=status))
    save(out / "visit_targets.json", targets)

    seeds = {"near_bare": seed(1e-6, 0.03), "low_low": seed(1, 0.005),
             "high_low": seed(100, 0.005), "high_high": seed(100, 0.15)}
    options = dict(years=cfg["spinup_years"], chunk_years=cfg["spinup_chunk_years"],
                   max_step=cfg["max_step_days"], rtol=cfg["rtol"], atol=cfg["atol"],
                   rhs_tolerance=cfg["rhs_tolerance"], drift_tolerance=cfg["drift_tolerance"])
    results = []
    for variant in cfg["variants"]:
        def evaluator(x, temperature, p, selected=variant):
            return evaluate_variant(x, temperature, p, visit, env, selected,
                                    cfg["reference_temperature_c"])
        for temperature in cfg["temperatures_c"]:
            for seed_name, initial in seeds.items():
                end, report = spinup(initial, temperature, base, evaluator=evaluator, **options)
                root = None if report["near_boundary"] else equilibrium(
                    end, temperature, base, evaluator=evaluator)
                row = dict(variant=variant, temperature_c=temperature, seed=seed_name,
                           converged=report["converged"], near_boundary=report["near_boundary"],
                           domain_limit=report["domain_limit"], endpoint=end.tolist(),
                           final_residual=report["records"][-1]["residual"],
                           annual_drift=report["records"][-1]["annual_drift"],
                           minimum_stock=min(r["min_stock"] for r in report["records"]),
                           max_integrated_budget_error=max(r["integrated_budget_error"]
                                                          for r in report["records"]),
                           root=(finite_linear_info(root, temperature, base, evaluator)
                                 if root is not None else None))
                results.append(row)
                save(out / "spinup_summary.json", results)
                print(variant, temperature, seed_name, row["converged"],
                      row["near_boundary"], flush=True)
    save(out / "summary.json", dict(runs=len(results),
        converged=sum(row["converged"] for row in results),
        near_boundary=sum(row["near_boundary"] for row in results),
        stable_root_solutions=sum(row["root"] is not None and row["root"]["max_real"] < 0
                                  for row in results),
        distinct_variant_temperature_conditions=len(cfg["variants"])*len(cfg["temperatures_c"]),
        scope="V0-V3 only; V4-V7 and rate experiments remain gated follow-up"))


if __name__ == "__main__":
    main()
