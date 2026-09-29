"""Reproducible analytic controls and native VISIT respiration audit.

No calibrated forest model, no forcing-rate experiment, no parameter fitting.
"""
import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from control_carbon.reserve_feedback_screen import (
    ReserveControl, positive_equilibrium, reserve_rhs, reserve_jacobian,
)
from control_carbon.temperature_nsc_visit import (
    load_tky_tree_parameters, size_dependent_respiration, SOURCE_COMMIT,
)
from control_carbon.visit_native import (
    build_native_visit_plant_respiration_bridge, run_native_visit_plant_respiration,
)
from control_carbon.visit_plant import VISITPlantStructuralState, VISITPlantRespirationParameters


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('/mnt/d/ct/VISIT-matrix/visit_local'))
    parser.add_argument('--output', type=Path,
                        default=ROOT/'artifacts/reserve_feedback_screen/evidence_20260929_v2')
    args = parser.parse_args()
    src, out = args.source.resolve(), args.output.resolve()
    commit = subprocess.check_output(['git', '-C', str(src), 'rev-parse', 'HEAD'], text=True).strip()
    if commit != SOURCE_COMMIT:
        raise ValueError('source revision mismatch')
    config = src/'INPUT/Config_TKY.xlsx'
    p = load_tky_tree_parameters(config)  # No fallback to dataclass defaults.
    out.mkdir(parents=True, exist_ok=False)
    binary = build_native_visit_plant_respiration_bridge(src, out/'respiration_bridge')
    native_p = VISITPlantRespirationParameters(
        p.rgf, p.rgc, p.rgr, p.rmf0, p.rmc_s, p.rmc_h, p.rmr_s, p.rmr_h,
        p.qtf0, p.qtc0, p.qtr0, p.f_size_stem, p.f_size_root)
    audit = []
    for t in (5., 15., 25., 35.):
        for mass in (0., 1e-8, .1, .999, 1., 1.001, 10., 40., 100., 300.):
            state = VISITPlantStructuralState(mass, mass, mass)
            native = run_native_visit_plant_respiration(
                binary, state, native_p, surface_temperature=t, upper_soil_temperature=t,
                allocation_fluxes=[0., 0., 0.], positive_epp=False).maintenance
            py = size_dependent_respiration(state.as_array(), t, p)
            np.testing.assert_allclose(py, native, rtol=2e-13, atol=1e-18)
            audit.append(dict(temperature_C=t, each_organ_MgC_ha=mass,
                              native_flux=native.tolist(), python_flux=py.tolist(),
                              max_abs_error=float(np.max(abs(native-py)))))
    slopes = []
    for mass in (.1, 2., 10., 40., 100., 300.):
        rates = size_dependent_respiration([mass]*3, 15., p)/mass
        estimates = []
        for h in (1e-4, 5e-5):
            plus, minus = mass*np.exp(h), mass*np.exp(-h)
            rp = size_dependent_respiration([plus]*3, 15., p)/plus
            rm = size_dependent_respiration([minus]*3, 15., p)/minus
            estimates.append((rp-rm)/(2*h))
        np.testing.assert_allclose(estimates[0], estimates[1], rtol=3e-5, atol=1e-13)
        if mass >= 2:
            assert np.all(estimates[1][1:] < 0)
        slopes.append(dict(mass_MgC_ha=mass, per_day=rates.tolist(),
                           derivative_wrt_log_mass=estimates[1].tolist(),
                           refinement_abs_error=float(np.max(abs(estimates[0]-estimates[1])))))
    controls = []
    for hill in (1., 2., 4., 8.):
        control = replace(ReserveControl(), hill=hill)
        eq = positive_equilibrium(control)
        j = reserve_jacobian(eq, control)
        controls.append(dict(parameters=asdict(control), equilibrium=eq.tolist(),
                             residual=reserve_rhs(eq, control).tolist(),
                             eigenvalues=np.linalg.eigvals(j).real.tolist(),
                             trace=float(np.trace(j)), determinant=float(np.linalg.det(j))))
    import openpyxl
    book = openpyxl.load_workbook(config, data_only=True, read_only=True)
    labels = {'rmf0', 'rmc_s', 'rmc_h', 'rmr_s', 'rmr_h', 'f_sz_s', 'f_sz_r',
              'qTf0', 'qTc0', 'qTr0', 'alloc_ass', 'alloc_abg'}
    cells = [{'label': str(row[3].value), 'cell': row[4].coordinate, 'value': row[4].value}
             for row in book['Parameter-daily'].iter_rows()
             if str(row[3].value) in labels]
    book.close()
    report = dict(approximation='model-agnostic derived analysis plus native-source respiration slice',
                  source_repository='Sachitama2001/VISIT-matrix', source_commit=commit,
                  scope='Synthetic component audit; NOT a TKY trajectory or R-tipping detection',
                  source_order='f_q10_ar -> f_spcfc_resp -> f_rfm/f_rcm/f_rrm; all use supplied masses',
                  flux_units='Mg C ha^-1 day^-1', stock_units='Mg C ha^-1',
                  config_cells=cells, native_comparisons=audit, size_slopes=slopes,
                  dimensionless_theory_controls=controls)
    (out/'verification.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    files = [Path(__file__), ROOT/'src/control_carbon/reserve_feedback_screen.py',
             ROOT/'tests/test_reserve_feedback_screen.py',
             ROOT/'docs/allocation_storage_respiration_evidence.md',
             ROOT/'docs/allocation_storage_respiration_validation.md',
             ROOT/'src/control_carbon/temperature_nsc_visit.py',
             ROOT/'src/control_carbon/visit_native.py', ROOT/'src/control_carbon/visit_plant.py',
             ROOT/'native/visit_plant_respiration_bridge.c', config,
             *[src/name for name in ('ecophysiology.c', 'respiration.c', 'allocation.c',
                                     'plant_proc.c', 'structure.h', 'prototype.h', 'definition.h', 'setting.h')]]
    snapshots = out/'source_snapshot'
    snapshots.mkdir()
    manifest = {}
    for i, path in enumerate(files):
        name = f'{i:02d}_{path.name}'
        shutil.copy2(path, snapshots/name)
        manifest[str(path)] = dict(sha256=hashlib.sha256(path.read_bytes()).hexdigest(), snapshot=name)
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps(dict(output=str(out), native_cases=len(audit),
                          max_native_error=max(r['max_abs_error'] for r in audit),
                          slopes=slopes, controls=controls, config_cells=cells), indent=2))


if __name__ == '__main__':
    main()
