"""Reproducible first autonomous NSC campaign; no native-model claims."""
import argparse
from dataclasses import asdict, replace
import hashlib
import itertools
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from control_carbon.temperature_nsc import Parameters, PLANT, NAMES, evaluate, concentrations
from control_carbon.nsc_analysis import (spinup, equilibrium, linear_info, arclength,
                                          integrate, smooth_temperature, bare_invasion)


def save(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')


def seed(b, chi):
    x = np.zeros(8)
    x[:3] = b*np.array([.03, .8, .17])
    x[4:7] = b*chi/(1-chi)*np.array([.1, .75, .15])
    x[[3, 7]] = [100, 1]
    return x


def diagnostics(states, temperatures, p):
    rows = []
    for x, t in zip(states, temperatures):
        _, d = evaluate(x, t, p)
        rows.append([t, d['gpp'], *d['demand'], *d['paid'], *d['unpaid'],
                     d['flush'], d['mortality'], d['growth_resp'], d['rh'],
                     d['export'], d['chi'], *d['organ_chi'], d['budget']])
    return np.array(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', default=str(ROOT/'configs/temperature_nsc_v1.json'))
    ap.add_argument('--output')
    args = ap.parse_args()
    cfg = json.loads(Path(args.config).read_text())
    out = Path(args.output or ROOT/'artifacts/temperature_nsc'/cfg['experiment_id'])
    out.mkdir(parents=True, exist_ok=False)
    p = Parameters(**cfg['parameters'])
    sources = [Path(__file__), ROOT/'src/control_carbon/temperature_nsc.py',
               ROOT/'src/control_carbon/nsc_analysis.py', Path(args.config)]
    (out/'source_snapshot').mkdir()
    for source in sources:
        shutil.copy2(source, out/'source_snapshot'/source.name)
    save(out/'manifest.json', dict(config=cfg, parameters=asdict(p),
         commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
         dirty=subprocess.check_output(['git', 'status', '--short'], cwd=ROOT, text=True),
         hashes={str(f): hashlib.sha256(f.read_bytes()).hexdigest() for f in sources},
         python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
         state_names=NAMES, unit='Mg C ha-1', time_unit='day',
         calibration='all numeric parameters assumed; not empirically constrained'))
    save(out/'bare_invasion.json', [dict(temperature=t, **bare_invasion(t,p))
                                   for t in cfg['temperatures']])
    options = dict(years=cfg['spinup_years'], chunk_years=cfg['spinup_chunk_years'],
                   max_step=cfg['max_step_days'], rtol=cfg['rtol'], atol=cfg['atol'],
                   rhs_tolerance=cfg['rhs_tolerance'], drift_tolerance=cfg['drift_tolerance'])
    seeds = {'near_bare': seed(1e-6, .03), 'low_low': seed(1, .005),
             'low_high': seed(1, .15), 'high_low': seed(100, .005),
             'high_high': seed(100, .15)}
    frozen = []
    q0 = equilibrium(seed(100, .1), cfg['ramp_temperatures'][0], p)
    if q0 is None:
        raise RuntimeError('No positive reference root')
    seeds['around_qse'] = q0*np.array([1.05,.95,1.05,.95,1.05,.95,1.05,.95])
    for temperature in cfg['temperatures']:
        for name, x in seeds.items():
            end, report = spinup(x, temperature, p, **options)
            report['seed_name'] = name
            root_x = equilibrium(end, temperature, p) if not report['near_boundary'] else None
            report['root'] = linear_info(root_x, temperature, p) if root_x is not None else None
            frozen.append(report)
            tt = np.concatenate([r['sample_days'] for r in report['records']])
            xx = np.concatenate([r['sample_states'] for r in report['records']])
            np.savez_compressed(out/f'spinup_T{temperature}_{name}.npz', times=tt, states=xx,
                                diagnostics=diagnostics(xx, np.full(len(xx), temperature), p))
            print('spinup', temperature, name, 'converged', report['converged'],
                  'near_boundary', report['near_boundary'], flush=True)
        save(out/'frozen.json', frozen)
    branches = {}
    for direction in (-1, 1):
        branches[str(direction)] = arclength(q0, cfg['ramp_temperatures'][0], p,
                                             direction=direction, steps=500)
    save(out/'branches.json', branches)
    print('continuation complete', flush=True)
    # Predeclared structural grid: roots from several seeds, NOT basin proof.
    scan = []
    keys = list(cfg['scan'])
    for values in itertools.product(*(cfg['scan'][k] for k in keys)):
        setting = dict(zip(keys, values)); changes = setting.copy()
        m = changes.pop('transfer_multiplier')
        pp = replace(p, **changes, transfer=tuple(m*v for v in p.transfer))
        roots = []
        for temperature in (20, 30):
            found = []
            for x in (seed(1,.005), seed(1,.15), seed(100,.005), seed(100,.15)):
                q = equilibrium(x, temperature, pp)
                if q is not None and not any(np.max(np.abs(np.log(q/r))) < 1e-4 for r in found):
                    found.append(q)
            roots.extend(linear_info(q, temperature, pp) for q in found)
        scan.append(dict(setting=setting, roots=roots))
        save(out/'structural_scan.json', scan)
        print('scan', setting, 'roots', len(roots), flush=True)
    # A cost-free flush violates the carbon budget: disable the entire flush
    # route instead. This is a mechanism removal, not a clean cost-only test.
    ablations = dict(no_leaf_dependence=replace(p, leaf_dependence=False),
                     no_high_T_decline=replace(p, high_temperature_decline=False),
                     no_q10=replace(p, q10=1),
                     no_flush_route=replace(p, flush_enabled=False),
                     no_nsc_mortality=replace(p, nsc_mortality=False),
                     no_transfer=replace(p, transfer=(0,0,0,0)))
    ablation_results = []
    for name, pp in ablations.items():
        end, report = spinup(seeds['high_high'], 25, pp, **options)
        report.update(name=name, parameters=asdict(pp))
        ablation_results.append(report)
        save(out/'ablations.json', ablation_results)
        print('ablation', name, flush=True)
    save(out/'ablations.json', ablation_results)
    # Resolve and check all sampled path equilibria before running rate family.
    t0, t1 = cfg['ramp_temperatures']; path = []
    q = q0
    for temperature in np.linspace(t0, t1, 101):
        q = equilibrium(q, temperature, p)
        if q is None:
            raise RuntimeError('Path lost positive root; rate experiment gated')
        info = linear_info(q, temperature, p); path.append(info)
        if info['max_real'] >= 0:
            raise RuntimeError('Path lost stability; rate experiment gated')
    save(out/'ramp_path_equilibria.json', path)
    # Start from integrated spinup, not from the root solver's returned state.
    start, spin = spinup(seeds['high_high'], t0, p, **options)
    if not spin['converged']:
        raise RuntimeError('Reference spinup did not converge')
    rates = []
    for duration in cfg['ramp_days']:
        ramp_t = np.unique(np.r_[np.arange(0, duration, 30.), duration])
        driver = lambda t: smooth_temperature(t, duration, t0, t1)
        tr = integrate(start, ramp_t, driver, p, cfg['max_step_days'], cfg['rtol'], cfg['atol'])
        hold_t = np.arange(0, cfg['hold_years']*365+1, 365.)
        hold = integrate(tr.states[-1], hold_t, lambda t: t1, p,
                         cfg['max_step_days'], cfg['rtol'], cfg['atol'])
        ref = integrate(start, ramp_t, driver, p, cfg['max_step_days']/2,
                        cfg['rtol']/10, cfg['atol']/10)
        refhold = integrate(ref.states[-1], hold_t, lambda t: t1, p,
                            cfg['max_step_days']/2, cfg['rtol']/10, cfg['atol']/10)
        diff = max(np.max(np.abs(tr.states-ref.states)/np.maximum(1,ref.states)),
                   np.max(np.abs(hold.states-refhold.states)/np.maximum(1,refhold.states)))
        xx = np.vstack([tr.states, hold.states[1:]])
        tt = np.r_[tr.times, duration+hold.times[1:]]
        dd = diagnostics(xx, [driver(t) for t in tt], p)
        np.savez_compressed(out/f'ramp_{duration}days.npz', times=tt, states=xx, diagnostics=dd)
        row = dict(duration_days=duration, final_state=hold.states[-1].tolist(),
                   scaled_distance_final_qse=float(np.max(np.abs(hold.states[-1]-q)/np.maximum(1,q))),
                   min_chi=float(dd[:,16].min()), refinement=float(diff),
                   algebraic_budget=float(np.abs(dd[:,-1]).max()),
                   integrated_budget_error=max(tr.integrated_budget_error, hold.integrated_budget_error),
                   minimum_stock=min(tr.min_stock, hold.min_stock),
                   hold_years=cfg['hold_years'],
                   hold_half_state=hold.states[len(hold.states)//2].tolist())
        if diff > cfg['refinement_tolerance'] or row['algebraic_budget'] > cfg['budget_tolerance']:
            raise RuntimeError('Rate experiment failed numerical audit')
        rates.append(row)
        print('ramp', duration, 'final distance', row['scaled_distance_final_qse'], flush=True)
        save(out/'rates.json', rates)
    save(out/'summary.json', dict(frozen_runs=len(frozen),
         converged=sum(r['converged'] for r in frozen),
         near_boundary=sum(r['near_boundary'] for r in frozen),
         scan_settings=len(scan), rate_runs=len(rates),
         conclusion='See report: finite exploration; no automatic R-tipping classification'))


if __name__ == '__main__':
    main()
