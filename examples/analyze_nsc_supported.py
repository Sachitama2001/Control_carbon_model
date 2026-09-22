"""Follow-up prompted by v1 roots: supported canopy, branches and basins.

Exploratory follow-up, not part of the initial preregistered rate experiment.
All numerical parameters unchanged except target_mode. Saves independent data.
"""
import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import sys
import shutil
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from control_carbon.temperature_nsc import Parameters, PLANT
from control_carbon.nsc_analysis import (equilibrium, linear_info, arclength, spinup,
                                        integrate, smooth_temperature, bare_invasion)
from run_temperature_nsc import seed, save, diagnostics


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', default=str(ROOT/'artifacts/temperature_nsc/supported_followup_v1'))
    args = ap.parse_args()
    out = Path(args.output); out.mkdir(parents=True, exist_ok=False)
    p = replace(Parameters(), target_mode='supported')
    sources = [Path(__file__), ROOT/'src/control_carbon/nsc_analysis.py',
               ROOT/'src/control_carbon/temperature_nsc.py', ROOT/'examples/run_temperature_nsc.py']
    (out/'source_snapshot').mkdir()
    for s in sources:
        shutil.copy2(s, out/'source_snapshot'/s.name)
    save(out/'manifest.json', dict(parameters=asdict(p),
        status='exploratory after initial scan, not confirmatory preregistration',
        source_hashes={str(s):hashlib.sha256(s.read_bytes()).hexdigest() for s in sources},
        basin_b=[.001,.01,.1,1,10,100], basin_chi=[.005,.15], temperatures=[20,30],
        ramp_t0=20, ramp_t1=30, durations=[30,365,3650,36500], hold_years=500))
    q = equilibrium(seed(100,.15),20,p)
    saddle = equilibrium(seed(1,.005),20,p)
    branches = {name:arclength(x,20,p,direction=direction,steps=650)
                for name,x,direction in [('stable_up',q,1),('saddle_up',saddle,1)]}
    save(out/'branches.json', branches)
    print('supported continuation done',flush=True)
    results=[]
    for t in (20,30):
        stable=equilibrium(seed(100,.15),t,p)
        for b in (.001,.01,.1,1,10,100):
            for chi in (.005,.15):
                x, report=spinup(seed(b,chi),t,p)
                report.update(b=b,chi=chi,
                              distance_positive=float(np.max(np.abs(x-stable)/np.maximum(1,stable))))
                results.append(report)
                save(out/'basins.json',results)
                print('basin',t,b,chi,'nearbare',report['near_boundary'],flush=True)
    rates=[]
    q1=equilibrium(q,30,p)
    path=[]
    for t in np.linspace(20,30,101):
        qt=equilibrium(q,t,p)
        if qt is None or linear_info(qt,t,p)['max_real']>=0:
            raise RuntimeError('stable path gate failed')
        path.append(linear_info(qt,t,p))
    save(out/'path.json',path)
    start, report=spinup(seed(100,.15),20,p)
    save(out/'reference_spinup.json',report)
    for duration in (30,365,3650,36500):
        driver=lambda t:smooth_temperature(t,duration,20,30)
        tt=np.unique(np.r_[np.arange(0,duration,30.),duration])
        a=integrate(start,tt,driver,p)
        ht=np.arange(0,500*365+1,365.)
        h=integrate(a.states[-1],ht,lambda t:30,p)
        ar=integrate(start,tt,driver,p,max_step=15,rtol=1e-9,atol=1e-13)
        hr=integrate(ar.states[-1],ht,lambda t:30,p,max_step=15,rtol=1e-9,atol=1e-13)
        diff=max(np.max(np.abs(a.states-ar.states)/np.maximum(1,ar.states)),
                 np.max(np.abs(h.states-hr.states)/np.maximum(1,hr.states)))
        states=np.vstack([a.states,h.states[1:]])
        times=np.r_[a.times,duration+h.times[1:]]
        np.savez_compressed(out/f'ramp_{duration}days.npz',times=times,states=states,
                            diagnostics=diagnostics(states,[driver(t) for t in times],p))
        rates.append(dict(duration=duration,distance_positive=float(np.max(np.abs(h.states[-1]-q1)/np.maximum(1,q1))),
                          refinement=float(diff),budget=max(a.integrated_budget_error,h.integrated_budget_error),
                          final_state=h.states[-1].tolist()))
        save(out/'rates.json',rates)
        print('supported ramp',duration,rates[-1]['distance_positive'],flush=True)
    save(out/'bare_invasion.json',[dict(temperature=t,**bare_invasion(t,p)) for t in (20,30)])


if __name__=='__main__':
    main()
