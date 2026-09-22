"""Declared follow-up controls, saddle-sided tests, and two-state diagnostics."""
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import shutil
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from control_carbon.temperature_nsc import Parameters, PLANT
from control_carbon.nsc_analysis import (equilibrium, linear_info, spinup, integrate,
                                        smooth_temperature, bare_invasion)
from control_carbon.nsc_reduction import reduced_equilibrium, ReducedParameters
from run_temperature_nsc import seed, save, diagnostics


def main():
    out=ROOT/'artifacts/temperature_nsc/controls_v1'
    out.mkdir(parents=True,exist_ok=False)
    p=replace(Parameters(),target_mode='supported')
    # Fixed before these runs: endpoint grid below the observed fold (~37.507 C).
    endpoints=[30,35,37,37.4]
    protocols=['step','smooth_fast','smooth_slow','linear','pulse','slow_return']
    sources=[Path(__file__),ROOT/'src/control_carbon/nsc_analysis.py',
             ROOT/'src/control_carbon/temperature_nsc.py',ROOT/'src/control_carbon/nsc_reduction.py',
             ROOT/'examples/run_temperature_nsc.py']
    (out/'source_snapshot').mkdir()
    for s in sources:
        shutil.copy2(s,out/'source_snapshot'/s.name)
    save(out/'manifest.json',dict(parameters=asdict(p),endpoints=endpoints,protocols=protocols,
         status='exploratory follow-up; endpoints chosen after frozen-branch calculation',
         hashes={str(s):hashlib.sha256(s.read_bytes()).hexdigest() for s in sources},
         fast_days=30,slow_days=36500,pulse_plateau_days=365,hold_years=500))
    start, report=spinup(seed(100,.15),20,p)
    save(out/'reference_spinup.json',report)
    results=[]
    for target in endpoints:
        qt=equilibrium(start,target,p)
        if qt is None or linear_info(qt,target,p)['max_real']>=0:
            raise RuntimeError('unstable endpoint')
        for protocol in protocols:
            segments=[]
            if protocol=='step':
                segments=[(500*365,lambda t,tar=target:tar)]
            elif protocol=='pulse':
                segments=[(30,lambda t,tar=target:smooth_temperature(t,30,20,tar)),
                          (365,lambda t,tar=target:tar),
                          (30,lambda t,tar=target:smooth_temperature(t,30,tar,20)),
                          (500*365,lambda t:20)]
            elif protocol=='slow_return':
                segments=[(30,lambda t,tar=target:smooth_temperature(t,30,20,tar)),
                          (500*365,lambda t,tar=target:tar),
                          (36500,lambda t,tar=target:smooth_temperature(t,36500,tar,20)),
                          (500*365,lambda t:20)]
            else:
                duration=36500 if protocol=='smooth_slow' else 30
                if protocol=='linear':
                    forcing=lambda t,tar=target:20+(tar-20)*min(t/30,1)
                else:
                    forcing=lambda t,tar=target,dur=duration:smooth_temperature(t,dur,20,tar)
                segments=[(duration,forcing),(500*365,lambda t,tar=target:tar)]
            x=start.copy(); xx=[]; tt=[]; dd=[]; offset=0; audit=0
            for duration,driver in segments:
                times=np.unique(np.r_[np.arange(0,duration,365 if duration>36500 else 30),duration])
                tr=integrate(x,times,driver,p)
                xx.extend(tr.states);tt.extend(offset+tr.times)
                dd.extend(diagnostics(tr.states,[driver(t) for t in tr.times],p))
                x=tr.states[-1];offset+=tr.times[-1]
                audit=max(audit,tr.integrated_budget_error)
                if tr.endpoint_event or tr.domain_event:
                    break
            final_temp=20 if protocol in ('pulse','slow_return') else target
            qf=equilibrium(start,final_temp,p)
            row=dict(target=target,protocol=protocol,final_state=x.tolist(),
                     distance_qse=float(np.max(np.abs(x-qf)/np.maximum(1,qf))),
                     near_boundary=tr.endpoint_event,domain_limit=tr.domain_event,
                     integrated_budget_error=audit)
            results.append(row); save(out/'controls.json',results)
            np.savez_compressed(out/f'{protocol}_{target}.npz',times=tt,states=xx,diagnostics=dd)
            print('control',target,protocol,row['distance_qse'],flush=True)
    edges=[]
    for t in (20,30):
        saddle=equilibrium(seed(1,.005),t,p)
        info=linear_info(saddle,t,p)
        w,v=np.linalg.eig(np.asarray(info['jacobian']))
        ev=v[:,np.argmax(w.real)].real
        ev/=np.max(np.abs(ev)/saddle)
        for side in (-1,1):
            x,report=spinup(saddle+side*.001*ev,t,p)
            report.update(side=side,saddle=info)
            edges.append(report)
            print('edge',t,side,report['near_boundary'],flush=True)
    save(out/'saddle_sides.json',edges)
    reduced=[]
    for t in np.linspace(0,55,111):
        q=reduced_equilibrium(t)
        reduced.append(dict(temperature=float(t),state=q.tolist() if q is not None else None))
    save(out/'reduced.json',dict(parameters=asdict(ReducedParameters()),points=reduced))


if __name__=='__main__':
    main()
