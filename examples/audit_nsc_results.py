"""Post-campaign fold, long-hold and provenance audit; saves separate outputs."""
from dataclasses import replace
import hashlib
from pathlib import Path
import sys
import numpy as np
from scipy.optimize import root
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from control_carbon.temperature_nsc import Parameters
from control_carbon.nsc_analysis import (relative_field, equilibrium, linear_info,
                                         integrate, smooth_temperature)
from run_temperature_nsc import save, seed


def main():
    out=ROOT/'artifacts/temperature_nsc/audit_v2';out.mkdir(parents=True,exist_ok=False)
    p=replace(Parameters(),target_mode='supported')
    save(out/'manifest.json',dict(purpose='independent post-campaign numerical checks',
         hold_years_total=2000,refine_max_step=15,refine_rtol=1e-9,
         hashes={str(s):hashlib.sha256(s.read_bytes()).hexdigest() for s in
                 [Path(__file__),ROOT/'src/control_carbon/nsc_analysis.py',
                  ROOT/'src/control_carbon/temperature_nsc.py']}))
    def jlog(z,t,h=1e-4):
        return np.column_stack([(relative_field(z+np.eye(8)[i]*h,t,p)
                                 -relative_field(z-np.eye(8)[i]*h,t,p))/(2*h) for i in range(8)])
    def foldfun(v):
        eig=np.linalg.eigvals(jlog(v[:8],v[8]))
        return np.r_[relative_field(v[:8],v[8],p),eig[np.argmin(abs(eig))].real]
    guess=np.r_[np.log([.0397,.6477,.051,.259,.0158,.0829,.0223,.001324]),37.507]
    fit=root(foldfun,guess,tol=1e-10)
    if not fit.success or max(abs(foldfun(fit.x)))>1e-7:
        raise RuntimeError('fold root failed')
    z,t=fit.x[:8],fit.x[8]
    u,s,vh=np.linalg.svd(jlog(z,t));right=vh[-1];left=u[:,-1]
    dt=(relative_field(z,t+1e-4,p)-relative_field(z,t-1e-4,p))/2e-4
    h=1e-3
    second=(relative_field(z+h*right,t,p)-2*relative_field(z,t,p)
            +relative_field(z-h*right,t,p))/h**2
    fold=dict(temperature=float(t),state=np.exp(z).tolist(),residual=float(max(abs(foldfun(fit.x)))),
              singular_values=s.tolist(),left_FT=float(left@dt),left_Fxx_vv=float(left@second),
              linear_info=linear_info(np.exp(z),t,p))
    save(out/'fold.json',fold)
    nomort=replace(p,nsc_mortality=False)
    roots=[equilibrium(x,20,nomort) for x in
           [np.array([.03,.8,.17,100,.001,.003,.001,1]),
            np.array([3,80,17,100,1,10,2,1])]]
    if any(q is None for q in roots):
        raise RuntimeError('no-mortality root search incomplete')
    save(out/'supported_no_nsc_mortality.json',[linear_info(q,20,nomort) for q in roots])
    holds=[]
    for name in ('step','smooth_fast','smooth_slow','linear'):
        a=np.load(ROOT/f'artifacts/temperature_nsc/controls_v1/{name}_37.4.npz')
        tr=integrate(a['states'][-1],np.arange(0,1500*365+1,365.),lambda t:37.4,p,max_step=90)
        q=equilibrium(seed(100,.15),37.4,p)
        np.savez_compressed(out/f'extended_{name}.npz',times=tr.times,states=tr.states)
        holds.append(dict(protocol=name,distance_after_2000years=float(np.max(np.abs(tr.states[-1]-q)/np.maximum(1,q))),
                          budget=tr.integrated_budget_error))
        print('extended hold',name,holds[-1],flush=True)
    save(out/'extended_holds.json',holds)
    original=np.load(ROOT/'artifacts/temperature_nsc/controls_v1/smooth_fast_37.4.npz')
    start=original['states'][0]
    a=integrate(start,[0,30],lambda t:smooth_temperature(t,30,20,37.4),p,
                max_step=15,rtol=1e-9,atol=1e-13)
    h=integrate(a.states[-1],np.arange(0,500*365+1,365.),lambda t:37.4,p,
                max_step=15,rtol=1e-9,atol=1e-13)
    distance=float(np.max(abs(h.states[-1]-original['states'][-1])/np.maximum(1,h.states[-1])))
    save(out/'near_fold_refinement.json',dict(final_scaled_difference=distance,
                                            budget=max(a.integrated_budget_error,h.integrated_budget_error)))


if __name__=='__main__':
    main()
