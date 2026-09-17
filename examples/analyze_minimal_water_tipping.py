"""Reproduce the minimal R-tipping research note (not a calibrated forest run)."""
from pathlib import Path
import argparse
import json
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))
from control_carbon.minimal_water_tipping import (
    MinimalTipping, finite_critical_rate, simulate_ramp, simulate_drying, scalar_rhs,
)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('artifacts/minimal_water_tipping'))
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    model=MinimalTipping()
    rc=finite_critical_rate(model)
    cases={label:simulate_ramp(rate,model,hold_years=500.)
           for label,rate in [('slow',.003),('fast',.01)]}
    critical_updates={}
    for interval in [1.,1/12]:
        critical_updates[str(interval)]=brentq(
            lambda r:simulate_ramp(r,model,update_years=interval,hold_years=0)['margin'],
            .9*rc,1.1*rc,xtol=1e-12)
    dynamic={str(T):simulate_drying(T) for T in [100.,1000.]}
    smooth_dynamic={str(T):simulate_drying(T,profile="smoothstep") for T in [10.,100.]}
    smooth_scalar={}
    for r in [.003,.01]:
        duration=np.log(model.final_scale)/r
        def smooth_rhs(t,y):
            s=t/duration
            scale=np.exp(np.log(model.final_scale)*(3*s*s-2*s*s*s))
            return [scalar_rhs(y[0],scale,model)]
        result=solve_ivp(smooth_rhs,[0,duration],[model.roots[1]],
                         rtol=1e-11,atol=1e-13,first_step=.01,max_step=1.)
        if not result.success: raise RuntimeError(result.message)
        smooth_scalar[str(r)]=float(result.y[0,-1]/model.final_scale-model.roots[0])
    # Matched quasi-steady, rainfall-only control with identical rainfall endpoints.
    quasi={}
    for T in [100.,1000.]:
        m=1.8
        c0=(10+np.sqrt(96))/2
        def rhs(t,z):
            c=np.exp(z[0]); a=10+(2.002-10)*min(t/T,1)
            return [m*(a*c/(1+c*c)-1)]
        sol=solve_ivp(rhs,[0,T+150],[np.log(c0)],rtol=1e-10,atol=1e-12,
                       max_step=min(T/100,1.))
        if not sol.success: raise RuntimeError(sol.message)
        quasi[str(T)]=float(np.exp(sol.y[0,-1]))
    # Baseline uptake exponent one: no positive density facilitation or bistability.
    # This is an analytic control documented in the TeX, not a tuned simulation.
    summary={
        'model':'quasi-steady homogeneous Klausmeier-type analytical benchmark',
        'calibrated':False,
        'parameters':vars(model),
        'roots':list(model.roots),
        'comoving_fold_rate_per_year':model.comoving_fold_rate,
        'finite_critical_rate_per_year':rc,
        'critical_duration_years':float(np.log(model.final_scale)/rc),
        'monthly_annual_flow_critical_rates':critical_updates,
        'scalar_cases':{k:{'ramp_duration_years':v['stop'], 'end_ramp_basin_margin':v['margin'],
                           'final_carbon_normalized':float(v['carbon'][-1]),
                           'outcome':v['outcome'],
                           # Rate comparisons have matched endpoints, NOT equal exposure.
                           'ramp_integrated_P_over_P0_years':(model.final_scale-1)/r,
                           'ramp_integrated_e_over_e0_years':(model.final_scale**2-1)/(2*r)}
                        for (k,v),r in zip(cases.items(),[.003,.01])},
        'dynamic_rainfall_only':{k:{'final_state':v['state'].tolist(),
                                   'bare_basin_rectangle':v['certified_bare'],
                                   'quasi_steady_final_carbon':quasi[k]}
                                for k,v in dynamic.items()},
        'dynamic_parameters':{'rainfall_start':10.,'rainfall_end':2.002,
                              'mortality':1.8,'water_loss':1.,'hold':150.},
        'smooth_log_scale_scalar_basin_margins':smooth_scalar,
        'smoothstep_rainfall_dynamic':{k:{'final_state':v['state'].tolist(),
                                           'bare_basin_rectangle':v['certified_bare']}
                                       for k,v in smooth_dynamic.items()},
        'provenance':{'visita_commit':'5c513196f21c1b1efb9ead540e3d40865b15e07e',
                      'visita_coefficients_used':False,
                      'water_uptake_family':'https://doi.org/10.1086/701669 Appendix A1/A6',
                      'forcing_path':'proposed analytical test P=P0*lambda, e=e0*lambda**2'},
    }
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(10.2,3.6),constrained_layout=True)
    colors={'slow':'#1565c0','fast':'#c62828'}
    for label,result in cases.items():
        axes[0].plot(result['time'],result['carbon'],label=label,color=colors[label])
        mask=result['time']<=result['stop']+1e-9
        axes[1].plot(result['scale'][mask],result['carbon'][mask],label=label,color=colors[label])
        indices=np.flatnonzero(mask)
        j=indices[len(indices)//2]
        j2=indices[min(len(indices)//2+max(1,len(indices)//20),len(indices)-1)]
        axes[1].annotate('',xy=(result['scale'][j2],result['carbon'][j2]),
                         xytext=(result['scale'][j],result['carbon'][j]),
                         arrowprops={'arrowstyle':'->','color':colors[label]})
    lam=np.linspace(1,2,200)
    axes[1].plot(lam,model.roots[1]*lam,'k--',label='stable forest')
    axes[1].plot(lam,model.roots[0]*lam,'k:',label='basin boundary')
    axes[0].set(xlabel='Time (years)',ylabel='Live carbon C / C0',xlim=(0,500))
    axes[1].set(xlabel='Environmental path scale lambda',ylabel='Live carbon C / C0')
    for ax in axes: ax.legend(fontsize=8); ax.grid(alpha=.2)
    fig.savefig(args.output/'scalar_rate_tipping.pdf')
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10.2,3.5),constrained_layout=True)
    for (T,result),col in zip(dynamic.items(),['#c62828','#1565c0']):
        duration=float(T)
        t=np.linspace(duration,duration+25,500)
        z=result['after'].sol(t)
        c=np.exp(z[0]); w=z[1]
        axes[0].plot(t-duration,c,label=f'ramp duration {T}',color=col)
        axes[1].plot(c,w,color=col)
        z_marks=result['after'].sol(duration+np.array([0.,5.,10.]))
        axes[1].scatter(np.exp(z_marks[0]),z_marks[1],s=12,color=col)
        axes[1].annotate('',xy=(c[100],w[100]),xytext=(c[70],w[70]),
                         arrowprops={'arrowstyle':'->','color':col})
    lo=(2.002-np.sqrt(2.002**2-4))/2
    hi=1/lo
    axes[0].axhline(hi,color='k',ls='--',lw=1)
    axes[1].scatter([lo,hi],[1/lo,1/hi],c=['white','black'],edgecolors='black',zorder=4)
    axes[1].fill_between([0,1/2.002],0,2.002,color='0.8',alpha=.5,
                         label='sufficient bare-basin region (w < a)')
    axes[0].set(xlabel='Time since rainfall ramp ended (water-loss time units)',ylabel='Live carbon c')
    axes[1].set(xlabel='Live carbon c',ylabel='Normalized water w',xlim=(0,1.4),ylim=(.65,2.05))
    axes[0].legend(fontsize=8); axes[1].legend(fontsize=7)
    for ax in axes: ax.grid(alpha=.2)
    fig.savefig(args.output/'dynamic_water_tipping.pdf')
    plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__': main()
