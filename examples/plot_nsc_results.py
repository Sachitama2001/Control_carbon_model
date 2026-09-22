"""Figures and compact numeric inventory from saved NSC experiments."""
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from control_carbon.nsc_reduction import reduced_rhs, reduced_equilibrium
from run_temperature_nsc import save


def read(path):
    return json.loads(path.read_text())


def main():
    root=ROOT/'artifacts/temperature_nsc'
    out=root/'report_v1';figdir=out/'figures';figdir.mkdir(parents=True,exist_ok=True)
    base=root/'autonomous_v1_final';supported=root/'supported_followup_v1'
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,2,figsize=(10,4),sharey=True)
    for ax,path,title in zip(axs,[base,supported],['Constant leaf target','Support-dependent leaf target']):
        branches=read(path/'branches.json')
        for branch in branches.values():
            rows=branch['points'];t=np.array([r['temperature'] for r in rows])
            b=np.array([sum(r['state'][:3]) for r in rows]);stable=np.array([r['max_real']<0 for r in rows])
            ax.scatter(t[stable],b[stable],s=3,c='tab:blue')
            ax.scatter(t[~stable],b[~stable],s=3,c='tab:red')
        ax.set(xlabel='Temperature (C)',title=title,yscale='log',xlim=(5,45),ylim=(.005,300))
    axs[0].set_ylabel('Plant structural C (Mg C / ha)')
    fig.tight_layout();fig.savefig(figdir/'branches.pdf');fig.savefig(figdir/'branches.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10,4))
    for ax,path,title in zip(axs,[base,supported],['Constant target: 20 to 25 C','Supported target: 20 to 30 C']):
        for duration in [30,365,3650,36500]:
            a=np.load(path/f'ramp_{duration}days.npz');tt=a['times']/365
            b=a['states'][:,:3].sum(axis=1)
            ax.plot(tt,b,label=f'{duration/365:.3g} yr ramp')
        ax.set(xlabel='Time since ramp start (yr)',ylabel='Plant structural C (Mg C / ha)',title=title,xlim=(0,160))
        ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(figdir/'rates.pdf');fig.savefig(figdir/'rates.png',dpi=160);plt.close(fig)
    basins=read(supported/'basins.json')
    fig,axs=plt.subplots(1,2,figsize=(9,3.5),sharey=True)
    for ax,t in zip(axs,[20,30]):
        for r in basins:
            if r['temperature']==t:
                ax.scatter(r['b'],r['chi'],s=100,c='tab:red' if r['near_boundary'] else 'tab:blue',
                           marker='x' if r['near_boundary'] else 'o')
        ax.set(xscale='log',xlabel='Initial structural C (Mg C / ha)',title=f'{t} C')
    axs[0].set_ylabel('Initial plant NSC fraction')
    fig.suptitle('Supported canopy: blue = positive QSE, red = near-bare approach')
    fig.tight_layout();fig.savefig(figdir/'basins.pdf');plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,4))
    bb,ss=np.meshgrid(np.geomspace(.01,1000,160),np.geomspace(.0001,100,160))
    f=np.array([reduced_rhs([b,s],20) for b,s in zip(bb.ravel(),ss.ravel())]).reshape(160,160,2)
    ax.contour(bb,ss,f[:,:,0],levels=[0],colors=['tab:blue'])
    ax.contour(bb,ss,f[:,:,1],levels=[0],colors=['tab:red'])
    q=reduced_equilibrium(20);ax.scatter(*q,color='black')
    ax.set(xscale='log',yscale='log',xlabel='B (Mg C / ha)',ylabel='S (Mg C / ha)',
           title='Auxiliary system: dB/dt=0 (blue), dS/dt=0 (red)')
    fig.tight_layout();fig.savefig(figdir/'nullclines.pdf');plt.close(fig)
    frozen=read(base/'frozen.json');rates=read(base/'rates.json');srates=read(supported/'rates.json')
    scan=read(base/'structural_scan.json')
    summary=dict(frozen_runs=len(frozen),frozen_converged=sum(r['converged'] for r in frozen),
        spinup_year_range=[min(r['records'][-1]['year'] for r in frozen),max(r['records'][-1]['year'] for r in frozen)],
        max_spinup_budget=max(z['integrated_budget_error'] for r in frozen for z in r['records']),
        max_rate_refinement=max(r['refinement'] for r in rates+srates),
        max_rate_algebraic_budget=max(r['algebraic_budget'] for r in rates),
        basin_nearbare=sum(r['near_boundary'] for r in basins),basin_positive=sum(r['converged'] for r in basins),
        structural_settings=len(scan),
        scan_roots=sum(len(r['roots']) for r in scan),
        fold=read(root/'audit_v2/fold.json'),
        extended_holds=read(root/'audit_v2/extended_holds.json'),
        near_fold_refinement=read(root/'audit_v2/near_fold_refinement.json'),
        conclusion='No R-tipping detected in these runs; not a global nonexistence proof.')
    save(out/'summary.json',summary)
    save(out/'diagnostic_columns.json', ['temperature','gpp','demand_L','demand_S','demand_R',
        'paid_L','paid_S','paid_R','unpaid_L','unpaid_S','unpaid_R','flush','mortality',
        'growth_resp','rh','export','chi','chi_L','chi_S','chi_R','algebraic_budget'])
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
