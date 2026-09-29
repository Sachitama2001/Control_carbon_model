"""Render analytic certificates and functions, without integrating ODE paths."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from control_carbon.qse_rate_theory import (
    COMOVING_FOLD_RATE,critical_rate,crossing_displacement,
    bounded_input_coefficients,smooth_certificates,
)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,default=ROOT/"artifacts/qse_rate_theory/analytic_proofs_20260929")
    args=ap.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=False)
    r=critical_rate()
    summary=dict(displacement=1.5,critical_linear_rate=r,critical_linear_duration=1.5/r,
        moving_frame_fold=COMOVING_FOLD_RATE,
        cubic_smooth_bounds=smooth_certificates(),
        bounded_smooth_bounds=smooth_certificates(bounded=True),
        bounded_coefficients_at_lambda_zero=bounded_input_coefficients(0),
        result="analytic existence and finite-linear-ramp uniqueness; no ODE trajectories used",
        provenance="model-agnostic derived construction; not VISIT or calibrated forest rates")
    (out/"certificates.json").write_text(json.dumps(summary,indent=2)+"\n")
    sources=[Path(__file__),ROOT/"src/control_carbon/qse_rate_theory.py",
             ROOT/"tests/test_qse_rate_theory.py",ROOT/"docs/qse_rate_tipping_analytic_proofs.tex"]
    snapshot=out/"source_snapshot";snapshot.mkdir()
    for f in sources:
        shutil.copy2(f,snapshot/f.name)
    (out/"manifest.json").write_text(json.dumps({
        "source_sha256":{str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sources},
        "claims":"Proofs are in the TeX; quadrature checks their explicit constants."},indent=2)+"\n")
    fig,axes=plt.subplots(1,3,figsize=(14,4))
    lam=np.linspace(0,1.5,100)
    axes[0].plot(lam,lam+3,label="High stable QSE",color="tab:green")
    axes[0].plot(lam,lam+2,"--",label="Unstable basin boundary",color="k")
    axes[0].plot(lam,lam+1,label="Low stable QSE",color="tab:orange")
    axes[0].scatter([0,1.5,1.5],[3,4.5,2.5],c=["k","tab:green","tab:orange"],zorder=5)
    axes[0].set(xlabel="Environment (dimensionless)",ylabel="Stock (dimensionless)",
                title="Both stable branches persist")
    axes[0].legend(fontsize=8,loc="upper left")
    rates=np.linspace(COMOVING_FOLD_RATE+.02,2,150)
    axes[1].plot(rates,[crossing_displacement(v) for v in rates],label="H(r)")
    axes[1].axhline(1.5,color="k",ls="--",label="D = 1.5")
    axes[1].scatter([r],[1.5],color="tab:red",label=f"r* = {r:.6f}")
    axes[1].set(xlabel="Linear ramp rate r",ylabel="Displacement to cross boundary",
                title="Unique finite-ramp threshold")
    axes[1].legend(fontsize=8)
    x=np.linspace(0,4,401);p,a,h2=bounded_input_coefficients(0)
    h=p/(1+x)+a*x/(h2+x*x)-1
    axes[2].plot(x,h,label="P(x)/x - 1")
    axes[2].axhline(0,color="k",lw=.8)
    axes[2].scatter([1,2,3],[0,0,0],c=["tab:green","tab:red","tab:green"],zorder=5)
    axes[2].set(xlabel="Stock (dimensionless)",ylabel="Per-capita net growth",
                title="Bounded production, linear loss")
    axes[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out/"proof_diagram.pdf")
    fig.savefig(out/"proof_diagram.png",dpi=160)
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
