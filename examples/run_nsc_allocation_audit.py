"""V4/V5 branch audit; no hard-switch trajectory is advanced across a surface."""
import argparse
from dataclasses import asdict,replace
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from control_carbon.temperature_nsc import Parameters
from control_carbon.temperature_nsc_visit import load_tky_tree_parameters,FixedEnvironment
from control_carbon.nsc_allocation import Mobilization
from control_carbon.nsc_allocation_analysis import frozen_audit,exclusion_bounds,bare_invasion,integrate_far_branch,rescue_audit


def save(path,obj):
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--config",type=Path,default=ROOT/"configs/temperature_nsc_allocation_v2.json")
    ap.add_argument("--visit-config",type=Path,default=ROOT.parent/"VISIT-matrix/visit_local/INPUT/Config_TKY.xlsx")
    ap.add_argument("--output",type=Path)
    ap.add_argument("--roots-only",action="store_true")
    args=ap.parse_args(); cfg=json.loads(args.config.read_text())
    out=args.output or ROOT/"artifacts/temperature_nsc_visit_allocation"/cfg["experiment_id"]
    out.mkdir(parents=True,exist_ok=False)
    visit=load_tky_tree_parameters(args.visit_config)
    env=FixedEnvironment(**cfg["fixed_environment"])
    original=replace(Parameters(),nsc_mortality=False,mu_max=0)
    sources=[Path(__file__),args.config,ROOT/"src/control_carbon/nsc_allocation.py",
             ROOT/"src/control_carbon/nsc_allocation_analysis.py",
             ROOT/"src/control_carbon/nsc_allocation_followups.py",
             ROOT/"src/control_carbon/temperature_nsc.py",
             ROOT/"src/control_carbon/temperature_nsc_visit.py",
             args.visit_config,args.visit_config.parent.parent/"allocation.c"]
    snapshot=out/"source_snapshot";snapshot.mkdir()
    for source in sources:
        shutil.copy2(source,snapshot/source.name)
    save(out/"manifest.json",dict(config=cfg,roots_only=args.roots_only,
        base=asdict(original),visit=asdict(visit),environment=asdict(env),
        source_commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=args.visit_config.parent,text=True).strip(),
        code_commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        hashes={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in sources},
        source_units="Mg C ha-1; day",status="exploratory; conditional fixed forcing; not calibrated TKY prediction"))
    results=[]
    for turnover in cfg["turnover_cases"]:
        p=(original if turnover=="inherited" else replace(original,turnover=(
            visit.leaf_turnover,visit.stem_turnover,visit.root_turnover)))
        for variant in cfg["variants"]:
            for multiplier in cfg["multipliers"]:
                m=Mobilization(tuple(multiplier*np.asarray(cfg["mobilization_rates_per_day"])),cfg["reference_days"])
                for t in cfg["temperatures_c"]:
                    row=frozen_audit(t,p,visit,env,m,variant,cfg["root_scan_points"])
                    row.update(turnover_case=turnover,multiplier=multiplier,
                        bounds=exclusion_bounds(t,p,visit,env,m.reference_days),
                        bare=bare_invasion(t,p,visit,env,m))
                    results.append(row)
                    save(out/"frozen.json",results)
                    print(turnover,variant,multiplier,t,
                          "admissible",sum(q["admissible"] for q in row["roots"]),
                          "surfaces",len(row["boundaries"]),flush=True)
    if args.roots_only:
        return
    if cfg.get("rescue_audit",False):
        rescue=[]
        for t in cfg["temperatures_c"]:
            rescue.append(rescue_audit(t,original,visit,env,
                Mobilization(tuple(cfg["mobilization_rates_per_day"]),cfg["reference_days"]),
                points=cfg["root_scan_points"]))
            save(out/"rescue.json",rescue)
    # Initial experiment retains prior turnover; sensitivity roots are separate.
    # Sample the endpoints and central temperature, full predeclared seed grid.
    m=Mobilization(tuple(cfg["mobilization_rates_per_day"]),cfg["reference_days"])
    trajectories=[]
    for variant in cfg["variants"]:
        for t in (10,20,30):
            for b in cfg["seed_biomass"]:
                for chi in cfg["seed_nsc_fraction"]:
                    x=np.zeros(8);x[:3]=b*np.array([.03,.8,.17])
                    x[4:7]=b*chi/(1-chi)*np.array([.1,.75,.15]);x[[3,7]]=[100,1]
                    opts=dict(years=cfg["spinup_years"],max_step=cfg["max_step_days"],rtol=cfg["rtol"],atol=cfg["atol"])
                    tr=integrate_far_branch(x,t,original,visit,env,m,variant,**opts)
                    name=f"{variant}_T{t}_B{b}_chi{chi}"
                    np.savez_compressed(out/(name+".npz"),times=tr["times"],states=tr["states"],external_integral=tr["external_integral"])
                    row={k:v for k,v in tr.items() if k not in ("times","states","external_integral")}
                    row.update(name=name,variant=variant,temperature=t,biomass=b,chi=chi,
                               end_day=float(tr["times"][-1]),endpoint=tr["states"][-1].tolist())
                    # One large initial condition per temperature/variant: tighter solve.
                    if b==100 and chi==.15:
                        refined=integrate_far_branch(x,t,original,visit,env,m,variant,
                            years=cfg["spinup_years"],max_step=cfg["max_step_days"]/2,
                            rtol=cfg["rtol"]/10,atol=cfg["atol"]/10)
                        # Compare only common scheduled times, not differing event endpoints.
                        common,ia,ib=np.intersect1d(tr["times"],refined["times"],return_indices=True)
                        row["refinement_scaled_error"]=float(np.max(np.abs(tr["states"][ia]-refined["states"][ib])/np.maximum(1,refined["states"][ib])))
                        row["refinement_common_samples"]=len(common)
                    trajectories.append(row);save(out/"trajectories.json",trajectories)
                    print(name,"near_bare",tr["near_bare_event"],"boundary",tr["boundary_event"],flush=True)
    save(out/"summary.json",dict(frozen_conditions=len(results),trajectories=len(trajectories),
        near_bare_events=sum(r["near_bare_event"] for r in trajectories),
        boundary_events=sum(r["boundary_event"] for r in trajectories),
        minimum_stock=min(r["minimum_stock"] for r in trajectories),
        budget_error=max(r["integrated_budget_error"] for r in trajectories)))


if __name__=="__main__":
    main()
