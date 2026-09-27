"""Summarize completed allocation experiments, with auditable scalar bounds."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("directory",type=Path)
    args=ap.parse_args(); out=args.directory
    rows=json.loads((out/"frozen.json").read_text())
    trajectories=json.loads((out/"trajectories.json").read_text())
    rescue=json.loads((out/"rescue.json").read_text())
    inherited=[q for q in rows if q["turnover_case"]=="inherited"]
    sensitivity=[q for q in rows if q["turnover_case"]!="inherited"]
    refinement=[q["refinement_scaled_error"] for q in trajectories if "refinement_scaled_error" in q]
    summary=dict(frozen_conditions=len(rows),trajectories=len(trajectories),
        inherited_positive_equilibria_excluded=all(q["bounds"]["positive_equilibria_excluded"] for q in inherited),
        tky_attracting_sliding=sum(b["classification"]=="attracting_sliding" and b["tangent_max_real"]<0 for q in sensitivity for b in q["boundaries"]),
        rescue_attracting_sliding=sum(b["attracting"] and b["tangent_max_real"]<0 for q in rescue for b in q["boundaries"]),
        near_bare_events=sum(q["near_bare_event"] for q in trajectories),
        boundary_events=sum(q["boundary_event"] for q in trajectories),
        domain_events=sum(q["domain_event"] for q in trajectories),
        minimum_stock=min(q["minimum_stock"] for q in trajectories),
        maximum_budget_error=max(q["integrated_budget_error"] for q in trajectories),
        maximum_refinement_error=max(refinement),refinement_runs=len(refinement),
        minimum_final_plant_stock=min(sum(np.array(q["endpoint"])[[0,1,2,4,5,6]]) for q in trajectories),
        maximum_final_plant_stock=max(sum(np.array(q["endpoint"])[[0,1,2,4,5,6]]) for q in trajectories),
        interpretation="No multistability or R-tipping established; sliding candidates conditional on continuous-time interpretation.")
    (out/"audited_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    if summary["maximum_refinement_error"]>1e-5:
        raise RuntimeError("refinement audit failed")
    def sci(v):
        a,b=f"{v:.3e}".split("e")
        return rf"{a}\times10^{{{int(b)}}}"
    paragraph=(
        f"凍結環境の解析は{len(rows)}条件、長期積分は{len(trajectories)}本である。"
        f"切替面への到達は{summary['boundary_events']}件、数値領域境界への到達は"
        f"{summary['domain_events']}件だった。"
        f"裸地近傍（植物総量$10^{{-12}}$ Mg C ha$^{{-1}}$）への到達は"
        f"{summary['near_bare_events']}件である。"
        f"数値領域境界では、いずれかの対数在庫が$-600$に達した時点で安全停止した。"
        "この停止を完全収束に数えない。\n"
        rf"監査点の最小在庫は$ {sci(summary['minimum_stock'])}$ Mg C ha$^{{-1}}$、"
        rf"累積収支誤差の最大値は$ {sci(summary['maximum_budget_error'])}$ Mg C ha$^{{-1}}$である。"
        f"刻み半減と許容誤差強化の比較は{len(refinement)}条件で行い、"
        rf"共通の年次標本での最大尺度化差は$ {sci(max(refinement))}$だった。"
        "境界イベントの到達は消滅を有限時刻に証明したものではない。\n")
    (out/"numerical_summary.tex").write_text(paragraph)
    fig,axes=plt.subplots(1,3,figsize=(14,4))
    selected=[q for q in inherited if q["variant"]=="V4" and q["multiplier"]==1]
    t=[q["temperature"] for q in selected]
    axes[0].plot(t,[q["bounds"]["far_gain_ratio"] for q in selected],"o-",label="Far: y c G'(0) / loss")
    axes[0].plot(t,[q["bounds"]["gpp_upper_bound"]/q["bounds"]["near_min_mobilization"] for q in selected],"s-",label="Near: Gmax / Vmin")
    axes[0].axhline(1,color="k",ls="--")
    axes[0].set(xlabel="Temperature (C)",ylabel="Necessary gain ratio",title="Inherited turnover: both ratios < 1")
    axes[0].legend(fontsize=8)
    q=[q for q in sensitivity if q["variant"]=="V4" and q["multiplier"]==1]
    axes[1].plot([r["temperature"] for r in q],[r["boundaries"][0]["state"][0] for r in q],"o-",label="V4, TKY basal turnover")
    axes[1].plot([r["temperature"] for r in rescue],[r["critical_leaf"] for r in rescue],"s-",label="V7, inherited turnover")
    axes[1].set(yscale="log",xlabel="Temperature (C)",ylabel="Leaf carbon (Mg C/ha)",title="Different models: sliding candidates")
    axes[1].legend(fontsize=8)
    for variant in ("V4","V5"):
        for t in (10,20,30):
            q=next(q for q in trajectories if q["variant"]==variant and q["temperature"]==t and q["biomass"]==100 and q["chi"]==.15)
            data=np.load(out/(q["name"]+".npz"))
            plant=data["states"][:,[0,1,2,4,5,6]].sum(axis=1)
            axes[2].plot(data["times"]/365,plant,label=f"{variant}, {t} C",ls="-" if variant=="V4" else "--")
    axes[2].set(yscale="log",xlabel="Years (fixed temperature)",ylabel="Plant C (Mg C/ha)",title="Inherited turnover: long relaxation")
    axes[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out/"comparison.pdf")
    fig.savefig(out/"comparison.png",dpi=160)
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
