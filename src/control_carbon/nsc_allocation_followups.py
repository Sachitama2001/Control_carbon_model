"""Source audit for V6/V7; no new ODE or calibrated parameters implied.

Sachitama2001/VISIT-matrix@3285bd8e131a932e338b59892751648fd9edcc7b:
ecophysiology.c:77-87, plant_proc.c:209-233, allocation.c:125-154.
Units: stocks/one-day increments Mg C ha-1. Native order is sequential.
"""
import numpy as np
from .temperature_nsc_visit import living_tissue,leaf_carbon_from_lai


def storage_capacity(stem,root,visit):
    return (.1*living_tissue(stem,visit.f_size_stem)
            +.3*living_tissue(root,visit.f_size_root))


def native_positive_storage_step(stocks,allocation,growth_resp,capacity):
    """[leaf,stem,root,storage] sequential update for positive allocations.

    Intentionally returns negative structure if native source does: this is an
    audit, not a nonnegative ODE rule. Respiration is debited from structure
    even when positive stem/root allocation was redirected into storage.
    """
    x=np.asarray(stocks,float).copy()
    a=np.asarray(allocation,float); r=np.asarray(growth_resp,float)
    if x.shape!=(4,) or a.shape!=(3,) or r.shape!=(3,):
        raise ValueError("invalid shape")
    if not np.isfinite(np.r_[x,a,r,capacity]).all() or min(np.r_[x,a,r,capacity])<0:
        raise ValueError("finite nonnegative input required")
    x[0]+=a[0]-r[0]
    routed=[]
    for i in (1,2):
        if a[i]>0 and x[3]<capacity:
            x[3]+=a[i];routed.append("storage")
        else:
            x[i]+=a[i];routed.append("structure")
        x[i]-=r[i]
    return x,routed


def rescue_increment(leaf,stem,root,visit):
    """Tree daily structural redistribution, active below critical LAI=0.1."""
    if not np.isfinite([leaf,stem,root]).all() or min(leaf,stem,root)<0:
        raise ValueError("finite nonnegative structure required")
    critical=leaf_carbon_from_lai(.1,visit)
    if leaf>=critical:
        return np.zeros(3)
    hs=critical*visit.alloc_abg*(.05*stem)/(.5*critical+.05*stem)
    hr=critical*(1-visit.alloc_abg)*(.1*root)/(.5*critical+.1*root)
    return np.array([hs+hr,-hs,-hr])
