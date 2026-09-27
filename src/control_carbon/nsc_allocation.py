"""Branch-aware VISIT allocation adapter (source-grounded reduced model).

Source: Sachitama2001/VISIT-matrix, commit
3285bd8e131a932e338b59892751648fd9edcc7b, visit_local/allocation.c,
f_allocation:16-114. Only positive EPP, active non-crop C3 tree branch.
Stocks Mg C ha-1; ODE time days; source EPP is a one-day amount.
Mobilization replaces native EPP and native sequential updates; not exact VISIT.
"""
from dataclasses import dataclass
import numpy as np

from .temperature_nsc import evaluate, transfers
from .temperature_nsc_visit import (
    optimum_leaf_carbon, visit_canopy_gpp, size_dependent_respiration,
)
from .nsc_allocation_followups import rescue_increment


@dataclass(frozen=True)
class Mobilization:
    """Explicit reduced-model assumptions, not TKY parameters."""
    rates: tuple
    reference_days: float = 1.0

    def __post_init__(self):
        if len(self.rates) != 3 or not np.all(np.isfinite(self.rates)) or min(self.rates) <= 0:
            raise ValueError("three positive finite mobilization rates required")
        if not np.isfinite(self.reference_days) or self.reference_days <= 0:
            raise ValueError("positive finite reference_days required")


def source_allocation(lai, target_lai, epp, visit):
    """Native positive-EPP algebra; input EPP is Mg C ha-1 per day step.

    Equality follows the C branch. Zero EPP has no positive-EPP ratio here.
    """
    values = [lai, target_lai, epp, visit.sla, visit.alloc_ass, visit.alloc_abg]
    if not np.all(np.isfinite(values)) or min(lai, target_lai) < 0 or epp <= 0:
        raise ValueError("finite nonnegative LAI and positive EPP required")
    if visit.sla <= 0 or not 0 < visit.alloc_ass <= 1 or not 0 <= visit.alloc_abg <= 1:
        raise ValueError("invalid source allocation parameters")
    if lai > target_lai:
        cf, branch = 0.0, "above"
    else:
        deficit = (target_lai-lai)*200.0/(2.2*visit.sla)
        if deficit <= visit.alloc_ass*epp:
            cf, branch = visit.alloc_ass, "near"
        else:
            cf, branch = min(deficit/epp, 0.05), "far"
    return np.array([cf, (1-cf)*visit.alloc_abg, (1-cf)*(1-visit.alloc_abg)]), branch


def branch_conditions(x, target, mobilization, visit):
    """h0=0: target; h1=0: far/near boundary, both in Mg C ha-1."""
    v = np.dot(mobilization.rates, x[4:7])
    h0 = x[0]-target
    h1 = x[0]-target+visit.alloc_ass*v*mobilization.reference_days
    return float(h0), float(h1)


def evaluate_allocation(x, temperature, base, visit, env, mobilization,
                        variant="V4", forced_fraction=None, forced_rescue=None):
    """Conservative donor-funded allocation; forced_fraction is branch analysis.

    V4 uses baseline demand; V5 uses source size-dependent demand. Common yield
    is kept unchanged to isolate allocation from growth-respiration parameters.
    Hard-switch trajectories need event handling; use branches for roots first.
    """
    if variant not in ("V4", "V5", "V7_V4"):
        raise ValueError("supported variants: V4,V5,V7_V4")
    if forced_rescue is not None and (variant != "V7_V4" or not 0 <= forced_rescue <= 1):
        raise ValueError("rescue override requires V7_V4 and fraction in [0,1]")
    if base.nsc_mortality or base.mu_max != 0:
        raise ValueError("NSC mortality must be disabled")
    x = np.asarray(x, dtype=float)
    # Reuse validation and background/litter budget; replace ALL old investment.
    dx, diag = evaluate(x, temperature, base)
    target, target_lai, status = optimum_leaf_carbon(temperature, visit, env)
    if not np.isfinite(target):
        raise ValueError(f"invalid target: {status}")
    n = x[4:7]
    donor = np.asarray(mobilization.rates)*n
    total = donor.sum()
    if forced_fraction is not None:
        if not np.isfinite(forced_fraction) or not 0 <= forced_fraction <= 1:
            raise ValueError("invalid forced allocation fraction")
        cf, branch = forced_fraction, "forced"
        fractions = np.array([cf, (1-cf)*visit.alloc_abg, (1-cf)*(1-visit.alloc_abg)])
    elif total == 0:
        fractions, branch = np.zeros(3), "no_mobilization"
    else:
        gamma = visit.sla*2.2/200.0
        fractions, branch = source_allocation(gamma*x[0], target_lai,
                                               total*mobilization.reference_days, visit)
    investment = fractions*total
    loss = np.asarray(base.turnover)+base.mu0
    demand = (size_dependent_respiration(x[:3], temperature, visit)
              if variant == "V5" else diag["demand"])
    paid = demand*n/(np.asarray(base.substrate_half)+n)
    gpp = visit_canopy_gpp(x[0], temperature, visit, env)
    j = transfers(n, base)
    dx[:3] = base.yield_c*investment-loss*x[:3]
    rescue=np.zeros(3)
    if variant=="V7_V4":
        rescue=(rescue_increment(*x[:3],visit) if forced_rescue is None else
                forced_rescue*rescue_increment(0,x[1],x[2],visit))
        rescue=rescue/mobilization.reference_days
        dx[:3]+=rescue
    dx[4:7] = -paid-loss*n-donor
    dx[4] += gpp-j[0]-j[1]+j[2]+j[3]
    dx[5] += j[0]-j[2]
    dx[6] += j[1]-j[3]
    rg = (1-base.yield_c)*total
    external = gpp-paid.sum()-rg-diag["rh"]-diag["export"]
    diag.update(gpp=gpp, demand=demand, paid=paid, unpaid=demand-paid,
                growth_resp=rg, investment=investment, donor=donor,
                fractions=fractions, branch=branch, target=target, rescue=rescue,
                target_lai=target_lai, external=external,
                budget=float(dx.sum()-external), variant=variant)
    # Old flush/growth diagnostics do not describe this new investment rule.
    diag.pop("flush")
    diag.pop("growth")
    return dx, diag
