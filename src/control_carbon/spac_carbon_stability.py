"""Outward-rounded B2 certificates for conducting E0--E3 frozen branches.

Decimal endpoints enclose real arithmetic with the supplied binary64
constants interpreted as exact real numbers. This is NOT a certificate
for nonautonomous tracking, mortality, or E4 wood feedback.
"""

from dataclasses import dataclass
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING

import numpy as np


PRECISION = 40


def _decimal(value):
    if isinstance(value, Decimal):
        return value
    if isinstance(value, float):
        return Decimal.from_float(value)
    return Decimal(value)


def _rounded(operation, rounding):
    with localcontext() as ctx:
        ctx.prec = PRECISION
        ctx.rounding = rounding
        return operation()


@dataclass(frozen=True)
class _Interval:
    lo: Decimal
    hi: Decimal

    def __init__(self, lo, hi=None):
        low = _decimal(lo)
        high = low if hi is None else _decimal(hi)
        if not low.is_finite() or not high.is_finite() or low > high:
            raise ValueError("interval endpoints must be finite and ordered")
        object.__setattr__(self, "lo", low)
        object.__setattr__(self, "hi", high)

    def __add__(self, other):
        other = other if isinstance(other, _Interval) else _Interval(other)
        return _Interval(_rounded(lambda: self.lo+other.lo, ROUND_FLOOR),
                         _rounded(lambda: self.hi+other.hi, ROUND_CEILING))

    __radd__ = __add__

    def __neg__(self):
        return _Interval(self.hi.copy_negate(), self.lo.copy_negate())

    def __sub__(self, other):
        other = other if isinstance(other, _Interval) else _Interval(other)
        return self+-other

    def __rsub__(self, other):
        return -self+other

    def __mul__(self, other):
        other = other if isinstance(other, _Interval) else _Interval(other)
        pairs = [(a, b) for a in [self.lo, self.hi] for b in [other.lo, other.hi]]
        return _Interval(min(_rounded(lambda: a*b, ROUND_FLOOR) for a, b in pairs),
                         max(_rounded(lambda: a*b, ROUND_CEILING) for a, b in pairs))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = other if isinstance(other, _Interval) else _Interval(other)
        if other.lo <= 0 <= other.hi:
            raise ValueError("interval denominator contains zero")
        inverse = _Interval(_rounded(lambda: Decimal(1)/other.hi, ROUND_FLOOR),
                            _rounded(lambda: Decimal(1)/other.lo, ROUND_CEILING))
        return self*inverse

    def __rtruediv__(self, other):
        return _Interval(other)/self

    def exp(self):
        # Decimal.exp is correctly rounded HALF_EVEN, regardless of context
        # rounding. One adjacent 40-digit value on each side encloses it.
        with localcontext() as ctx:
            ctx.prec = PRECISION
            return _Interval(self.lo.exp().next_minus(), self.hi.exp().next_plus())

    def strings(self):
        return [str(self.lo), str(self.hi)]


def _branch_intervals(leaf, root_potential, hydraulic, carbon, feedback):
    if feedback.wood_supply != 0:
        raise ValueError("cubic certificate requires zero wood feedback")
    constants = {name: _Interval(value) for name, value in {
        "gxmax": hydraulic.gx_max, "gsmax": hydraulic.gs_max,
        "ax": -2*hydraulic.psi50_x, "bs": -2*hydraulic.psi50_s,
        "cx": hydraulic.cx, "cl": hydraulic.cl, "root": root_potential,
        "alpha": carbon.allocation[0], "p": carbon.assimilation_max,
        "k": carbon.losses[0], "q": carbon.canopy_coefficient,
        "lf": feedback.leaf_floor, "lr": feedback.leaf_reference,
        "al": feedback.leaf_supply, "ae": feedback.leaf_demand,
    }.items()}
    g = constants
    decay = (-g["q"]*leaf).exp()
    canopy = 1-decay
    alpha = g["alpha"]*g["p"]
    h = g["k"]*leaf/(alpha*canopy)
    psi = g["bs"]*(h-1)
    gx = g["gxmax"]*(1+psi/g["ax"])
    gs = g["gsmax"]*h
    dfl = g["al"]/(g["lr"]+g["lf"])
    dfe = g["ae"]/g["lr"]
    fl = 1+g["al"]*((leaf+g["lf"])/(g["lr"]+g["lf"])-1)
    fe = 1+g["ae"]*(leaf/g["lr"]-1)
    x = (g["root"]+fl*psi)/(1+fl)
    geff = 2*gx*fl/(1+fl)
    drop = g["root"]-psi
    d = geff*drop/(gs*fe)
    dgx, dgs = g["gxmax"]/g["ax"], g["gsmax"]/g["bs"]
    psiprime = g["bs"]*g["k"]/alpha*(canopy-leaf*g["q"]*decay)/(canopy*canopy)
    geffprime = 2*dgx*psiprime*fl/(1+fl)+2*gx*dfl/((1+fl)*(1+fl))
    dprime = (geffprime*drop-geff*psiprime)/(gs*fe)-d*(dgs*psiprime/gs+dfe/fe)
    gr, gl = 2*gx, 2*gx*fl
    a = [
        [-(gr+gl)/g["cx"],
         (2*dgx*(g["root"]-x)-2*dgx*fl*(x-psi)+gl)/g["cx"],
         -2*gx*dfl*(x-psi)/g["cx"]],
        [gl/g["cl"], (2*dgx*fl*(x-psi)-gl-dgs*fe*d)/g["cl"],
         (2*gx*dfl*(x-psi)-gs*dfe*d)/g["cl"]],
        [_Interval(0), alpha/g["bs"]*canopy, alpha*h*g["q"]*decay-g["k"]],
    ]
    a1 = -(a[0][0]+a[1][1]+a[2][2])
    a2 = sum(a[i][i]*a[j][j]-a[i][j]*a[j][i] for i, j in [(0, 1), (0, 2), (1, 2)])
    a3 = -(a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
           -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
           +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))
    return {"h": h, "psi": psi, "gx": gx, "D": d, "Dprime_L": dprime,
            "a1": a1, "a2": a2, "a3": a3, "delta": a1*a2-a3}


def certify_stable_path(leaf_bounds, vpd_bounds, root_potential, hydraulic, carbon, feedback, *,
                        max_depth=18, max_cells=50000):
    """Cover a continuous leaf interval, certify Hurwitz and monotone D(L).

    Endpoint D enclosures must cover the requested whole forcing interval.
    Any undecided/unstable cell fails explicitly; no grid-only fallback.
    """
    lo, hi = map(_decimal, leaf_bounds)
    dlow, dhigh = map(_decimal, vpd_bounds)
    if not lo.is_finite() or not hi.is_finite() or not Decimal(0) < lo < hi:
        raise ValueError("positive finite ordered leaf bounds required")
    if not dlow.is_finite() or not dhigh.is_finite() or not Decimal(0) <= dlow < dhigh:
        raise ValueError("ordered nonnegative finite VPD interval required")
    if not np.isfinite(root_potential) or root_potential > 0 or lo <= _decimal(feedback.minimum_leaf):
        raise ValueError("root/leaf domain invalid for continuous-path certificate")
    if max_depth < 1 or max_cells < 1:
        raise ValueError("positive certificate subdivision limits required")
    low_data = _branch_intervals(_Interval(lo), root_potential, hydraulic, carbon, feedback)
    high_data = _branch_intervals(_Interval(hi), root_potential, hydraulic, carbon, feedback)
    if low_data["D"].lo < dhigh or high_data["D"].hi > dlow:
        raise ValueError("endpoint D enclosures do not cover requested forcing path")
    stack, accepted = [(lo, hi, 0)], []
    while stack:
        left, right, depth = stack.pop()
        data = _branch_intervals(_Interval(left, right), root_potential, hydraulic, carbon, feedback)
        physical = (data["h"].lo > 0 and data["h"].hi < 1 and data["gx"].lo > 0
                    and data["psi"].hi < _decimal(root_potential))
        stable = all(data[key].lo > 0 for key in ["a1", "a2", "a3", "delta"])
        if physical and stable and data["Dprime_L"].hi < 0:
            accepted.append({"leaf_interval": [str(left), str(right)],
                             **{key: data[key].strings() for key in ["D", "Dprime_L", "a1", "a2", "a3", "delta"]}})
            continue
        if depth >= max_depth or len(accepted)+len(stack)+2 > max_cells:
            raise RuntimeError(f"certificate unresolved/nonstable on leaf [{left},{right}], depth={depth}")
        with localcontext() as ctx:
            ctx.prec = PRECISION
            midpoint = (left+right)/2
        if not left < midpoint < right:
            raise RuntimeError("certificate subdivision exhausted representable resolution")
        stack.extend([(midpoint, right, depth+1), (left, midpoint, depth+1)])
    return {
        "precision_decimal_digits": PRECISION,
        "parameter_semantics": "input binary64 constants interpreted as exact real numbers",
        "leaf_bounds": [str(lo), str(hi)], "requested_D_bounds": [str(dlow), str(dhigh)],
        "endpoint_D": {"low_leaf": low_data["D"].strings(), "high_leaf": high_data["D"].strings()},
        "cell_count": len(accepted), "cells": accepted,
        "strict_lower_bounds": {key: str(min(_decimal(cell[key][0]) for cell in accepted))
                                for key in ["a1", "a2", "a3", "delta"]},
        "Dprime_L_upper_bound": str(max(_decimal(cell["Dprime_L"][1]) for cell in accepted)),
        "claim": "unique stable positive equilibrium along the certified branch for every D in the requested interval",
        "not_claimed": "global absence of other branches; nonautonomous tracking; no-R; physiological mortality",
    }
