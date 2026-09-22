"""Analytical two-state auxiliary system, NOT an exact 8-state reduction.

Fixed leaf fraction, equal structural/NSC turnover and one linear growth sink
are additional closure assumptions. Positive equilibria have unique chi.
"""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq
from .temperature_nsc import Parameters, mortality


@dataclass(frozen=True)
class ReducedParameters:
    growth: float = .0015
    turnover: float = .0001
    respiration: float = .0001
    substrate_half: float = .5
    leaf_fraction: float = .03


def reduced_rhs(x, temperature, p=Parameters(), a=ReducedParameters()):
    b, s = x
    if b < 0 or s < 0:
        raise ValueError('nonnegative stocks required')
    if b+s == 0:
        return np.zeros(2)
    ell = mortality(s/(b+s), p)+a.turnover
    thermal = np.exp(-.5*((temperature-p.t_opt)/p.sigma_g)**2)
    g = p.gmax*thermal*(-np.expm1(-p.leaf_extinction*a.leaf_fraction*b))
    r = a.respiration*p.q10**((temperature-p.t_ref)/10)*b*s/(a.substrate_half+s)
    u = a.growth*s
    return np.array([p.yield_c*u-ell*b, g-r-u-ell*s])


def reduced_equilibrium(temperature, p=Parameters(), a=ReducedParameters()):
    chi = brentq(lambda c: p.yield_c*a.growth*c/(1-c)
                 -mortality(c,p)-a.turnover, 0, 1-1e-12)
    z = chi/(1-chi)
    ell = mortality(chi,p)+a.turnover
    g0 = p.gmax*np.exp(-.5*((temperature-p.t_opt)/p.sigma_g)**2)
    k = p.leaf_extinction*a.leaf_fraction
    r = a.respiration*p.q10**((temperature-p.t_ref)/10)
    def balance(b):
        assimilated_per_b = g0*k if b == 0 else g0*(-np.expm1(-k*b))/b
        return assimilated_per_b-r*z*b/(a.substrate_half+z*b)-(a.growth+ell)*z
    if balance(0) <= 0:
        return None
    upper = 1.0
    while balance(upper) > 0:
        upper *= 2
    b = brentq(balance, 0, upper)
    return np.array([b,z*b])
