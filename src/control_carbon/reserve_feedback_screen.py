"""Model-agnostic two-pool screening; NOT a fitted forest or native VISIT ODE.

See docs/allocation_storage_respiration_evidence.md for assumptions and proofs.
The default numbers below are dimensionless verification fixtures only.
"""
from dataclasses import dataclass
import math
import numpy as np


@dataclass(frozen=True)
class ReserveControl:
    pmax: float = 3.0
    kp: float = 1.0
    vmax: float = 1.0
    kv: float = 0.2
    hill: float = 1.0
    rmax: float = 0.2
    kr: float = 0.1
    yield_: float = 0.7
    turnover: float = 0.1
    reserve_loss: float = 0.02

    def __post_init__(self):
        vals = tuple(self.__dict__.values())
        if not all(math.isfinite(v) for v in vals):
            raise ValueError("finite parameters required")
        if min(self.pmax, self.kp, self.vmax, self.kv, self.kr,
               self.yield_, self.turnover) <= 0:
            raise ValueError("positive rates, scales and yield required")
        if self.hill < 1 or self.yield_ > 1 or min(self.rmax, self.reserve_loss) < 0:
            raise ValueError("Hill exponent >=1, yield <=1 and nonnegative losses required")


def concentration_rates(n, p):
    """Growth/respiration per structural carbon and concentration derivatives."""
    if not math.isfinite(n) or n < 0:
        raise ValueError("nonnegative finite reserve concentration required")
    if n == 0:
        v, vp = 0.0, (p.vmax/p.kv if p.hill == 1 else 0.0)
    else:
        log_ratio = p.hill*math.log(n/p.kv)
        # Stable logistic evaluation of the Hill fraction.
        small = math.exp(-abs(log_ratio))
        z = 1/(1+small) if log_ratio >= 0 else small/(1+small)
        v, vp = p.vmax*z, p.vmax*p.hill*z*(1-z)/n
    r = p.rmax*n/(p.kr+n)
    rp = p.rmax*p.kr/(p.kr+n)**2
    return v, vp, r, rp


def reserve_fluxes(state, p):
    c, nsc = map(float, state)
    if not all(math.isfinite(z) and z >= 0 for z in (c, nsc)):
        raise ValueError("finite nonnegative carbon stocks required")
    if c == 0:
        v = r = 0.0
    else:
        v, _, r, _ = concentration_rates(nsc/c, p)
    production = p.pmax*c/(p.kp+c)
    mobilization = c*v
    maintenance = c*r
    return dict(production=production, mobilization=mobilization,
                maintenance=maintenance,
                growth_respiration=(1-p.yield_)*mobilization,
                structural_export=p.turnover*c,
                reserve_export=p.reserve_loss*nsc)


def reserve_rhs(state, p):
    f = reserve_fluxes(state, p)
    return np.array([p.yield_*f['mobilization']-f['structural_export'],
                     f['production']-f['mobilization']-f['maintenance']-f['reserve_export']])


def positive_equilibrium(p):
    """Exact elimination, not fast-NSC approximation. None means no positive root."""
    if p.yield_*p.vmax <= p.turnover:
        return None
    n = p.kv*(p.turnover/(p.yield_*p.vmax-p.turnover))**(1/p.hill)
    v, _, r, _ = concentration_rates(n, p)
    cost = v+r+p.reserve_loss*n
    c = p.pmax/cost-p.kp
    return None if c <= 0 else np.array([c, n*c])


def reserve_jacobian(state, p):
    c, nsc = map(float, state)
    if c <= 0 or nsc < 0:
        raise ValueError("positive structural carbon and nonnegative NSC required")
    n = nsc/c
    v, vp, r, rp = concentration_rates(n, p)
    pp = p.pmax*p.kp/(p.kp+c)**2
    return np.array([[p.yield_*(v-n*vp)-p.turnover, p.yield_*vp],
                     [pp-v+n*vp-r+n*rp, -vp-rp-p.reserve_loss]])


def dulac_divergence(state, p):
    """div(F/C), strictly negative in the positive quadrant for this control."""
    c, nsc = map(float, state)
    if c <= 0 or nsc <= 0:
        raise ValueError("strictly positive stocks required")
    n = nsc/c
    _, vp, _, rp = concentration_rates(n, p)
    return -((1+p.yield_*n)*vp+rp+p.reserve_loss)/c


def feedback_jacobian(c, n, *, yield_, turnover, reserve_loss,
                      v, vp, r, r_n, r_c, production_prime):
    """General r(C,n): exact full Jacobian, derivatives r_c at FIXED n.

    At equilibrium, det(J)=-yield_*vp*C*H'(C), with H=P/C-v-r-delta*n.
    No VISIT-derived coefficients are introduced here.
    """
    return np.array([[yield_*(v-n*vp)-turnover, yield_*vp],
                     [production_prime-v+n*vp-r-c*r_c+n*r_n,
                      -vp-r_n-reserve_loss]])
