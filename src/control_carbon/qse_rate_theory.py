"""Algebra and quadrature supporting the QSE/rate-tipping proofs.

Model-agnostic analytic constructions, NOT VISIT closures or forest calibration.
States and time are nondimensional. No ODE simulation is needed for certificates.
See docs/qse_rate_tipping_analytic_proofs.tex for proofs and assumptions.
"""
import math
from scipy.integrate import quad
from scipy.optimize import brentq


COMOVING_FOLD_RATE = 2/(3*math.sqrt(3))


def cubic_field(x, lam):
    return -(x-lam-1)*(x-lam-2)*(x-lam-3)


def crossing_displacement(rate):
    """H(r)=r*first crossing time, r>2/(3 sqrt(3)).

    H is strictly decreasing from infinity to one. H(r)=D is the UNIQUE
    critical finite-linear-ramp rate when total displacement D>1.
    """
    if not math.isfinite(rate) or rate <= COMOVING_FOLD_RATE:
        raise ValueError("crossing formula requires rate above comoving fold")
    return quad(lambda y: rate/(rate-y+y**3),0,1,
                points=[1/math.sqrt(3)],epsabs=2e-11,epsrel=2e-11)[0]


def critical_rate(displacement=1.5):
    """Evaluate an analytically unique threshold, not a simulation classifier."""
    if not math.isfinite(displacement) or displacement <= 1:
        raise ValueError("finite critical rate exists only for displacement > 1")
    lo=COMOVING_FOLD_RATE*(1+1e-6)
    hi=1.0
    while crossing_displacement(hi)>displacement:
        hi*=2
    # Avoid silently giving a bracket for extreme D that was not certified.
    if crossing_displacement(lo)<displacement:
        raise ValueError("requested displacement exceeds quadrature bracket")
    return brentq(lambda r:crossing_displacement(r)-displacement,lo,hi,xtol=1e-12)


def bounded_input_coefficients(lam):
    """Two bounded production terms, linear loss, and zero production at x=0.

    P=p*x/(1+x)+a*x^2/(h2+x^2), R=x.
    F=-x*(x-L)*(x-U)*(x-H)/[(1+x)*(h2+x^2)].
    These coefficients were constructed algebraically, not fitted to data.
    """
    if not math.isfinite(lam) or lam<0:
        raise ValueError("finite nonnegative environmental parameter required")
    s=lam+2
    sigma1=3*s
    sigma2=3*s*s-1
    sigma3=s**3-s
    a=(sigma1-sigma2+math.sqrt((sigma1+sigma2)**2-4*sigma3))/2
    p=sigma1+1-a
    return p,a,a+sigma2


def bounded_production(x,lam):
    p,a,h2=bounded_input_coefficients(lam)
    return p*x/(1+x)+a*x*x/(h2+x*x)


def bounded_field(x,lam):
    """Factored net balance, evaluated without subtracting near-equal fluxes."""
    _,_,h2=bounded_input_coefficients(lam)
    return x*cubic_field(x,lam)/((1+x)*(h2+x*x))


def smooth_certificates(displacement=1.5,bounded=False):
    """Conservative sufficient durations for a quintic monotone finite ramp.

    s(u)=10u^3-15u^4+6u^5 has max s'=15/8.
    No claim of a unique smooth-ramp threshold is made.
    """
    d=displacement
    if not math.isfinite(d) or not 1<d<2:
        raise ValueError("source/sink reversal construction requires 1<D<2")
    if bounded:
        # On the tracking tube x>=2.5; x<=D+3; h2<sigma1+sigma2.
        smax=d+2
        h2_upper=3*smax+3*smax*smax-1
        qmin=2.5/((d+4)*(h2_upper+(d+3)**2))
        qmax=1/11  # h2>=11 and x/(1+x)<=1 on the invariant rectangle.
    else:
        qmin=qmax=1.0
    max_cubic=(d+1)**3-(d+1)
    return dict(q_lower=qmin,q_upper=qmax,
                fast_duration_strict_upper=(d-1)/(qmax*max_cubic),
                slow_duration_lower=5*d/qmin,
                net_stock_change_fast=d-2,net_stock_change_slow=d)
