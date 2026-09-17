"""Analytical, uncalibrated carbon-water R-tipping benchmarks.

Not a VISIT adapter. Water uptake proportional to W*C**2 follows the
spatially homogeneous, ungrazed Klausmeier family (Siero et al. 2019,
doi:10.1086/701669, A6). Quasi-steady elimination, the paired environmental
path, and the finite-ramp criterion are derived in water_tipping_minimal.tex.
Carbon c is C/C0; time is years when mortality is per year.
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq


@dataclass(frozen=True)
class MinimalTipping:
    productivity: float = 2.05  # a > 2, dimensionless
    mortality: float = 0.1  # per year; illustrative, NOT calibrated
    final_scale: float = 2.0  # Lambda > 1

    def __post_init__(self):
        if not all(np.isfinite(v) for v in (
            self.productivity, self.mortality, self.final_scale
        )):
            raise ValueError("parameters must be finite")
        if self.productivity <= 2 or self.mortality <= 0 or self.final_scale <= 1:
            raise ValueError("require a > 2, mortality > 0, final_scale > 1")

    @property
    def roots(self):
        upper = (self.productivity + np.sqrt(self.productivity**2 - 4)) / 2
        return 1 / upper, upper

    @property
    def comoving_fold_rate(self):
        """Infinite-ramp lower bound, NOT the finite-shift critical rate."""
        return self.mortality * (self.productivity / 2 - 1)


def scalar_rhs(carbon, scale, model=MinimalTipping()):
    """dc/dt after quasi-steady water elimination, c >= 0, scale > 0."""
    return model.mortality * carbon * (
        model.productivity * scale * carbon / (scale**2 + carbon**2) - 1
    )


def moving_rhs(relative_carbon, rate, model=MinimalTipping()):
    """x=c/lambda for lambda=exp(rate*t)."""
    x = relative_carbon
    return x * (
        model.mortality * (model.productivity * x / (1 + x*x) - 1) - rate
    )


def crossing_log_shift(rate, model=MinimalTipping()):
    """Log shift needed to cross the final frozen-system separatrix.

    Infinity below the moving-frame saddle-node; above it integrate the
    exact separation-of-variables expression, not a fitted rate formula.
    """
    if not np.isfinite(rate) or rate <= 0:
        raise ValueError("rate must be positive and finite")
    if rate <= model.comoving_fold_rate:
        return np.inf
    lower, upper = model.roots
    def integrand(x):
        intrinsic = model.mortality * (model.productivity*x/(1+x*x)-1)
        return rate / (x * (rate-intrinsic))
    return quad(integrand, lower, upper, points=[1.0], epsabs=1e-10,
                epsrel=1e-10)[0]


def finite_critical_rate(model=MinimalTipping()):
    """Unique finite-shift critical rate, or inf if the shift is too small."""
    lower, upper = model.roots
    if model.final_scale <= upper/lower:
        return np.inf
    log_shift = np.log(model.final_scale)
    left = model.comoving_fold_rate * (1 + 1e-4)
    # A shift extremely large may require a bracket closer to the fold.
    while crossing_log_shift(left, model) <= log_shift:
        left = (left + model.comoving_fold_rate) / 2
        if left == model.comoving_fold_rate:
            raise ValueError("critical rate too close to fold for float precision")
    right = max(2*left, model.mortality)
    while crossing_log_shift(right, model) > log_shift:
        right *= 2
    return brentq(lambda r: crossing_log_shift(r, model)-log_shift,
                  left, right, xtol=1e-13)


def advance_relative(x, duration, rate, model=MinimalTipping(), *, rtol=1e-10):
    """One month/year (or remainder) flow map via adaptive scalar integration.

    This is an accuracy-controlled approximation to the interval flow map,
    not one forward-Euler update and not a daily physiological simulation.
    """
    if x < 0 or duration < 0 or rate < 0 or not np.all(np.isfinite([x, duration, rate])):
        raise ValueError("finite nonnegative x, duration and rate required")
    if duration == 0 or x == 0:
        return float(x)
    result = solve_ivp(lambda t, z: [moving_rhs(z[0], rate, model)],
                       (0, duration), [x], method="DOP853", rtol=rtol,
                       atol=rtol*1e-3, max_step=min(duration, 1/model.mortality))
    if not result.success or result.y[0, -1] < 0:
        raise RuntimeError("flow integration failed or violated positivity")
    return float(result.y[0, -1])


def simulate_ramp(rate, model=MinimalTipping(), *, update_years=1.,
                  hold_years=400., rtol=1e-10):
    """Monthly/yearly endpoint states; exactly split at the ramp endpoint."""
    if rate <= 0 or update_years <= 0 or hold_years < 0 or not np.all(
        np.isfinite([rate, update_years, hold_years])
    ):
        raise ValueError("invalid rate, update interval or hold")
    stop = np.log(model.final_scale)/rate
    times, xs, scales = [0.], [model.roots[1]], [1.]
    x = xs[0]
    for begin, end, velocity in [(0., stop, rate), (stop, stop+hold_years, 0.)]:
        t = begin
        while t < end - 1e-11:
            step = min(update_years, end-t)
            x = advance_relative(x, step, velocity, model, rtol=rtol)
            t += step
            times.append(t)
            xs.append(x)
            scales.append(np.exp(rate*min(t, stop)))
    times, xs, scales = np.array(times), np.array(xs), np.array(scales)
    end_index = int(np.argmin(abs(times-stop)))
    margin = xs[end_index]-model.roots[0]
    return dict(time=times, carbon=xs*scales, relative=xs, scale=scales,
                stop=stop, margin=float(margin),
                outcome="tipped" if margin < 0 else "tracking")


def dynamic_rhs(state, rainfall, mortality=1.8, water_loss=1.):
    """Two-state homogeneous Klausmeier benchmark, normalized water.

    c'=m(w*c*c-c); w'=ell(a-(1+c*c)*w). No explicit seasonal cycle.
    Rainfall a > 2 supports a vegetated equilibrium if trace < 0.
    """
    c, w = state
    return np.array([mortality*(w*c*c-c), water_loss*(rainfall-(1+c*c)*w)])


def dynamic_jacobian(state, mortality=1.8, water_loss=1.):
    c, w = state
    return np.array([[mortality*(2*w*c-1), mortality*c*c],
                     [-2*water_loss*w*c, -water_loss*(1+c*c)]])


def simulate_drying(duration, *, rainfall_start=10., rainfall_end=2.002,
                    mortality=1.8, water_loss=1., hold=150., rtol=1e-9,
                    profile="linear"):
    """Illustrative finite linear rainfall ramp; not forest calibration.

    Log-carbon integration enforces positive live carbon without clipping.
    The final invariant rectangle c < 1/a, 0 < w < a certifies the bare basin.
    """
    if (not np.all(np.isfinite([duration, rainfall_start, rainfall_end,
                               mortality, water_loss, hold, rtol]))
            or duration <= 0 or hold <= 0 or mortality <= 0 or water_loss <= 0
            or not rainfall_start > rainfall_end > 2 or rtol <= 0):
        raise ValueError("require finite positive controls and start > end > 2")
    if profile not in ("linear", "smoothstep"):
        raise ValueError("profile must be linear or smoothstep")
    c0 = (rainfall_start + np.sqrt(rainfall_start**2-4))/2
    def rhs(t, z, a):
        c = np.exp(z[0]); w = z[1]
        return [mortality*(w*c-1), water_loss*(a-(1+c*c)*w)]
    def ramp_rainfall(t):
        s=t/duration
        fraction=s if profile == "linear" else 3*s*s-2*s*s*s
        return rainfall_start+(rainfall_end-rainfall_start)*fraction
    ramp = solve_ivp(lambda t,z: rhs(t,z,ramp_rainfall(t)),
                     [0,duration], [np.log(c0),1/c0], method="Radau",
                     rtol=rtol, atol=rtol*1e-2, dense_output=True,
                     max_step=max(duration/200, .01))
    after = solve_ivp(lambda t,z: rhs(t,z,rainfall_end), [duration,duration+hold],
                      ramp.y[:,-1], method="Radau", rtol=rtol, atol=rtol*1e-2,
                      max_step=1., dense_output=True)
    if not ramp.success or not after.success:
        raise RuntimeError("dynamic water integration failed")
    state = np.array([np.exp(after.y[0,-1]), after.y[1,-1]])
    certified_bare = bool(state[0] < 1/rainfall_end and 0 < state[1] <= rainfall_end*(1+1e-10))
    return dict(state=state, ramp=ramp, after=after, certified_bare=certified_bare)
