"""Analytic checks for a published grass--forest cover model, not native VISIT.

Kumar K & Dutta (2026), doi:10.1098/rspa.2025.0803, (2.3), (4.1).
The constructed environmental path and reflected noise are our assumptions.
Time is in model units; no conversion to years or temperature is calibrated.
"""
from dataclasses import dataclass, replace
from functools import lru_cache

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import expit


@dataclass(frozen=True)
class Parameters:
    alpha: float = 1.2
    phi0: float = 0.1709
    phi1: float = 0.9
    theta: float = 0.432
    width: float = 0.082

    def __post_init__(self):
        if not (self.alpha > 0 and self.phi1 > self.phi0 > 0
                and self.width > 0 and 0 < self.theta < 1):
            raise ValueError("Require positive rates, phi1 > phi0, and 0 < theta < 1.")


def logistic(g, p=Parameters()):
    return expit((np.asarray(g) - p.theta) / p.width)


def phi(g, p=Parameters()):
    return p.phi0 + (p.phi1 - p.phi0) * logistic(g, p)


def phi_prime(g, p=Parameters()):
    z = logistic(g, p)
    return (p.phi1 - p.phi0) * z * (1 - z) / p.width


def h(g, p=Parameters()):
    return phi(g, p) - p.alpha * np.asarray(g)


def drift(g, p=Parameters()):
    return (1 - np.asarray(g)) * h(g, p)


def drift_prime(g, p=Parameters()):
    return -h(g, p) + (1 - np.asarray(g)) * (phi_prime(g, p) - p.alpha)


def turning_points(p=Parameters()):
    """All interior critical points of h, analytically, rather than grid search."""
    ratio = p.alpha * p.width / (p.phi1 - p.phi0)
    if ratio >= 0.25:
        return []
    z = (1 + np.array([-1., 1.]) * np.sqrt(1 - 4 * ratio)) / 2
    g = p.theta + p.width * np.log(z / (1 - z))
    return [float(x) for x in g if 0 < x < 1]


def interior_roots(p=Parameters()):
    """Find roots of h in (0,1); the separate boundary equilibrium g=1 is omitted."""
    edges = [0., *turning_points(p), 1.]
    roots = []
    for x in edges[1:-1]:
        if abs(h(x, p)) < 1e-13:
            roots.append(x)
    for a, b in zip(edges[:-1], edges[1:]):
        if h(a, p) * h(b, p) < 0:
            roots.append(brentq(lambda x: h(x, p), a, b, xtol=5e-15))
    return np.array(sorted(set(roots)))


def fold_points(p=Parameters()):
    """The two folds for the baseline family varying alpha alone."""
    equation = lambda g: g * phi_prime(g, p) - phi(g, p)
    roots = [brentq(equation, 1e-12, p.theta, xtol=5e-15),
             brentq(equation, p.theta, 1 - 1e-12, xtol=5e-15)]
    return [(g, float(phi(g, p) / g)) for g in roots]


def path_parameters(lam):
    """Analytically designed single-driver path; NOT an empirical response law.

    Exact low stable / unstable roots are a=.4-.1*lam and b=a+.02.
    No arbitrary multiplication or modification of the published vector field.
    """
    if not 0 <= lam <= 1:
        raise ValueError("The proved environmental path is restricted to [0,1].")
    a = 0.4 - 0.1 * lam
    b = a + 0.02
    la, lb = logistic(a), logistic(b)
    denominator = a * lb - b * la
    p = Parameters()
    delta = p.phi0 * (b - a) / denominator
    alpha = p.phi0 * (lb - la) / denominator
    return replace(p, alpha=float(alpha), phi1=float(p.phi0 + delta))


def smooth_step(u):
    u = np.clip(u, 0, 1)
    # Bound floating-point endpoint roundoff; do not clip the ecosystem state.
    return np.clip(u**3 * (10 - 15 * u + 6 * u**2), 0, 1)


def ramp_solution(duration, *, rtol=2e-10, atol=2e-12, max_step=0.01,
                  sensitivity=False):
    """Integrate in u=t/duration, starting at the exact past equilibrium."""
    if duration <= 0:
        raise ValueError("duration must be positive")

    def rhs(u, y):
        p = path_parameters(float(smooth_step(u)))
        f = float(drift(y[0], p))
        if sensitivity:
            return [duration * f,
                    f + duration * float(drift_prime(y[0], p)) * y[1]]
        return [duration * f]

    sol = solve_ivp(rhs, (0, 1), [0.4, 0.] if sensitivity else [0.4],
                    method="DOP853", rtol=rtol, atol=atol, max_step=max_step,
                    dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def analytic_bounds():
    """Conservative explicit proof bounds, NOT measured critical durations."""
    p = Parameters()
    l0, l1 = float(logistic(.3)), float(logistic(.42))
    dmin = .02 * (.3 * l0 * (1 - l0) / p.width - float(logistic(.4)))
    dmax = .4 * l1
    delta_min = p.phi0 * .02 / dmax
    delta_max = p.phi0 * .02 / dmin
    lpp_min = l0 * (1 - l1) * (1 - 2 * l1) / p.width**2
    epsilon = .005
    m = (1 - .405) * delta_min * lpp_min * epsilon * (.02 - epsilon) / 2
    speed_bound = p.phi0 + delta_max + delta_max / (4 * p.width)
    return dict(denominator_min=dmin, denominator_max=dmax,
                delta_min=delta_min, delta_max=delta_max,
                logistic_second_derivative_min=lpp_min, tube_epsilon=epsilon,
                inward_drift_bound=m, absolute_drift_bound=speed_bound,
                fast_duration_upper=.08 / speed_bound,
                slow_duration_lower=.1875 / m)


def matched_solutions(duration, *, match=.5, rtol=2e-10, atol=2e-12,
                      max_step=.01, sensitivity=False):
    """Shoot forward from low past QSE and backward from future basin boundary.

    Direct forward endpoint shooting is exponentially ill-conditioned at threshold.
    Match in the interior instead. The analytic threshold remains g(1;tau)=.32.
    """
    if duration <= 0 or not 0 < match < 1:
        raise ValueError("Require duration > 0 and 0 < match < 1.")

    def rhs(u, y):
        p = path_parameters(float(smooth_step(u)))
        f = float(drift(y[0], p))
        return ([duration * f,
                 f + duration * float(drift_prime(y[0], p)) * y[1]]
                if sensitivity else [duration * f])

    options = dict(method="DOP853", rtol=rtol, atol=atol, max_step=max_step,
                   dense_output=True)
    forward = solve_ivp(rhs, [0, match], [.4, 0.] if sensitivity else [.4], **options)
    backward = solve_ivp(rhs, [1, match], [.32, 0.] if sensitivity else [.32], **options)
    if not (forward.success and backward.success):
        raise RuntimeError("Forward/backward matching failed.")
    return forward, backward


def critical_duration(**options):
    """Numerical value of the analytically unique endpoint separatrix crossing."""
    bounds = analytic_bounds()
    def residual(tau):
        fw, bw = matched_solutions(tau, **options)
        return fw.y[0, -1] - bw.y[0, -1]
    return brentq(residual, bounds["fast_duration_upper"] / 2,
                  bounds["slow_duration_lower"] * 1.1, xtol=1e-9, rtol=1e-13)


@lru_cache(maxsize=None)
def _gauss(n):
    return np.polynomial.legendre.leggauss(n)


def potential(g, p=Parameters(), order=80):
    """V(0)=0, V'=-f, using Gaussian quadrature."""
    z, weights = _gauss(order)
    g = np.asarray(g)
    points = g[..., None] * (z + 1) / 2
    return -g / 2 * np.sum(weights * drift(points, p), axis=-1)


def mean_first_passage(x, target, sigma, p=Parameters(), order=100):
    """Exact integral formula evaluated by quadrature: reflect 0, absorb target.

    Continued stochastic dynamics can also return from the second basin.
    This is not a deterministic or irreversible transition threshold.
    """
    if not (0 <= x < target <= 1 and sigma > 0):
        raise ValueError("Require 0 <= x < target <= 1 and sigma > 0.")
    nodes, weights = _gauss(order)
    y = x + (target - x) * (nodes + 1) / 2
    z = y[:, None] * (nodes[None, :] + 1) / 2
    # Combine exponentials before integration to cancel arbitrary V offsets.
    exponent = 2 * (potential(y, p)[:, None] - potential(z, p)) / sigma**2
    inner = y / 2 * np.sum(weights[None, :] * np.exp(exponent), axis=1)
    return float((target - x) / sigma**2 * np.sum(weights * inner))
