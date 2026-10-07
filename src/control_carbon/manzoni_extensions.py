"""Manzoni (2014) Appendix A and soil/root-continuity extensions.

Same units as manzoni_hydraulics. Soil coefficients are required explicitly:
the supplied article gives a scaling law, not a complete numerical calibration.
"""

from dataclasses import dataclass, replace
from math import fsum

import numpy as np
from numpy.polynomial import polynomial as poly
from scipy.optimize import brentq
from scipy.special import expit

from .manzoni_hydraulics import (
    HydraulicParameters, _fluxes, _forcing, _potential, _state,
    _state_jacobian, _storage_rhs, conductances, hydraulic_jacobian, hydraulic_rhs,
)


@dataclass(frozen=True)
class SigmoidShape:
    alpha_x: float
    alpha_s: float

    def __post_init__(self) -> None:
        if not np.all(np.isfinite([self.alpha_x, self.alpha_s])):
            raise ValueError("sigmoid exponents must be finite")
        if min(self.alpha_x, self.alpha_s) <= 0:
            raise ValueError("sigmoid exponents must be positive")


def sigmoid_conductances(
    psi_leaf: float, p: HydraulicParameters, shape: SigmoidShape,
) -> tuple[float, float]:
    """Equations (A1),(A2), evaluated without overflowing powers."""
    _potential(psi_leaf)
    if psi_leaf == 0:
        return p.gx_max, p.gs_max
    log_u = np.log(-psi_leaf)
    return (
        p.gx_max * expit(-shape.alpha_x * (log_u-np.log(-p.psi50_x))),
        p.gs_max * expit(-shape.alpha_s * (log_u-np.log(-p.psi50_s))),
    )


def sigmoid_rhs(
    state: np.ndarray, psi_root: float, vpd: float,
    p: HydraulicParameters, shape: SigmoidShape,
) -> np.ndarray:
    state = _state(state)
    _forcing(psi_root, vpd)
    gx, gs = sigmoid_conductances(state[1], p, shape)
    return _storage_rhs(_fluxes(state, psi_root, vpd, gx, gs), p)


def sigmoid_jacobian(
    state: np.ndarray, psi_root: float, vpd: float,
    p: HydraulicParameters, shape: SigmoidShape,
) -> np.ndarray:
    state = _state(state)
    _forcing(psi_root, vpd)
    leaf = state[1]
    gx, gs = sigmoid_conductances(leaf, p, shape)
    derivatives = []
    for maximum, value, half, exponent in [
        (p.gx_max, gx, -p.psi50_x, shape.alpha_x),
        (p.gs_max, gs, -p.psi50_s, shape.alpha_s),
    ]:
        if leaf == 0:
            if exponent < 1:
                raise ValueError("sigmoid derivative is singular at zero potential for alpha < 1")
            derivative = maximum / half if exponent == 1 else 0.
        else:
            weight = expit(exponent * (np.log(-leaf)-np.log(half)))
            derivative = value * exponent * weight / (-leaf)
        derivatives.append(derivative)
    return _state_jacobian(state, psi_root, vpd, p, gx, *derivatives)


def sigmoid_steady_fluxes(
    psi_leaf: float, psi_root: float, vpd: float,
    p: HydraulicParameters, shape: SigmoidShape,
) -> tuple[float, float]:
    _forcing(psi_root, vpd)
    gx, gs = sigmoid_conductances(psi_leaf, p, shape)
    return gx*(psi_root-psi_leaf), gs*vpd


def sigmoid_log_flux_ratio(
    psi_leaf: float, psi_root: float, vpd: float,
    p: HydraulicParameters, shape: SigmoidShape,
) -> float:
    """Log(supply/demand), also when both physical fluxes underflow."""
    _potential(psi_leaf)
    _forcing(psi_root, vpd)
    if psi_leaf >= psi_root or vpd == 0:
        raise ValueError("log flux ratio requires positive supply and demand")
    log_u = np.log(-psi_leaf)
    return float(
        np.log(p.gx_max/p.gs_max)+np.log(psi_root-psi_leaf)-np.log(vpd)
        + np.logaddexp(0., shape.alpha_s*(log_u-np.log(-p.psi50_s)))
        - np.logaddexp(0., shape.alpha_x*(log_u-np.log(-p.psi50_x)))
    )


def sigmoid_equilibrium_jacobian_scaled(
    psi_leaf: float, psi_root: float, vpd: float,
    p: HydraulicParameters, shape: SigmoidShape,
) -> np.ndarray:
    """Full equilibrium Jacobian divided by positive gx; preserves signs.

    Unlike absolute eigenvalues, this remains representable at extremely dry
    finite roots. Only valid at a checked positive-demand equilibrium.
    """
    residual = sigmoid_log_flux_ratio(psi_leaf, psi_root, vpd, p, shape)
    if abs(residual) > 1e-8:
        raise ValueError("scaled equilibrium Jacobian requires flux balance")
    u = -psi_leaf
    wx = expit(shape.alpha_x*(np.log(u)-np.log(-p.psi50_x)))
    ws = expit(shape.alpha_s*(np.log(u)-np.log(-p.psi50_s)))
    return np.array([
        [-4/p.cx, 2/p.cx],
        [2/p.cl, (shape.alpha_x*wx*(psi_root-psi_leaf)/u-2
                   - shape.alpha_s*ws*(psi_root-psi_leaf)/u)/p.cl],
    ])


def _exponential_roots(terms: dict[float, float]) -> tuple[float, ...]:
    """Isolate real roots of a finite exponential sum, including tangencies.

    Divide out its least exponential; derivative recursion removes one term.
    Dominant-term bounds bracket both infinite tails, so no dry-end cutoff
    or sampling grid determines the root count. Near-zero stationary values
    use a relative roundoff tolerance; this is floating-point, not interval
    arithmetic certification.
    """
    items = sorted((power, value) for power, value in terms.items() if value != 0)
    if len(items) < 2:
        return ()
    powers = np.array([power-items[0][0] for power, value in items])
    coefficients = np.array([value for power, value in items])
    logs = np.log(abs(coefficients))
    signs = np.sign(coefficients)
    count = len(items)
    low = float(np.min((logs[0]-np.log(2*(count-1))-logs[1:])/powers[1:])) - 1.
    high = float(np.max(
        (np.log(2*(count-1))+logs[:-1]-logs[-1])/(powers[-1]-powers[:-1]),
    )) + 1.

    def evaluate(z):
        exponents = logs + powers*z
        return fsum(signs*np.exp(exponents-np.max(exponents)))

    critical = _exponential_roots({
        float(power): float(coefficient*power)
        for power, coefficient in zip(powers[1:], coefficients[1:])
    })
    nodes = [low, *[z for z in critical if low < z < high], high]
    values = [evaluate(z) for z in nodes]
    zeros = [abs(value) <= 128*np.finfo(float).eps for value in values]
    roots = [z for z, zero in zip(nodes, zeros) if zero]
    for left, right, fl, fr, zl, zr in zip(
        nodes[:-1], nodes[1:], values[:-1], values[1:], zeros[:-1], zeros[1:],
    ):
        if not zl and not zr and fl*fr < 0:
            roots.append(brentq(evaluate, left, right, xtol=1e-12, rtol=1e-14))
    return tuple(sorted(roots))


def sigmoid_equilibria(
    psi_root: float, vpd: float, p: HydraulicParameters, shape: SigmoidShape,
) -> tuple[float, ...]:
    """All representable finite leaf equilibria, wettest first.

    For D>0, u=-psi_leaf>R=-psi_root solves
    (u-R)(1+(u/|psi50_s|)^alpha_s)
        - t(1+(u/|psi50_x|)^alpha_x) = 0.
    A log-domain exponential sum isolates roots on the entire positive axis.
    """
    _forcing(psi_root, vpd)
    if vpd == 0:
        return (float(psi_root),)
    r, t = -psi_root, p.gs_max*vpd/p.gx_max
    x, s = -p.psi50_x, -p.psi50_s
    terms: dict[float, float] = {}
    for power, coefficient in [
        (0., -r-t), (1., 1.),
        (shape.alpha_s, -r/s**shape.alpha_s),
        (shape.alpha_s+1., 1./s**shape.alpha_s),
        (shape.alpha_x, -t/x**shape.alpha_x),
    ]:
        terms[power] = terms.get(power, 0.) + coefficient
    logs = _exponential_roots(terms)
    if any(z > np.log(np.finfo(float).max) for z in logs):
        raise OverflowError("finite equilibrium potential exceeds floating-point range")
    return tuple(-float(np.exp(z)) for z in logs if np.exp(z) > r)


@dataclass(frozen=True)
class SoilRootParameters:
    gsr_sat: float
    psi_sat: float

    def __post_init__(self) -> None:
        if not np.all(np.isfinite([self.gsr_sat, self.psi_sat])):
            raise ValueError("soil-root parameters must be finite")
        if self.gsr_sat <= 0 or self.psi_sat >= 0:
            raise ValueError("gsr_sat must be positive and psi_sat negative")


def soil_root_conductance(psi_soil: float, soil: SoilRootParameters) -> float:
    """Section 2 scaling: gsr = gsr_sat*(psi_sat/psi_soil)^2."""
    _potential(psi_soil)
    if psi_soil > soil.psi_sat:
        raise ValueError("soil potential must not exceed the prescribed saturation potential")
    return soil.gsr_sat*(soil.psi_sat/psi_soil)**2


def root_node(
    state: np.ndarray, psi_soil: float, p: HydraulicParameters,
    soil: SoilRootParameters,
) -> float:
    """Instantaneous Esr=Erx; root node has no storage."""
    state = _state(state)
    gsr = soil_root_conductance(psi_soil, soil)
    gx, gs = conductances(state[1], p)
    return (gsr*psi_soil+2*gx*state[0])/(gsr+2*gx)


def soil_hydraulic_rhs(
    state: np.ndarray, psi_soil: float, vpd: float,
    p: HydraulicParameters, soil: SoilRootParameters,
) -> np.ndarray:
    return hydraulic_rhs(state, root_node(state, psi_soil, p, soil), vpd, p)


def soil_hydraulic_jacobian(
    state: np.ndarray, psi_soil: float, vpd: float,
    p: HydraulicParameters, soil: SoilRootParameters,
) -> np.ndarray:
    state = _state(state)
    root = root_node(state, psi_soil, p, soil)
    jacobian = hydraulic_jacobian(state, root, vpd, p)
    gx, gs = conductances(state[1], p)
    gsr = soil_root_conductance(psi_soil, soil)
    dgx = -p.gx_max/(2*p.psi50_x) if gx > 0 else 0.
    jacobian[0, 0] += 2*gx/p.cx * (2*gx/(gsr+2*gx))
    jacobian[0, 1] += 2*gx/p.cx * (2*dgx*(state[0]-root)/(gsr+2*gx))
    return jacobian


def soil_steady_fluxes(
    psi_leaf: float, psi_soil: float, vpd: float,
    p: HydraulicParameters, soil: SoilRootParameters,
) -> tuple[float, float]:
    _forcing(psi_soil, vpd)
    gx, gs = conductances(psi_leaf, p)
    gsr = soil_root_conductance(psi_soil, soil)
    return gx*gsr/(gx+gsr)*(psi_soil-psi_leaf), gs*vpd


def soil_equilibria(
    psi_soil: float, vpd: float, p: HydraulicParameters, soil: SoilRootParameters,
) -> tuple[float, ...]:
    """Active leaf equilibria of the series soil+whole-xylem pathway."""
    _forcing(psi_soil, vpd)
    gsr = soil_root_conductance(psi_soil, soil)
    a, b, r = -2*p.psi50_x, -2*p.psi50_s, -psi_soil
    t, ratio = p.gs_max*vpd/p.gx_max, p.gx_max/gsr
    supply = poly.polymul([1., -1/a], [-r, 1.])
    demand = t*poly.polymul([1., -1/b], [1+ratio, -ratio/a])
    coefficients = poly.polysub(supply, demand)
    c, linear, quadratic = coefficients
    discriminant = linear**2-4*quadratic*c
    tolerance = 64*np.finfo(float).eps*max(linear**2, abs(4*quadratic*c), 1.)
    if abs(discriminant) <= tolerance:
        roots = (-linear/(2*quadratic),)
    elif discriminant < 0:
        roots = ()
    else:
        q = -.5*(linear+np.copysign(np.sqrt(discriminant), linear))
        roots = (q/quadratic, c/q)
    return tuple(sorted(
        (-float(z) for z in roots if r <= z < min(a, b)),
        reverse=True,
    ))


def soil_equilibrium_state(
    psi_leaf: float, psi_soil: float, p: HydraulicParameters, soil: SoilRootParameters,
) -> np.ndarray:
    gx, gs = conductances(psi_leaf, p)
    gsr = soil_root_conductance(psi_soil, soil)
    flux = gx*gsr/(gx+gsr)*(psi_soil-psi_leaf)
    if gx == 0:
        raise ValueError("nonconducting equilibrium has no unique xylem state")
    root = psi_soil-flux/gsr
    return np.array([(root+psi_leaf)/2, psi_leaf])


def soil_fold(
    psi_soil: float, p: HydraulicParameters, soil: SoilRootParameters,
) -> tuple[float, float] | None:
    """Active fold with prescribed soil, including the series resistance."""
    gsr = soil_root_conductance(psi_soil, soil)
    a, b, r = -2*p.psi50_x, -2*p.psi50_s, -psi_soil
    if b <= a or r >= a:
        return None
    ratio = p.gx_max/gsr
    numerator = poly.polymul([1., -1/a], [-r, 1.])
    denominator = poly.polymul([1., -1/b], [1+ratio, -ratio/a])
    derivative = poly.polysub(
        poly.polymul(poly.polyder(numerator), denominator),
        poly.polymul(numerator, poly.polyder(denominator)),
    )
    u = brentq(lambda value: poly.polyval(value, derivative), r, a, xtol=1e-13)
    critical = p.gx_max/p.gs_max*poly.polyval(u, numerator)/poly.polyval(u, denominator)
    return float(critical), -float(u)


def deteriorated_traits(p: HydraulicParameters, fraction: float) -> HydraulicParameters:
    """Section 3.3: fixed damage factor, not an extra dynamic state."""
    if not np.isfinite(fraction) or not 0 < fraction <= 1:
        raise ValueError("deterioration fraction must lie in (0, 1]")
    return replace(p, gx_max=fraction*p.gx_max, psi50_x=fraction*p.psi50_x)
