"""Stage B1 synthetic structural flux feedback; Stage A remains unchanged.

MPa, day, kg C/m2 ground, D mol/mol. Fixed capacities; no Leuning,
soil state, damage dynamics or mortality regime. Analysis is restricted
to positive viable leaf carbon, not to a physiological dead-plant state.
"""

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root

from .manzoni_hydraulics import HydraulicParameters, _forcing, _state, conductances
from .spac_carbon import (
    CarbonParameters, combined_jacobian, carbon_rhs, hydraulic_factor,
    positive_carbon_equilibrium,
)


@dataclass(frozen=True)
class StructuralFeedback:
    leaf_reference: float
    wood_reference: float
    leaf_floor: float
    minimum_leaf: float
    leaf_supply: float = 0.
    leaf_demand: float = 0.
    wood_supply: float = 0.
    wood_exponent: float = 1.

    def __post_init__(self):
        values = tuple(self.__dict__.values())
        if not np.all(np.isfinite(values)):
            raise ValueError("structural feedback parameters must be finite")
        if min(self.leaf_reference, self.wood_reference, self.minimum_leaf, self.wood_exponent) <= 0 or self.leaf_floor < 0:
            raise ValueError("references, minimum leaf and exponent must be positive; floor nonnegative")
        if self.minimum_leaf >= self.leaf_reference:
            raise ValueError("minimum viable leaf must be below fixed leaf reference")
        if min(self.leaf_supply, self.leaf_demand, self.wood_supply) < 0 or max(
            self.leaf_supply, self.leaf_demand, self.wood_supply,
        ) > 1:
            raise ValueError("feedback strengths must lie in [0,1]")


def feedback_factors(leaf: float, wood: float, p: StructuralFeedback):
    """Return FL, FE, FW, dFL/dL, dFE/dL, dFW/dW.

    Strength zero means no feedback. Leaf floor is in supply only; no
    extra LAI multiplier enters the unchanged Stage A photosynthesis.
    """
    if not np.all(np.isfinite([leaf, wood])) or leaf <= 0 or wood < 0:
        raise ValueError("feedback evaluation requires positive leaf and nonnegative wood carbon")
    fl = 1+p.leaf_supply*((leaf+p.leaf_floor)/(p.leaf_reference+p.leaf_floor)-1)
    fe = 1+p.leaf_demand*(leaf/p.leaf_reference-1)
    if p.wood_supply == 0:
        fw, dfw = 1., 0.
    else:
        if wood == 0 and p.wood_exponent < 1:
            raise ValueError("wood feedback derivative is singular at zero for exponent < 1")
        fw = 1+p.wood_supply*((wood/p.wood_reference)**p.wood_exponent-1)
        dfw = p.wood_supply*p.wood_exponent/p.wood_reference*(wood/p.wood_reference)**(p.wood_exponent-1)
    return fl, fe, fw, p.leaf_supply/(p.leaf_reference+p.leaf_floor), p.leaf_demand/p.leaf_reference, dfw


def effective_conductance(g_rx: float, g_xl: float) -> float:
    if not np.all(np.isfinite([g_rx, g_xl])) or min(g_rx, g_xl) < 0:
        raise ValueError("series conductances must be finite nonnegative")
    return g_rx*g_xl/(g_rx+g_xl) if g_rx+g_xl > 0 else 0.


def _combined(state):
    state = np.asarray(state, dtype=float)
    if state.shape != (6,) or not np.all(np.isfinite(state)) or np.any(state[2:] < 0) or state[2] <= 0:
        raise ValueError("six-state vector requires finite nonnegative carbon and positive leaf")
    _state(state[:2])
    return state


def structural_fluxes(state, root_potential, vpd, hydraulic, feedback):
    state = _combined(state)
    _forcing(root_potential, vpd)
    gx, gs = conductances(state[1], hydraulic)
    fl, fe, fw, *derivatives = feedback_factors(state[2], state[3], feedback)
    return np.array([
        2*gx*fw*(root_potential-state[0]),
        2*gx*fl*(state[0]-state[1]),
        gs*fe*vpd,
    ])


def two_way_rhs(state, root_potential, vpd, hydraulic, carbon, feedback):
    """Formal positive-leaf vector field; integrator stops at L_min."""
    state = _combined(state)
    erx, exl, eatm = structural_fluxes(state, root_potential, vpd, hydraulic, feedback)
    return np.r_[
        (erx-exl)/hydraulic.cx, (exl-eatm)/hydraulic.cl,
        carbon_rhs(state[2:], hydraulic_factor(state[1], hydraulic), carbon),
    ]


def two_way_jacobian(state, root_potential, vpd, hydraulic, carbon, feedback):
    state = _combined(state)
    # Reuse the checked carbon blocks and the source kink guards.
    jac = combined_jacobian(state, root_potential, vpd, hydraulic, carbon)
    gx, gs = conductances(state[1], hydraulic)
    dgx = -hydraulic.gx_max/(2*hydraulic.psi50_x) if gx > 0 else 0.
    dgs = -hydraulic.gs_max/(2*hydraulic.psi50_s) if gs > 0 else 0.
    fl, fe, fw, dfl, dfe, dfw = feedback_factors(state[2], state[3], feedback)
    x, leaf = state[:2]
    gr, gl = 2*gx*fw, 2*gx*fl
    jac[:2, :] = 0.
    jac[0, 0] = -(gr+gl)/hydraulic.cx
    jac[0, 1] = (2*dgx*fw*(root_potential-x)-2*dgx*fl*(x-leaf)+gl)/hydraulic.cx
    jac[0, 2] = -2*gx*dfl*(x-leaf)/hydraulic.cx
    jac[0, 3] = 2*gx*dfw*(root_potential-x)/hydraulic.cx
    jac[1, 0] = gl/hydraulic.cl
    jac[1, 1] = (2*dgx*fl*(x-leaf)-gl-dgs*fe*vpd)/hydraulic.cl
    jac[1, 2] = (2*gx*dfl*(x-leaf)-gs*dfe*vpd)/hydraulic.cl
    return jac


def equilibrium_curve(leaf_potential, root_potential, hydraulic, carbon, feedback):
    """Return (full equilibrium state, D, dD/dpsi_leaf) on positive branch.

    Carbon roots are eliminated exactly. This parameterization traverses
    D-folds and contains every conducting positive-carbon equilibrium
    within the chosen leaf-potential domain.
    """
    _forcing(root_potential, 0.)
    h = hydraulic_factor(leaf_potential, hydraulic)
    stocks = positive_carbon_equilibrium(h, carbon)
    if stocks is None or stocks[0] <= feedback.minimum_leaf:
        raise ValueError("equilibrium is outside the viable positive-carbon domain")
    gx, gs = conductances(leaf_potential, hydraulic)
    if gx <= 0 or gs <= 0 or leaf_potential >= root_potential:
        raise ValueError("curve requires positive conductances and root-to-leaf potential drop")
    fl, fe, fw, dfl, dfe, dfw = feedback_factors(stocks[0], stocks[1], feedback)
    gr, gl = 2*gx*fw, 2*gx*fl
    geff = effective_conductance(gr, gl)
    state = np.r_[(gr*root_potential+gl*leaf_potential)/(gr+gl), leaf_potential, stocks]
    drop = root_potential-leaf_potential
    d = geff*drop/(gs*fe)
    dgx, dgs = -hydraulic.gx_max/(2*hydraulic.psi50_x), -hydraulic.gs_max/(2*hydraulic.psi50_s)
    canopy = -np.expm1(-carbon.canopy_coefficient*stocks[0])
    carbon_slope = carbon.allocation[0]*carbon.assimilation_max*h*carbon.canopy_coefficient*np.exp(
        -carbon.canopy_coefficient*stocks[0],
    )-carbon.losses[0]
    dl = -carbon.allocation[0]*carbon.assimilation_max*(dgs/hydraulic.gs_max)*canopy/carbon_slope
    dw = stocks[1]/stocks[0]*dl
    dgr = 2*(dgx*fw+gx*dfw*dw)
    dgl = 2*(dgx*fl+gx*dfl*dl)
    dgeff = (gl**2*dgr+gr**2*dgl)/(gr+gl)**2
    derivative = d*(dgeff/geff-1/drop-dgs/gs-dfe*dl/fe)
    return state, float(d), float(derivative)


def stationary_points(bounds, root_potential, hydraulic, carbon, feedback, *, samples=300):
    """Bracket all detected simple stationary points; not a global certificate."""
    low, high = bounds
    if not np.all(np.isfinite(bounds)) or low >= high or samples < 3:
        raise ValueError("ordered finite curve bounds and at least three samples required")
    grid = np.linspace(low, high, samples)
    slopes = [equilibrium_curve(z, root_potential, hydraulic, carbon, feedback)[2] for z in grid]
    points = []
    for left, right, fl, fr in zip(grid[:-1], grid[1:], slopes[:-1], slopes[1:]):
        if fl == 0:
            points.append(left)
        elif fl*fr < 0:
            points.append(brentq(
                lambda z: equilibrium_curve(z, root_potential, hydraulic, carbon, feedback)[2],
                left, right, xtol=1e-12,
            ))
    if slopes[-1] == 0:
        points.append(high)
    return tuple(points)


def equilibria_at_vpd(vpd, bounds, stationary, root_potential, hydraulic, carbon, feedback):
    """Bracket positive roots on intervals partitioned by stationary points."""
    _forcing(root_potential, vpd)
    nodes = sorted([bounds[0], *stationary, bounds[1]])
    values = [equilibrium_curve(z, root_potential, hydraulic, carbon, feedback)[1]-vpd for z in nodes]
    roots = []
    for z, value in zip(nodes, values):
        if abs(value) < 1e-12:
            roots.append(z)
    for left, right, fl, fr in zip(nodes[:-1], nodes[1:], values[:-1], values[1:]):
        if abs(fl) >= 1e-12 and abs(fr) >= 1e-12 and fl*fr < 0:
            roots.append(brentq(
                lambda z: equilibrium_curve(z, root_potential, hydraulic, carbon, feedback)[1]-vpd,
                left, right, xtol=1e-12,
            ))
    return tuple(equilibrium_curve(z, root_potential, hydraulic, carbon, feedback)[0] for z in sorted(roots, reverse=True))


def continue_curve(start_leaf, step, count, bounds, root_potential, hydraulic, carbon, feedback, *, d_scale=.02):
    """Reduced pseudo-arclength predictor/corrector in (psi_leaf,D/d_scale).

    Positive carbon is uniquely eliminated. Correct F=D_curve(psi)/scale-d=0
    plus the tangent hyperplane, not a fixed-D forward root solve.
    """
    if not np.all(np.isfinite([step, d_scale])) or step <= 0 or d_scale <= 0 or count < 1:
        raise ValueError("positive finite arclength step/scale and count required")
    if not bounds[0] < start_leaf < bounds[1]:
        raise ValueError("continuation start must lie inside bounds")
    state, d, derivative = equilibrium_curve(start_leaf, root_potential, hydraulic, carbon, feedback)
    point = np.array([start_leaf, d/d_scale])
    tangent = np.array([-1., -derivative/d_scale])
    tangent /= np.linalg.norm(tangent)
    records = [(state, d)]
    for _ in range(count):
        predictor = point+step*tangent
        if not bounds[0] < predictor[0] < bounds[1]:
            return records, "domain_boundary"
        def residual(z):
            if not bounds[0] < z[0] < bounds[1]:
                raise ValueError("continuation corrector left prescribed curve domain")
            dcurve = equilibrium_curve(z[0], root_potential, hydraulic, carbon, feedback)[1]
            return np.array([dcurve/d_scale-z[1], np.dot(z-predictor, tangent)])
        def jacobian(z):
            slope = equilibrium_curve(z[0], root_potential, hydraulic, carbon, feedback)[2]
            return np.array([[slope/d_scale, -1.], tangent])
        corrected = root(residual, predictor, jac=jacobian, tol=1e-10)
        if not corrected.success or np.max(abs(residual(corrected.x))) > 1e-9:
            raise RuntimeError(f"pseudo-arclength corrector failed: {corrected.message}")
        point = corrected.x
        state, d, derivative = equilibrium_curve(point[0], root_potential, hydraulic, carbon, feedback)
        next_tangent = np.array([1., derivative/d_scale])
        next_tangent /= np.linalg.norm(next_tangent)
        if np.dot(next_tangent, tangent) < 0:
            next_tangent *= -1
        tangent = next_tangent
        records.append((state, d))
    return records, "step_limit"


def simulate_two_way(initial, times, forcing: Callable[[float], float], root_potential,
                     hydraulic, carbon, feedback, *, rtol=1e-9):
    """Log-carbon Radau integration terminating at L_min, not at 'death'."""
    initial = _combined(initial)
    times = np.asarray(times, dtype=float)
    if np.any(initial[2:] <= 0) or initial[2] <= feedback.minimum_leaf:
        raise ValueError("simulation requires strictly positive stocks and L>L_min")
    if times.ndim != 1 or len(times) < 2 or not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0):
        raise ValueError("sample times must be finite and strictly increasing")
    if not np.isfinite(rtol) or not 0 < rtol < 1:
        raise ValueError("rtol must lie in (0,1)")
    def rhs(t, z):
        stocks = np.exp(z[2:])
        if not np.all(np.isfinite(stocks)) or np.any(stocks == 0):
            raise FloatingPointError("carbon exceeds representable positive domain")
        values = two_way_rhs(np.r_[z[:2], stocks], root_potential, forcing(t), hydraulic, carbon, feedback)
        return np.r_[values[:2], values[2:]/stocks]
    def boundary(t, z):
        return z[2]-np.log(feedback.minimum_leaf)
    boundary.terminal = True
    boundary.direction = -1
    result = solve_ivp(
        rhs, (times[0], times[-1]), np.r_[initial[:2], np.log(initial[2:])],
        method="Radau", rtol=rtol, atol=rtol*.01, dense_output=True, events=boundary,
    )
    if not result.success:
        raise RuntimeError(result.message)
    status = "vegetation_domain_exit" if result.status == 1 else "completed"
    end = float(result.t[-1])
    sampled_times = np.unique(np.r_[times[times <= end], end])
    z = result.sol(sampled_times)
    states = np.vstack([z[:2], np.exp(z[2:])]).T
    if not np.all(np.isfinite(states)) or np.any(states[:, :2] > 0) or np.any(result.y[:2] > 0):
        raise RuntimeError("accepted/sample trajectory left physical domain")
    if np.any(states[:, 2:] <= 0) or np.any(states[:, 2] < feedback.minimum_leaf*(1-1e-8)):
        raise RuntimeError("vegetation event failed to stop in the declared domain")
    return {"status": status, "time_days": sampled_times, "state": states, "nfev": result.nfev}
