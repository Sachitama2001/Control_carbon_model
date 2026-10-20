"""Local ISS of the unchanged six-state E4 system with VPD as its only input.

Decimal intervals enclose real arithmetic with binary64 coefficients taken as
exact constants. A certified ellipsoid is local, not a global basin estimate.
"""

from dataclasses import dataclass
from decimal import localcontext
import math

import numpy as np
from scipy.linalg import solve_continuous_lyapunov

from .spac_carbon_rational import RationalSPAC
from .spac_carbon_stability import _Interval, PRECISION


def _sqrt(value):
    if value.lo < 0:
        raise ValueError("Interval square root requires nonnegative endpoints")
    with localcontext() as context:
        context.prec = PRECISION
        low = value.lo.sqrt().next_minus() if value.lo else value.lo
        return _Interval(low, value.hi.sqrt().next_plus())


def _lower(value):
    return math.nextafter(float(value.lo), -math.inf)


def _upper(value):
    return math.nextafter(float(value.hi), math.inf)


def _magnitude(value):
    return max(value.lo.copy_abs(), value.hi.copy_abs())


def _point_matrix(matrix):
    return [[_Interval(float(value)) for value in row] for row in matrix]


def _transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def _matmul(left, right):
    return [[sum((a*b for a, b in zip(row, col)), _Interval(0))
             for col in zip(*right)] for row in left]


def _gershgorin(matrix):
    lower, upper = [], []
    for i, row in enumerate(matrix):
        radius = sum((_Interval(_magnitude(x)) for j, x in enumerate(row) if i != j),
                     _Interval(0))
        lower.append((row[i]-radius).lo)
        upper.append((row[i]+radius).hi)
    return _Interval(min(lower), max(upper))


def _state_on_leaf_interval(model, leaf):
    c, f, h = model.base.carbon, model.base.feedback, model.base.hydraulic
    losses = [_Interval(x) for x in c.losses]
    allocation = [_Interval(x) for x in c.allocation]
    opening = losses[0]*(1+_Interval(c.canopy_coefficient)*leaf)/(
        allocation[0]*_Interval(c.assimilation_max)*_Interval(c.canopy_coefficient))
    psi = _Interval(-2*h.psi50_s)*(opening-1)
    wood = allocation[1]*losses[0]/(allocation[0]*losses[1])*leaf
    root = allocation[2]*losses[0]/(allocation[0]*losses[2])*leaf
    soil = sum((_Interval(eta)*loss*stock
                for eta, loss, stock in zip(c.litter_fractions, losses, [leaf, wood, root])),
               _Interval(0))/losses[3]
    fl = (leaf+_Interval(f.leaf_floor))/(
        _Interval(f.leaf_reference)+_Interval(f.leaf_floor))
    fw, fe = wood/_Interval(f.wood_reference), leaf/_Interval(f.leaf_reference)
    xylem = fl*psi/(fl+fw)
    gx = _Interval(h.gx_max)*(1-psi/_Interval(2*h.psi50_x))
    gs = _Interval(h.gs_max)*opening
    vpd = (2*gx*fl*fw*(-psi)/((fl+fw)*gs*fe)
           *_Interval(model.base.atmospheric_pressure_mpa))
    return [xylem, psi, leaf, wood, root, soil], vpd


def equilibrium_enclosure(model: RationalSPAC, vpd_mpa: float):
    """Enclose the unique positive E4 root; endpoints are verified, not sampled."""
    if model.label != "E4" or not np.isfinite(vpd_mpa) or vpd_mpa <= 0:
        raise ValueError("A positive VPD and rational E4 are required")
    shape = model.shape_bounds()
    roots = model.equilibria(vpd_mpa)
    if len(roots) != 1:
        raise ValueError("A unique positive E4 equilibrium is required")
    leaf = float(roots[0][2])
    width = max(leaf*1e-12, 1e-15)
    bracket = _Interval(leaf-width, leaf+width)
    if bracket.lo <= 0:
        raise ValueError("Positive root is too close to zero for this enclosing bracket")
    _, low_vpd = _state_on_leaf_interval(model, _Interval(bracket.lo))
    _, high_vpd = _state_on_leaf_interval(model, _Interval(bracket.hi))
    if not low_vpd.lo > _Interval(vpd_mpa).hi or not high_vpd.hi < _Interval(vpd_mpa).lo:
        raise RuntimeError("Interval endpoint signs did not enclose the E4 equilibrium")
    state, _ = _state_on_leaf_interval(model, bracket)
    return roots[0], state, {
        "leaf_bracket": bracket.strings(), "low_leaf_vpd": low_vpd.strings(),
        "high_leaf_vpd": high_vpd.strings(), "state_intervals": [x.strings() for x in state],
        "monotonicity": shape,
    }


def conducting_jacobian_enclosure(model: RationalSPAC, state, vpd_mpa: float):
    """Full six-state Jacobian on a strictly conducting interval box."""
    h, c, f = model.base.hydraulic, model.base.carbon, model.base.feedback
    if (model.label != "E4" or not np.isfinite(vpd_mpa) or vpd_mpa < 0
            or len(state) != 6 or state[1].lo <= 2*h.psi50_x
            or any(x.hi > 0 for x in state[:2]) or any(x.lo < 0 for x in state[2:])):
        raise ValueError("Jacobian enclosure needs a nonnegative-carbon conducting physical box")
    x, psi, leaf, wood, _, _ = state
    gx = _Interval(h.gx_max)*(1-psi/_Interval(2*h.psi50_x))
    gs = _Interval(h.gs_max)*(1-psi/_Interval(2*h.psi50_s))
    dgx = -_Interval(h.gx_max)/_Interval(2*h.psi50_x)
    dgs = -_Interval(h.gs_max)/_Interval(2*h.psi50_s)
    fl = (leaf+_Interval(f.leaf_floor))/(
        _Interval(f.leaf_reference)+_Interval(f.leaf_floor))
    fe, fw = leaf/_Interval(f.leaf_reference), wood/_Interval(f.wood_reference)
    dfl = 1/(_Interval(f.leaf_reference)+_Interval(f.leaf_floor))
    dfe, dfw = 1/_Interval(f.leaf_reference), 1/_Interval(f.wood_reference)
    d = _Interval(vpd_mpa)/_Interval(model.base.atmospheric_pressure_mpa)
    cx, cl = _Interval(h.cx), _Interval(h.cl)
    jac = [[_Interval(0) for _ in range(6)] for _ in range(6)]
    jac[0][0] = -2*gx*(fw+fl)/cx
    jac[0][1] = (2*dgx*(fw*(-x)-fl*(x-psi))+2*gx*fl)/cx
    jac[0][2] = -2*gx*dfl*(x-psi)/cx
    jac[0][3] = 2*gx*dfw*(-x)/cx
    jac[1][0] = 2*gx*fl/cl
    jac[1][1] = (2*dgx*fl*(x-psi)-2*gx*fl-dgs*fe*d)/cl
    jac[1][2] = (2*gx*dfl*(x-psi)-gs*dfe*d)/cl
    losses = [_Interval(value) for value in c.losses]
    for i, loss in enumerate(losses):
        jac[i+2][i+2] = -loss
    for i, eta in enumerate(c.litter_fractions):
        jac[5][i+2] += _Interval(eta)*losses[i]
    q_leaf = _Interval(c.canopy_coefficient)*leaf
    for i, allocation in enumerate([*c.allocation, 0.]):
        jac[i+2][1] = (_Interval(allocation)*_Interval(c.assimilation_max)
                        *dgs/_Interval(h.gs_max)*q_leaf/(1+q_leaf))
        jac[i+2][2] += (_Interval(allocation)*_Interval(c.assimilation_max)
                        *gs/_Interval(h.gs_max)*_Interval(c.canopy_coefficient)
                        /((1+q_leaf)*(1+q_leaf)))
    return jac


@dataclass(frozen=True)
class LocalVPDISS:
    vpd_mpa: float
    equilibrium: np.ndarray
    scales: np.ndarray
    metric: np.ndarray
    state_box_radius_scaled: float
    metric_lower: float
    metric_upper: float
    dissipation_lower_per_day: float
    input_cross_gain_per_mpa: float
    energy_limit: float
    max_vpd_perturbation_mpa: float
    energy_decay_per_day: float
    energy_input_gain: float
    capacity_error_lipschitz: float
    evidence: dict


def certify_local_vpd_iss(model: RationalSPAC, vpd_mpa: float, scales,
                          state_box_radius_scaled: float) -> LocalVPDISS:
    """Verify a nonlinear local VPD-input ISS region using interval arithmetic."""
    scale = np.asarray(scales, dtype=float)
    radius = state_box_radius_scaled
    if (scale.shape != (6,) or not np.all(np.isfinite(scale)) or np.any(scale <= 0)
            or not np.isfinite(radius) or radius <= 0):
        raise ValueError("Local ISS requires six positive scales and a positive finite radius")
    equilibrium, root_box, root_evidence = equilibrium_enclosure(model, vpd_mpa)
    jacobian = model.jacobian(equilibrium, vpd_mpa)*scale[None, :]/scale[:, None]
    if np.linalg.eigvals(jacobian).real.max() >= 0:
        raise ValueError("Numerical frozen Jacobian is not Hurwitz; no local certificate attempted")
    metric = solve_continuous_lyapunov(jacobian.T, -np.eye(6))
    metric = (metric+metric.T)/2
    if not np.all(np.isfinite(metric)):
        raise RuntimeError("Nonfinite candidate Lyapunov metric")
    normalizer = np.linalg.inv(np.linalg.cholesky(metric))
    interval_metric, interval_normalizer = _point_matrix(metric), _point_matrix(normalizer)
    normalized = _matmul(_matmul(interval_normalizer, interval_metric),
                         _transpose(interval_normalizer))
    metric_check = _gershgorin(normalized)
    if metric_check.lo <= 0:
        raise RuntimeError("Interval congruence did not certify a positive metric")
    norm_squared = sum((_Interval(float(t))*_Interval(float(t))
                        for t in normalizer.ravel()), _Interval(0))
    p_lower = _lower(_Interval(metric_check.lo)/norm_squared)
    p_upper = _upper(_Interval(max(
        sum((_Interval(abs(float(t))) for t in row), _Interval(0)).hi for row in metric)))
    level = _lower(_Interval(p_lower)*_Interval(radius)*_Interval(radius))
    axis_bounds = [
        _sqrt(_Interval(level)*sum(
            (row[i]*row[i] for row in interval_normalizer), _Interval(0))
              /_Interval(metric_check.lo)) for i in range(6)]
    box = []
    for value, s, axis in zip(root_box, scale, axis_bounds):
        extent = _Interval(float(s))*_Interval(axis.hi)
        box.append(value+_Interval(extent.hi.copy_negate(), extent.hi))
    if (box[1].lo <= 2*model.base.hydraulic.psi50_x
            or any(x.hi >= 0 for x in box[:2]) or any(x.lo <= 0 for x in box[2:])):
        raise ValueError("Requested region reaches a physical or conductance boundary")

    physical_jacobian = conducting_jacobian_enclosure(model, box, vpd_mpa)
    scaled = [[value*_Interval(float(scale[j]))/_Interval(float(scale[i]))
               for j, value in enumerate(row)] for i, row in enumerate(physical_jacobian)]
    pa = _matmul(interval_metric, scaled)
    negative_symmetric = [[-pa[i][j]-pa[j][i] for j in range(6)] for i in range(6)]
    dissipation_check = _gershgorin(negative_symmetric)
    if dissipation_check.lo <= 0:
        raise ValueError("Interval nonlinear dissipation is unresolved on the requested box")
    q_lower = _lower(_Interval(dissipation_check.lo))
    h, c, f = model.base.hydraulic, model.base.carbon, model.base.feedback
    opening = 1-box[1]/_Interval(2*h.psi50_s)
    forcing_coefficient = (-_Interval(h.gs_max)*opening*box[2]/_Interval(f.leaf_reference)
                           /(_Interval(h.cl)*_Interval(model.base.atmospheric_pressure_mpa)
                             *_Interval(float(scale[1]))))
    p_column_norm = _sqrt(sum((row[1]*row[1] for row in interval_metric), _Interval(0)))
    b_upper = _upper(_Interval(_magnitude(forcing_coefficient))*p_column_norm)
    delta = _lower(_Interval(q_lower)/(4*_Interval(b_upper))
                    *_sqrt(_Interval(level)/_Interval(p_upper)))
    delta = min(delta, vpd_mpa/2)
    decay = _lower(_Interval(q_lower)/(2*_Interval(p_upper)))
    input_gain = _upper(2*_Interval(b_upper)*_Interval(b_upper)/_Interval(q_lower))
    if min(p_lower, p_upper, q_lower, b_upper, level, delta, decay, input_gain) <= 0:
        raise RuntimeError("Local ISS constants lost their positive margins")

    zero = _Interval(0)
    # Compute -M^-1 b directly from loss/transfer coefficients, not a rounded inverse.
    direction = [_Interval(c.allocation[i])/_Interval(c.losses[i]) for i in range(3)]
    direction.append(sum((_Interval(eta)*_Interval(allocation)
                          for eta, allocation in zip(c.litter_fractions, c.allocation)),
                         zero)/_Interval(c.losses[3]))
    q_leaf = _Interval(c.canopy_coefficient)*box[2]
    a_leaf = (_Interval(c.assimilation_max)*opening*_Interval(c.canopy_coefficient)
              /((1+q_leaf)*(1+q_leaf)))
    a_psi = -_Interval(c.assimilation_max)/_Interval(2*h.psi50_s)*q_leaf/(1+q_leaf)
    output_gradient = []
    for i in range(4):
        row = []
        for j in range(6):
            value = _Interval(1 if j == i+2 else 0)
            if j in [1, 2]:
                value -= direction[i]*(a_psi if j == 1 else a_leaf)*_Interval(
                    float(scale[j]))/_Interval(float(scale[i+2]))
            row.append(value)
        output_gradient.append(row)
    lipschitz = _upper(_sqrt(sum(
        (_Interval(_magnitude(x))*_Interval(_magnitude(x))
         for row in output_gradient for x in row), zero)))
    evidence = {
        "arithmetic": "40-digit outward Decimal; binary64 coefficients are exact real constants",
        "root": root_evidence, "state_box": [x.strings() for x in box],
        "ellipsoid_axis_bounds_scaled": [x.strings() for x in axis_bounds],
        "normalized_metric_eigenvalue_enclosure": metric_check.strings(),
        "dissipation_gershgorin_enclosure": dissipation_check.strings(),
        "forcing_coefficient_interval_per_mpa": forcing_coefficient.strings(),
        "output_gradient_intervals": [[x.strings() for x in row] for row in output_gradient],
        "interpretation": "Certified ellipsoid inside a conducting box, NOT the entire box or a global basin",
    }
    return LocalVPDISS(
        vpd_mpa, equilibrium, scale.copy(), metric, radius, p_lower, p_upper, q_lower,
        b_upper, level, delta, decay, input_gain, lipschitz, evidence)


def initial_energy_upper(certificate: LocalVPDISS, initial_state):
    initial = np.asarray(initial_state, dtype=float)
    if initial.shape != (6,) or not np.all(np.isfinite(initial)):
        raise ValueError("Initial state must have six finite entries")
    roots = certificate.evidence["root"]["state_intervals"]
    error = [(_Interval(float(y))-_Interval(*bounds))/_Interval(float(scale))
             for y, bounds, scale in zip(initial, roots, certificate.scales)]
    quadratic = sum((error[i]*_Interval(float(certificate.metric[i, j]))*error[j]
                     for i in range(6) for j in range(6)), _Interval(0))
    return _upper(_Interval(quadratic.hi))


def local_iss_bounds(certificate: LocalVPDISS, elapsed_days, initial_energy: float,
                     input_sup_mpa: float):
    times = np.asarray(elapsed_days, dtype=float)
    if (np.any(~np.isfinite(times)) or np.any(times < 0)
            or not np.isfinite(initial_energy) or not 0 <= initial_energy <= certificate.energy_limit
            or not np.isfinite(input_sup_mpa)
            or not 0 <= input_sup_mpa <= certificate.max_vpd_perturbation_mpa):
        raise ValueError("Times, initial energy or VPD bound are outside the local certificate")
    decay = np.exp(-certificate.energy_decay_per_day*times)
    energy = (decay*initial_energy+certificate.energy_input_gain
              /certificate.energy_decay_per_day*(-np.expm1(
                  -certificate.energy_decay_per_day*times))*input_sup_mpa**2)
    state = np.sqrt(energy/certificate.metric_lower)
    return {"energy": energy, "state_norm_scaled": state,
            "capacity_error_norm_scaled": certificate.capacity_error_lipschitz*state}


def quartic_instability_obstruction(model: RationalSPAC, vpd_mpa: float):
    """A strictly negative third Hurwitz determinant rules out local ISS."""
    from itertools import combinations, permutations
    _, box, root = equilibrium_enclosure(model, vpd_mpa)
    jac = conducting_jacobian_enclosure(model, box, vpd_mpa)
    indices = [0, 1, 2, 3]
    active = [[jac[i][j] for j in indices] for i in indices]
    def determinant(matrix):
        result = _Interval(0)
        for permutation in permutations(range(len(matrix))):
            inversions = sum(permutation[i] > permutation[j] for i in range(len(matrix))
                             for j in range(i+1, len(matrix)))
            term = _Interval(-1 if inversions % 2 else 1)
            for i, j in enumerate(permutation):
                term *= matrix[i][j]
            result += term
        return result
    coefficients = []
    for degree in range(1, 5):
        total = _Interval(0)
        for subset in combinations(range(4), degree):
            total += determinant([[active[i][j] for j in subset] for i in subset])
        coefficients.append(total*(-1 if degree % 2 else 1))
    a1, a2, a3, a4 = coefficients
    delta3 = a1*a2*a3-a3*a3-a1*a1*a4
    if delta3.hi >= 0:
        raise ValueError("This Hurwitz-determinant obstruction did not prove instability")
    return {
        "root": root, "characteristic_coefficients": [a.strings() for a in coefficients],
        "third_hurwitz_determinant": delta3.strings(),
        "conclusion": "A positive-real-part eigenvalue exists; no local ISS about this equilibrium",
    }
