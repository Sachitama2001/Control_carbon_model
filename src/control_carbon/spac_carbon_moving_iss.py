"""Verified moving-equilibrium tubes for rational E4; VPD is the only forcing."""

from dataclasses import dataclass

import numpy as np

from .spac_carbon_rational import RationalSPAC
from .spac_carbon_stability import _Interval
from .spac_carbon_tracking import SmoothVPDRamp
from .spac_carbon_vpd_iss import (
    _capacity_error_lipschitz, _gershgorin, _lower, _matmul, _point_matrix,
    _quadratic_error_upper, _sqrt, _state_on_leaf_interval, _upper,
    certify_local_vpd_iss, conducting_jacobian_enclosure, equilibrium_enclosure,
)


def equilibrium_tangent_enclosure(model: RationalSPAC, leaf: _Interval):
    """Analytic dy*/dv enclosure, independently of a rounded matrix inverse."""
    state, vpd = _state_on_leaf_interval(model, leaf)
    h, c, f = model.base.hydraulic, model.base.carbon, model.base.feedback
    losses = [_Interval(x) for x in c.losses]
    allocation = [_Interval(x) for x in c.allocation]
    slope = _Interval(-2*h.psi50_s)*losses[0]/(allocation[0]*_Interval(c.assimilation_max))
    fl = (leaf+_Interval(f.leaf_floor))/(
        _Interval(f.leaf_reference)+_Interval(f.leaf_floor))
    dfl = 1/(_Interval(f.leaf_reference)+_Interval(f.leaf_floor))
    ratio_w = allocation[1]*losses[0]/(allocation[0]*losses[1])
    fw, dfw = state[3]/_Interval(f.wood_reference), ratio_w/_Interval(f.wood_reference)
    gx = _Interval(h.gx_max)*(1-state[1]/_Interval(2*h.psi50_x))
    gs = _Interval(h.gs_max)*(1-state[1]/_Interval(2*h.psi50_s))
    dv_dleaf = vpd*(
        -_Interval(h.gx_max)*slope/(_Interval(2*h.psi50_x)*gx)+dfl/fl
        +slope/state[1]-(dfl+dfw)/(fl+fw)
        +_Interval(h.gs_max)*slope/(_Interval(2*h.psi50_s)*gs))
    if dv_dleaf.hi >= 0:
        raise ValueError("Equilibrium tangent enclosure did not resolve a strictly negative slope")
    dx_dleaf = ((dfl*state[1]+fl*slope)/(fl+fw)
                -fl*state[1]*(dfl+dfw)/((fl+fw)*(fl+fw)))
    ratio_r = allocation[2]*losses[0]/(allocation[0]*losses[2])
    ratio_s = sum((_Interval(eta)*loss*ratio
                   for eta, loss, ratio in zip(c.litter_fractions, losses,
                                               [_Interval(1), ratio_w, ratio_r])),
                  _Interval(0))/losses[3]
    tangent = [x/dv_dleaf for x in [dx_dleaf, slope, _Interval(1), ratio_w, ratio_r, ratio_s]]
    return tangent, dv_dleaf


@dataclass(frozen=True)
class MovingEquilibriumISS:
    vpd_band_mpa: tuple[float, float]
    scales: np.ndarray
    metric: np.ndarray
    metric_lower: float
    metric_upper: float
    energy_limit: float
    dissipation_lower_per_day: float
    decay_per_day: float
    equilibrium_speed_gain_per_mpa: float
    max_vpd_rate_mpa_per_day: float
    max_total_vpd_variation_mpa: float
    capacity_error_lipschitz: float
    evidence: dict


def certify_moving_equilibrium(model: RationalSPAC, vpd_band_mpa, scales,
                               state_box_radius_scaled=.0015) -> MovingEquilibriumISS:
    band = np.asarray(vpd_band_mpa, dtype=float)
    if (band.shape != (2,) or np.any(~np.isfinite(band))
            or not 0 < band[0] < band[1]):
        raise ValueError("Moving equilibrium needs an ordered positive VPD band")
    midpoint = float((band[0]+band[1])/2)
    fixed = certify_local_vpd_iss(model, midpoint, scales, state_box_radius_scaled)
    _, _, left = equilibrium_enclosure(model, float(band[0]))
    _, _, right = equilibrium_enclosure(model, float(band[1]))
    leaf = _Interval(right["leaf_bracket"][0], left["leaf_bracket"][1])
    branch_box, _ = _state_on_leaf_interval(model, leaf)
    box = []
    for value, scale, pair in zip(branch_box, fixed.scales,
                                  fixed.evidence["ellipsoid_axis_bounds_scaled"]):
        extent = _Interval(float(scale))*_Interval(pair[1])
        box.append(value+_Interval(extent.hi.copy_negate(), extent.hi))
    if any(x.lo <= 0 for x in box[2:]):
        raise ValueError("Moving tube reaches the zero-carbon boundary")
    jac = conducting_jacobian_enclosure(model, box, _Interval(float(band[0]), float(band[1])))
    scaled = [[x*_Interval(float(fixed.scales[j]))/_Interval(float(fixed.scales[i]))
               for j, x in enumerate(row)] for i, row in enumerate(jac)]
    pa = _matmul(_point_matrix(fixed.metric), scaled)
    check = _gershgorin([[-pa[i][j]-pa[j][i] for j in range(6)] for i in range(6)])
    if check.lo <= 0:
        raise ValueError("No common nonlinear dissipation certificate on the VPD band")
    q = _lower(_Interval(check.lo))
    decay = _lower(_Interval(q)/(2*_Interval(fixed.metric_upper)))
    tangent, slope = equilibrium_tangent_enclosure(model, leaf)
    scaled_tangent = [x/_Interval(float(s)) for x, s in zip(tangent, fixed.scales)]
    square_gain = sum(
        (scaled_tangent[i]*_Interval(float(fixed.metric[i, j]))*scaled_tangent[j]
         for i in range(6) for j in range(6)), _Interval(0))
    if square_gain.hi <= 0:
        raise RuntimeError("Equilibrium-speed gain enclosure is not positive")
    gain = _upper(_sqrt(_Interval(square_gain.hi)))
    speed = _lower(_Interval(decay)*_sqrt(_Interval(fixed.energy_limit))/(2*_Interval(gain)))
    variation = _lower(_sqrt(_Interval(fixed.energy_limit))/(2*_Interval(gain)))
    lipschitz, _ = _capacity_error_lipschitz(model, box, fixed.scales)
    return MovingEquilibriumISS(
        (float(band[0]), float(band[1])), fixed.scales, fixed.metric, fixed.metric_lower,
        fixed.metric_upper, fixed.energy_limit, q, decay, gain, speed, variation, lipschitz, {
            "left_equilibrium": left, "right_equilibrium": right,
            "metric_certificate": fixed.evidence,
            "state_box": [x.strings() for x in box],
            "leaf_band": leaf.strings(), "dv_dleaf_interval": slope.strings(),
            "equilibrium_tangent_intervals_per_mpa": [x.strings() for x in tangent],
            "nonlinear_dissipation_enclosure": check.strings(),
            "interpretation": "Common constant metric; moving ellipsoid, not the entire enclosing box",
        })


def moving_initial_energy_upper(model, certificate, initial_state, vpd_mpa):
    if not certificate.vpd_band_mpa[0] <= vpd_mpa <= certificate.vpd_band_mpa[1]:
        raise ValueError("Initial VPD lies outside the certified band")
    _, _, evidence = equilibrium_enclosure(model, vpd_mpa)
    return _quadratic_error_upper(certificate.metric, certificate.scales,
                                  evidence["state_intervals"], initial_state)


def ramp_speed_upper(ramp: SmoothVPDRamp):
    delta = _Interval(ramp.end_mpa)-_Interval(ramp.start_mpa)
    magnitude = max(delta.lo.copy_abs(), delta.hi.copy_abs())
    return _upper(_Interval("1.5")*_Interval(magnitude)/_Interval(ramp.duration_days))


def moving_ramp_bounds(certificate: MovingEquilibriumISS, ramp: SmoothVPDRamp,
                       elapsed_days, initial_energy: float):
    times = np.asarray(elapsed_days, dtype=float)
    speed = ramp_speed_upper(ramp)
    change = _Interval(ramp.end_mpa)-_Interval(ramp.start_mpa)
    variation = _upper(_Interval(max(change.lo.copy_abs(), change.hi.copy_abs())))
    if not np.isfinite(initial_energy) or not 0 <= initial_energy <= certificate.energy_limit:
        raise ValueError("Initial energy violates the moving-tube certificate")
    variation_budget = (_sqrt(_Interval(initial_energy))
                        +_Interval(certificate.equilibrium_speed_gain_per_mpa)*_Interval(variation))
    variation_ok = variation_budget.hi <= (_sqrt(_Interval(certificate.energy_limit))/2).lo
    if (np.any(~np.isfinite(times)) or np.any(times < 0)
            or min(ramp.start_mpa, ramp.end_mpa) < certificate.vpd_band_mpa[0]
            or max(ramp.start_mpa, ramp.end_mpa) > certificate.vpd_band_mpa[1]
            or (speed > certificate.max_vpd_rate_mpa_per_day and not variation_ok)):
        raise ValueError("Ramp, initial state or time violates the moving-tube certificate")
    during = np.minimum(times, ramp.duration_days)
    after = np.maximum(times-ramp.duration_days, 0)
    rate_envelope = (certificate.equilibrium_speed_gain_per_mpa*speed/certificate.decay_per_day
                     *(-np.expm1(-certificate.decay_per_day*during)))
    variation_envelope = certificate.equilibrium_speed_gain_per_mpa*variation
    radius = (np.exp(-certificate.decay_per_day*during)*np.sqrt(initial_energy)
              +np.minimum(rate_envelope, variation_envelope))
    radius *= np.exp(-certificate.decay_per_day*after)
    state = radius/np.sqrt(certificate.metric_lower)
    return {"metric_radius": radius, "moving_norm_scaled": state,
            "capacity_error_norm_scaled": certificate.capacity_error_lipschitz*state}
