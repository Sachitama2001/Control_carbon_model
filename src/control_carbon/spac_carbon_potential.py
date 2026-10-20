"""Signed carbon potential diagnostics and explicitly scoped lower bounds."""

from dataclasses import dataclass

import numpy as np

from .spac_carbon import carbon_matrix
from .spac_carbon_iss import capacity_diagnostics
from .spac_carbon_rational import RationalSPAC
from .spac_carbon_stability import _Interval
from .spac_carbon_vpd_iss import _lower, _magnitude, _upper, equilibrium_enclosure


@dataclass(frozen=True)
class CarbonPotentialStructure:
    matrix: np.ndarray
    allocation: np.ndarray
    capacity_direction_days: np.ndarray


def potential_structure(model: RationalSPAC):
    matrix = carbon_matrix(model.base.carbon)
    allocation = np.r_[model.base.carbon.allocation, 0.]
    return CarbonPotentialStructure(matrix, allocation, np.linalg.solve(-matrix, allocation))


def potential_diagnostics(model: RationalSPAC, state, vpd_mpa):
    diag = capacity_diagnostics(model, state, vpd_mpa)
    structure = potential_structure(model)
    potential = -diag["error"]
    potential_rate = -diag["error_rate"]
    return {
        "assimilation": float(model.assimilation(state)),
        "assimilation_rate": diag["assimilation_rate"], "carbon": np.asarray(state)[2:].copy(),
        "capacity": diag["capacity"], "potential": potential,
        "capacity_rate": diag["capacity_rate"], "potential_rate": potential_rate,
        "rhs": diag["rhs"],
        "potential_identity_residual": float(np.max(np.abs(
            potential_rate-structure.matrix@potential-diag["capacity_rate"]))),
        "total_carbon_rate": float(sum(diag["rhs"][2:])),
    }


def endpoint_peak_lower(matrix, duration_days, capacity_change, initial_potential=None):
    """L1 peak lower bound; endpoint inputs must be bounded/known independently.

    Passing ordinary ODE endpoints makes the result endpoint-informed numerical
    evidence, not a new VPD-only interval guarantee.
    """
    matrix, delta = np.asarray(matrix, dtype=float), np.asarray(capacity_change, dtype=float)
    if (matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or delta.shape != (len(matrix),)
            or np.any(~np.isfinite(matrix)) or np.any(~np.isfinite(delta))
            or not np.isfinite(duration_days) or duration_days < 0):
        raise ValueError("Peak bound needs compatible finite data and nonnegative duration")
    if duration_days == 0 and np.any(delta != 0):
        raise ValueError("Zero duration requires zero capacity change")
    norm = float(np.linalg.norm(matrix, 1))
    change = float(np.linalg.norm(delta, 1))
    if initial_potential is None:
        return change/(2+norm*duration_days)
    initial = np.asarray(initial_potential, dtype=float)
    if initial.shape != delta.shape or np.any(~np.isfinite(initial)):
        raise ValueError("Initial potential must have one finite entry per pool")
    initial_norm = float(np.linalg.norm(initial, 1))
    return max(initial_norm, max(0., change-initial_norm)/(1+norm*duration_days))


def conducting_rhs_enclosure(model: RationalSPAC, state, vpd_mpa):
    if (model.label != "E4" or len(state) != 6 or state[1].lo <= -3
            or any(x.hi > 0 for x in state[:2]) or any(x.lo < 0 for x in state[2:])):
        raise ValueError("Potential certificate requires a physical conducting E4 box")
    vpd = vpd_mpa if isinstance(vpd_mpa, _Interval) else _Interval(vpd_mpa)
    if vpd.lo < 0:
        raise ValueError("VPD must be nonnegative")
    h, c, f = model.base.hydraulic, model.base.carbon, model.base.feedback
    x, psi, leaf, wood, root, soil = state
    gx = _Interval(h.gx_max)*(1-psi/_Interval(2*h.psi50_x))
    gs = _Interval(h.gs_max)*(1-psi/_Interval(2*h.psi50_s))
    fl = (leaf+_Interval(f.leaf_floor))/(
        _Interval(f.leaf_reference)+_Interval(f.leaf_floor))
    fe, fw = leaf/_Interval(f.leaf_reference), wood/_Interval(f.wood_reference)
    erx, exl = 2*gx*fw*(-x), 2*gx*fl*(x-psi)
    evaporation = gs*fe*vpd/_Interval(model.base.atmospheric_pressure_mpa)
    q_leaf = _Interval(c.canopy_coefficient)*leaf
    assimilation = _Interval(c.assimilation_max)*gs/_Interval(h.gs_max)*q_leaf/(1+q_leaf)
    carbon = [-_Interval(loss)*stock+_Interval(allocation)*assimilation
              for loss, stock, allocation in zip(c.losses, [leaf, wood, root], c.allocation)]
    carbon.append(sum((_Interval(eta)*_Interval(loss)*stock
                       for eta, loss, stock in zip(c.litter_fractions, c.losses, [leaf, wood, root])),
                      _Interval(0))-_Interval(c.losses[3])*soil)
    return [(erx-exl)/_Interval(h.cx), (exl-evaporation)/_Interval(h.cl), *carbon]


def assimilation_rate_enclosure(model, state, rhs):
    c, h = model.base.carbon, model.base.hydraulic
    q_leaf = _Interval(c.canopy_coefficient)*state[2]
    a_leaf = (_Interval(c.assimilation_max)*(1-state[1]/_Interval(2*h.psi50_s))
              *_Interval(c.canopy_coefficient)/((1+q_leaf)*(1+q_leaf)))
    a_psi = -_Interval(c.assimilation_max)/_Interval(2*h.psi50_s)*q_leaf/(1+q_leaf)
    return a_leaf*rhs[2]+a_psi*rhs[1]


def certify_drying_step_potential(model: RationalSPAC, initial_vpd_mpa: float,
                                  step_vpd_mpa: float, duration_days: float):
    """A priori short-time enclosure and strictly negative potential certificate.

    No six-state time integration is used. The initial state is the exact
    positive equilibrium, enclosed by verified endpoint signs.
    """
    if (not np.isfinite(step_vpd_mpa) or step_vpd_mpa <= initial_vpd_mpa
            or not np.isfinite(duration_days) or duration_days <= 0):
        raise ValueError("Certificate requires a drying step and a positive finite horizon")
    equilibrium, initial, evidence = equilibrium_enclosure(model, initial_vpd_mpa)
    initial_rhs = conducting_rhs_enclosure(model, initial, step_vpd_mpa)
    duration = _Interval(duration_days)
    radius = [
        _Interval(max(2*(duration*_Interval(_magnitude(rate))).hi,
                      (value.hi-value.lo)*16, _Interval("1e-14").hi)).hi
        for value, rate in zip(initial, initial_rhs)]
    for iteration in range(20):
        box = [value+_Interval(r.copy_negate(), r) for value, r in zip(initial, radius)]
        rhs = conducting_rhs_enclosure(model, box, step_vpd_mpa)
        reach = [value+_Interval(0, duration.hi)*rate for value, rate in zip(initial, rhs)]
        if all(candidate.lo < reached.lo and reached.hi < candidate.hi
               for candidate, reached in zip(box, reach)):
            break
        radius = [
            max(r, (_Interval("1.25")*_Interval(max(
                (reached.lo-value.lo).copy_abs(), (reached.hi-value.hi).copy_abs()))).hi)
            for r, value, reached in zip(radius, initial, reach)]
    else:
        raise ValueError("A priori Picard inclusion did not close at this horizon")
    rate = assimilation_rate_enclosure(model, box, rhs)
    if rate.hi >= 0:
        raise ValueError("A priori box did not certify strictly decreasing assimilation")
    c = model.base.carbon
    losses = [_Interval(value) for value in c.losses]
    direction = [_Interval(b)/k for b, k in zip(c.allocation, losses)]
    direction.append(sum((_Interval(eta)*_Interval(b)
                          for eta, b in zip(c.litter_fractions, c.allocation)),
                         _Interval(0))/losses[3])
    minimum, maximum = _Interval(rate.hi.copy_negate()), _Interval(rate.lo.copy_negate())
    phi = [(1-(-k*duration).exp())/k for k in losses]
    lower = [minimum*v*p for v, p in zip(direction, phi)]
    upper = [maximum*v*p for v, p in zip(direction[:3], phi[:3])]
    transfer = sum((_Interval(eta)*k*v
                    for eta, k, v in zip(c.litter_fractions, losses, direction)), _Interval(0))
    upper.append(maximum*(direction[3]*phi[3]+transfer*duration*duration/2))
    change_lower = minimum*duration*sum(direction, _Interval(0))
    matrix_norm_upper = _Interval(max(
        [((_Interval(1)+_Interval(eta))*k).hi
         for eta, k in zip(c.litter_fractions, losses)] + [losses[3].hi]))
    endpoint_lower = change_lower/(1+matrix_norm_upper*duration)
    cumulative = sum((minimum*v*(duration-p)/k
                      for v, p, k in zip(direction, phi, losses)), _Interval(0))
    return {
        "initial_vpd_mpa": initial_vpd_mpa, "step_vpd_mpa": step_vpd_mpa,
        "duration_days": duration_days, "initial_equilibrium_approximation": equilibrium,
        "initial_equilibrium_enclosure": evidence,
        "state_box": [x.strings() for x in box], "reach_box": [x.strings() for x in reach],
        "rhs_intervals": [x.strings() for x in rhs], "picard_iterations": iteration+1,
        "assimilation_rate_interval": rate.strings(),
        "absolute_potential_component_lower": [_lower(x) for x in lower],
        "absolute_potential_component_upper": [_upper(x) for x in upper],
        "capacity_drop_l1_lower": _lower(change_lower),
        "peak_potential_l1_endpoint_lower": _lower(endpoint_lower),
        "final_potential_l1_lower": _lower(sum(lower, _Interval(0))),
        "cumulative_potential_l1_lower": _lower(cumulative),
        "scope": "B: initial exact equilibrium, constant post-step VPD, a priori short horizon only",
    }
