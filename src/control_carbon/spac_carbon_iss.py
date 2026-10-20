"""Capacity-error identities and ISS certificates for unchanged rational E4.

Rates are per day. Exogenous capacity-rate ISS does not certify stability of
the six-state feedback system. See docs/e4_rational_capacity_iss_analysis_ja.md.
"""

from dataclasses import dataclass

import numpy as np
from scipy.linalg import expm, solve_continuous_lyapunov

from .manzoni_hydraulics import conductances
from .spac_carbon import CarbonParameters, carbon_matrix
from .spac_carbon_equilibria import stock_ratios
from .spac_carbon_rational import RationalSPAC
from .spac_carbon_tracking import forcing_derivative


@dataclass(frozen=True)
class CarbonISSCertificate:
    matrix: np.ndarray
    l1_decay_per_day: float
    lyapunov_matrix_days: np.ndarray
    euclidean_prefactor: float
    euclidean_decay_per_day: float
    lyapunov_residual: float


def carbon_iss_certificate(carbon: CarbonParameters) -> CarbonISSCertificate:
    matrix = carbon_matrix(carbon)
    decay = float(-matrix.sum(axis=0).max())
    if decay <= 0:
        raise ValueError("This l1 certificate requires strictly negative column sums")
    metric = solve_continuous_lyapunov(matrix.T, -np.eye(4))
    metric = (metric+metric.T)/2
    eigenvalues = np.linalg.eigvalsh(metric)
    if not np.all(np.isfinite(metric)) or eigenvalues[0] <= 0:
        raise RuntimeError("Carbon Lyapunov solve did not produce a positive metric")
    residual = float(np.max(np.abs(matrix.T@metric+metric@matrix+np.eye(4))))
    if residual > 1e-10:
        raise RuntimeError("Carbon Lyapunov residual exceeded tolerance")
    return CarbonISSCertificate(
        matrix, decay, metric, float(np.sqrt(eigenvalues[-1]/eigenvalues[0])),
        float(1/(2*eigenvalues[-1])), residual)


def componentwise_iss_bound(carbon: CarbonParameters, elapsed_days: float,
                            initial_error, capacity_rate_bound):
    """Bound |e| using the componentwise essential supremum of |Kdot|."""
    initial = np.asarray(initial_error, dtype=float)
    rate = np.asarray(capacity_rate_bound, dtype=float)
    if (not np.isfinite(elapsed_days) or elapsed_days < 0 or initial.shape != (4,)
            or rate.shape != (4,) or not np.all(np.isfinite(initial))
            or not np.all(np.isfinite(rate)) or np.any(rate < 0)):
        raise ValueError("ISS bound needs nonnegative time, four finite errors and rate bounds")
    matrix = carbon_matrix(carbon)
    transition = expm(matrix*elapsed_days)
    return transition@np.abs(initial)+np.linalg.solve(
        -matrix, (np.eye(4)-transition)@rate)


def assimilation_gradient(model: RationalSPAC, state):
    """Classical six-state gradient; no hydraulic Jacobian is needed at -3 MPa."""
    model.assimilation(state)
    state = np.asarray(state, dtype=float)
    hydraulic, carbon = model.base.hydraulic, model.base.carbon
    threshold = 2*hydraulic.psi50_s
    if state[1] == threshold:
        raise ValueError("Classical assimilation gradient is not assigned at the stomatal kink")
    _, gs = conductances(state[1], hydraulic)
    opening_slope = -1/threshold if state[1] > threshold else 0.
    q_leaf = carbon.canopy_coefficient*state[2]
    gradient = np.zeros(6)
    gradient[1] = carbon.assimilation_max*opening_slope*q_leaf/(1+q_leaf)
    gradient[2] = (carbon.assimilation_max*(gs/hydraulic.gs_max)
                   *carbon.canopy_coefficient/(1+q_leaf)**2)
    return gradient


def capacity_diagnostics(model: RationalSPAC, state, vpd_mpa: float):
    """Exact along-path rates, including a directional value at psi_l=-4."""
    derivative = model.rhs(state, vpd_mpa)
    state = np.asarray(state, dtype=float)
    carbon, hydraulic = model.base.carbon, model.base.hydraulic
    matrix = carbon_matrix(carbon)
    allocation = np.r_[carbon.allocation, 0.]
    capacity_direction = np.linalg.solve(-matrix, allocation)
    threshold = 2*hydraulic.psi50_s
    at_kink = state[1] == threshold
    if at_kink:
        q_leaf = carbon.canopy_coefficient*state[2]
        gradient = np.zeros(6)
        gradient[1] = (carbon.assimilation_max*(-1/threshold)*q_leaf/(1+q_leaf)
                       if derivative[1] > 0 else 0.)
    else:
        gradient = assimilation_gradient(model, state)
    capacity = model.capacity(state)
    error = state[2:]-capacity
    capacity_rate = capacity_direction*float(gradient@derivative)
    closed_matrix = matrix.copy()
    closed_matrix[:, 0] += capacity_direction*carbon.losses[0]*gradient[2]
    water_input = -capacity_direction*gradient[1]*derivative[1]
    exact_error_rate = derivative[2:]-capacity_rate
    return {
        "rhs": derivative, "capacity": capacity, "error": error,
        "assimilation_rate": float(gradient@derivative),
        "capacity_rate": capacity_rate, "error_rate": exact_error_rate,
        "carbon_feedback_matrix": closed_matrix,
        "water_rate_input": water_input,
        "leaf_feedback_ratio": float(capacity_direction[0]*gradient[2]),
        "carbon_identity_residual": float(np.max(np.abs(derivative[2:]-matrix@error))),
        "closed_error_identity_residual": float(np.max(np.abs(
            exact_error_rate-closed_matrix@error-water_input))),
        "capacity_derivative_kind": "directional_at_stomatal_kink" if at_kink else "classical",
    }


def carbon_feedback_certificate(model: RationalSPAC, feedback_upper: float):
    """Uniform carbon-error ISS while v_K,L*dA/dC_L <= feedback_upper < 1.

    Leaf-water-potential rate is the input; invariance of this region and the
    hydraulic feedback gain are additional, unproved whole-system conditions.
    """
    if not np.isfinite(feedback_upper) or not 0 <= feedback_upper < 1:
        raise ValueError("Carbon feedback upper bound must lie in [0, 1)")
    carbon = model.base.carbon
    matrix = carbon_matrix(carbon)
    direction = np.linalg.solve(-matrix, np.r_[carbon.allocation, 0.])
    comparison = matrix.copy()
    comparison[:, 0] += direction*carbon.losses[0]*feedback_upper/direction[0]
    unnormalized = np.linalg.solve(-comparison, np.ones(4))
    if not np.all(np.isfinite(unnormalized)) or np.any(unnormalized <= 0):
        raise RuntimeError("Feedback comparison matrix did not produce positive weights")
    weights = unnormalized/unnormalized[0]
    alpha = float(np.min(1/unnormalized))
    input_slope = carbon.assimilation_max/(-2*model.base.hydraulic.psi50_s)
    gain = float(np.max(direction*input_slope/weights)/alpha)
    return {
        "feedback_upper": feedback_upper, "comparison_matrix": comparison,
        "weights": weights, "decay_per_day": alpha,
        "water_rate_gain": gain,
        "leaf_decay_per_day": carbon.losses[0]*(1-feedback_upper),
    }


def equilibrium_tracking_diagnostics(model: RationalSPAC, state, vpd_mpa: float):
    """Local error coordinates and analytic dy*/dVPD, only at a positive E4 root."""
    if model.label != "E4":
        raise ValueError("This full-equilibrium tangent implements only rational E4")
    state = np.asarray(state, dtype=float)
    derivative = model.rhs(state, vpd_mpa)
    if state[2] <= 0 or np.max(np.abs(derivative)) > 1e-10:
        raise ValueError("A positive full equilibrium is required")
    jacobian = model.jacobian(state, vpd_mpa)
    carbon = model.base.carbon
    direction = np.linalg.solve(-carbon_matrix(carbon), np.r_[carbon.allocation, 0.])
    gradient = assimilation_gradient(model, state)
    transform = np.eye(6)
    transform[2:] -= direction[:, None]*gradient[None, :]
    ratio = direction[0]*gradient[2]
    if abs(1-ratio) <= 1e-10:
        raise ValueError("Capacity error coordinates are singular near feedback ratio one")
    transformed = np.linalg.solve(transform.T, (transform@jacobian).T).T
    tangent = np.linalg.solve(jacobian, -forcing_derivative(state, vpd_mpa, model.base))
    feedback = model.base.feedback
    ratios = stock_ratios(model.base)
    fl, _, fw = model.factors(state)
    dfl = 1/(feedback.leaf_reference+feedback.leaf_floor)
    dfw = ratios[1]/feedback.wood_reference
    dpsi = (-2*model.base.hydraulic.psi50_s)*carbon.losses[0]/(
        carbon.allocation[0]*carbon.assimilation_max)
    dx = ((dfl*state[1]+fl*dpsi)/(fl+fw)
          -fl*state[1]*(dfl+dfw)/(fl+fw)**2)
    curve_tangent = np.r_[dx, dpsi, ratios]/model.curve(state[2])[2]
    return {
        "state": state, "jacobian": jacobian, "capacity_coordinate_transform": transform,
        "transformed_jacobian": transformed, "coordinate_determinant": float(1-ratio),
        "tangent_per_mpa": tangent, "curve_tangent_per_mpa": curve_tangent,
        "full_eigenvalues_per_day": np.linalg.eigvals(jacobian),
        "hydraulic_frozen_carbon_eigenvalues_per_day": np.linalg.eigvals(jacobian[:2, :2]),
        "carbon_feedback_eigenvalues_per_day": np.linalg.eigvals(
            capacity_diagnostics(model, state, vpd_mpa)["carbon_feedback_matrix"]),
        "leaf_feedback_ratio": float(ratio),
    }
