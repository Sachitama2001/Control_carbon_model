"""Local phase-fixed E4 saddle-cycle geometry, not a global basin boundary."""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import null_space

from .spac_carbon_e4 import core_rhs, core_jacobian


def local_cycle_section(cycle, model, scales, *, rtol=1e-10, steps_per_period=160,
                        equilibrium=None):
    """Remove the neutral phase and compute the three-dimensional return derivative.

    The zero set of the left unstable covector is the stable-manifold tangent
    plane on this section. It is not the nonlinear/global stable manifold.
    """
    if not np.isfinite(rtol) or not 0 < rtol < 1 or steps_per_period < 20:
        raise ValueError("section needs positive tolerance and at least 20 steps per period")
    scale = np.asarray(scales.state[:4])
    base = np.asarray(cycle["state"][0, :4])
    period, vpd = cycle["period_days"], cycle["vpd_mpa"]
    if not np.isfinite(period) or period <= 0:
        raise ValueError("cycle period must be finite positive")
    tangent = core_rhs(base, vpd, model)/scale
    speed = np.linalg.norm(tangent)
    if speed <= 1e-12:
        raise ValueError("section requires a nonstationary cycle base point")
    normal = tangent/speed
    basis = null_space(normal[None, :])
    def variation(t, augmented):
        y = augmented[:4]
        a = core_jacobian(y, vpd, model)*scale[None, :]/scale[:, None]
        return np.r_[core_rhs(y, vpd, model), (a@augmented[4:].reshape(4, 4)).ravel()]
    flow = solve_ivp(variation, (0., period), np.r_[base, np.eye(4).ravel()],
                     method="Radau", rtol=rtol, atol=rtol*.01, max_step=period/steps_per_period)
    if not flow.success:
        raise RuntimeError(f"E4 section variational IVP failed: {flow.message}")
    closure = float(np.linalg.norm((flow.y[:4, -1]-base)/scale))
    if closure > 1e-7:
        raise RuntimeError("E4 section cycle failed independent closure check")
    monodromy = flow.y[4:, -1].reshape(4, 4)
    end_tangent = core_rhs(flow.y[:4, -1], vpd, model)/scale
    projection = np.eye(4)-np.outer(end_tangent, normal)/np.dot(normal, end_tangent)
    derivative = basis.T@projection@monodromy@basis
    values, vectors = np.linalg.eig(derivative)
    if np.max(abs(values.imag)) > 1e-7 or np.max(abs(vectors.imag)) > 1e-7:
        raise RuntimeError("local E4 section requires resolved real transverse modes")
    order = np.argsort(values.real)[::-1]
    values, vectors = values.real[order], vectors.real[:, order]
    if not (values[0] > 1+1e-6 and np.all(abs(values[1:]) < 1-1e-6)):
        raise RuntimeError("section does not resolve one unstable and two stable transverse modes")
    vectors /= np.linalg.norm(vectors, axis=0)
    if equilibrium is not None:
        eq = np.asarray(equilibrium, dtype=float)
        if np.max(abs(model.rhs(eq, vpd))) > 1e-8:
            raise ValueError("section equilibrium reference must be a full frozen root")
        outward = basis.T@((base-eq[:4])/scale)
        if np.dot(vectors[:, 0], outward) < 0:
            vectors[:, 0] *= -1
    for j in range(1, 3):
        if vectors[np.argmax(abs(vectors[:, j])), j] < 0:
            vectors[:, j] *= -1
    covectors = np.linalg.inv(vectors)
    residual = float(np.max(abs(derivative@vectors-vectors*values)))
    if residual > 1e-7:
        raise RuntimeError("section eigendirections failed residual check")
    return {
        "base": base, "scale": scale, "vpd_mpa": vpd, "period_days": period,
        "normal_scaled": normal, "section_basis_scaled": basis,
        "monodromy_scaled": monodromy, "return_derivative": derivative,
        "multipliers": values, "right_vectors": vectors, "left_covectors": covectors,
        "eigenvector_residual": residual, "cycle_closure_scaled": closure,
        "unstable_growth_per_day": float(np.log(values[0])/period),
        "radial_e_folding_days": float(period/np.log(values[0])),
        "stable_tangent_covector_scaled": covectors[0]@basis.T,
        "interpretation": "linear_tangent_to_local_stable_manifold_not_nonlinear_global_basin",
    }


def section_return(section, coordinates, model, *, rtol=1e-10, steps_per_period=160):
    """One nonlinear positive-oriented return near one period, independent of DP."""
    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.shape != (3,) or not np.all(np.isfinite(coordinates)):
        raise ValueError("section coordinates must contain three finite scaled entries")
    if not np.isfinite(rtol) or not 0 < rtol < 1 or steps_per_period < 20:
        raise ValueError("section return needs positive tolerance and at least 20 steps per period")
    base, scale, basis = section["base"], section["scale"], section["section_basis_scaled"]
    normal, period, vpd = section["normal_scaled"], section["period_days"], section["vpd_mpa"]
    initial = base+scale*(basis@coordinates)
    if initial[2] <= model.feedback.minimum_leaf or initial[3] <= 0 or np.any(initial[:2] > 0):
        raise ValueError("section initial point lies outside the positive vegetation domain")
    def crossing(t, y):
        return float(np.dot(normal, (y-base)/scale))
    crossing.direction = 1
    def boundary(t, y):
        return y[2]-model.feedback.minimum_leaf
    boundary.terminal = True
    boundary.direction = -1
    flow = solve_ivp(lambda t, y: core_rhs(y, vpd, model), (0., 1.5*period), initial,
                     method="Radau", jac=lambda t, y: core_jacobian(y, vpd, model),
                     rtol=rtol, atol=rtol*.01, max_step=period/steps_per_period,
                     events=[crossing, boundary])
    if not flow.success:
        raise RuntimeError(f"nonlinear E4 section return failed: {flow.message}")
    if len(flow.t_events[1]):
        raise RuntimeError("local section trajectory exited the vegetation domain before return verification")
    indices = np.flatnonzero(flow.t_events[0] > .5*period)
    if len(indices) != 1:
        raise RuntimeError("expected exactly one local positive-oriented return near one period")
    i = indices[0]
    day, state = flow.t_events[0][i], flow.y_events[0][i]
    if not .75*period < day < 1.25*period:
        raise RuntimeError("return time left the declared local one-period neighborhood")
    if np.any(flow.y[:2] > 0) or np.any(flow.y[2:] <= 0):
        raise RuntimeError("accepted section trajectory left the physical domain")
    returned = basis.T@((state-base)/scale)
    return {
        "time_days": float(day), "state": state, "initial_state": initial,
        "coordinates": returned, "mode_coordinates": section["left_covectors"]@returned,
        "section_residual": float(abs(crossing(day, state))),
        "status": "finite_time_local_return_not_recovery_exit_classification",
    }
