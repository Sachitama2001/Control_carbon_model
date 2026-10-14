"""Moving full-equilibrium targets, not frozen-input K targets (MPa/day)."""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_continuous_lyapunov
from scipy.optimize import brentq

from .manzoni_hydraulics import conductances
from .spac_carbon_mpa import PressureSPAC
from .spac_carbon_capacity import CapacityScales
from .spac_carbon_two_way import feedback_factors


@dataclass(frozen=True)
class SmoothVPDRamp:
    start_mpa: float
    end_mpa: float
    duration_days: float

    def __post_init__(self):
        values = [self.start_mpa, self.end_mpa, self.duration_days]
        if not np.all(np.isfinite(values)) or min(values[:2]) < 0 or self.duration_days <= 0:
            raise ValueError("ramp requires nonnegative finite MPa endpoints and positive duration")

    def value_rate(self, day):
        if not np.isfinite(day):
            raise ValueError("ramp time must be finite (day)")
        if day <= 0:
            return self.start_mpa, 0.
        if day >= self.duration_days:
            return self.end_mpa, 0.
        s = day/self.duration_days
        delta = self.end_mpa-self.start_mpa
        return (self.start_mpa+delta*s*s*(3-2*s),
                delta*6*s*(1-s)/self.duration_days)


def forcing_derivative(state, vpd_mpa, model: PressureSPAC):
    """Exact partial F/partial v at fixed state; only leaf-water RHS depends on v."""
    model.rhs(state, vpd_mpa)
    _, gs = conductances(state[1], model.hydraulic)
    _, fe, *_ = feedback_factors(state[2], state[3], model.feedback)
    derivative = np.zeros(6)
    derivative[1] = -model.gs_max_per_mpa*(gs/model.hydraulic.gs_max)*fe/model.hydraulic.cl
    return derivative


@dataclass(frozen=True)
class FrozenBranch:
    model: PressureSPAC
    leaf_bounds_mpa: tuple[float, float]

    def __post_init__(self):
        bounds = np.asarray(self.leaf_bounds_mpa, dtype=float)
        if bounds.shape != (2,) or not np.all(np.isfinite(bounds)) or bounds[0] >= bounds[1]:
            raise ValueError("branch needs two ordered finite leaf-potential bounds")
        slopes = np.array([self.model.equilibrium_curve(z)[2] for z in np.linspace(*bounds, 101)])
        folds = self.model.stationary_points(bounds)
        if folds or not (np.all(slopes > 0) or np.all(slopes < 0)):
            raise ValueError("branch bracket must exclude folds and sampled slope reversals")

    def state(self, vpd_mpa):
        self.model.mole_fraction(vpd_mpa)
        endpoints = [self.model.equilibrium_curve(z)[1] for z in self.leaf_bounds_mpa]
        if not min(endpoints) <= vpd_mpa <= max(endpoints):
            raise ValueError("VPD lies outside this declared monotone frozen branch")
        psi = brentq(lambda z: self.model.equilibrium_curve(z)[1]-vpd_mpa,
                     *self.leaf_bounds_mpa, xtol=1e-13)
        y = self.model.equilibrium_curve(psi)[0]
        if np.max(abs(self.model.rhs(y, vpd_mpa))) > 1e-9:
            raise RuntimeError("moving equilibrium failed full frozen residual")
        return y

    def tangent(self, vpd_mpa, scales: CapacityScales):
        y = self.state(vpd_mpa)
        scale = np.asarray(scales.state)
        jac = self.model.jacobian(y, vpd_mpa)
        a = jac*scale[None, :]/scale[:, None]
        condition = np.linalg.cond(a)
        if not np.isfinite(condition) or condition > 1e12:
            raise ValueError("equilibrium Jacobian is singular or ill-conditioned near a fold")
        tangent = np.linalg.solve(a, -forcing_derivative(y, vpd_mpa, self.model)/scale)*scale
        return y, tangent, float(condition)


def tracking_diagnostics(state, vpd_mpa, vpd_rate_mpa_per_day, branch: FrozenBranch,
                         scales: CapacityScales):
    if not np.isfinite(vpd_rate_mpa_per_day):
        raise ValueError("VPD rate must be finite MPa/day")
    ystar, tangent, condition = branch.tangent(vpd_mpa, scales)
    scale = np.asarray(scales.state)
    rhs = branch.model.rhs(state, vpd_mpa)
    equilibrium_rhs = branch.model.rhs(ystar, vpd_mpa)
    error = (np.asarray(state)-ystar)/scale
    a = branch.model.jacobian(ystar, vpd_mpa)*scale[None, :]/scale[:, None]
    eigenvalues, vectors = np.linalg.eig(a)
    spectral_abscissa = float(max(eigenvalues.real))
    numerical_abscissa = float(max(np.linalg.eigvalsh((a+a.T)/2)))
    motion = tangent*vpd_rate_mpa_per_day/scale
    exact = (rhs-equilibrium_rhs)/scale-motion
    linear = a@error-motion
    return {
        "equilibrium": ystar, "equilibrium_tangent_per_mpa": tangent,
        "error_scaled": error, "error_norm_scaled": float(np.linalg.norm(error)),
        "branch_speed_norm_scaled_per_day": float(np.linalg.norm(motion)),
        "exact_error_rate_scaled": exact, "linearized_error_rate_scaled": linear,
        "nonlinear_remainder_norm_scaled_per_day": float(np.linalg.norm(exact-linear)),
        "equilibrium_residual_norm_scaled_per_day": float(np.linalg.norm(equilibrium_rhs/scale)),
        "scaled_jacobian": a, "eigenvalues_per_day": eigenvalues,
        "spectral_abscissa_per_day": spectral_abscissa,
        "numerical_abscissa_per_day": numerical_abscissa,
        "eigenvector_condition_2": float(np.linalg.cond(vectors)),
        "scaled_jacobian_condition_2": condition,
        "slowest_decay_days": -1/spectral_abscissa if spectral_abscissa < 0 else None,
    }


def linearized_tracking(times, ramp: SmoothVPDRamp, branch: FrozenBranch, scales: CapacityScales,
                        *, initial_error_scaled=None, rtol=1e-9):
    """Integrate e'=J_*(v(t))e-y*'(v(t))v' in fixed scaled coordinates."""
    times = np.asarray(times, dtype=float)
    if (times.ndim != 1 or len(times) < 2 or not np.all(np.isfinite(times))
            or times[0] != 0 or np.any(np.diff(times) <= 0)):
        raise ValueError("tracking times must start at zero and increase strictly")
    if not np.isfinite(rtol) or not 0 < rtol < 1:
        raise ValueError("rtol must lie in (0,1)")
    error0 = np.zeros(6) if initial_error_scaled is None else np.asarray(initial_error_scaled, dtype=float)
    if error0.shape != (6,) or not np.all(np.isfinite(error0)):
        raise ValueError("linearized initial error must have six finite scaled entries")
    scale = np.asarray(scales.state)
    def coefficients(t):
        v, rate = ramp.value_rate(t)
        y, tangent, _ = branch.tangent(v, scales)
        a = branch.model.jacobian(y, v)*scale[None, :]/scale[:, None]
        return a, tangent*rate/scale
    def rhs(t, e):
        a, motion = coefficients(t)
        return a@e-motion
    solution = solve_ivp(rhs, (0., times[-1]), error0, method="Radau", rtol=rtol,
                         atol=rtol*.01, t_eval=times, jac=lambda t, e: coefficients(t)[0])
    if not solution.success:
        raise RuntimeError(f"linearized tracking IVP failed: {solution.message}")
    if not np.all(np.isfinite(solution.y)):
        raise RuntimeError("linearized tracking produced nonfinite errors")
    return {"time_days": times, "error_scaled": solution.y.T, "nfev": solution.nfev}


def frozen_lyapunov_metric(vpd_mpa, branch: FrozenBranch, scales: CapacityScales, *, step_mpa=1e-7):
    """Frozen A.T Q+Q A=-I and its parameter derivative, not a uniform certificate.

    Q_v uses a central difference of A(v) and a differentiated Lyapunov solve.
    For moving v the linear metric drift is -I+Q_v*v', not merely -I.
    """
    if not np.isfinite(step_mpa) or step_mpa <= 0:
        raise ValueError("metric derivative step must be finite positive MPa")
    scale = np.asarray(scales.state)
    def matrix(v):
        y = branch.state(v)
        return branch.model.jacobian(y, v)*scale[None, :]/scale[:, None]
    a = matrix(vpd_mpa)
    if max(np.linalg.eigvals(a).real) >= -1e-10:
        raise ValueError("frozen Lyapunov metric requires a resolved stable equilibrium")
    q = solve_continuous_lyapunov(a.T, -np.eye(6))
    q = (q+q.T)/2
    eigenvalues = np.linalg.eigvalsh(q)
    if eigenvalues[0] <= 0 or not np.all(np.isfinite(q)):
        raise RuntimeError("frozen Lyapunov solve did not yield a positive-definite metric")
    av = (matrix(vpd_mpa+step_mpa)-matrix(vpd_mpa-step_mpa))/(2*step_mpa)
    forcing = -(av.T@q+q@av)
    qv = solve_continuous_lyapunov(a.T, forcing)
    qv = (qv+qv.T)/2
    residual = float(np.max(abs(a.T@q+q@a+np.eye(6))))
    derivative_residual = float(np.max(abs(a.T@qv+qv@a-forcing)))
    if residual > 1e-6 or derivative_residual > 1e-6*max(1., np.max(abs(forcing))):
        raise RuntimeError("frozen Lyapunov metric failed its equation residual")
    return {"metric": q, "metric_derivative_per_mpa": qv,
            "minimum_metric_eigenvalue_days": float(eigenvalues[0]),
            "maximum_metric_eigenvalue_days": float(eigenvalues[-1]),
            "lyapunov_residual": residual, "derivative_equation_residual": derivative_residual,
            "parameter_difference_step_mpa": step_mpa}
