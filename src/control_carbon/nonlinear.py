"""Model-agnostic deterministic continuous dynamics with one external driver."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike
from scipy.integrate import solve_ivp
from scipy.optimize import root


@dataclass(frozen=True)
class ContinuousNonlinear:
    """dx/dt = rhs(x, driver(t)); all other parameters are captured and fixed.

    Jacobians are supplied explicitly. This smooth continuous research API is
    not a representation of native VISIT's sequential daily update.
    """

    rhs: Callable[[np.ndarray, float], np.ndarray]
    jacobian: Callable[[np.ndarray, float], np.ndarray]
    state_names: tuple[str, ...]
    state_units: tuple[str, ...]
    driver_name: str
    driver_unit: str
    time_unit: str
    approximation: str

    def __post_init__(self) -> None:
        if not self.state_names or len(set(self.state_names)) != len(self.state_names):
            raise ValueError("state_names must be nonempty and unique")
        if len(self.state_units) != len(self.state_names):
            raise ValueError("state_units must match state_names")
        if not all((*self.state_units, self.driver_name, self.driver_unit,
                    self.time_unit, self.approximation)):
            raise ValueError("units, driver name and approximation must be specified")

    def evaluate(self, state: ArrayLike, driver: float) -> np.ndarray:
        state = _vector(state, len(self.state_names), "state")
        if not np.isfinite(driver):
            raise ValueError("driver must be finite")
        return _vector(self.rhs(state, driver), state.size, "rhs")

    def linearize(self, state: ArrayLike, driver: float) -> np.ndarray:
        self.evaluate(state, driver)
        state = np.asarray(state, dtype=float)
        matrix = np.asarray(self.jacobian(state, driver), dtype=float)
        if matrix.shape != (state.size, state.size) or not np.all(np.isfinite(matrix)):
            raise ValueError("jacobian must be finite with shape (n_state, n_state)")
        return matrix


def _vector(values: ArrayLike, size: int, name: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.shape != (size,) or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite with shape ({size},)")
    return array


@dataclass(frozen=True)
class NonlinearTrajectory:
    system: ContinuousNonlinear
    times: np.ndarray
    states: np.ndarray
    drivers: np.ndarray
    rtol: float
    atol: float
    max_step: float
    function_evaluations: int
    method: str = "DOP853"


def simulate_continuous(
    system: ContinuousNonlinear,
    initial_state: ArrayLike,
    times: ArrayLike,
    driver: Callable[[float], float],
    *,
    rtol: float = 1e-9,
    atol: float = 1e-11,
    max_step: float = np.inf,
) -> NonlinearTrajectory:
    """Integrate with SciPy DOP853; no state clipping or stochastic terms.

    Output times do not constrain internal solver steps. Set max_step to
    resolve fast forcing; discontinuous protocols require split integrations.
    """
    times = np.asarray(times, dtype=float)
    if (times.ndim != 1 or times.size < 2 or not np.all(np.isfinite(times))
            or np.any(np.diff(times) <= 0)):
        raise ValueError("times must be finite and strictly increasing, with >= 2 entries")
    if not np.isfinite(rtol) or not np.isfinite(atol) or rtol <= 0 or atol <= 0:
        raise ValueError("rtol and atol must be finite and positive")
    if np.isnan(max_step) or max_step <= 0:
        raise ValueError("max_step must be positive")
    initial_state = _vector(initial_state, len(system.state_names), "initial_state")
    solution = solve_ivp(
        lambda time, state: system.evaluate(state, driver(time)),
        (times[0], times[-1]), initial_state, t_eval=times,
        method="DOP853", rtol=rtol, atol=atol, max_step=max_step,
    )
    if not solution.success:
        raise RuntimeError(f"ODE integration failed: {solution.message}")
    drivers = _vector([driver(time) for time in times], times.size, "drivers")
    if not np.all(np.isfinite(solution.y)):
        raise RuntimeError("ODE integration returned nonfinite states")
    return NonlinearTrajectory(
        system, times.copy(), solution.y.T.copy(), drivers,
        rtol, atol, max_step, solution.nfev,
    )


@dataclass(frozen=True)
class ContinuousEquilibrium:
    """Local root diagnostics, not a certificate that all equilibria were found."""

    state: np.ndarray
    driver: float
    residual_norm: float
    converged: bool
    jacobian: np.ndarray
    eigenvalues: np.ndarray
    stability: str
    condition_number: float
    function_evaluations: int
    message: str


def continuous_equilibrium(
    system: ContinuousNonlinear,
    initial_guess: ArrayLike,
    driver: float,
    *,
    residual_tolerance: float = 1e-9,
    stability_tolerance: float = 1e-8,
) -> ContinuousEquilibrium:
    """Solve a frozen-driver root using an analytic state Jacobian.

    Zero/near-zero real parts are nonhyperbolic, not classified as stable.
    A moving sequence of these roots is a QSE, not an ODE solution.
    """
    if (not np.isfinite(residual_tolerance) or residual_tolerance <= 0
            or not np.isfinite(stability_tolerance) or stability_tolerance <= 0):
        raise ValueError("tolerances must be finite and positive")
    initial_guess = _vector(initial_guess, len(system.state_names), "initial_guess")
    solution = root(
        lambda state: system.evaluate(state, driver), initial_guess,
        jac=lambda state: system.linearize(state, driver), method="hybr",
    )
    residual = float(np.linalg.norm(system.evaluate(solution.x, driver), ord=np.inf))
    matrix = system.linearize(solution.x, driver)
    eigenvalues = np.linalg.eigvals(matrix)
    converged = bool(solution.success and residual <= residual_tolerance)
    stability = "unresolved"
    if converged:
        if np.any(eigenvalues.real > stability_tolerance):
            stability = "unstable"
        elif np.all(eigenvalues.real < -stability_tolerance):
            stability = "stable"
        else:
            stability = "nonhyperbolic"
    return ContinuousEquilibrium(
        solution.x.copy(), float(driver), residual, converged, matrix,
        eigenvalues, stability, float(np.linalg.cond(matrix)),
        solution.nfev, str(solution.message),
    )


def compartment_rhs(
    state: ArrayLike, transfer: ArrayLike, turnover: ArrayLike,
    environment: ArrayLike, allocation: ArrayLike, inputs: ArrayLike,
) -> np.ndarray:
    """A diag(xi) diag(k) x + B mu, with signed donor-column A.

    A has diagonal -1, nonnegative transfers and nonpositive column sums.
    Missing column mass exits the system. B columns sum to one, so mu
    represents total external input flux per channel. Arrays may be evaluated
    from state/forcing by the caller; this function validates each evaluation.
    """
    state = np.asarray(state, dtype=float)
    if state.ndim != 1 or state.size == 0 or not np.all(np.isfinite(state)):
        raise ValueError("state must be a nonempty finite vector")
    size = state.size
    transfer = np.asarray(transfer, dtype=float)
    if transfer.shape != (size, size) or not np.all(np.isfinite(transfer)):
        raise ValueError("transfer must be a finite square matrix matching state")
    off_diagonal = transfer.copy()
    np.fill_diagonal(off_diagonal, 0.0)
    if (not np.allclose(np.diag(transfer), -1.0, rtol=0, atol=1e-12)
            or np.any(off_diagonal < 0) or np.any(transfer.sum(axis=0) > 1e-12)):
        raise ValueError("transfer violates signed compartment convention")
    turnover = _vector(turnover, size, "turnover")
    environment = _vector(environment, size, "environment")
    inputs = np.asarray(inputs, dtype=float)
    if inputs.ndim != 1 or not np.all(np.isfinite(inputs)):
        raise ValueError("inputs must be a finite vector")
    allocation = np.asarray(allocation, dtype=float)
    if (allocation.shape != (size, inputs.size)
            or not np.all(np.isfinite(allocation)) or np.any(allocation < 0)
            or not np.allclose(allocation.sum(axis=0), 1.0, rtol=0, atol=1e-12)):
        raise ValueError("allocation columns must be nonnegative and sum to one")
    if np.any(turnover < 0) or np.any(environment < 0) or np.any(inputs < 0):
        raise ValueError("rates, environment and inputs must be nonnegative")
    return transfer @ (environment * turnover * state) + allocation @ inputs