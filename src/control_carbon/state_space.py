"""Generic continuous-time state-space and IRF utilities.

The module is intentionally model-agnostic. Model adapters (e.g. VISIT) should
construct A, B, C, D from source-grounded equations and pass them here.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike
from scipy.linalg import expm


@dataclass(frozen=True)
class ContinuousLTI:
    """Continuous-time linear time-invariant state-space system.

    dx/dt = A x + B u
        y = C x + D u
    """

    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    D: np.ndarray

    def __post_init__(self) -> None:
        A = np.asarray(self.A, dtype=float)
        B = np.asarray(self.B, dtype=float)
        C = np.asarray(self.C, dtype=float)
        D = np.asarray(self.D, dtype=float)
        n = A.shape[0]
        if A.shape != (n, n):
            raise ValueError("A must be square")
        if B.ndim != 2 or B.shape[0] != n:
            raise ValueError("B must have shape (n_state, n_input)")
        if C.ndim != 2 or C.shape[1] != n:
            raise ValueError("C must have shape (n_output, n_state)")
        if D.shape != (C.shape[0], B.shape[1]):
            raise ValueError("D must have shape (n_output, n_input)")
        object.__setattr__(self, "A", A)
        object.__setattr__(self, "B", B)
        object.__setattr__(self, "C", C)
        object.__setattr__(self, "D", D)

    @property
    def n_state(self) -> int:
        return self.A.shape[0]

    @property
    def n_input(self) -> int:
        return self.B.shape[1]

    @property
    def n_output(self) -> int:
        return self.C.shape[0]

    def poles(self) -> np.ndarray:
        """Return eigenvalues of A."""
        return np.linalg.eigvals(self.A)

    def frequency_response(self, omega: ArrayLike) -> np.ndarray:
        """Evaluate G(i omega) for each angular frequency.

        Returns shape (n_omega, n_output, n_input).
        """
        omega = np.atleast_1d(np.asarray(omega, dtype=float))
        eye = np.eye(self.n_state)
        out = np.empty((omega.size, self.n_output, self.n_input), dtype=complex)
        for k, w in enumerate(omega):
            out[k] = self.C @ np.linalg.solve(1j * w * eye - self.A, self.B) + self.D
        return out


@dataclass(frozen=True)
class DiscreteLTI:
    """Discrete-time linear time-invariant state-space system.

    x[k+1] = A x[k] + B u[k]
       y[k] = C x[k] + D u[k]
    """

    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    D: np.ndarray
    dt: float = 1.0

    def __post_init__(self) -> None:
        A = np.asarray(self.A, dtype=float)
        B = np.asarray(self.B, dtype=float)
        C = np.asarray(self.C, dtype=float)
        D = np.asarray(self.D, dtype=float)
        n = A.shape[0]
        if A.shape != (n, n):
            raise ValueError("A must be square")
        if B.ndim != 2 or B.shape[0] != n:
            raise ValueError("B must have shape (n_state, n_input)")
        if C.ndim != 2 or C.shape[1] != n:
            raise ValueError("C must have shape (n_output, n_state)")
        if D.shape != (C.shape[0], B.shape[1]):
            raise ValueError("D must have shape (n_output, n_input)")
        if not np.isfinite(self.dt) or self.dt <= 0:
            raise ValueError("dt must be finite and positive")
        object.__setattr__(self, "A", A)
        object.__setattr__(self, "B", B)
        object.__setattr__(self, "C", C)
        object.__setattr__(self, "D", D)
        object.__setattr__(self, "dt", float(self.dt))

    @property
    def n_state(self) -> int:
        return self.A.shape[0]

    @property
    def n_input(self) -> int:
        return self.B.shape[1]

    @property
    def n_output(self) -> int:
        return self.C.shape[0]

    def poles(self) -> np.ndarray:
        """Return eigenvalues of A."""
        return np.linalg.eigvals(self.A)

    def is_stable(self) -> bool:
        """Return whether every pole lies strictly inside the unit circle."""
        return bool(np.all(np.abs(self.poles()) < 1.0))

    def frequency_response(self, omega: ArrayLike) -> np.ndarray:
        """Evaluate G(exp(i omega dt)) for angular frequencies."""
        omega = np.atleast_1d(np.asarray(omega, dtype=float))
        eye = np.eye(self.n_state)
        out = np.empty((omega.size, self.n_output, self.n_input), dtype=complex)
        for k, w in enumerate(omega):
            z = np.exp(1j * w * self.dt)
            out[k] = self.C @ np.linalg.solve(z * eye - self.A, self.B) + self.D
        return out

    def impulse_response(self, n_steps: int) -> np.ndarray:
        """Return h[0]=D and h[k]=C A^(k-1) B for k >= 1."""
        if n_steps < 1:
            raise ValueError("n_steps must be at least 1")
        out = np.empty((n_steps, self.n_output, self.n_input), dtype=float)
        out[0] = self.D
        state_response = self.B.copy()
        for k in range(1, n_steps):
            out[k] = self.C @ state_response
            state_response = self.A @ state_response
        return out

    def step_response(self, n_steps: int) -> np.ndarray:
        """Return the response to a unit step in every input channel."""
        return np.cumsum(self.impulse_response(n_steps), axis=0)

    def forced_response(
        self, inputs: ArrayLike, initial_state: ArrayLike | None = None
    ) -> tuple[np.ndarray, np.ndarray]:
        """Simulate input samples, returning states and same-index outputs.

        The returned states include the initial state and therefore have one
        more time sample than the outputs.
        """
        inputs = np.asarray(inputs, dtype=float)
        if inputs.ndim != 2 or inputs.shape[1] != self.n_input:
            raise ValueError("inputs must have shape (n_steps, n_input)")
        if initial_state is None:
            initial_state_array = np.zeros(self.n_state, dtype=float)
        else:
            initial_state_array = np.asarray(initial_state, dtype=float)
            if initial_state_array.shape != (self.n_state,):
                raise ValueError("initial_state must have shape (n_state,)")

        states = np.empty((inputs.shape[0] + 1, self.n_state), dtype=float)
        outputs = np.empty((inputs.shape[0], self.n_output), dtype=float)
        states[0] = initial_state_array
        for k, input_sample in enumerate(inputs):
            outputs[k] = self.C @ states[k] + self.D @ input_sample
            states[k + 1] = self.A @ states[k] + self.B @ input_sample
        return states, outputs


@dataclass(frozen=True)
class DiscreteFixedPointResult:
    """Diagnostics for a fixed point under a constant discrete-time input."""

    state: np.ndarray
    residual_norm: float
    condition_number: float
    stable: bool
    converged: bool = True
    method: str = "linear solve"


@dataclass(frozen=True)
class DiscreteModalAnalysis:
    """Grouped spectral decomposition of a diagonalizable discrete system."""

    poles: np.ndarray
    multiplicities: np.ndarray
    decay_times: np.ndarray
    projectors: np.ndarray
    state_participation: np.ndarray
    residues: np.ndarray
    eigenvector_condition_number: float
    reconstruction_error: float
    direct_term: np.ndarray

    def impulse_response(self, n_steps: int) -> np.ndarray:
        """Reconstruct ``h[0]=D`` and dynamic samples from grouped residues."""
        if n_steps < 1:
            raise ValueError("n_steps must be at least 1")
        response = np.empty(
            (n_steps, self.direct_term.shape[0], self.direct_term.shape[1]),
            dtype=complex,
        )
        response[0] = self.direct_term
        for step in range(1, n_steps):
            response[step] = np.sum(
                self.poles[:, np.newaxis, np.newaxis] ** (step - 1)
                * self.residues,
                axis=0,
            )
        return np.real_if_close(response)


def discrete_fixed_point(
    system: DiscreteLTI, constant_input: ArrayLike
) -> DiscreteFixedPointResult:
    """Solve ``x = A x + B u`` for a constant input and report diagnostics."""
    constant_input_array = np.asarray(constant_input, dtype=float)
    if constant_input_array.shape != (system.n_input,):
        raise ValueError("constant_input must have shape (n_input,)")
    if not np.all(np.isfinite(constant_input_array)):
        raise ValueError("constant_input must be finite")

    fixed_point_matrix = np.eye(system.n_state) - system.A
    condition_number = float(np.linalg.cond(fixed_point_matrix))
    try:
        state = np.linalg.solve(fixed_point_matrix, system.B @ constant_input_array)
    except np.linalg.LinAlgError as error:
        raise ValueError("system does not have a unique finite fixed point") from error
    residual = system.A @ state + system.B @ constant_input_array - state
    return DiscreteFixedPointResult(
        state=state,
        residual_norm=float(np.linalg.norm(residual)),
        condition_number=condition_number,
        stable=system.is_stable(),
    )


def discrete_modal_analysis(
    system: DiscreteLTI, *, pole_tolerance: float = 1e-9
) -> DiscreteModalAnalysis:
    """Group repeated poles and compute spectral projectors and I/O residues.

    The diagonal of each grouped spectral projector is a basis-invariant state
    participation measure for that pole group. Non-diagonalizable systems are
    rejected because a simple residue sum cannot represent Jordan terms.
    """
    if not np.isfinite(pole_tolerance) or pole_tolerance <= 0:
        raise ValueError("pole_tolerance must be finite and positive")
    eigenvalues, eigenvectors = np.linalg.eig(system.A)
    grouped: list[list[complex]] = []
    for eigenvalue in sorted(eigenvalues, key=lambda value: (value.real, value.imag)):
        for group in grouped:
            representative = sum(group) / len(group)
            if abs(eigenvalue - representative) <= pole_tolerance * max(
                1.0, abs(eigenvalue), abs(representative)
            ):
                group.append(complex(eigenvalue))
                break
        else:
            grouped.append([complex(eigenvalue)])
    poles = np.asarray([sum(group) / len(group) for group in grouped], dtype=complex)
    multiplicities = np.asarray([len(group) for group in grouped], dtype=int)
    identity = np.eye(system.n_state, dtype=complex)
    matrix_a = system.A.astype(complex)
    projectors = []
    for pole_index, pole in enumerate(poles):
        projector = identity.copy()
        for other_index, other_pole in enumerate(poles):
            if other_index != pole_index:
                projector = projector @ (
                    (matrix_a - other_pole * identity) / (pole - other_pole)
                )
        projectors.append(projector)
    projector_array = np.asarray(projectors)
    reconstructed_identity = projector_array.sum(axis=0)
    reconstructed_a = np.sum(
        poles[:, np.newaxis, np.newaxis] * projector_array, axis=0
    )
    reconstruction_error = float(
        max(
            np.linalg.norm(reconstructed_identity - identity),
            np.linalg.norm(reconstructed_a - matrix_a),
        )
    )
    scale = max(1.0, np.linalg.norm(matrix_a))
    if reconstruction_error > 100.0 * pole_tolerance * scale:
        raise ValueError(
            "system is not diagonalizable into simple grouped spectral projectors"
        )
    residues = np.asarray(
        [system.C @ projector @ system.B for projector in projector_array]
    )
    magnitudes = np.abs(poles)
    decay_times = np.full(poles.shape, np.nan, dtype=float)
    zero = magnitudes == 0.0
    stable = (magnitudes > 0.0) & (magnitudes < 1.0)
    decay_times[zero] = 0.0
    decay_times[stable] = -system.dt / np.log(magnitudes[stable])
    decay_times[magnitudes == 1.0] = np.inf
    return DiscreteModalAnalysis(
        poles=np.real_if_close(poles),
        multiplicities=multiplicities,
        decay_times=decay_times,
        projectors=np.real_if_close(projector_array),
        state_participation=np.real_if_close(
            np.diagonal(projector_array, axis1=1, axis2=2)
        ),
        residues=np.real_if_close(residues),
        eigenvector_condition_number=float(np.linalg.cond(eigenvectors)),
        reconstruction_error=reconstruction_error,
        direct_term=system.D.copy(),
    )


def discrete_transfer_function_matrix(system: DiscreteLTI, z: complex) -> np.ndarray:
    """Evaluate G(z) = C (zI-A)^-1 B + D."""
    eye = np.eye(system.n_state)
    return system.C @ np.linalg.solve(z * eye - system.A, system.B) + system.D


def transfer_function_matrix(system: ContinuousLTI, s: complex) -> np.ndarray:
    """Evaluate G(s) = C (sI-A)^-1 B + D."""
    eye = np.eye(system.n_state)
    return system.C @ np.linalg.solve(s * eye - system.A, system.B) + system.D


def impulse_response(system: ContinuousLTI, times: ArrayLike) -> np.ndarray:
    """Continuous-time impulse response for t > 0.

    Returns C exp(A t) B with shape (n_time, n_output, n_input).
    The Dirac feedthrough contribution D delta(t) is not represented in the
    returned finite-valued array and must be handled separately if D != 0.
    """
    times = np.atleast_1d(np.asarray(times, dtype=float))
    if np.any(times < 0):
        raise ValueError("times must be non-negative")
    out = np.empty((times.size, system.n_output, system.n_input), dtype=float)
    for k, t in enumerate(times):
        out[k] = system.C @ expm(system.A * t) @ system.B
    return out


def state_impulse_response(A: ArrayLike, B: ArrayLike, times: ArrayLike) -> np.ndarray:
    """Return exp(A t) B, useful for carbon-pool IRFs."""
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    times = np.atleast_1d(np.asarray(times, dtype=float))
    out = np.empty((times.size, A.shape[0], B.shape[1]), dtype=float)
    for k, t in enumerate(times):
        out[k] = expm(A * t) @ B
    return out
