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
