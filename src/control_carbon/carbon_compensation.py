"""Minimum terminal compensation for a prescribed, exogenous capacity path.

These are ordinary floating-point control calculations, not interval
certificates or a controlled nonlinear E4 model.
"""

from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import expm, solve_triangular


@dataclass(frozen=True)
class CompensationGramian:
    matrix: np.ndarray
    input_matrix: np.ndarray
    weight: np.ndarray
    duration_days: float
    gramian: np.ndarray
    left_vectors: np.ndarray
    singular_values: np.ndarray
    resolved_rank: int
    rank_tolerance: float
    quadrature_order: int

    def minimum_action(self, natural_terminal):
        target = np.asarray(natural_terminal, dtype=float)
        n = len(self.matrix)
        if target.shape != (n,) or np.any(~np.isfinite(target)):
            raise ValueError("Terminal potential must have one finite entry per pool")
        coordinates = self.left_vectors.T@target
        unresolved = np.linalg.norm(coordinates[self.resolved_rank:])
        if np.any(coordinates[self.resolved_rank:] != 0):
            raise ValueError(
                "Terminal target is outside the numerically resolved reachable subspace; "
                "it is unreachable or reachability needs higher precision")
        values = self.singular_values[:self.resolved_rank]
        dual = self.left_vectors[:, :self.resolved_rank]@(coordinates[:self.resolved_rank]/values**2)
        energy = .5*np.sum((coordinates[:self.resolved_rank]/values)**2)
        if not np.isfinite(energy) or np.any(~np.isfinite(dual)):
            raise RuntimeError("Minimum compensation is not resolved at working precision")
        return MinimumCompensation(self, target.copy(), dual, float(energy), float(unresolved))

    def output_budget_interval(self, natural_terminal, output_weights, energy_budget):
        target, weights = np.asarray(natural_terminal, dtype=float), np.asarray(output_weights, dtype=float)
        n = len(self.matrix)
        if (target.shape != (n,) or weights.shape != (n,)
                or np.any(~np.isfinite(target)) or np.any(~np.isfinite(weights))
                or not np.isfinite(energy_budget) or energy_budget < 0):
            raise ValueError("Output budget needs finite vectors and a nonnegative energy")
        # Retain every quadrature singular value, not only resolved control directions.
        variance = np.sum((self.singular_values*(self.left_vectors.T@weights))**2)
        radius = float(np.sqrt(2*energy_budget*variance))
        center = float(weights@target)
        if not np.isfinite(radius) or not np.isfinite(center):
            raise RuntimeError("Output compensation bound overflowed at working precision")
        return {"signed_lower": center-radius, "signed_upper": center+radius,
                "unavoidable_absolute_lower": max(0., abs(center)-radius),
                "compensation_radius": radius}


@dataclass(frozen=True)
class MinimumCompensation:
    problem: CompensationGramian
    natural_terminal: np.ndarray
    multiplier: np.ndarray
    action: float
    unresolved_projection_norm: float

    def control(self, day):
        if not np.isfinite(day) or not 0 <= day <= self.problem.duration_days:
            raise ValueError("Control time must lie in the declared horizon")
        rhs = (self.problem.input_matrix.T
               @expm(self.problem.matrix.T*(self.problem.duration_days-day))@self.multiplier)
        value = -np.linalg.solve(self.problem.weight, rhs)
        if np.any(~np.isfinite(value)):
            raise RuntimeError("Nonfinite optimal compensation flux")
        return value


def compensation_gramian(matrix, input_matrix, weight, duration_days, *, quadrature_order=64):
    m, b, r = [np.asarray(x, dtype=float) for x in [matrix, input_matrix, weight]]
    if (m.ndim != 2 or len(m) == 0 or m.shape[0] != m.shape[1] or b.ndim != 2 or b.shape[0] != len(m)
            or b.shape[1] == 0 or r.shape != (b.shape[1], b.shape[1])
            or any(np.any(~np.isfinite(x)) for x in [m, b, r])
            or not np.array_equal(r, r.T) or not np.isfinite(duration_days) or duration_days <= 0
            or not isinstance(quadrature_order, int) or quadrature_order < 4):
        raise ValueError("Gramian needs finite compatible matrices, symmetric R and positive horizon")
    try:
        factor_r = np.linalg.cholesky(r)
    except np.linalg.LinAlgError as error:
        raise ValueError("Compensation weight R must be positive definite") from error
    whitened = solve_triangular(factor_r, b.T, lower=True).T
    nodes, weights = leggauss(quadrature_order)
    times = (nodes+1)*duration_days/2
    factor = np.concatenate([
        np.sqrt(w*duration_days/2)*expm(m*t)@whitened for t, w in zip(times, weights)], axis=1)
    left, values, _ = np.linalg.svd(factor, full_matrices=False)
    if len(values) != len(m):
        raise ValueError("Quadrature factor has fewer columns than states; increase its order")
    cutoff = np.finfo(float).eps*max(factor.shape)*(values[0] if len(values) else 0.)
    rank = int(np.count_nonzero(values > cutoff))
    gramian = factor@factor.T
    if np.any(~np.isfinite(gramian)):
        raise RuntimeError("Nonfinite compensation Gramian")
    return CompensationGramian(
        m.copy(), b.copy(), r.copy(), float(duration_days), gramian, left, values,
        rank, float(cutoff), quadrature_order)
