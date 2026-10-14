"""Capacity-first diagnostics for unchanged E0-E4 equations, MPa and day.

K is a pointwise frozen-input map, P=K-C is potential, and neither is
automatically a coupled equilibrium or an attracting target. M is constant
in these models; the derivative identities below rely on that restriction.
"""

from dataclasses import dataclass

import numpy as np

from .spac_carbon import (
    assimilation, carbon_budget, carbon_jacobian, carbon_matrix,
    frozen_coefficient_capacity, hydraulic_factor, positive_carbon_equilibrium,
    reproduction_number,
)
from .spac_carbon_mpa import PressureSPAC


@dataclass(frozen=True)
class CapacityScales:
    """Fixed diagonal state scales: two MPa, then four kg C/m2 stocks."""

    state: tuple[float, float, float, float, float, float]

    def __post_init__(self):
        values = np.asarray(self.state, dtype=float)
        if values.shape != (6,) or not np.all(np.isfinite(values)) or np.any(values <= 0):
            raise ValueError("capacity scales require six finite positive state scales")


@dataclass(frozen=True)
class ConditionalCarbonRoot:
    state: np.ndarray
    kind: str
    stability: str
    eigenvalues: np.ndarray
    in_positive_leaf_domain: bool


def reference_scales(model: PressureSPAC) -> CapacityScales:
    """Keep structural references fixed, not the current moving equilibrium."""
    c, f = model.carbon, model.feedback
    root = f.leaf_reference*c.allocation[2]*c.losses[0]/(c.allocation[0]*c.losses[2])
    soil = np.dot(np.array(c.litter_fractions)*np.array(c.losses[:3]),
                  [f.leaf_reference, f.wood_reference, root])/c.losses[3]
    return CapacityScales((1., 1., f.leaf_reference, f.wood_reference, float(root), float(soil)))


def conditional_carbon_roots(leaf_potential_mpa: float, model: PressureSPAC):
    """All carbon roots at fixed water, including the formal zero root."""
    c = model.carbon
    h = hydraulic_factor(leaf_potential_mpa, model.hydraulic)
    r0 = reproduction_number(h, c)
    stability = ("unstable" if r0 > 1 else "attracting_nonhyperbolic_nonnegative_side"
                 if r0 == 1 else "attracting")
    zero = np.zeros(4)
    roots = [ConditionalCarbonRoot(zero, "formal_zero", stability,
                                   np.linalg.eigvals(carbon_jacobian(zero, h, c)), False)]
    positive = positive_carbon_equilibrium(h, c)
    if positive is not None:
        roots.append(ConditionalCarbonRoot(
            positive, "positive", "attracting",
            np.linalg.eigvals(carbon_jacobian(positive, h, c)),
            bool(positive[0] > model.feedback.minimum_leaf),
        ))
    return tuple(roots)


def capacity_map(state, vpd_mpa: float, model: PressureSPAC):
    model.rhs(state, vpd_mpa)
    y = np.asarray(state, dtype=float)
    return frozen_coefficient_capacity(y[2:], hydraulic_factor(y[1], model.hydraulic), model.carbon)


def capacity_jacobian(state, vpd_mpa: float, model: PressureSPAC):
    """D_y K = E_C - M^-1 J_Cy, exact for the fixed-M E0-E4 family."""
    jac = model.jacobian(state, vpd_mpa)
    return np.eye(6)[2:] - np.linalg.solve(carbon_matrix(model.carbon), jac[2:])


def spectrum_status(eigenvalues, *, tolerance=1e-10):
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("spectral tolerance must be finite positive (day^-1)")
    values = np.asarray(eigenvalues, dtype=complex)
    if values.ndim != 1 or not len(values) or not np.all(np.isfinite(values)):
        raise ValueError("a finite nonempty eigenvalue vector is required")
    largest = max(values.real)
    return ("unstable" if largest > tolerance else "stable" if largest < -tolerance
            else "nonhyperbolic_within_tolerance")


def diagnose_capacity(state, vpd_mpa: float, model: PressureSPAC, scales: CapacityScales,
                      *, equilibrium=None):
    """Pure diagnostics; an optional reference must be a full frozen root.

    Norms use the declared fixed diagonal scales, never mixed raw MPa/carbon.
    Rates are fixed-v derivatives; K has no direct v dependence in E0-E4.
    Local spectra away from equilibrium are not stability certificates.
    """
    rhs = model.rhs(state, vpd_mpa)
    y = np.asarray(state, dtype=float)
    scale = np.asarray(scales.state)
    c, h = model.carbon, hydraulic_factor(y[1], model.hydraulic)
    matrix = carbon_matrix(c)
    k = frozen_coefficient_capacity(y[2:], h, c)
    potential = k-y[2:]
    dk = capacity_jacobian(y, vpd_mpa, model)
    kdot = dk@rhs
    pdot = kdot-rhs[2:]
    roots = conditional_carbon_roots(y[1], model)
    target = roots[-1]
    jac = model.jacobian(y, vpd_mpa)
    carbon_eigenvalues = np.linalg.eigvals(jac[2:, 2:])
    full_eigenvalues = np.linalg.eigvals(jac)
    gpp = assimilation(y[2], h, c)
    input_vector = np.r_[c.allocation, 0.]*gpp
    erx, _, eatm = model.fluxes(y, vpd_mpa)
    norm = lambda value: float(np.linalg.norm(value))
    result = {
        "input": input_vector, "matrix": matrix, "matrix_condition_2": float(np.linalg.cond(matrix)),
        "matrix_eigenvalues": np.linalg.eigvals(matrix),
        "capacity": k, "potential": potential, "potential_sum": float(sum(potential)),
        "potential_norm_scaled": norm(potential/scale[2:]),
        "capacity_jacobian": dk, "capacity_rate": kdot, "potential_rate": pdot,
        "capacity_rate_norm_scaled_per_day": norm(kdot/scale[2:]),
        "potential_energy_rate_scaled_per_day": float(np.dot(potential/scale[2:], pdot/scale[2:])),
        "rhs": rhs, "full_residual_norm_scaled_per_day": norm(rhs/scale),
        "carbon_residual_norm_scaled_per_day": norm(rhs[2:]/scale[2:]),
        "water_residual_norm_scaled_per_day": norm(rhs[:2]/scale[:2]),
        "identity_residual_max": float(max(abs(matrix@(y[2:]-k)-rhs[2:]))),
        "potential_identity_residual_max": float(max(abs(potential+np.linalg.solve(matrix, rhs[2:])))),
        "carbon_budget_residual": carbon_budget(y[2:], h, c)["residual"],
        "water_budget_residual": float(model.hydraulic.cx*rhs[0]+model.hydraulic.cl*rhs[1]-erx+eatm),
        "reproduction_number": float(reproduction_number(h, c)),
        "conditional_roots": roots, "conditional_target_kind": target.kind,
        "conditional_target_in_positive_leaf_domain": target.in_positive_leaf_domain,
        "conditional_distance_norm_scaled": norm((y[2:]-target.state)/scale[2:]),
        "local_carbon_eigenvalues": carbon_eigenvalues,
        "local_full_eigenvalues": full_eigenvalues,
        "state_above_leaf_exit_boundary": bool(y[2] > model.feedback.minimum_leaf),
        "equilibrium_reference_status": "not_supplied",
        "equilibrium_error_norm_scaled": None,
        "equilibrium_eigenvalues": None, "equilibrium_stability": None,
        "equilibrium_unstable_dimensions": None, "equilibrium_center_dimensions": None,
        "equilibrium_type": None,
    }
    if equilibrium is not None:
        residual = model.rhs(equilibrium, vpd_mpa)
        if norm(residual/scale) > 1e-8:
            raise ValueError("equilibrium reference must satisfy the full frozen equations")
        eq = np.asarray(equilibrium, dtype=float)
        eigenvalues = np.linalg.eigvals(model.jacobian(eq, vpd_mpa))
        unstable = int(sum(eigenvalues.real > 1e-10))
        center = int(sum(abs(eigenvalues.real) <= 1e-10))
        kind = ("nonhyperbolic" if center else "hyperbolic_sink" if unstable == 0
                else "hyperbolic_source" if unstable == 6 else "hyperbolic_saddle")
        result.update({
            "equilibrium_reference_status": "verified_frozen_root",
            "equilibrium_error_norm_scaled": norm((y-eq)/scale),
            "equilibrium_eigenvalues": eigenvalues,
            "equilibrium_stability": spectrum_status(eigenvalues),
            "equilibrium_unstable_dimensions": unstable, "equilibrium_center_dimensions": center,
            "equilibrium_type": kind,
        })
    return result


def frozen_equilibria(vpd_mpa: float, leaf_bounds, model: PressureSPAC, scales: CapacityScales,
                      *, samples=300):
    """Enumerate conducting positive roots in a stated, sampled branch domain.

    The carbon elimination is exact; detection of stationary points is not
    interval-certified. No zero-leaf full root or post-exit attractor is added.
    """
    model.mole_fraction(vpd_mpa)
    stationary = model.stationary_points(leaf_bounds, samples=samples)
    states = model.equilibria(vpd_mpa, leaf_bounds, stationary)
    records = []
    for y in states:
        data = diagnose_capacity(y, vpd_mpa, model, scales, equilibrium=y)
        if data["full_residual_norm_scaled_per_day"] > 1e-8:
            raise RuntimeError("frozen root failed the full scaled residual check")
        records.append({"state": y, "stability": data["equilibrium_stability"], "diagnostics": data})
    return {
        "roots": tuple(records), "leaf_bounds_mpa": tuple(leaf_bounds),
        "detected_stationary_leaf_mpa": stationary, "stationary_scan_samples": samples,
        "coverage": "conducting_positive_roots_in_sampled_partition_not_global_certificate",
        "formal_zero_full_roots": "not_defined_by_positive_leaf_model",
    }
