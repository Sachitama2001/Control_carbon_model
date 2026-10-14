from dataclasses import replace

import numpy as np
import pytest

from control_carbon.capacity_feedback import CubicCapacityFeedback
from control_carbon.manzoni_hydraulics import HydraulicParameters, transpiring_equilibria
from control_carbon.spac_carbon import (
    YEAR_DAYS, CarbonParameters, carbon_rhs, combined_rhs,
    frozen_coefficient_capacity, hydraulic_factor, positive_carbon_equilibrium,
)
from control_carbon.spac_carbon_two_way import StructuralFeedback
from control_carbon.spac_carbon_mpa import PressureSPAC
from control_carbon.spac_carbon_capacity import (
    CapacityScales, capacity_map, capacity_jacobian, conditional_carbon_roots,
    diagnose_capacity, frozen_equilibria, reference_scales, spectrum_status,
)

H = HydraulicParameters(.0086, .47, -1.5, -2., .001, .001)
C = CarbonParameters(1.2/YEAR_DAYS, (.3, .4, .3), tuple(k/YEAR_DAYS for k in [.8, .05, .2, .04]),
                     (.7, .9, .8), 10., .5, .25)
LEAF = transpiring_equilibria(0., .01, H)[0]
REF = positive_carbon_equilibrium(hydraulic_factor(LEAF, H), C)
F = StructuralFeedback(REF[0], REF[1], .01, 1e-6)
M = PressureSPAC(H, C, F, 0., .101325)
VARIANTS = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1)]


def model(strengths):
    return replace(M, feedback=replace(F, leaf_supply=strengths[0], leaf_demand=strengths[1],
                                      wood_supply=strengths[2]))


@pytest.mark.parametrize("strengths", VARIANTS)
def test_capacity_identities_jacobian_budgets_and_purity(strengths):
    m = model(strengths)
    y = np.array([-.8, -1.2, .2, 2., .4, 5.])
    saved = y.copy()
    scales = reference_scales(m)
    data = diagnose_capacity(y, .002, m, scales)
    np.testing.assert_array_equal(y, saved)
    np.testing.assert_allclose(data["potential"], data["capacity"]-y[2:], atol=1e-14)
    assert data["identity_residual_max"] < 1e-15
    assert data["potential_identity_residual_max"] < 1e-14
    assert abs(data["carbon_budget_residual"]) < 1e-15
    assert abs(data["water_budget_residual"]) < 1e-15
    assert np.all(data["capacity"] >= 0)
    assert np.isfinite(data["matrix_condition_2"])
    np.testing.assert_allclose(capacity_map(y, .002, m),
                               frozen_coefficient_capacity(y[2:], hydraulic_factor(y[1], H), C))
    step = 1e-6
    fd = np.column_stack([(capacity_map(y+step*e, .002, m)-capacity_map(y-step*e, .002, m))/(2*step)
                         for e in np.eye(6)])
    np.testing.assert_allclose(capacity_jacobian(y, .002, m), fd, rtol=1e-8, atol=1e-8)
    np.testing.assert_allclose(data["capacity_rate"], fd@m.rhs(y, .002), rtol=1e-8, atol=1e-8)
    if strengths == VARIANTS[0]:
        np.testing.assert_allclose(data["rhs"], combined_rhs(y, 0., .002/.101325, H, C), atol=1e-14)
    assert data["equilibrium_reference_status"] == "not_supplied"
    assert data["equilibrium_error_norm_scaled"] is None


@pytest.mark.parametrize("strengths", VARIANTS)
def test_frozen_roots_conditional_vs_full_stability(strengths):
    m = model(strengths)
    v = .00244 if strengths == VARIANTS[0] else .0013
    roots = frozen_equilibria(v, (-2.222, -.001), m, reference_scales(m))
    assert len(roots["roots"]) == (2 if strengths[1] == 0 else 1)
    for root in roots["roots"]:
        data = root["diagnostics"]
        assert data["potential_norm_scaled"] < 1e-11
        assert data["full_residual_norm_scaled_per_day"] < 1e-10
        assert max(data["local_carbon_eigenvalues"].real) < 0
        dk = data["capacity_jacobian"]
        np.testing.assert_allclose(m.jacobian(root["state"], v)[2:, 2:],
                                   data["matrix"]@(np.eye(4)-dk[:, 2:]), atol=1e-15)
    assert roots["formal_zero_full_roots"] == "not_defined_by_positive_leaf_model"
    if strengths[1] == 0:
        assert {r["stability"] for r in roots["roots"]} == {"stable", "unstable"}
        assert {r["diagnostics"]["equilibrium_unstable_dimensions"] for r in roots["roots"]} == {0, 1}
        assert {r["diagnostics"]["equilibrium_type"] for r in roots["roots"]} == {
            "hyperbolic_sink", "hyperbolic_saddle",
        }
    if strengths == VARIANTS[-1]:
        high = frozen_equilibria(.0028, (-2.222, -.001), m, reference_scales(m))
        assert len(high["roots"]) == 1
        assert high["roots"][0]["stability"] == "unstable"
        assert high["roots"][0]["diagnostics"]["equilibrium_unstable_dimensions"] == 2
        assert max(high["roots"][0]["diagnostics"]["local_carbon_eigenvalues"].real) < 0


@pytest.mark.parametrize("leaf,count,zero_status", [
    (-1., 2, "unstable"), (-3., 1, "attracting"), (-4.2, 1, "attracting"),
])
def test_complete_conditional_carbon_roots(leaf, count, zero_status):
    roots = conditional_carbon_roots(leaf, M)
    assert len(roots) == count
    assert roots[0].kind == "formal_zero"
    assert not roots[0].in_positive_leaf_domain
    assert roots[0].stability == zero_status
    for root in roots:
        np.testing.assert_allclose(carbon_rhs(root.state, hydraulic_factor(leaf, H), C), 0., atol=1e-15)


def test_nonhyperbolic_threshold_and_positive_root_below_leaf_boundary():
    c = replace(C, assimilation_max=C.losses[0]/(C.allocation[0]*C.canopy_coefficient))
    m = replace(M, carbon=c)
    roots = conditional_carbon_roots(0., m)
    assert len(roots) == 1
    assert roots[0].stability == "attracting_nonhyperbolic_nonnegative_side"
    m = replace(M, feedback=replace(F, minimum_leaf=.2))
    roots = conditional_carbon_roots(-2., m)
    assert roots[1].state[0] < .2
    assert not roots[1].in_positive_leaf_domain


def test_zero_potential_is_not_full_equilibrium_and_sum_cancels():
    m = model(VARIANTS[-1])
    stocks = positive_carbon_equilibrium(hydraulic_factor(-1.2, H), C)
    y = np.r_[-.8, -1.2, stocks]
    data = diagnose_capacity(y, .002, m, reference_scales(m))
    assert data["potential_norm_scaled"] < 1e-13
    assert data["water_residual_norm_scaled_per_day"] > .1
    y[3] += .05
    y[4] -= .05
    data = diagnose_capacity(y, .002, m, reference_scales(m))
    assert abs(data["potential_sum"]) < 1e-12
    assert data["potential_norm_scaled"] > .01
    assert data["carbon_residual_norm_scaled_per_day"] > 1e-6
    with pytest.raises(ValueError, match="full frozen"):
        diagnose_capacity(y, .002, m, reference_scales(m), equilibrium=y)


def test_fold_and_absent_roots_are_bounded_not_zero_attractors():
    folds = M.stationary_points((-2.222, -.001))
    assert len(folds) == 1
    y, v, _ = M.equilibrium_curve(folds[0])
    record = frozen_equilibria(v, (-2.222, -.001), M, reference_scales(M))
    assert len(record["roots"]) == 1
    assert record["roots"][0]["stability"] == "nonhyperbolic_within_tolerance"
    assert record["roots"][0]["diagnostics"]["equilibrium_center_dimensions"] == 1
    assert not frozen_equilibria(.0028, (-2.222, -.001), M, reference_scales(M))["roots"]


def test_scalar_potential_growth_reuses_existing_counterexample():
    m = CubicCapacityFeedback()
    assert m.capacity_derivative(2.8, 0.) > 1
    assert m.potential_energy_derivative(2.8, 0.) > 0
    step = 1e-5
    fd = (m.potential(2.8+step, 0.)-m.potential(2.8-step, 0.))/(2*step)
    assert m.potential_energy_derivative(2.8, 0.) == pytest.approx(
        m.potential(2.8, 0.)*fd*m.rhs(2.8, 0.), rel=1e-8,
    )


@pytest.mark.parametrize("strengths", VARIANTS)
def test_diagnostics_preserve_positive_short_trajectories(strengths):
    m = model(strengths)
    initial = np.r_[LEAF/2, LEAF, REF*np.array([.9, .95, .9, .95])]
    solution = m.simulate(initial, np.linspace(0., 100., 21), lambda t: .0013)
    assert solution["status"] == "completed"
    assert np.all(solution["state"][:, 2:] > 0)
    assert np.all(solution["state"][:, 2] > F.minimum_leaf)
    for y in solution["state"][::5]:
        data = diagnose_capacity(y, .0013, m, reference_scales(m))
        assert data["identity_residual_max"] < 1e-15
        assert abs(data["carbon_budget_residual"]) < 1e-15


def test_fixed_scaling_and_directional_potential_derivative():
    scales = reference_scales(M)
    np.testing.assert_allclose(scales.state[2:], REF)
    y = np.array([-.8, -1.2, .2, 2., .4, 5.])
    data = diagnose_capacity(y, .002, M, scales)
    step = 1e-5
    direction = data["rhs"]
    def energy(state):
        p = capacity_map(state, .002, M)-state[2:]
        return np.dot(p/np.array(scales.state[2:]), p/np.array(scales.state[2:]))/2
    fd = (energy(y+step*direction)-energy(y-step*direction))/(2*step)
    assert data["potential_energy_rate_scaled_per_day"] == pytest.approx(fd, rel=1e-7)


@pytest.mark.parametrize("bad", [(-1., 1., 1., 1., 1., 1.), (1.,)*5, (float("nan"),)*6])
def test_invalid_scales(bad):
    with pytest.raises(ValueError, match="scales"):
        CapacityScales(bad)


def test_invalid_states_forcing_spectra_and_kinks_fail_explicitly():
    y = np.array([-.8, -1.2, .2, 2., .4, 5.])
    with pytest.raises(ValueError):
        capacity_map(y, -.001, M)
    with pytest.raises(ValueError):
        diagnose_capacity(np.r_[y[:2], 0., y[3:]], .002, M, reference_scales(M))
    with pytest.raises(ValueError):
        conditional_carbon_roots(.1, M)
    y[1] = -3.
    with pytest.raises(ValueError):
        capacity_jacobian(y, .002, M)
    with pytest.raises(ValueError):
        spectrum_status([np.nan])
    with pytest.raises(ValueError):
        spectrum_status([-1], tolerance=0)
