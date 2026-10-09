from dataclasses import replace

import numpy as np
import pytest
from scipy.optimize import brentq

from control_carbon.manzoni_hydraulics import HydraulicParameters, hydraulic_fold, transpiring_equilibria
from control_carbon.spac_carbon import (
    YEAR_DAYS, CarbonParameters, carbon_budget, carbon_matrix, combined_rhs,
    combined_jacobian, hydraulic_factor, positive_carbon_equilibrium,
    frozen_coefficient_capacity, simulate,
)
from control_carbon.spac_carbon_two_way import (
    StructuralFeedback, feedback_factors, effective_conductance, structural_fluxes,
    two_way_rhs, two_way_jacobian, equilibrium_curve, stationary_points,
    equilibria_at_vpd, continue_curve, simulate_two_way,
)


H = HydraulicParameters(.0086, .47, -1.5, -2., .001, .001)
C = CarbonParameters(1.2/YEAR_DAYS, (.3, .4, .3), tuple(k/YEAR_DAYS for k in [.8, .05, .2, .04]),
                     (.7, .9, .8), 10., .5, .25)
LEAF = transpiring_equilibria(0., .01, H)[0]
STOCKS = positive_carbon_equilibrium(hydraulic_factor(LEAF, H), C)
REF = np.r_[LEAF/2, LEAF, STOCKS]
F = StructuralFeedback(STOCKS[0], STOCKS[1], .01, 1e-6)
BOUNDS = (-2.222, -.001)


def variant(supply, demand, wood=0., **kwargs):
    return replace(F, leaf_supply=supply, leaf_demand=demand, wood_supply=wood, **kwargs)


def test_no_feedback_matches_stage_a_rhs_jacobian_capacity_and_trajectory():
    y = np.array([-.8, -1.2, .2, 2., .4, 5.])
    np.testing.assert_allclose(two_way_rhs(y, -.5, .01, H, C, F),
                               combined_rhs(y, -.5, .01, H, C), atol=1e-14)
    np.testing.assert_allclose(two_way_jacobian(y, -.5, .01, H, C, F),
                               combined_jacobian(y, -.5, .01, H, C), atol=1e-14)
    np.testing.assert_allclose(two_way_rhs(REF, 0., .01, H, C, F), 0., atol=1e-12)
    times = np.linspace(0, 100, 31)
    a = simulate(y, times, lambda t: (-.5, .01, 1., 1.), H, C)
    b = simulate_two_way(y, times, lambda t: .01, -.5, H, C, F)
    assert b["status"] == "completed"
    np.testing.assert_allclose(a["state"], b["state"], atol=1e-9)


@pytest.mark.parametrize("feedback", [F, variant(1, 0), variant(0, 1), variant(1, 1), variant(1, 1, 1, wood_exponent=.7)])
def test_full_jacobian_against_central_differences_and_budgets(feedback):
    y = np.array([-.8, -1.2, .2, 2., .4, 5.])
    jac = two_way_jacobian(y, -.5, .01, H, C, feedback)
    step = 1e-6
    fd = np.column_stack([
        (two_way_rhs(y+step*e, -.5, .01, H, C, feedback)
         -two_way_rhs(y-step*e, -.5, .01, H, C, feedback))/(2*step)
        for e in np.eye(6)
    ])
    np.testing.assert_allclose(jac, fd, rtol=1e-8, atol=1e-8)
    if feedback != F:
        assert np.any(jac[:2, 2:] != 0)
    erx, exl, e = structural_fluxes(y, -.5, .01, H, feedback)
    f = two_way_rhs(y, -.5, .01, H, C, feedback)
    assert H.cx*f[0]+H.cl*f[1] == pytest.approx(erx-e, abs=1e-15)
    h = hydraulic_factor(y[1], H)
    assert sum(f[2:]) == pytest.approx(carbon_budget(y[2:], h, C)["nep"], abs=1e-15)
    capacity = frozen_coefficient_capacity(y[2:], h, C)
    np.testing.assert_allclose(f[2:], carbon_matrix(C)@(y[2:]-capacity), atol=1e-15)
    np.testing.assert_allclose(f[2:], combined_rhs(y, -.5, .01, H, C)[2:], atol=1e-15)
    assert np.linalg.norm(capacity-y[2:]) > .1


@pytest.mark.parametrize("feedback", [F, variant(1, 0), variant(0, 1), variant(1, 1), variant(1, 1, 1)])
def test_series_equilibrium_reduction_and_derivative(feedback):
    leaf = -1.5
    y, d, slope = equilibrium_curve(leaf, 0., H, C, feedback)
    np.testing.assert_allclose(two_way_rhs(y, 0., d, H, C, feedback), 0., atol=1e-12)
    fluxes = structural_fluxes(y, 0., d, H, feedback)
    np.testing.assert_allclose(fluxes, fluxes[0], atol=1e-15)
    np.testing.assert_allclose(frozen_coefficient_capacity(y[2:], hydraulic_factor(leaf, H), C), y[2:], atol=1e-12)
    step = 1e-5
    fd = (equilibrium_curve(leaf+step, 0., H, C, feedback)[1]
          -equilibrium_curve(leaf-step, 0., H, C, feedback)[1])/(2*step)
    assert slope == pytest.approx(fd, rel=1e-8, abs=1e-10)
    if feedback.wood_supply == 0:
        fl = feedback_factors(y[2], y[3], feedback)[0]
        gx = H.gx_max*(1+leaf/3)
        assert fluxes[0]/(-leaf) == pytest.approx(2*gx*fl/(1+fl))
        if feedback.leaf_supply == 1:
            assert not np.isclose(fl, 1)
            assert abs(fl*gx-fluxes[0]/(-leaf)) > 1e-4


def test_reference_is_fixed_and_all_variants_retain_it():
    for f in [F, variant(1, 0), variant(0, 1), variant(1, 1), variant(1, 1, 1)]:
        fl, fe, fw, *_ = feedback_factors(*STOCKS[:2], f)
        np.testing.assert_allclose([fl, fe, fw], 1., atol=1e-15)
        np.testing.assert_allclose(two_way_rhs(REF, 0., .01, H, C, f), 0., atol=1e-12)


def test_stage_a_fold_schur_and_root_counts():
    folds = stationary_points(BOUNDS, 0., H, C, F)
    assert len(folds) == 1
    expected_d, expected_leaf = hydraulic_fold(0., H)
    assert folds[0] == pytest.approx(expected_leaf, abs=1e-10)
    y, d, _ = equilibrium_curve(folds[0], 0., H, C, F)
    assert d == pytest.approx(expected_d, abs=1e-12)
    jac = two_way_jacobian(y, 0., d, H, C, F)
    schur = jac[:2, :2]-jac[:2, 2:]@np.linalg.solve(jac[2:, 2:], jac[2:, :2])
    assert min(abs(np.linalg.eigvals(schur))) < 1e-10
    eigen = np.linalg.eigvals(jac)
    assert np.count_nonzero(abs(eigen) < 1e-9) == 1
    for vpd in [.01, .0242]:
        equilibria = equilibria_at_vpd(vpd, BOUNDS, folds, 0., H, C, F)
        expected = [leaf for leaf in transpiring_equilibria(0., vpd, H)
                    if BOUNDS[0] <= leaf <= BOUNDS[1]]
        np.testing.assert_allclose([y[1] for y in equilibria], expected, atol=1e-10)


def test_pseudo_arclength_crosses_fold_and_refines():
    curves = []
    for step in [.02, .01]:
        records, status = continue_curve(-.1, step, 1000, (-2.21, -.001), 0., H, C, F)
        assert status == "domain_boundary"
        states = np.array([r[0] for r in records])
        ds = np.array([r[1] for r in records])
        assert states[-1, 1] < -2.1
        assert np.any(np.diff(ds) < 0)  # continuation passes the D fold
        for y, d in records[::15]:
            np.testing.assert_allclose(two_way_rhs(y, 0., d, H, C, F), 0., atol=1e-11)
        curves.append(ds.max())
    assert abs(curves[0]-curves[1]) < 2e-6


def test_formal_closed_transport_is_not_an_isolated_sink():
    y = np.r_[-1., -5., STOCKS]
    f = variant(1, 1)
    assert two_way_rhs(y, 0., .03, H, C, f)[0] == 0
    assert two_way_rhs(y, 0., .03, H, C, f)[1] == 0
    np.testing.assert_array_equal(two_way_jacobian(y, 0., .03, H, C, f)[:2], 0)


def test_vegetation_boundary_event_not_mortality():
    f = variant(1, 0, minimum_leaf=.25)
    result = simulate_two_way(REF, [0., 2*YEAR_DAYS], lambda t: .05, 0., H, C, f)
    assert result["status"] == "vegetation_domain_exit"
    assert result["state"][-1, 2] == pytest.approx(.25, abs=1e-10)
    assert result["time_days"][-1] < 2*YEAR_DAYS
    with pytest.raises(ValueError):
        simulate_two_way(np.r_[REF[:2], np.zeros(4)], [0, 10], lambda t: .01, 0., H, C, F)


def test_supply_only_fold_and_demand_combination_full_stability_not_curve_slope():
    e1 = variant(1, 0)
    folds = stationary_points(BOUNDS, 0., H, C, e1)
    assert len(folds) == 1
    y, d, _ = equilibrium_curve(folds[0], 0., H, C, e1)
    assert d < hydraulic_fold(0., H)[0]
    jac = two_way_jacobian(y, 0., d, H, C, e1)
    schur = jac[:2, :2]-jac[:2, 2:]@np.linalg.solve(jac[2:, 2:], jac[2:, :2])
    assert min(abs(np.linalg.eigvals(schur))) < 1e-10
    e3 = variant(1, 1)
    assert stationary_points(BOUNDS, 0., H, C, e3) == ()
    def alpha(leaf):
        y, d, _ = equilibrium_curve(leaf, 0., H, C, e3)
        return max(np.linalg.eigvals(two_way_jacobian(y, 0., d, H, C, e3)).real)
    leaf = brentq(alpha, -2.1, -1.9, xtol=1e-12)
    y, d, slope = equilibrium_curve(leaf, 0., H, C, e3)
    eig = np.linalg.eigvals(two_way_jacobian(y, 0., d, H, C, e3))
    assert slope != pytest.approx(0., abs=1e-4)
    assert abs(eig[np.argmax(eig.real)].imag) > .02
    assert max(eig.real) == pytest.approx(0., abs=1e-10)


def test_rate_comparison_on_persistent_branch_without_a_mortality_claim():
    e3 = variant(1, 1)
    for d in np.linspace(.01, .03, 11):
        roots = equilibria_at_vpd(d, BOUNDS, (), 0., H, C, e3)
        assert len(roots) == 1
        assert max(np.linalg.eigvals(two_way_jacobian(roots[0], 0., d, H, C, e3)).real) < 0
    final = equilibria_at_vpd(.03, BOUNDS, (), 0., H, C, e3)[0]
    for duration in [10., 10*YEAR_DAYS]:
        def forcing(t):
            s = min(t/duration, 1)
            return .01+.02*(3*s*s-2*s*s*s)
        result = simulate_two_way(REF, [0., duration+500*YEAR_DAYS], forcing, 0., H, C, e3)
        if duration == 10.:
            assert result["status"] == "vegetation_domain_exit"
        else:
            assert result["status"] == "completed"
            np.testing.assert_allclose(result["state"][-1], final, rtol=1e-6)


def test_e2_rate_independent_invariant_box_vector_field_faces():
    e2 = variant(0, 1)
    leaf_barrier, x_barrier = -1.8, -.9
    leaf_lower = positive_carbon_equilibrium(hydraulic_factor(leaf_barrier, H), C)[0]
    leaf_upper = positive_carbon_equilibrium(1., C)[0]
    dmax = H.gx_max*(1+leaf_barrier/3)*(-leaf_barrier)/(H.gs_max*(1+leaf_barrier/4)*(leaf_upper/F.leaf_reference))
    assert dmax > .018
    for x in [x_barrier, 0.]:
        for psi in [leaf_barrier, 0.]:
            for leaf in [leaf_lower, leaf_upper]:
                for d in [0., dmax]:
                    f = two_way_rhs(np.r_[x, psi, leaf, STOCKS[1:]], 0., d, H, C, e2)
                    assert f[0]*(1 if x == x_barrier else -1) >= -1e-12
                    assert f[1]*(1 if psi == leaf_barrier else -1) >= -1e-12
                    assert f[2]*(1 if leaf == leaf_lower else -1) >= -1e-12


@pytest.mark.parametrize("change", [
    {"leaf_reference": 0}, {"wood_reference": np.nan}, {"leaf_floor": -1},
    {"minimum_leaf": 1}, {"leaf_supply": -1}, {"leaf_demand": 2},
    {"wood_exponent": 0},
])
def test_invalid_feedback(change):
    with pytest.raises(ValueError):
        replace(F, **change)


def test_invalid_states_curves_and_series_conductance():
    with pytest.raises(ValueError):
        two_way_rhs(np.r_[REF[:2], np.zeros(4)], 0., .01, H, C, F)
    with pytest.raises(ValueError):
        equilibrium_curve(-2.3, 0., H, C, F)
    with pytest.raises(ValueError):
        two_way_jacobian(np.r_[REF[0], -3., STOCKS], 0., .01, H, C, F)
    with pytest.raises(ValueError):
        effective_conductance(-1., 1.)
    assert effective_conductance(0., 0.) == 0.
    with pytest.raises(ValueError):
        feedback_factors(.1, 0., variant(1, 1, 1, wood_exponent=.5))
