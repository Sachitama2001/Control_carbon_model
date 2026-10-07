from dataclasses import replace

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

from control_carbon.manzoni_hydraulics import (
    HydraulicParameters, conductances, hydraulic_fold, hydraulic_jacobian,
    hydraulic_rhs, steady_fluxes, transpiring_equilibria, vpd_mole_fraction,
    water_fluxes,
)


def parameters(**kwargs):
    return replace(
        HydraulicParameters(.0086, .47, -1.5, -2., .001, .001), **kwargs,
    )


def test_conductance_algebra_and_units():
    p = parameters()
    assert conductances(0., p) == pytest.approx((.0086, .47))
    assert conductances(-1.5, p)[0] == pytest.approx(.0086 / 2)
    assert conductances(-2., p)[1] == pytest.approx(.47 / 2)
    assert conductances(-5., p) == (0., 0.)
    assert vpd_mole_fraction(2.0265, 101.325) == pytest.approx(.02)


def test_hand_calculated_fluxes_and_water_budget():
    p = parameters(gx_max=1., gs_max=2., cx=3., cl=4.)
    state = np.array([-1., -2.])
    erx, exl, evaporation = water_fluxes(state, 0., .1, p)
    assert (erx, exl, evaporation) == pytest.approx((2/3, 2/3, .1))
    rhs = hydraulic_rhs(state, 0., .1, p)
    assert rhs == pytest.approx((0., (2/3 - .1)/4))
    assert p.cx * rhs[0] + p.cl * rhs[1] == pytest.approx(erx - evaporation)


@pytest.mark.parametrize("root", [0., -1.])
def test_active_roots_against_independent_bracketed_flux_balance(root):
    p = parameters()
    fold_d, fold_leaf = hydraulic_fold(root, p)
    d = .7 * fold_d
    leaves = transpiring_equilibria(root, d, p)
    def balance(leaf):
        supply, demand = steady_fluxes(leaf, root, d, p)
        return supply - demand
    expected = (
        brentq(balance, fold_leaf, root),
        brentq(balance, 2*p.psi50_x, fold_leaf),
    )
    assert leaves == pytest.approx(expected, abs=1e-12)
    for i, leaf in enumerate(leaves):
        state = np.array([(root + leaf)/2, leaf])
        np.testing.assert_allclose(hydraulic_rhs(state, root, d, p), 0., atol=1e-12)
        eig = np.linalg.eigvals(hydraulic_jacobian(state, root, d, p))
        assert np.max(eig.real) < 0 if i == 0 else np.prod(eig.real) < 0


@pytest.mark.parametrize("root", [0., -1., -2.9])
def test_fold_location_residual_tangency_and_full_eigenvalues(root):
    p = parameters()
    d, leaf = hydraulic_fold(root, p)
    assert transpiring_equilibria(root, d, p) == pytest.approx((leaf,), abs=1e-12)
    state = np.array([(root + leaf)/2, leaf])
    np.testing.assert_allclose(hydraulic_rhs(state, root, d, p), 0., atol=1e-12)
    eig = np.linalg.eigvals(hydraulic_jacobian(state, root, d, p))
    assert min(abs(eig)) < 1e-10
    assert min(eig.real) < 0
    h = 1e-6
    def balance(value):
        supply, demand = steady_fluxes(value, root, d, p)
        return supply - demand
    assert (balance(leaf+h) - balance(leaf-h))/(2*h) == pytest.approx(0., abs=1e-10)
    assert len(transpiring_equilibria(root, .999*d, p)) == 2
    assert transpiring_equilibria(root, 1.001*d, p) == ()
    def demand_limit(value):
        gx, gs = conductances(value, p)
        return -gx*(root-value)/gs
    optimum = minimize_scalar(
        demand_limit, bounds=(2*p.psi50_x, root), method="bounded",
        options={"xatol": 1e-12},
    )
    assert optimum.success
    assert -optimum.fun == pytest.approx(d, abs=1e-12)


def test_fig2_wet_root_fold_hand_value():
    d, leaf = hydraulic_fold(0., parameters())
    assert d == pytest.approx((.0086/.47)*(4/3))
    assert leaf == pytest.approx(-2.)


def test_safe_ordering_no_air_drying_fold_and_no_spurious_quadratic_roots():
    p = parameters(psi50_x=-2., psi50_s=-1.5)
    assert hydraulic_fold(0., p) is None
    for d in [0., .01, .04, .1, 1., 10.]:
        leaves = transpiring_equilibria(0., d, p)
        assert len(leaves) == 1
        assert leaves[0] > 2*p.psi50_s
    assert hydraulic_fold(-3., parameters()) is None
    assert transpiring_equilibria(-3.5, .01, parameters()) == ()


@pytest.mark.parametrize("leaf", [-1.2, -3.3, -4.5])
def test_full_jacobian_finite_difference_in_all_piecewise_regions(leaf):
    p = parameters()
    state = np.array([-.8, leaf])
    h = 1e-6
    fd = np.column_stack([
        (hydraulic_rhs(state+h*e, -.2, .01, p)
         - hydraulic_rhs(state-h*e, -.2, .01, p))/(2*h)
        for e in np.eye(2)
    ])
    np.testing.assert_allclose(hydraulic_jacobian(state, -.2, .01, p), fd, atol=1e-8)


def test_desiccated_full_system_is_a_continuum_not_a_hyperbolic_sink():
    p = parameters()
    for x in [-.2, -1., -2.]:
        for leaf in [-4., -5.]:
            np.testing.assert_array_equal(hydraulic_rhs([x, leaf], 0., .03, p), [0., 0.])
    np.testing.assert_array_equal(hydraulic_jacobian([-.5, -5.], 0., .03, p), np.zeros((2, 2)))
    with pytest.raises(ValueError, match="kink"):
        hydraulic_jacobian([-.5, -4.], 0., .03, p)


@pytest.mark.parametrize("d,initial,dry", [
    (.01, [-.2, -.4], False),
    (.01, [-1.7, -3.2], True),
    (.03, [-.2, -.4], True),
])
def test_full_ode_outcomes_and_tolerance_refinement(d, initial, dry):
    p = parameters()
    trajectories = []
    for rtol, atol in [(1e-10, 1e-12), (1e-12, 1e-14)]:
        sol = solve_ivp(
            lambda t, y: hydraulic_rhs(y, 0., d, p),
            (0., 100.), initial, rtol=rtol, atol=atol, max_step=.1,
        )
        assert sol.success, sol.message
        assert np.all(sol.y <= 0)
        trajectories.append(sol.y[:, -1])
    np.testing.assert_allclose(*trajectories, atol=2e-8, rtol=0.)
    target = -4. if dry else transpiring_equilibria(0., d, p)[0]
    assert trajectories[-1][1] == pytest.approx(target, abs=1e-7)


@pytest.mark.parametrize("kwargs", [
    {"cx": 0.}, {"cl": -1.}, {"psi50_x": 0.}, {"psi50_s": np.nan},
    {"gx_max": -1.}, {"gs_max": np.inf},
])
def test_invalid_parameters(kwargs):
    with pytest.raises(ValueError):
        parameters(**kwargs)


@pytest.mark.parametrize("root,d", [(1., .01), (0., -.1), (np.nan, .01), (0., np.inf)])
def test_invalid_forcing(root, d):
    with pytest.raises(ValueError):
        hydraulic_rhs([-.1, -.2], root, d, parameters())
    with pytest.raises(ValueError):
        transpiring_equilibria(root, d, parameters())


def test_invalid_state_and_pressure_conversion():
    with pytest.raises(ValueError):
        hydraulic_rhs([-.1], 0., .01, parameters())
    with pytest.raises(ValueError):
        hydraulic_rhs([-.1, np.nan], 0., .01, parameters())
    with pytest.raises(ValueError):
        vpd_mole_fraction(1., 0.)
    with pytest.raises(ValueError):
        vpd_mole_fraction(-1., 101.325)
