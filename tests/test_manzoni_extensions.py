from dataclasses import replace

import numpy as np
from numpy.polynomial import polynomial as poly
import pytest
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from control_carbon.manzoni_hydraulics import (
    HydraulicParameters, conductances, hydraulic_fold, hydraulic_rhs,
    transpiring_equilibria, water_fluxes,
)
from control_carbon.manzoni_extensions import (
    SigmoidShape, SoilRootParameters, _exponential_roots, deteriorated_traits,
    root_node, sigmoid_conductances, sigmoid_equilibria, sigmoid_jacobian,
    sigmoid_equilibrium_jacobian_scaled, sigmoid_log_flux_ratio,
    sigmoid_rhs, sigmoid_steady_fluxes, soil_equilibria, soil_equilibrium_state,
    soil_fold, soil_hydraulic_jacobian, soil_hydraulic_rhs,
    soil_root_conductance, soil_steady_fluxes,
)


P = HydraulicParameters(.0086, .47, -3., -2., .001, .001)
SOIL = SoilRootParameters(.1, -.1)


def test_sigmoid_hand_values_budget_and_stable_extreme_evaluation():
    shape = SigmoidShape(4., 2.)
    assert sigmoid_conductances(-3., P, shape)[0] == pytest.approx(P.gx_max/2)
    assert sigmoid_conductances(-2., P, shape)[1] == pytest.approx(P.gs_max/2)
    assert sigmoid_conductances(0., P, shape) == (P.gx_max, P.gs_max)
    assert np.all(np.isfinite(sigmoid_conductances(-1e100, P, shape)))
    state = np.array([-1., -2.])
    gx, gs = sigmoid_conductances(state[1], P, shape)
    rhs = sigmoid_rhs(state, -.5, .02, P, shape)
    assert P.cx*rhs[0]+P.cl*rhs[1] == pytest.approx(2*gx*(-.5-state[0])-gs*.02)


@pytest.mark.parametrize("shape", [SigmoidShape(4., 2.), SigmoidShape(1., 1.), SigmoidShape(.5, .8)])
def test_sigmoid_full_jacobian_finite_difference(shape):
    state = np.array([-1., -2.3])
    h = 1e-6
    fd = np.column_stack([
        (sigmoid_rhs(state+h*e, -.5, .02, P, shape)
         - sigmoid_rhs(state-h*e, -.5, .02, P, shape))/(2*h)
        for e in np.eye(2)
    ])
    np.testing.assert_allclose(sigmoid_jacobian(state, -.5, .02, P, shape), fd, rtol=1e-8, atol=1e-8)


def test_exponential_root_isolation_tangencies_cancellation_and_far_tail():
    assert _exponential_roots({0.: 4., 1.: -4., 2.: 1.}) == pytest.approx((np.log(2.),))
    assert _exponential_roots({0.: -6., 1.: 11., 2.: -6., 3.: 1.}) == pytest.approx(
        (0., np.log(2.), np.log(3.)), abs=1e-11,
    )
    assert _exponential_roots({0.: -1., .01: 1e-10}) == pytest.approx((1000*np.log(10),))
    assert _exponential_roots({0.: 1., 2.: 1.}) == ()
    assert _exponential_roots({0.: 0., 1.: 2.}) == ()


@pytest.mark.parametrize("root,d,x,s,ax,ass,count", [
    (-.5, .01, -3., -2., 1, 1, 1),
    (-.5, .01, -3., -2., 3, 1, 2),
    (-.5, .06, -3., -2., 3, 1, 0),
    (-.5, .01, -1.5, -2., 6, 6, 3),
    (0., .02, -3., -2., 4, 12, 1),
    (0., .02, -3., -2., 12, 4, 2),
])
def test_sigmoid_count_and_potentials_against_independent_integer_polynomial(root, d, x, s, ax, ass, count):
    p = replace(P, psi50_x=x, psi50_s=s)
    leaves = sigmoid_equilibria(root, d, p, SigmoidShape(ax, ass))
    coefficients = np.zeros(max(ax, ass+1)+1)
    coefficients[0] = root-p.gs_max*d/p.gx_max
    coefficients[1] += 1.
    coefficients[ass] += root/(-s)**ass
    coefficients[ass+1] += 1./(-s)**ass
    coefficients[ax] -= p.gs_max*d/p.gx_max/(-x)**ax
    polynomial_roots = poly.polyroots(coefficients)
    expected = sorted(
        [-z.real for z in polynomial_roots if abs(z.imag) < 1e-7 and z.real > -root],
        reverse=True,
    )
    assert len(leaves) == count == len(expected)
    assert leaves == pytest.approx(expected, rel=1e-9, abs=1e-9)
    for i, leaf in enumerate(leaves):
        state = np.array([(root+leaf)/2, leaf])
        np.testing.assert_allclose(sigmoid_rhs(state, root, d, p, SigmoidShape(ax, ass)), 0., atol=1e-10)
        eig = np.linalg.eigvals(sigmoid_jacobian(state, root, d, p, SigmoidShape(ax, ass)))
        assert np.max(eig.real) < 0 if i % 2 == 0 else np.prod(eig.real) < 0
        scaled = sigmoid_equilibrium_jacobian_scaled(leaf, root, d, p, SigmoidShape(ax, ass))
        gx, gs = sigmoid_conductances(leaf, p, SigmoidShape(ax, ass))
        np.testing.assert_allclose(gx*scaled, sigmoid_jacobian(state, root, d, p, SigmoidShape(ax, ass)), atol=1e-10)


def test_noninteger_far_equilibrium_and_equal_tail_exponents():
    shape = SigmoidShape(4.01, 3.)
    leaves = sigmoid_equilibria(0., .01, P, shape)
    assert len(leaves) == 2
    assert abs(leaves[1]) > 1e40
    leaf = leaves[1]
    assert sigmoid_log_flux_ratio(leaf, 0., .01, P, shape) == pytest.approx(0., abs=1e-10)
    eig = np.linalg.eigvals(sigmoid_equilibrium_jacobian_scaled(leaf, 0., .01, P, shape))
    assert np.prod(eig.real) < 0
    for shape in [SigmoidShape(4., 3.), SigmoidShape(4., 4.)]:
        assert len(sigmoid_equilibria(0., .01, P, shape)) == 1
    assert sigmoid_equilibria(-.5, 0., P, shape) == (-.5,)


def test_sigmoid_double_equilibrium_at_independent_fold():
    shape = SigmoidShape(4., 2.)
    def derivative(leaf):
        gx, gs = sigmoid_conductances(leaf, P, shape)
        state = np.array([leaf/2, leaf])
        d = gx*(-leaf)/gs
        j = sigmoid_jacobian(state, 0., d, P, shape)
        return np.linalg.det(j)
    leaf = brentq(derivative, -7., -2.)
    gx, gs = sigmoid_conductances(leaf, P, shape)
    d = gx*(-leaf)/gs
    leaves = sigmoid_equilibria(0., d, P, shape)
    assert leaves == pytest.approx((leaf,), abs=1e-8)
    assert len(sigmoid_equilibria(0., .999*d, P, shape)) == 2
    assert sigmoid_equilibria(0., 1.001*d, P, shape) == ()


def test_soil_root_hand_continuity_and_budget():
    p = replace(P, gx_max=1., gs_max=2., psi50_x=-2., psi50_s=-2., cx=3., cl=4.)
    soil = SoilRootParameters(2., -1.)
    state = np.array([-2., -3.])
    assert soil_root_conductance(-2., soil) == .5
    root = root_node(state, -1., p, soil)
    assert root == pytest.approx(-1.2)
    erx, exl, e = water_fluxes(state, root, .1, p)
    assert (erx, exl, e) == pytest.approx((.4, .5, .05))
    assert soil_root_conductance(-1., soil)*(-1.-root) == pytest.approx(erx)
    rhs = soil_hydraulic_rhs(state, -1., .1, p, soil)
    assert p.cx*rhs[0]+p.cl*rhs[1] == pytest.approx(erx-e)


def test_soil_jacobian_with_both_root_chain_rule_terms():
    state = np.array([-1., -2.3])
    h = 1e-6
    fd = np.column_stack([
        (soil_hydraulic_rhs(state+h*e, -.5, .02, P, SOIL)
         - soil_hydraulic_rhs(state-h*e, -.5, .02, P, SOIL))/(2*h)
        for e in np.eye(2)
    ])
    np.testing.assert_allclose(soil_hydraulic_jacobian(state, -.5, .02, P, SOIL), fd, rtol=1e-8, atol=1e-8)


@pytest.mark.parametrize("gsr", [.01, .1, 1.])
def test_soil_fold_tangency_equilibria_full_jacobian_and_deterioration(gsr):
    damaged = deteriorated_traits(P, .6)
    soil = SoilRootParameters(gsr, -.1)
    root = -.5
    d, leaf = soil_fold(root, damaged, soil)
    assert soil_equilibria(root, d, damaged, soil) == pytest.approx((leaf,), abs=1e-9)
    state = soil_equilibrium_state(leaf, root, damaged, soil)
    np.testing.assert_allclose(soil_hydraulic_rhs(state, root, d, damaged, soil), 0., atol=1e-12)
    eig = np.linalg.eigvals(soil_hydraulic_jacobian(state, root, d, damaged, soil))
    assert min(abs(eig)) < 1e-9
    assert min(eig.real) < 0
    for factor, count in [(.999, 2), (1.001, 0)]:
        leaves = soil_equilibria(root, factor*d, damaged, soil)
        assert len(leaves) == count
        for leaf in leaves:
            supply, demand = soil_steady_fluxes(leaf, root, factor*d, damaged, soil)
            assert supply == pytest.approx(demand, abs=1e-14)
    assert soil_fold(root, P, soil) is None
    assert soil_fold(-4., damaged, soil) is None
    assert soil_fold(root, damaged, soil)[0] < hydraulic_fold(root, damaged)[0]


def test_soil_large_conductance_limit_and_fixed_damage_mapping():
    soil = SoilRootParameters(1e12, -.1)
    damaged = deteriorated_traits(P, .6)
    assert damaged.gx_max == pytest.approx(.6*P.gx_max)
    assert damaged.psi50_x == pytest.approx(.6*P.psi50_x)
    assert damaged.gs_max == P.gs_max and damaged.psi50_s == P.psi50_s
    assert deteriorated_traits(P, 1.) == P
    assert soil_equilibria(-.5, .01, damaged, soil) == pytest.approx(
        transpiring_equilibria(-.5, .01, damaged), rel=1e-9,
    )
    state = [-.5, -1.]
    np.testing.assert_allclose(soil_hydraulic_rhs(state, -.1, .01, P, soil),
                               hydraulic_rhs(state, -.1, .01, P), rtol=1e-12)


def test_sigmoid_three_branch_destinations_and_refinement():
    p = replace(P, psi50_x=-1.5)
    shape = SigmoidShape(6., 6.)
    leaves = sigmoid_equilibria(-.5, .01, p, shape)
    for initial, target in [([-.7, -1.], leaves[0]), ([-1.8, -3.5], leaves[-1])]:
        endings = []
        for rtol in [1e-9, 1e-11]:
            sol = solve_ivp(lambda t, y: sigmoid_rhs(y, -.5, .01, p, shape),
                            (0., 500.), initial, rtol=rtol, atol=rtol*.01, max_step=.2)
            assert sol.success, sol.message
            endings.append(sol.y[:, -1])
        np.testing.assert_allclose(*endings, atol=1e-8, rtol=0.)
        assert endings[-1][1] == pytest.approx(target, abs=1e-7)


@pytest.mark.parametrize("factor", [.7, 1.3])
def test_soil_full_ode_destinations_and_refinement(factor):
    p = deteriorated_traits(P, .6)
    d = factor*soil_fold(-.5, p, SOIL)[0]
    leaves = soil_equilibria(-.5, d, p, SOIL)
    target = leaves[0] if leaves else 2*p.psi50_s
    endings = []
    for rtol in [1e-10, 1e-12]:
        sol = solve_ivp(
            lambda t, y: soil_hydraulic_rhs(y, -.5, d, p, SOIL),
            (0., 200.), [-.6, -.8], rtol=rtol, atol=rtol*.01, max_step=.1,
        )
        assert sol.success, sol.message
        assert np.all(sol.y <= 0)
        endings.append(sol.y[:, -1])
    np.testing.assert_allclose(*endings, atol=1e-8, rtol=0.)
    assert endings[-1][1] == pytest.approx(target, abs=1e-7)


@pytest.mark.parametrize("factory,args", [
    (SigmoidShape, (0., 1.)), (SigmoidShape, (1., np.nan)),
    (SoilRootParameters, (0., -.1)), (SoilRootParameters, (.1, 0.)),
    (SoilRootParameters, (.1, np.inf)),
])
def test_invalid_extension_parameters(factory, args):
    with pytest.raises(ValueError):
        factory(*args)


def test_invalid_extension_inputs():
    with pytest.raises(ValueError):
        soil_root_conductance(0., SOIL)
    with pytest.raises(ValueError):
        deteriorated_traits(P, 0.)
    with pytest.raises(ValueError):
        deteriorated_traits(P, 1.1)
    with pytest.raises(ValueError):
        sigmoid_equilibria(-.5, -.01, P, SigmoidShape(4., 4.))
    with pytest.raises(ValueError):
        sigmoid_jacobian([0., 0.], 0., .01, P, SigmoidShape(.5, .5))
