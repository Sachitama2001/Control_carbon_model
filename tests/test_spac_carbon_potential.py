from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm

from control_carbon.carbon_compensation import compensation_gramian
from control_carbon.spac_carbon_potential import (
    assimilation_rate_enclosure, certify_drying_step_potential,
    conducting_rhs_enclosure, endpoint_peak_lower, potential_diagnostics, potential_structure,
)
from control_carbon.spac_carbon_rational import RationalSPAC
from control_carbon.spac_carbon_stability import _Interval


@pytest.fixture
def model(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]/"examples"))
    from prepare_manuscript_phase_portraits import model_family
    return RationalSPAC(model_family()["E4"])


def test_potential_sign_and_both_balance_identities(model):
    state = model.equilibria(.00244)[0]
    state[1] -= .01
    diag = potential_diagnostics(model, state, .00254)
    structure = potential_structure(model)
    np.testing.assert_allclose(diag["capacity"], structure.capacity_direction_days*diag["assimilation"])
    np.testing.assert_allclose(diag["potential"], diag["capacity"]-state[2:])
    np.testing.assert_allclose(diag["rhs"][2:], -structure.matrix@diag["potential"], atol=1e-17)
    assert diag["potential_identity_residual"] < 1e-17
    assert diag["total_carbon_rate"] == pytest.approx(-np.ones(4)@structure.matrix@diag["potential"])
    np.testing.assert_allclose(structure.capacity_direction_days/365.25, [.375, 8, 1.5, 20.25])


@pytest.mark.parametrize("initial", [None, np.zeros(4), np.array([.2, -.1, .1, -.3])])
def test_endpoint_bound_uses_correct_initial_condition_and_l1_norm(model, initial):
    matrix = potential_structure(model).matrix
    p0 = np.zeros(4) if initial is None else initial
    rate = np.array([-.001, .002, -.003, .004])
    duration = 50.
    solve = solve_ivp(lambda t, p: matrix@p+rate, (0, duration), p0,
                      dense_output=True, rtol=1e-11, atol=1e-13)
    assert solve.success
    peak = np.abs(solve.sol(np.linspace(0, duration, 2001))).sum(axis=0).max()
    lower = endpoint_peak_lower(matrix, duration, rate*duration, initial)
    assert 0 < lower <= peak
    if initial is None:
        assert lower == pytest.approx(np.linalg.norm(rate*duration, 1)/(2+np.linalg.norm(matrix, 1)*duration))
    else:
        assert lower >= np.linalg.norm(initial, 1)


def test_endpoint_bound_does_not_mistake_closed_capacity_path_for_no_transient():
    matrix = -np.eye(4)
    assert endpoint_peak_lower(matrix, 1, np.zeros(4), np.zeros(4)) == 0
    with pytest.raises(ValueError, match="Zero duration"):
        endpoint_peak_lower(matrix, 0, np.ones(4))
    with pytest.raises(ValueError, match="Initial potential"):
        endpoint_peak_lower(matrix, 1, np.ones(4), [1])


@pytest.mark.parametrize("weight", [
    np.eye(3), np.diag([1., 4., .25]), np.array([[2., .2, .1], [.2, 1., .3], [.1, .3, 2.]])])
def test_gramian_and_action_match_hand_computable_integrators(weight):
    duration, terminal = 7., np.array([.3, -.4, .8])
    problem = compensation_gramian(np.zeros((3, 3)), np.eye(3), weight, duration)
    np.testing.assert_allclose(problem.gramian, duration*np.linalg.inv(weight), atol=1e-13)
    minimum = problem.minimum_action(terminal)
    assert minimum.action == pytest.approx(.5*terminal@weight@terminal/duration)
    for day in [0., 2., 7.]:
        np.testing.assert_allclose(minimum.control(day), -terminal/duration, atol=1e-15)
    integral = quad(lambda t: .5*minimum.control(t)@weight@minimum.control(t), 0, duration)[0]
    assert integral == pytest.approx(minimum.action)


def test_singular_gramian_checks_reachability_before_using_resolved_inverse():
    matrix, inputs = -np.eye(2), np.array([[1.], [0.]])
    problem = compensation_gramian(matrix, inputs, np.eye(1), 2.)
    assert problem.resolved_rank == 1
    minimum = problem.minimum_action([.5, 0])
    assert minimum.action == pytest.approx(.5*.5**2/((1-np.exp(-4))/2))
    with pytest.raises(ValueError, match="reachable subspace"):
        problem.minimum_action([0., .5])
    with pytest.raises(ValueError, match="reachable subspace"):
        problem.minimum_action([.5, 1e-12])
    zero = compensation_gramian(matrix, np.zeros((2, 1)), np.eye(1), 2)
    assert zero.minimum_action([0., 0.]).action == 0
    with pytest.raises(ValueError, match="reachable subspace"):
        zero.minimum_action([1., 0.])


def test_output_budget_is_sharp_and_does_not_require_full_reachability():
    matrix, inputs = -np.eye(2), np.array([[1.], [0.]])
    problem = compensation_gramian(matrix, inputs, np.eye(1), 2)
    bound = problem.output_budget_interval([1., 3.], [0., 1.], 100.)
    assert bound["signed_lower"] == bound["signed_upper"] == 3
    bound = problem.output_budget_interval([1., 3.], [1., 0.], .1)
    variance = problem.gramian[0, 0]
    amplitude = np.sqrt(.2/variance)
    shift = quad(lambda t: np.exp(-(2-t))*(-amplitude*np.exp(-(2-t))), 0, 2)[0]
    assert 1+shift == pytest.approx(bound["signed_lower"])
    assert quad(lambda t: .5*(amplitude*np.exp(-(2-t)))**2, 0, 2)[0] == pytest.approx(.1)
    with pytest.raises(ValueError, match="nonnegative energy"):
        problem.output_budget_interval([1, 3], [1, 0], -1)


@pytest.mark.parametrize("weight", [np.zeros((2, 2)), np.diag([-1., 1.]), np.array([[1., 1.], [0., 1.]])])
def test_invalid_compensation_weight_is_explicit(weight):
    with pytest.raises(ValueError):
        compensation_gramian(-np.eye(2), np.eye(2), weight, 1)


def test_gramian_quadrature_and_independent_direct_terminal_integration(model):
    structure, duration = potential_structure(model), 3652.5
    inputs = structure.allocation[:, None]
    terminal = np.array([.01, .2, .04, .5])
    problem = compensation_gramian(structure.matrix, inputs, np.eye(1), duration)
    refined = compensation_gramian(structure.matrix, inputs, np.eye(1), duration, quadrature_order=128)
    minimum = problem.minimum_action(terminal)
    assert minimum.action == pytest.approx(refined.minimum_action(terminal).action, rel=2e-11)
    # Homogeneous natural trajectory has this terminal value.
    initial = np.linalg.solve(expm(structure.matrix*duration), terminal)
    run = solve_ivp(lambda t, p: structure.matrix@p+inputs@minimum.control(t),
                    (0, duration), initial, method="DOP853", rtol=1e-11, atol=1e-12)
    assert run.success
    assert np.linalg.norm(run.y[:, -1]) < 1e-8
    with pytest.raises(ValueError, match="horizon"):
        minimum.control(-1)


def test_interval_rhs_and_assimilation_rate_contain_independent_point_evaluations(model):
    state = model.equilibria(.00244)[0]
    widths = np.array([1e-4, 2e-4, 1e-8, 1e-8, 1e-8, 1e-8])
    box = [_Interval(float(x-r), float(x+r)) for x, r in zip(state, widths)]
    rhs = conducting_rhs_enclosure(model, box, _Interval(.0025, .00254))
    rate = assimilation_rate_enclosure(model, box, rhs)
    generator = np.random.default_rng(414)
    for _ in range(50):
        point = state+widths*generator.uniform(-1, 1, 6)
        vpd = generator.uniform(.0025, .00254)
        diag = potential_diagnostics(model, point, vpd)
        for value, enclosing in zip(diag["rhs"], rhs):
            assert float(enclosing.lo) <= value <= float(enclosing.hi)
        assert float(rate.lo) <= diag["assimilation_rate"] <= float(rate.hi)


def test_same_vpd_and_weather_rate_do_not_determine_capacity_rate_sign(model):
    equilibrium = model.equilibria(.00244)[0]
    for displacement, sign in [(-.01, -1), (.01, 1)]:
        point = equilibrium.copy()
        point[0] += displacement
        box = [_Interval(float(x)) for x in point]
        rate = assimilation_rate_enclosure(model, box, conducting_rhs_enclosure(model, box, .00244))
        assert (rate.lo > 0) if sign > 0 else (rate.hi < 0)


def test_short_step_proves_membership_sign_and_nonzero_component_lower_bounds(model):
    cert = certify_drying_step_potential(model, .00244, .00254, .02)
    assert float(cert["assimilation_rate_interval"][1]) < 0
    for box, reach in zip(cert["state_box"], cert["reach_box"]):
        assert float(box[0]) < float(reach[0]) < float(reach[1]) < float(box[1])
    lower, upper = np.array(cert["absolute_potential_component_lower"]), np.array(cert["absolute_potential_component_upper"])
    assert np.all(0 < lower)
    assert np.all(lower < upper)
    assert cert["final_potential_l1_lower"] > .002
    initial = model.equilibria(.00244)[0]
    run = solve_ivp(lambda t, y: model.rhs(y, .00254), (0, .02), initial,
                    method="Radau", rtol=2e-12, atol=1e-14, dense_output=True)
    assert run.success
    final = potential_diagnostics(model, run.y[:, -1], .00254)["potential"]
    assert np.all(final < 0)
    assert np.all(lower <= -final)
    assert np.all(-final <= upper)
    sample = run.sol(np.linspace(0, .02, 101)).T
    for state in sample:
        assert all(float(box[0]) <= y <= float(box[1]) for y, box in zip(state, cert["state_box"]))
    cumulative = quad(lambda t: np.abs(potential_diagnostics(
        model, run.sol(t), .00254)["potential"]).sum(), 0, .02, epsabs=1e-13)[0]
    assert cumulative >= cert["cumulative_potential_l1_lower"] > 0


@pytest.mark.parametrize("vpd,duration", [(.00244, .02), (.0024, .02), (.00254, -1), (.00254, 10.)])
def test_failed_step_conditions_are_not_silent_zero_guarantees(model, vpd, duration):
    with pytest.raises(ValueError):
        certify_drying_step_potential(model, .00244, vpd, duration)
