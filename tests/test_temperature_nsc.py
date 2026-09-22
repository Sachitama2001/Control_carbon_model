from dataclasses import replace
import numpy as np
import pytest
from control_carbon.temperature_nsc import (
    Parameters, concentrations, dry_to_carbon_fraction, photosynthesis,
    maintenance, flush_draw, evaluate, rhs, PLANT)
from control_carbon.nsc_analysis import equilibrium, linear_info, integrate, smooth_temperature

P = Parameters()
X = np.array([2., 80., 15., 100., .2, 5., 1., 1.])


def test_leaf_zero_and_temperature_respiration():
    assert photosynthesis(0, 20, P) == 0
    d0, _ = maintenance(X[:3], X[4:7], 20, P)
    d1, _ = maintenance(X[:3], X[4:7], 30, P)
    np.testing.assert_allclose(d1, P.q10*d0)


def test_flush_hand_budget_and_empty_donor():
    draw = flush_draw(X[:3], X[4:7], P)
    np.testing.assert_allclose(draw, .03*np.array([5., 1.])/6.5)
    n = np.array([.2, 0, 1])
    assert flush_draw(X[:3], n, P)[0] == 0
    dx, d = evaluate(X, 20, P)
    dx0, d0 = evaluate(X, 20, replace(P, flush_enabled=False))
    assert dx[0]-dx0[0] == pytest.approx(P.yield_c*draw.sum())
    assert (dx-dx0).sum() == pytest.approx(-(d['growth_resp']-d0['growth_resp']))


@pytest.mark.parametrize('retention', [0., .5, 1.])
def test_budget_and_boundary_inward(retention):
    rng = np.random.default_rng(10)
    p = replace(P, retention=retention)
    for _ in range(20):
        x = np.exp(rng.uniform(-8, 5, 8))
        dx, d = evaluate(x, 32, p)
        assert abs(d['budget']) < 1e-13
        for i in range(8):
            b = x.copy(); b[i] = 0
            assert rhs(b, 32, p)[i] >= 0


def test_mortality_is_internal_and_preserves_concentration():
    d = rhs(X, 20, P)-rhs(X, 20, replace(P, nsc_mortality=False))
    assert abs(d.sum()) < 1e-14
    np.testing.assert_allclose(d[:3]/X[:3], d[4:7]/X[4:7])


def test_soil_not_plant_and_bare_no_epsilon():
    x = X.copy(); x[[3, 7]] *= 1000
    assert concentrations(x)[0] == concentrations(X)[0]
    bare = np.zeros(8); bare[[3, 7]] = [10, 2]
    assert np.isnan(concentrations(bare)[0])
    np.testing.assert_array_equal(rhs(bare, 20, P)[PLANT], 0)
    assert dry_to_carbon_fraction(.1, .4, .5) == pytest.approx(.04/.49)


def test_positive_equilibrium_and_step_refinement():
    q = equilibrium(X, 20, P)
    assert q is not None
    info = linear_info(q, 20, P)
    assert info['residual'] < 1e-9
    assert info['max_real'] < 0
    times = np.linspace(0, 365, 30)
    a = integrate(X, times, lambda t: 25, P, max_step=30)
    b = integrate(X, times, lambda t: 25, P, max_step=15, rtol=1e-9)
    assert a.min_stock > 0
    assert a.integrated_budget_error < 1e-6
    assert np.max(np.abs(a.states-b.states)/np.maximum(1, b.states)) < 1e-6


def test_temperature_same_endpoints():
    for duration in (30, 365, 3650):
        assert smooth_temperature(0, duration, 20, 25) == 20
        assert smooth_temperature(duration, duration, 20, 25) == 25
        assert smooth_temperature(2*duration, duration, 20, 25) == 25


def test_invalid_coefficients_and_states():
    with pytest.raises(ValueError):
        Parameters(flush_half=0)
    with pytest.raises(ValueError):
        Parameters(yield_c=2)
    with pytest.raises(ValueError):
        rhs(-X, 20, P)


def test_bare_invasion_limit_and_reduced_root():
    from control_carbon.nsc_analysis import bare_invasion
    from control_carbon.nsc_reduction import reduced_equilibrium, reduced_rhs
    info = bare_invasion(20, P)
    v = np.array(info['direction'])
    x = np.zeros(8); x[PLANT] = 1e-8*v
    np.testing.assert_allclose(rhs(x,20,P)[PLANT]/1e-8,
                               info['invasion_exponent']*v, atol=1e-9)
    assert info['invasion_exponent'] > 0
    q = reduced_equilibrium(20)
    np.testing.assert_allclose(reduced_rhs(q,20), 0, atol=1e-12)


def test_slow_tracking_and_hold_invariance():
    q0 = equilibrium(X, 20, P)
    q1 = equilibrium(q0, 25, P)
    errors = []
    for duration in (365, 36500):
        tr = integrate(q0, [0, duration],
                       lambda t: smooth_temperature(t, duration, 20, 25), P,
                       max_step=180)
        errors.append(np.linalg.norm((tr.states[-1]-q1)/np.maximum(1,q1)))
    assert errors[1] < .05*errors[0]
    hold = integrate(tr.states[-1], [0, 250*365, 500*365], lambda t: 25, P,
                     max_step=180)
    np.testing.assert_allclose(hold.states[-2:], np.tile(q1,(2,1)), rtol=1e-6,atol=1e-8)


def test_transfer_removal_reaches_domain_limit_without_clipping():
    # Regression: default adaptive finite-difference Jacobian previously
    # overflowed on insensitive log-state/budget columns in this ablation.
    from control_carbon.nsc_analysis import spinup
    _, report = spinup(X,25,replace(P,transfer=(0,0,0,0)),max_step=90)
    assert report['domain_limit'] or report['near_boundary']
    assert not report['converged']


def test_reduced_negative_divergence():
    from control_carbon.nsc_reduction import reduced_rhs, ReducedParameters
    from control_carbon.temperature_nsc import mortality
    a=ReducedParameters()
    for x in (np.array([.1,.001]),np.array([100.,10.])):
        h=1e-6
        j=np.column_stack([(reduced_rhs(x+np.eye(2)[i]*h,25)-
                            reduced_rhs(x-np.eye(2)[i]*h,25))/(2*h) for i in range(2)])
        b,s=x
        expected=(-a.growth-2*(mortality(s/(b+s),P)+a.turnover)
                  -a.respiration*P.q10**.5*b*a.substrate_half/(a.substrate_half+s)**2)
        assert np.trace(j)==pytest.approx(expected,rel=1e-6)
        assert expected<0
