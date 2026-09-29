from dataclasses import replace
import numpy as np
import pytest
from scipy.integrate import solve_ivp
from control_carbon.reserve_feedback_screen import (
    ReserveControl, concentration_rates, reserve_fluxes, reserve_rhs,
    positive_equilibrium, reserve_jacobian, dulac_divergence, feedback_jacobian,
)


def test_hand_computable_rates_and_budget():
    p = ReserveControl()
    v, vp, r, rp = concentration_rates(.2, p)
    assert v == pytest.approx(.5)
    assert vp == pytest.approx(1.25)
    assert r == pytest.approx(2/15)
    x = [2., .4]
    f = reserve_fluxes(x, p)
    net = f['production']-f['maintenance']-f['growth_respiration']-f['structural_export']-f['reserve_export']
    assert sum(reserve_rhs(x, p)) == pytest.approx(net)


@pytest.mark.parametrize('x', [[0, 0], [0, 2], [2, 0]])
def test_nonnegative_boundary(x):
    dx = reserve_rhs(x, ReserveControl())
    assert all(dx[i] >= 0 for i in range(2) if x[i] == 0)


@pytest.mark.parametrize('hill', [1., 2., 4., 8.])
def test_equilibrium_and_full_stability_even_with_steep_storage_priority(hill):
    p = replace(ReserveControl(), hill=hill)
    x = positive_equilibrium(p)
    np.testing.assert_allclose(reserve_rhs(x, p), 0, atol=2e-14)
    j = reserve_jacobian(x, p)
    assert np.trace(j) < 0 and np.linalg.det(j) > 0
    assert np.max(np.linalg.eigvals(j).real) < 0
    c, nsc = x
    v, vp, r, rp = concentration_rates(nsc/c, p)
    cost = v+r+p.reserve_loss*nsc/c
    expected_det = p.yield_*vp*(cost-p.pmax*p.kp/(p.kp+c)**2)
    assert np.linalg.det(j) == pytest.approx(expected_det)
    for xx in ([.1, 3], [4, .02], x):
        assert dulac_divergence(xx, p) < 0


def test_no_positive_equilibrium_conditions():
    assert positive_equilibrium(replace(ReserveControl(), turnover=1)) is None
    assert positive_equilibrium(replace(ReserveControl(), pmax=.001)) is None


def test_analytic_jacobian_and_both_couplings():
    p = ReserveControl()
    x = np.array([3., .2])
    eps = 1e-5
    numeric = np.column_stack([(reserve_rhs(x+eps*e, p)-reserve_rhs(x-eps*e, p))/(2*eps)
                               for e in np.eye(2)])
    j = reserve_jacobian(x, p)
    np.testing.assert_allclose(j, numeric, rtol=1e-7, atol=1e-8)
    assert j[0, 1] != 0 and j[1, 0] != 0


def test_size_feedback_changes_determinant_not_trace():
    p = ReserveControl()
    c, nsc = positive_equilibrium(p)
    n = nsc/c
    v, vp, r, rp = concentration_rates(n, p)
    pp = p.pmax*p.kp/(p.kp+c)**2
    cost = v+r+p.reserve_loss*n
    for rc in (0., -.1):
        j = feedback_jacobian(c, n, yield_=p.yield_, turnover=p.turnover,
                              reserve_loss=p.reserve_loss, v=v, vp=vp, r=r,
                              r_n=rp, r_c=rc, production_prime=pp)
        hprime = (pp-cost)/c-rc
        assert np.linalg.det(j) == pytest.approx(-p.yield_*vp*c*hprime)
        assert np.trace(j) < 0
        assert (np.linalg.det(j) > 0) == (hprime < 0)


def test_reference_trajectories_budget_positivity_and_refinement():
    p = ReserveControl()
    eq = positive_equilibrium(p)
    def augmented(t, z):
        f = reserve_fluxes(z[:2], p)
        net = f['production']-f['maintenance']-f['growth_respiration']-f['structural_export']-f['reserve_export']
        return np.r_[reserve_rhs(z[:2], p), net]
    for initial in ([.1, .001], [1, 3], [30, .1]):
        runs = [solve_ivp(augmented, [0, 400], [*initial, 0],
                          rtol=tol, atol=tol*.01, max_step=step,
                          t_eval=np.linspace(0, 400, 401))
                for tol, step in ((1e-8, 1), (1e-10, .5))]
        for sol in runs:
            assert sol.success
            assert np.min(sol.y[:2]) >= 0
            np.testing.assert_allclose(sol.y[:2, -1], eq, rtol=1e-7)
            np.testing.assert_allclose(sol.y[:2].sum(axis=0)-sum(initial), sol.y[2], atol=1e-10)
        np.testing.assert_allclose(runs[0].y, runs[1].y, rtol=2e-7, atol=2e-8)
