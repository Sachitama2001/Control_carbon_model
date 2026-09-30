import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import solve_ivp
from scipy.linalg import solve_banded
from dataclasses import replace
from control_carbon.forest_grass_tipping import (
    Parameters, phi, h, drift, drift_prime, interior_roots, fold_points,
    turning_points, path_parameters, smooth_step, analytic_bounds,
    ramp_solution, matched_solutions, critical_duration, potential, mean_first_passage,
)


def test_published_three_roots_and_stability():
    roots = interior_roots()
    assert_allclose(roots, [.164949831889184, .415204066760038, .735328312252550],
                    atol=1e-12)
    assert_allclose(h(roots), 0, atol=1e-13)
    assert np.array_equal(np.sign(drift_prime(roots)), [-1, 1, -1])
    assert drift(0) > 0 and drift(1) == 0
    assert drift_prime(1) > 0


def test_folds_and_bistability_disappearance():
    (g0, a0), (g1, a1) = fold_points()
    assert_allclose([a0, a1], [.961809985648403, 1.38224561592477])
    for g, alpha in [(g0, a0), (g1, a1)]:
        p = replace(Parameters(), alpha=alpha)
        assert abs(h(g, p)) < 1e-13
        assert abs(drift_prime(g, p)) < 1e-13
    assert len(interior_roots(replace(Parameters(), alpha=a0 - .01))) == 1
    assert len(interior_roots(replace(Parameters(), alpha=a1 + .01))) == 1


def test_constructed_path_and_analytic_certificates():
    bounds = analytic_bounds()
    assert min(bounds.values()) > 0
    for lam in np.linspace(0, 1, 201):
        p = path_parameters(float(lam))
        roots = interior_roots(p)
        a = .4 - .1 * lam
        assert_allclose(roots[:2], [a, a + .02], atol=5e-12)
        assert len(roots) == 3 and roots[2] > p.theta
        assert p.alpha > p.phi1
        assert np.array_equal(np.sign(drift_prime(roots, p)), [-1, 1, -1])
        assert -drift(a + .005, p) >= bounds["inward_drift_bound"]
        assert drift(a - .005, p) > 0
        assert np.max(abs(drift(np.linspace(0, 1, 101), p))) < bounds["absolute_drift_bound"]
        c0, c1 = turning_points(p)
        assert h(c0, p) < 0 < h(c1, p)


def test_ramp_fixed_path_not_rate_dependent_amplitude():
    assert smooth_step(0) == 0 and smooth_step(1) == 1
    assert_allclose(smooth_step(.5), .5)
    u = np.linspace(0, 1, 501)
    assert np.all(np.diff(smooth_step(u)) >= 0)


def test_r_tipping_and_threshold_refinement():
    tau = critical_duration()
    refined = critical_duration(rtol=2e-12, atol=2e-14, max_step=.004)
    assert_allclose(tau, refined, rtol=1e-7)
    fw, bw = matched_solutions(refined, sensitivity=True, rtol=2e-12, atol=2e-14)
    assert abs(fw.y[0, -1] - bw.y[0, -1]) < 1e-10
    assert fw.y[1, -1] - bw.y[1, -1] < 0
    assert_allclose(refined, critical_duration(match=.4, rtol=2e-12, atol=2e-14),
                    rtol=1e-8)
    for factor, tipped in [(.5, True), (2., False)]:
        sol = ramp_solution(factor * tau)
        assert (sol.y[0, -1] > .32) == tipped
        assert np.all((sol.y[0] >= 0) & (sol.y[0] <= 1))
        g = sol.y[0]
        assert_allclose(g + (1 - g), 1, atol=0)
    bounds = analytic_bounds()
    assert ramp_solution(.5 * bounds["fast_duration_upper"]).y[0, -1] > .32
    assert ramp_solution(1.1 * bounds["slow_duration_lower"]).y[0, -1] < .305


def test_alpha_only_invariant_barrier():
    minimum, maximum = 1.02, 1.32
    barrier = interior_roots(replace(Parameters(), alpha=minimum))[1]
    start = interior_roots(replace(Parameters(), alpha=maximum))[0]
    for alpha in np.linspace(minimum, maximum, 51):
        assert drift(barrier, replace(Parameters(), alpha=alpha)) <= 1e-13
        assert interior_roots(replace(Parameters(), alpha=alpha))[1] >= barrier - 1e-13
    for duration in [.01, 10., 1000.]:
        rhs = lambda u, y: duration * drift(y, replace(
            Parameters(), alpha=maximum - (maximum - minimum) * smooth_step(u)))
        sol = solve_ivp(rhs, [0, 1], [start], rtol=1e-10, atol=1e-12)
        assert np.max(sol.y) < barrier


def test_potential_and_mfpt_against_independent_backward_equation():
    low, saddle, high = interior_roots()
    eps = 1e-5
    for g in [.1, .3, .6, .9]:
        derivative = (potential(g + eps) - potential(g - eps)) / (2 * eps)
        assert_allclose(derivative, -drift(g), atol=1e-9)
    assert potential(saddle) > potential(low)
    assert potential(saddle) > potential(high)
    target = (saddle + high) / 2
    sigma = .05
    exact = mean_first_passage(low, target, sigma)
    assert exact > 0
    assert_allclose(exact, mean_first_passage(low, target, sigma, order=140), rtol=1e-10)
    # Independent finite-difference backward equation; ghost-point reflecting wall.
    n = 4001
    x = np.linspace(0, target, n)
    dx = x[1]
    diffusion = sigma**2 / 2
    bands = np.zeros((3, n))
    bands[1, 1:-1] = -2 * diffusion / dx**2
    bands[0, 2:] = diffusion / dx**2 + drift(x[1:-1]) / (2 * dx)
    bands[2, :-2] = diffusion / dx**2 - drift(x[1:-1]) / (2 * dx)
    bands[1, 0] = -2 * diffusion / dx**2
    bands[0, 1] = 2 * diffusion / dx**2
    bands[1, -1] = 1
    rhs = -np.ones(n)
    rhs[-1] = 0
    numeric = solve_banded((1, 1), bands, rhs)
    assert_allclose(np.interp(low, x, numeric), exact, rtol=2e-5)
    assert mean_first_passage(low, target, .04) > exact
