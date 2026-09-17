import numpy as np
import pytest
from control_carbon.minimal_nitrogen_tipping import (
    Parameters, growth, rhs, jacobian, equilibria, simulate, classify,
    influx_at, fast_step_certificate,
)


@pytest.mark.parametrize("k", [None, .02, .1, 1.])
def test_budget_and_jacobian(k):
    p = Parameters(decomposition=k)
    y = np.array([12., 2.]) if k is None else np.array([12., 6., 2.])
    f = rhs(y, 5., p)
    export = (1-p.recycle)*(p.m*y[0] if k is None else k*y[1])
    assert sum(f) == pytest.approx(5-p.loss*y[-1]-export)
    if k is not None:
        assert sum(f[:2]) == pytest.approx(growth(y[-1], p)*y[0]-k*y[1])
    eps = 1e-5
    numerical = np.column_stack([(rhs(y+eps*e, 5., p)-rhs(y-eps*e, 5., p))/(2*eps)
                                  for e in np.eye(len(y))])
    np.testing.assert_allclose(jacobian(y, p), numerical, atol=1e-8)
    for j in range(len(y)):
        boundary = y.copy()
        boundary[j] = 0
        assert rhs(boundary, 5., p)[j] >= 0


@pytest.mark.parametrize("k", [None, .02, .1, 1., 10.])
def test_frozen_branches(k):
    p = Parameters(decomposition=k)
    for i in np.linspace(3., 10., 31):
        for name, y in equilibria(i, p).items():
            np.testing.assert_allclose(rhs(y, i, p), 0., atol=1e-12)
            eig = np.linalg.eigvals(jacobian(y, p)).real
            assert (np.max(eig) < 0) if name != "saddle" else (np.max(eig) > 0)


def test_fast_basin_certificate():
    a = fast_step_certificate()
    assert a["c_upper"] < a["c_boundary"]
    assert min(a["inward"], a["decay"], a["n_speed_lower"]) > 0


@pytest.mark.parametrize("inhibition", [False, True])
def test_rate_and_negative_control(inhibition):
    p = Parameters(inhibition=inhibition)
    outcomes = []
    for duration in (10., 40.):
        _, y = simulate(duration, p)
        assert np.all(y >= 0)
        outcomes.append(classify(y[-1], p=p))
    assert outcomes == (["bare", "forest"] if inhibition else ["forest", "forest"])


def test_same_dose():
    # Centered ramps have equal integrals, including for smoothstep.
    for smooth in (False, True):
        for duration in (10., 40., 100.):
            local = np.linspace(0, duration, 1001)
            ramp = np.trapz(influx_at(local, duration, smooth=smooth), local)
            dose = 3*(100-duration/2)+ramp+10*(1000-100-duration/2)
            assert dose == pytest.approx(9300.)


def test_monthly_annual_and_solver_refinement():
    _, annual = simulate(30., hold=100., sample_step=1.)
    _, monthly = simulate(30., hold=100., sample_step=1/12,
                          max_step=1/12, rtol=1e-10)
    np.testing.assert_allclose(annual[-1], monthly[-1], rtol=1e-7)


def test_parameter_guards():
    with pytest.raises(ValueError):
        Parameters(recycle=1.)
