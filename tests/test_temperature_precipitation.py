from dataclasses import replace

import numpy as np
import pytest
from numpy.testing import assert_allclose

from control_carbon.temperature_precipitation import (
    Climate, Controls, Parameters, WATER_INCIDENCE, budget_residuals, capacity,
    domain_violation, equilibrium, fluxes, frozen_carbon_capacity, integrate,
    integrate_months, jacobian, parameter_manifest, positive_flow, rhs,
)


@pytest.fixture(scope="module")
def baseline():
    return equilibrium()


def test_month_lengths_and_parameter_units():
    assert Climate.from_month(27, 290, 2024, 2).precipitation == 10
    assert Climate.from_month(27, 280, 2023, 2).precipitation == 10
    assert all(v["unit"] and v["source"] and v["status"] == "uncalibrated assumption"
               for v in parameter_manifest().values())
    assert_allclose(WATER_INCIDENCE.sum(axis=0), 0)
    with pytest.raises(ValueError):
        Parameters(eta=(1.1, 1., 1.))
    with pytest.raises(ValueError):
        Parameters(epsilon=(0., .1, .1))
    with pytest.raises(ValueError):
        Climate(27, -1)


def test_hand_computable_carbon_limit():
    p = replace(Parameters(), allocation=(.25, .5, .25), turnover=(.1, .2, .3),
                respiration=(.01, .02, .03), soil_respiration=.04,
                gpp_max=2., k_leaf=np.log(2))
    x = np.array([1., 2., 3., 4., .1, .1, .1, 100.])
    control = Controls(temperature_response=False, water_response=False)
    # G=2*(1-1/2)=1; internal litter=0.1+0.4+0.9=1.4.
    expected = [.25-.11, .5-.44, .25-.99, 1.4-.16]
    assert_allclose(rhs(0, x, Climate(80, 6), p, control)[:4], expected)
    assert_allclose(sum(expected), 1-(.01+.04+.09)-.16)


def test_flow_and_hazard_hand_values():
    assert positive_flow(-1., .01) == 0.
    assert positive_flow(1., 1.) == .5
    p = Parameters()
    c = np.array([3., 50., 10., 100.])
    x = np.r_[c, p.drought_half*capacity(c, p), 200.]
    f = fluxes(x, Climate(p.heat_half, 6), p, Controls(mortality="additive"))
    assert_allclose(f.heat_hazard, p.heat_max/2)
    assert_allclose(f.drought_hazard, p.drought_max/2)
    assert_allclose(f.mortality_rates, (p.heat_max+p.drought_max)/2)


@pytest.mark.parametrize("mortality", ["none", "heat", "water", "additive"])
def test_budgets_and_invariant_boundaries(mortality):
    rng = np.random.default_rng(583)
    p = replace(Parameters(), eta=(.8, .6, .7), zeta=(.7, .9, .5))
    control = Controls(mortality=mortality)
    for _ in range(15):
        c = rng.uniform(.1, 150, 4)
        x = np.r_[c, capacity(c, p)*rng.uniform(0, 1, 3), rng.uniform(0, p.soil_sat)]
        u = Climate(rng.uniform(20, 40), rng.uniform(0, 15))
        assert_allclose(budget_residuals(x, u, p, control), 0, atol=1e-10)
        for k in range(8):
            z = x.copy()
            z[k] = 0
            assert rhs(0, z, u, p, control)[k] >= -1e-12
        for k in range(3):
            z = x.copy()
            z[k+4] = capacity(c, p)[k]
            dx = rhs(0, z, u, p, control)
            assert p.omega[k]*dx[k]-dx[k+4] >= -1e-12
        z = x.copy()
        z[7] = p.soil_sat
        assert rhs(0, z, u, p, control)[7] <= 1e-12


def test_equilibrium_capacity_and_both_coupling_directions(baseline):
    x = baseline.state
    assert x.min() > 0
    assert baseline.residual < 1e-8
    assert baseline.max_real_part < -1e-8
    assert domain_violation(x) == 0
    assert_allclose(frozen_carbon_capacity(x, Climate()), x[:4], atol=1e-6)
    j = jacobian(x, Climate())
    j2 = jacobian(x, Climate(), step=5e-6)
    assert np.linalg.norm(j[:4, 4:]) > 1e-6
    assert np.linalg.norm(j[4:, :4]) > 1e-6
    assert_allclose(j, j2, rtol=2e-5, atol=1e-8)
    d = np.array([.1, .2, .1, .3, .03, .03, .03, .1])
    errors = [np.linalg.norm(rhs(0, x+eps*d, Climate())-rhs(0, x, Climate())-j@(eps*d))
              for eps in (1e-2, 5e-3)]
    assert 3.8 < errors[0]/errors[1] < 4.2


def test_matrix_budget_and_capacity_is_not_nonlinear_equilibrium(baseline):
    x = baseline.state.copy()
    x[0] *= .6
    f = fluxes(x, Climate())
    assert_allclose(f.carbon_matrix.sum(axis=0),
                    -np.r_[np.asarray(Parameters().respiration),
                            -f.carbon_matrix[3, 3]])
    assert_allclose(f.carbon_input+f.carbon_matrix@x[:4], rhs(0, x, Climate())[:4])
    cap = frozen_carbon_capacity(x, Climate())
    assert np.linalg.norm(rhs(0, np.r_[cap, x[4:]], Climate())[:4]) > 1e-6


def test_refinement_cross_solver_and_drydown(baseline):
    times = np.linspace(0, 3650, 122)
    forcing = lambda t: Climate(27+4*min(t/365, 1), 6-2*min(t/365, 1))
    a = integrate(baseline.state, times, forcing)
    b = integrate(baseline.state, times, forcing, method="BDF", max_step=5, rtol=2e-9, atol=2e-11)
    assert np.max(np.abs(a-b)/np.maximum(1, np.abs(b))) < 1e-5
    assert a.min() >= -1e-8
    dry = integrate(baseline.state, [0., 30., 100.], Climate(40, 0))
    assert max(domain_violation(x) for x in dry) < 1e-6


def test_monthly_segments_match_explicit_restart(baseline):
    climates = [Climate.from_month(27, 180, 2024, 1), Climate.from_month(31, 60, 2024, 2)]
    times, x = integrate_months(baseline.state, climates, [31, 29])
    assert_allclose(times, [0, 31, 60])
    y = integrate(baseline.state, [0., 31.], climates[0])[-1]
    z = integrate(y, [31., 60.], climates[1])[-1]
    assert_allclose(x[-1], z, rtol=0, atol=1e-12)


def test_clamps_are_process_specific(baseline):
    a = fluxes(baseline.state, Climate(27, 6), control=Controls(mortality="additive"))
    b = fluxes(baseline.state, Climate(32, 6), control=Controls(
        mortality="additive", heat_temperature=27, photo_temperature=27, respiration_temperature=27))
    assert_allclose(b.gpp, a.gpp)
    assert_allclose(b.respiration, a.respiration)
    assert_allclose(b.mortality, a.mortality)
    assert b.transpiration > a.transpiration


def test_rainfall_plateau_is_separate_from_water_limited_response(baseline):
    p = Parameters()
    mild = equilibrium(Climate(27, 4), guess=baseline.state)
    dry = equilibrium(Climate(27, 2), guess=baseline.state)
    assert mild.state[7] > p.soil_fc
    assert dry.state[7] < p.soil_fc
    # At r_O=1, additional rain changes overflow/drainage, not plant supply.
    assert_allclose(mild.state[:4], baseline.state[:4], rtol=1e-8)
    assert dry.state[:4].sum() < .9*baseline.state[:4].sum()
