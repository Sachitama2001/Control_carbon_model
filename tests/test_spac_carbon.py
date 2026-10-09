from dataclasses import replace

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from control_carbon.manzoni_hydraulics import HydraulicParameters, hydraulic_fold, transpiring_equilibria
from control_carbon.manzoni_extensions import SigmoidShape, sigmoid_equilibria
from control_carbon.spac_carbon import (
    YEAR_DAYS, CarbonParameters, TemperatureResponse, assimilation, carbon_matrix,
    carbon_rhs, carbon_budget, frozen_coefficient_capacity, reproduction_number,
    positive_carbon_equilibrium, carbon_jacobian, hydraulic_factor,
    combined_rhs, combined_jacobian, simulate,
)


C = CarbonParameters(1.2/YEAR_DAYS, (.3, .4, .3), tuple(k/YEAR_DAYS for k in [.8, .05, .2, .04]),
                     (.7, .9, .8), 10., .5, .25)
H = HydraulicParameters(.0086, .47, -1.5, -2., .001, .001)


def test_units_canopy_and_hand_budget():
    assert C.canopy_coefficient == 5.
    assert assimilation(.2, .5, C) == pytest.approx(.6/YEAR_DAYS*(1-np.exp(-1)))
    p = CarbonParameters(2., (.5, .3, .2), (1., 2., 3., 4.), (.5, .75, 1.), 1., .5, .5)
    x = np.array([1., 2., 3., 4.])
    budget = carbon_budget(x, .5, p)
    assert budget["autotrophic_respiration"] == pytest.approx(1.5)
    assert budget["heterotrophic_respiration"] == 16.
    assert budget["residual"] == pytest.approx(0., abs=1e-14)
    assert budget["npp"] == pytest.approx(budget["gpp"]-1.5)


def test_compartment_structure_positivity_and_capacity_identity():
    matrix = carbon_matrix(C)
    assert np.all(matrix-np.diag(np.diag(matrix)) >= 0)
    assert np.all(matrix.sum(axis=0) <= 0)
    assert np.max(np.linalg.eigvals(matrix).real) < 0
    x = np.array([.1, 2., .3, 4.])
    cap = frozen_coefficient_capacity(x, .8, C)
    np.testing.assert_allclose(carbon_rhs(x, .8, C), matrix@(x-cap), atol=1e-15)
    for j in range(4):
        boundary = x.copy()
        boundary[j] = 0
        assert carbon_rhs(boundary, .8, C)[j] >= 0


@pytest.mark.parametrize("h", [0., .2, 1/2.25, .8, 1.])
def test_threshold_zero_and_unique_positive_equilibria(h):
    r0 = reproduction_number(h, C)
    zero_jac = carbon_jacobian(np.zeros(4), h, C)
    assert zero_jac[0, 0] == pytest.approx(C.losses[0]*(r0-1))
    np.testing.assert_array_equal(carbon_rhs(np.zeros(4), h, C), np.zeros(4))
    positive = positive_carbon_equilibrium(h, C)
    if r0 <= 1:
        assert positive is None
        assert zero_jac[0, 0] <= 1e-18
    else:
        assert positive is not None and np.all(positive > 0)
        np.testing.assert_allclose(carbon_rhs(positive, h, C), 0., atol=1e-15)
        np.testing.assert_allclose(frozen_coefficient_capacity(positive, h, C), positive, atol=1e-12)
        assert np.max(np.linalg.eigvals(carbon_jacobian(positive, h, C)).real) < 0
        # Per-capita canopy assimilation is strictly decreasing.
        leaves = np.geomspace(1e-8, 100, 100)
        values = np.array([C.allocation[0]*assimilation(l, h, C)/l-C.losses[0] for l in leaves])
        assert np.all(np.diff(values) < 0)


def test_input_off_exact_decay_and_zero_capacity():
    x = np.array([.3, 2., .4, 5.])
    np.testing.assert_array_equal(frozen_coefficient_capacity(x, 0., C), np.zeros(4))
    end = expm(carbon_matrix(C)*500*YEAR_DAYS)@x
    assert np.max(end) < 1e-6
    assert end[3] > end[0]


@pytest.mark.parametrize("shape", [None, SigmoidShape(6., 6.)])
def test_full_jacobian_coupling_and_finite_differences(shape):
    y = np.array([-.8, -1.2, .2, 2., .4, 5.])
    jac = combined_jacobian(y, -.5, .01, H, C, shape=shape)
    step = 1e-6
    fd = np.column_stack([
        (combined_rhs(y+step*e, -.5, .01, H, C, shape=shape)
         -combined_rhs(y-step*e, -.5, .01, H, C, shape=shape))/(2*step)
        for e in np.eye(6)
    ])
    np.testing.assert_allclose(jac, fd, rtol=1e-8, atol=1e-8)
    np.testing.assert_array_equal(jac[:2, 2:], 0)
    assert np.any(jac[2:, 1] > 0)
    eig = np.sort(np.linalg.eigvals(jac).real)
    expected = np.sort(np.r_[np.linalg.eigvals(jac[:2, :2]).real, np.linalg.eigvals(jac[2:, 2:]).real])
    np.testing.assert_allclose(eig, expected, atol=1e-12)


def test_hydraulic_three_roots_lift_to_positive_saddle_and_zero_carbon():
    shape = SigmoidShape(6., 6.)
    leaves = sigmoid_equilibria(-.5, .01, H, shape)
    assert len(leaves) == 3
    labels = []
    for leaf in leaves:
        h = hydraulic_factor(leaf, H, shape)
        positive = positive_carbon_equilibrium(h, C)
        carbon = np.zeros(4) if positive is None else positive
        y = np.r_[(-.5+leaf)/2, leaf, carbon]
        np.testing.assert_allclose(combined_rhs(y, -.5, .01, H, C, shape=shape), 0., atol=1e-10)
        eig = np.linalg.eigvals(combined_jacobian(y, -.5, .01, H, C, shape=shape))
        labels.append("stable" if np.max(eig.real) < 0 else "unstable")
    assert labels == ["stable", "unstable", "stable"]
    assert positive_carbon_equilibrium(hydraulic_factor(leaves[0], H, shape), C) is not None
    assert positive_carbon_equilibrium(hydraulic_factor(leaves[2], H, shape), C) is None
    assert hydraulic_factor(leaves[2], H, shape) > 0  # zero biomass is NOT exact stomatal closure


def test_linear_dry_state_remains_nonisolated_after_carbon_lifting():
    for x in [-.5, -1.5]:
        y = np.r_[x, -5., np.zeros(4)]
        np.testing.assert_array_equal(combined_rhs(y, 0., .03, H, C), 0.)
        eig = np.linalg.eigvals(combined_jacobian(y, 0., .03, H, C))
        assert np.count_nonzero(eig == 0) == 2


def test_temperature_threshold_and_safe_final_recovery_without_rate_tipping():
    response = TemperatureResponse(25., 25., 15., 2.)
    assert response.factors(25.) == (1., 1.)
    leaves = transpiring_equilibria(0., .01, H)
    h = hydraulic_factor(leaves[0], H)
    initial = positive_carbon_equilibrium(h, C)
    assert initial is not None
    assert reproduction_number(h, C, *response.factors(45.)) < 1
    # Same 25 -> 45 -> 25 temperature pulse, different durations.
    for duration in [1., 10.]:
        def log_rhs(t, z):
            temperature = 25+20*np.sin(np.pi*t/(duration*YEAR_DAYS))**2
            x = np.exp(z)
            return carbon_rhs(x, h, C, *response.factors(temperature))/x
        pulse = solve_ivp(log_rhs, (0., duration*YEAR_DAYS), np.log(initial), rtol=1e-9, atol=1e-11)
        assert pulse.success
        hold = solve_ivp(lambda t, z: carbon_rhs(np.exp(z), h, C)/np.exp(z),
                         (0., 500*YEAR_DAYS), pulse.y[:, -1], rtol=1e-10, atol=1e-12)
        assert hold.success
        np.testing.assert_allclose(np.exp(hold.y[:, -1]), initial, rtol=1e-7)


def test_fold_independent_of_carbon_and_near_threshold_positive_root():
    fold = hydraulic_fold(0., H)
    for p in [C, replace(C, assimilation_max=2*C.assimilation_max)]:
        assert hydraulic_fold(0., H) == fold
        h = (1+1e-7)/reproduction_number(1., p)
        positive = positive_carbon_equilibrium(h, p)
        assert positive[0] > 0
        np.testing.assert_allclose(carbon_rhs(positive, h, p), 0., atol=1e-15)


def test_simulator_constant_equilibrium_and_positive_log_stocks():
    leaf = transpiring_equilibria(0., .01, H)[0]
    stocks = positive_carbon_equilibrium(hydraulic_factor(leaf, H), C)
    y = np.r_[leaf/2, leaf, stocks]
    result = simulate(y, np.linspace(0., 30., 10), lambda t: (0., .01, 1., 1.), H, C)
    np.testing.assert_allclose(result["state"], np.broadcast_to(y, result["state"].shape), atol=1e-9)
    with pytest.raises(ValueError):
        simulate(np.r_[y[:2], np.zeros(4)], [0., 10.], lambda t: (0., .01, 1., 1.), H, C)


@pytest.mark.parametrize("change", [
    {"allocation": (.3, .3, .3)}, {"losses": (0., .1, .2, .3)},
    {"litter_fractions": (1.1, .3, .3)}, {"specific_leaf_area": np.nan},
    {"assimilation_max": -1.}, {"carbon_dry_fraction": 2.},
])
def test_invalid_parameters(change):
    with pytest.raises(ValueError):
        replace(C, **change)


def test_invalid_states_environment_and_temperature():
    with pytest.raises(ValueError):
        carbon_rhs([-1., 0., 0., 0.], .5, C)
    with pytest.raises(ValueError):
        reproduction_number(2., C)
    with pytest.raises(ValueError):
        combined_rhs([0., 0.], 0., .01, H, C)
    with pytest.raises(ValueError):
        TemperatureResponse(25., 25., 0., 2.)
    with pytest.raises(ValueError):
        TemperatureResponse(25., 25., 15., 2.).factors(np.nan)
