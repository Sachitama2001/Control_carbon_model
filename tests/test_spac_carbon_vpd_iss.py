from decimal import Decimal
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from control_carbon.spac_carbon_capacity import reference_scales
from control_carbon.spac_carbon_rational import RationalSPAC
from control_carbon.spac_carbon_stability import _Interval
from control_carbon.spac_carbon_vpd_iss import (
    _magnitude, _sqrt, certify_local_vpd_iss, conducting_jacobian_enclosure,
    equilibrium_enclosure, initial_energy_upper, local_iss_bounds,
    quartic_instability_obstruction,
)


@pytest.fixture
def model(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]/"examples"))
    from prepare_manuscript_phase_portraits import model_family
    return RationalSPAC(model_family()["E4"])


@pytest.fixture
def certificate(model):
    return certify_local_vpd_iss(model, .00244, reference_scales(model.base).state, .0015)


def test_interval_helpers_do_not_round_endpoints_to_default_decimal_precision():
    value = Decimal("-1.000000000000000000000000000000000000001")
    assert _magnitude(_Interval(value)) == value.copy_abs()
    root = _sqrt(_Interval(2))
    assert (root*root).lo <= 2 <= (root*root).hi
    with pytest.raises(ValueError):
        _sqrt(_Interval(-1, 1))


@pytest.mark.parametrize("vpd", [.00244, .003])
def test_positive_root_is_enclosed_with_verified_endpoint_signs(model, vpd):
    equilibrium, intervals, evidence = equilibrium_enclosure(model, vpd)
    for value, interval in zip(equilibrium, intervals):
        assert interval.lo < Decimal.from_float(float(value)) < interval.hi
    assert Decimal(evidence["low_leaf_vpd"][0]) > Decimal.from_float(vpd)
    assert Decimal(evidence["high_leaf_vpd"][1]) < Decimal.from_float(vpd)
    np.testing.assert_allclose(model.rhs(equilibrium, vpd), 0, atol=2e-14)


def test_metric_nonlinear_box_and_input_gain_are_independently_consistent(model, certificate):
    cert = certificate
    eigenvalues = np.linalg.eigvalsh(cert.metric)
    assert 0 < cert.metric_lower < eigenvalues.min()
    assert eigenvalues.max() < cert.metric_upper
    assert 0 < cert.max_vpd_perturbation_mpa < 1e-9
    assert cert.dissipation_lower_per_day > .5
    box = [_Interval(*pair) for pair in cert.evidence["state_box"]]
    lower = np.array([float(x.lo) for x in box])
    upper = np.array([float(x.hi) for x in box])
    jacobian_box = conducting_jacobian_enclosure(model, box, cert.vpd_mpa)
    for point in np.random.default_rng(1020).uniform(lower, upper, (20, 6)):
        jacobian = model.jacobian(point, cert.vpd_mpa)
        for i in range(6):
            for j in range(6):
                # The interval encloses real arithmetic; the comparison uses rounded numpy arithmetic.
                tolerance = 1e-13*max(1., abs(jacobian[i, j]))
                assert float(jacobian_box[i][j].lo)-tolerance <= jacobian[i, j]
                assert jacobian[i, j] <= float(jacobian_box[i][j].hi)+tolerance
        scaled = jacobian*cert.scales[None, :]/cert.scales[:, None]
        dissipation = -(scaled.T@cert.metric+cert.metric@scaled)
        assert np.linalg.eigvalsh(dissipation).min() >= cert.dissipation_lower_per_day
        forcing = ((model.rhs(point, cert.vpd_mpa+1e-7)
                    -model.rhs(point, cert.vpd_mpa-1e-7))/2e-7)/cert.scales
        assert np.linalg.norm(cert.metric@forcing) <= cert.input_cross_gain_per_mpa


def test_ellipsoid_projection_contains_its_boundary(certificate):
    c = certificate
    inverse = np.linalg.inv(c.metric)
    bounds = np.array([float(x[1]) for x in c.evidence["ellipsoid_axis_bounds_scaled"]])
    assert np.all(np.sqrt(c.energy_limit*np.diag(inverse)) <= bounds*(1+1e-13))
    for j in range(6):
        point = np.sqrt(c.energy_limit/inverse[j, j])*inverse[:, j]
        assert point@c.metric@point == pytest.approx(c.energy_limit)
        assert np.all(np.abs(point) <= bounds*(1+1e-12))


def test_nonlinear_full_state_step_obeys_local_bounds(model, certificate):
    c = certificate
    initial = c.equilibrium.copy()
    initial[0] += 1e-5
    initial_energy = initial_energy_upper(c, initial)
    assert 0 < initial_energy < c.energy_limit
    amplitude = .8*c.max_vpd_perturbation_mpa
    day = np.r_[0., np.geomspace(.001, 365.25, 70)]
    run = solve_ivp(lambda t, y: model.rhs(y, c.vpd_mpa+amplitude),
                    (0., day[-1]), initial, method="Radau", t_eval=day,
                    rtol=2e-12, atol=1e-14)
    assert run.success
    x = (run.y.T-c.equilibrium)/c.scales
    energy = np.einsum("ni,ij,nj->n", x, c.metric, x)
    error = np.array([(y[2:]-model.capacity(y))/c.scales[2:] for y in run.y.T])
    bound = local_iss_bounds(c, day, initial_energy, amplitude)
    assert np.max(energy-bound["energy"]) < 1e-15
    assert np.all(np.linalg.norm(x, axis=1) <= bound["state_norm_scaled"]+1e-12)
    assert np.all(np.linalg.norm(error, axis=1)
                  <= bound["capacity_error_norm_scaled"]+1e-12)
    assert np.min(run.y[1]) > -3.


def test_verified_instability_not_just_failed_sufficient_condition(model):
    obstruction = quartic_instability_obstruction(model, .003)
    assert Decimal(obstruction["third_hurwitz_determinant"][1]) < 0
    equilibrium = model.equilibria(.003)[0]
    active = model.jacobian(equilibrium, .003)[:4, :4]
    coefficients = np.poly(active)[1:]
    for coefficient, interval in zip(coefficients, obstruction["characteristic_coefficients"]):
        assert float(interval[0]) <= coefficient <= float(interval[1])
    assert np.max(np.linalg.eigvals(active).real) > .01
    with pytest.raises(ValueError, match="did not prove"):
        quartic_instability_obstruction(model, .00244)


@pytest.mark.parametrize("radius", [0., -1., np.inf, np.nan, .01])
def test_bad_or_unresolved_region_is_not_certified(model, radius):
    with pytest.raises(ValueError):
        certify_local_vpd_iss(model, .00244, reference_scales(model.base).state, radius)


@pytest.mark.parametrize("case", ["large_input", "large_initial", "negative_time", "nan_input"])
def test_bounds_refuse_extrapolation(certificate, case):
    energy, amplitude, time = 0., 0., 1.
    if case == "large_input":
        amplitude = 2*certificate.max_vpd_perturbation_mpa
    elif case == "large_initial":
        energy = 2*certificate.energy_limit
    elif case == "negative_time":
        time = -1.
    else:
        amplitude = np.nan
    with pytest.raises(ValueError):
        local_iss_bounds(certificate, time, energy, amplitude)


def test_zero_leaf_audit_preserves_transport_but_not_assimilation(model):
    state = np.r_[-.3, -1., 0., 1., .2, 2.]
    assert model.assimilation(state) == 0
    assert model.fluxes(state, .00244)[2] == 0
    np.testing.assert_array_equal(model.capacity(state), np.zeros(4))
    assert model.fluxes(state, .00244)[1] > 0
    assert model.rhs(state, .00244)[2] == 0
    assert np.linalg.norm(model.rhs(state, .00244)) > 0
    all_zero_carbon = np.r_[state[:2], np.zeros(4)]
    np.testing.assert_array_equal(model.capacity(all_zero_carbon), np.zeros(4))
    assert np.linalg.norm(model.rhs(all_zero_carbon, .00244)[:2]) > 0


def test_water_and_carbon_balance_audit(model):
    state = np.r_[-.85, -1.62, .068, 1.45, .273, 3.68]
    rhs = model.rhs(state, .00244)
    erx, _, evaporation = model.fluxes(state, .00244)
    h, c = model.base.hydraulic, model.base.carbon
    assert h.cx*rhs[0]+h.cl*rhs[1] == pytest.approx(erx-evaporation, abs=1e-18)
    output = sum((1-np.array(c.litter_fractions))*np.array(c.losses[:3])*state[2:5])
    output += c.losses[3]*state[5]
    assert sum(rhs[2:]) == pytest.approx(model.assimilation(state)-output, abs=1e-18)
