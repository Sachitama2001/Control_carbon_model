from dataclasses import replace
from decimal import Decimal, localcontext

import numpy as np
import pytest

from control_carbon.manzoni_hydraulics import HydraulicParameters, conductances, transpiring_equilibria
from control_carbon.spac_carbon import YEAR_DAYS, CarbonParameters, hydraulic_factor, positive_carbon_equilibrium
from control_carbon import spac_carbon_two_way as legacy
from control_carbon import spac_carbon_bifurcation as bifurcation
from control_carbon.spac_carbon_stability import certify_stable_path
from control_carbon.spac_carbon_mpa import PressureSPAC


H = HydraulicParameters(.0086, .47, -1.5, -2., .001, .001)
C = CarbonParameters(1.2/YEAR_DAYS, (.3, .4, .3), tuple(k/YEAR_DAYS for k in [.8, .05, .2, .04]),
                     (.7, .9, .8), 10., .5, .25)
LEAF = transpiring_equilibria(0., .01, H)[0]
STOCKS = positive_carbon_equilibrium(hydraulic_factor(LEAF, H), C)
F = legacy.StructuralFeedback(STOCKS[0], STOCKS[1], .01, 1e-6)
P = .101325
BOUNDS = (-2.222, -.001)


@pytest.mark.parametrize("pressure", [.08, P, .12])
@pytest.mark.parametrize("factors", [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1)])
def test_pressure_units_preserve_flux_rhs_jacobian_equilibria_and_curve_derivatives(pressure, factors):
    f = replace(F, leaf_supply=factors[0], leaf_demand=factors[1], wood_supply=factors[2])
    model = PressureSPAC(H, C, f, 0., pressure)
    params = 0., H, C, f
    y, d, slope = legacy.equilibrium_curve(-1.5, *params)
    yp, vpd, sp = model.equilibrium_curve(-1.5)
    np.testing.assert_array_equal(yp, y)
    assert vpd == pytest.approx(d*pressure)
    assert sp == pytest.approx(slope*pressure)
    fd = (model.equilibrium_curve(-1.5+1e-5)[1]-model.equilibrium_curve(-1.5-1e-5)[1])/2e-5
    assert sp == pytest.approx(fd, rel=1e-8)
    perturbed = y*np.array([1.01, 1.02, .97, 1.04, 1.02, 1.01])
    np.testing.assert_allclose(model.rhs(perturbed, vpd), legacy.two_way_rhs(perturbed, 0., d, H, C, f), atol=1e-14)
    np.testing.assert_allclose(model.jacobian(perturbed, vpd), legacy.two_way_jacobian(perturbed, 0., d, H, C, f), atol=1e-14)
    erx, exl, e = model.fluxes(perturbed, vpd)
    gs = conductances(perturbed[1], H)[1]
    fe = legacy.feedback_factors(perturbed[2], perturbed[3], f)[1]
    assert e == pytest.approx(gs/pressure*fe*vpd)
    assert H.cx*model.rhs(perturbed, vpd)[0]+H.cl*model.rhs(perturbed, vpd)[1] == pytest.approx(erx-e, abs=1e-15)
    critical = model.stationary_points(BOUNDS)
    np.testing.assert_allclose(model.equilibria(vpd, BOUNDS, critical),
                               legacy.equilibria_at_vpd(d, BOUNDS, critical, *params), atol=1e-12)


def test_continuation_and_nonconstant_forcing_trajectories_unchanged():
    f = replace(F, leaf_supply=1., leaf_demand=1.)
    model = PressureSPAC(H, C, f, 0., P)
    expected, status = legacy.continue_curve(-.5, .02, 15, BOUNDS, 0., H, C, f)
    records, got_status = model.continue_curve(-.5, .02, 15, BOUNDS, vpd_scale_mpa=.02*P)
    assert status == got_status
    for (y, d), (yp, vp) in zip(expected, records):
        np.testing.assert_allclose(yp, y, atol=1e-12)
        assert vp == pytest.approx(d*P)
    initial = np.r_[LEAF/2, LEAF, STOCKS]
    times = [0., 10., 50., 200.]
    def forcing(t):
        return .01+.02*min(t/10, 1.)
    expected = legacy.simulate_two_way(initial, times, forcing, 0., H, C, f)
    result = model.simulate(initial, times, lambda t: P*forcing(t))
    assert result["status"] == expected["status"]
    np.testing.assert_allclose(result["time_days"], expected["time_days"], atol=1e-6)
    np.testing.assert_allclose(result["state"], expected["state"], rtol=1e-7, atol=1e-10)


def test_hopf_chain_rule_tensors_and_shooting_shift_preserve_periodic_orbit():
    f = replace(F, leaf_supply=1., leaf_demand=1.)
    model = PressureSPAC(H, C, f, 0., P)
    old = bifurcation.hopf_point((-2.1, -1.9), 0., H, C, f)
    hopf = model.hopf_point((-2.1, -1.9))
    assert "D" not in hopf and "curve_slope" not in hopf
    assert hopf["vpd_mpa"] == pytest.approx(old["D"]*P)
    np.testing.assert_allclose(hopf["eigenvalue_slope_per_vpd_mpa_refinement"],
                               np.array(old["eigenvalue_slope_D_refinement"])/P)
    np.testing.assert_array_equal(hopf["eigenvalues"], old["eigenvalues"])
    assert hopf["l1"] == old["l1"] > 0
    y = hopf["state"][:3]
    np.testing.assert_allclose(model.core_rhs(y, hopf["vpd_mpa"]), bifurcation.core_rhs(y, 0., old["D"], H, C, f), atol=1e-14)
    np.testing.assert_allclose(model.core_jacobian(y, hopf["vpd_mpa"]), bifurcation.core_jacobian(y, 0., old["D"], H, C, f), atol=1e-14)
    for expected, got in zip(bifurcation.core_tensors(y, 0., old["D"], H, C, f), model.core_tensors(y, hopf["vpd_mpa"])):
        np.testing.assert_allclose(got, expected, atol=1e-12)
    amplitude = .001
    cycle = model.periodic_orbit(hopf, amplitude, leaf_bracket=(-2.1, -1.9))
    expected = bifurcation.periodic_orbit(old, amplitude, 0., H, C, f, leaf_bracket=(-2.1, -1.9))
    assert "D" not in cycle and "predicted_D_shift" not in cycle
    assert cycle["predicted_vpd_shift_mpa"] == pytest.approx(
        -hopf["omega"]*hopf["l1"]*amplitude**2/hopf["eigenvalue_slope_per_vpd_mpa_refinement"][-1])
    assert cycle["vpd_mpa"] == pytest.approx(expected["D"]*P, abs=1e-12)
    assert cycle["period_days"] == pytest.approx(expected["period_days"], rel=1e-7)
    np.testing.assert_allclose(cycle["full_floquet"], expected["full_floquet"], rtol=1e-7, atol=1e-10)
    np.testing.assert_allclose(cycle["state"], expected["state"], rtol=1e-7, atol=1e-9)


@pytest.mark.parametrize("supply,leaf_bounds", [(0., (.19, .30)), (1., (.145, .30))])
def test_certificate_outward_conversion_covers_exact_mpa_input(supply, leaf_bounds):
    f = replace(F, leaf_supply=supply, leaf_demand=1.)
    model = PressureSPAC(H, C, f, 0., P)
    bounds = (.00101325, .00303975)
    result = model.certify_stable_path(leaf_bounds, bounds)
    original = certify_stable_path(leaf_bounds, (.01, .03), 0., H, C, f)
    assert result["cell_count"] == original["cell_count"]
    with localcontext() as ctx:
        ctx.prec = 80
        p = Decimal.from_float(P)
        assert Decimal(result["endpoint_vpd_mpa"]["low_leaf"][0]) >= Decimal.from_float(bounds[1])
        assert Decimal(result["endpoint_vpd_mpa"]["high_leaf"][1]) <= Decimal.from_float(bounds[0])
        for old, new in zip(original["cells"], result["cells"]):
            assert "D" not in new and "Dprime_L" not in new
            for key, new_key in [("D", "vpd_mpa"), ("Dprime_L", "vpd_slope_per_leaf_carbon")]:
                assert Decimal(new[new_key][0]) <= Decimal(old[key][0])*p
                assert Decimal(new[new_key][1]) >= Decimal(old[key][1])*p
        assert Decimal(result["vpd_slope_per_leaf_carbon_upper_bound"]) < 0


@pytest.mark.parametrize("bad", [0., -1., np.inf, np.nan])
def test_invalid_pressure_rejected(bad):
    with pytest.raises(ValueError, match="atmospheric pressure"):
        PressureSPAC(H, C, F, 0., bad)


@pytest.mark.parametrize("bad", [-.001, np.inf, np.nan])
def test_invalid_vpd_rejected(bad):
    model = PressureSPAC(H, C, F, 0., P)
    with pytest.raises(ValueError, match="VPD"):
        model.rhs(np.r_[LEAF/2, LEAF, STOCKS], bad)


def test_temperature_labelled_vpd_range_and_no_clamping_of_formal_values():
    model = PressureSPAC(H, C, F, 0., P)
    safe = model.air_diagnostic(.003, 25.)
    assert safe["saturation_vapor_pressure_mpa"] == pytest.approx(.003168, rel=1e-3)
    assert 0 < safe["implied_relative_humidity"] < 1
    assert safe["admissible_at_diagnostic_temperature"]
    formal = model.air_diagnostic(.015338, 25.)
    assert not formal["admissible_at_diagnostic_temperature"]
    assert formal["implied_relative_humidity"] < 0
    assert model.pressure(model.mole_fraction(.015338)) == pytest.approx(.015338)
    with pytest.raises(ValueError):
        model.air_diagnostic(.003, -273.15)
