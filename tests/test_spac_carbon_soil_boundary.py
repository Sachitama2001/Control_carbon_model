from pathlib import Path

import numpy as np
import pytest

from control_carbon.manzoni_extensions import SoilRootParameters, soil_root_conductance
from control_carbon.manzoni_hydraulics import conductances
from control_carbon.spac_carbon_rational import RationalSPAC
from control_carbon.spac_carbon_soil_boundary import SoilBoundaryRationalSPAC


@pytest.fixture
def model(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]/"examples"))
    from prepare_manuscript_phase_portraits import model_family
    return SoilBoundaryRationalSPAC(
        RationalSPAC(model_family()["E4"]), SoilRootParameters(.1, -.1), -.3)


def test_root_node_continuity_and_plant_water_budget(model):
    state = np.array([-1.3, -2., .03, .7, .12, 1.5])
    root = model.root_potential(state)
    gx, _ = conductances(state[1], model.plant.base.hydraulic)
    _, _, fw = model.plant.factors(state)
    gsr = soil_root_conductance(model.psi_soil_mpa, model.soil)
    erx, exl, evaporation = model.fluxes(state, .00244)
    assert min(model.psi_soil_mpa, state[0]) <= root <= max(model.psi_soil_mpa, state[0])
    assert erx == pytest.approx(gsr*(model.psi_soil_mpa-root), abs=1e-17)
    assert erx == pytest.approx(2*gx*fw*(root-state[0]), abs=1e-17)
    rhs = model.rhs(state, .00244)
    h = model.plant.base.hydraulic
    assert h.cx*rhs[0]+h.cl*rhs[1] == pytest.approx(erx-evaporation, abs=1e-17)
    np.testing.assert_array_equal(rhs[2:], model.plant.rhs(state, .00244)[2:])
    np.testing.assert_array_equal(model.capacity(state), model.plant.capacity(state))


def test_full_analytic_jacobian_includes_root_chain_rule(model):
    state = np.array([-1.3, -2., .03, .7, .12, 1.5])
    step = 1e-6
    finite_difference = np.column_stack([
        (model.rhs(state+step*e, .00244)-model.rhs(state-step*e, .00244))/(2*step)
        for e in np.eye(6)])
    np.testing.assert_allclose(model.jacobian(state, .00244), finite_difference,
                               rtol=3e-8, atol=1e-9)


@pytest.mark.parametrize("psi", [-.1, -.3, -.5])
@pytest.mark.parametrize("vpd", [0., .00244, .00254])
def test_quartic_equilibria_close_all_six_balances(model, psi, vpd):
    coupled = SoilBoundaryRationalSPAC(model.plant, model.soil, psi)
    assert coupled.equilibrium_polynomial(.00244).degree() == 4
    roots = coupled.equilibria(vpd)
    assert len(roots) == 1
    for root in roots:
        np.testing.assert_allclose(coupled.rhs(root, vpd), 0, atol=2e-13)
        np.testing.assert_allclose(root[2:], coupled.capacity(root), atol=2e-13)
        assert coupled.root_potential(root) <= psi+1e-14


def test_nonconducting_node_and_zero_leaf_do_not_gain_a_storage_state(model):
    zero_leaf = np.array([-1., -2., 0., .7, .1, 1.])
    assert len(model.rhs(zero_leaf, .00244)) == 6
    assert model.fluxes(zero_leaf, .00244)[2] == 0
    np.testing.assert_array_equal(model.capacity(zero_leaf), np.zeros(4))
    closed = np.array([-1., -3.5, .02, .7, .1, 1.])
    assert model.root_potential(closed) == pytest.approx(model.psi_soil_mpa)
    assert model.fluxes(closed, .00244)[0] == 0
    assert model.rhs(closed, .00244)[0] == 0


def test_large_gsr_has_stable_flux_limit_and_soil_domain_is_explicit(model):
    state = np.array([-1.3, -2., .03, .7, .12, 1.5])
    coupled = SoilBoundaryRationalSPAC(
        model.plant, SoilRootParameters(1e14, -.1), -.3)
    gx, _ = conductances(state[1], model.plant.base.hydraulic)
    _, _, fw = model.plant.factors(state)
    assert coupled.fluxes(state, .00244)[0] == pytest.approx(
        2*gx*fw*(-.3-state[0]), rel=1e-12)
    with pytest.raises(ValueError):
        SoilBoundaryRationalSPAC(model.plant, model.soil, -.01)
