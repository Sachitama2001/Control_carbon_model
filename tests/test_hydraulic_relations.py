import numpy as np
import pytest

from control_carbon.coupled_carbon_water import (
    CoupledCarbonWaterForcing,
    CoupledCarbonWaterParameters,
    evaluate_coupled_carbon_water,
)
from control_carbon.hydraulic_relations import (
    PlantPressureVolumeParameters,
    VISITCSoilTexture,
    plant_pressure_volume_capacitance,
    plant_pressure_volume_potential,
    visitc_soil_matric_potential,
    visitc_soil_total_potential,
)


def test_visitc_soil_retention_matches_hand_calculation_for_all_textures():
    water = 50.0
    capacity = 100.0
    for texture, coefficient, exponent in (
        (VISITCSoilTexture.SANDY, 0.121, 4.05),
        (VISITCSoilTexture.MEDIUM, 0.478, 5.39),
        (VISITCSoilTexture.FINE, 0.405, 11.4),
    ):
        expected = -coefficient * (water / capacity) ** (-exponent)
        assert visitc_soil_matric_potential(
            water, capacity, texture
        ) == pytest.approx(expected)


def test_visitc_soil_retention_preserves_source_dry_floor_and_gravity():
    at_zero = visitc_soil_total_potential(
        0.0, 100.0, 0, gravitational_potential=-0.05
    )
    at_floor = visitc_soil_total_potential(
        0.2, 100.0, 0, gravitational_potential=-0.05
    )
    assert at_zero == pytest.approx(at_floor)


def test_pressure_volume_curve_is_continuous_at_turgor_loss():
    parameters = PlantPressureVolumeParameters(10.0, -1.5, 10.0)
    water_tlp = parameters.saturated_water * parameters.relative_water_at_turgor_loss
    epsilon = 1.0e-8
    left = plant_pressure_volume_potential(water_tlp - epsilon, parameters)
    right = plant_pressure_volume_potential(water_tlp + epsilon, parameters)
    assert left == pytest.approx(parameters.turgor_loss_potential, rel=1e-7)
    assert right == pytest.approx(parameters.turgor_loss_potential, rel=1e-7)
    assert plant_pressure_volume_potential(10.0, parameters) == pytest.approx(0.0)


def test_pressure_volume_capacitance_matches_finite_difference():
    parameters = PlantPressureVolumeParameters(10.0, -1.5, 10.0)
    water = 8.0
    step = 1.0e-5
    dpsi_dwater = (
        plant_pressure_volume_potential(water + step, parameters)
        - plant_pressure_volume_potential(water - step, parameters)
    ) / (2.0 * step)
    expected = 1.0 / dpsi_dwater
    assert plant_pressure_volume_capacitance(water, parameters) == pytest.approx(
        expected, rel=2e-6
    )


def test_reduced_transport_is_conductance_times_sourced_potential_gradient():
    parameters = CoupledCarbonWaterParameters()
    state = np.array([2.0, 30.0, 7.0, 120.0, 0.8, 8.0, 3.0, 150.0])
    evaluation = evaluate_coupled_carbon_water(
        state, CoupledCarbonWaterForcing(precipitation=1.5), parameters
    )
    psi = evaluation.fluxes.water_potential
    expected_gradient = np.array((psi[3] - psi[2], psi[2] - psi[1], psi[1] - psi[0]))
    np.testing.assert_allclose(
        evaluation.fluxes.water_internal,
        evaluation.fluxes.effective_hydraulic_conductance * expected_gradient,
    )


def test_invalid_hydraulic_parameters_are_rejected():
    with pytest.raises(ValueError, match="negative"):
        PlantPressureVolumeParameters(1.0, 1.0, 10.0)
    with pytest.raises(ValueError, match="selector"):
        visitc_soil_matric_potential(1.0, 2.0, 4)
