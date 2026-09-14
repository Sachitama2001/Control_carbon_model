from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from control_carbon.visit_decomposition import (
    VISITDecompositionEnvironment,
    VISITDecompositionParameters,
)
from control_carbon.visit_native import (
    DECOMPOSITION_ENVIRONMENT_NAMES,
    build_native_visit_soil_bridge,
    compare_native_visit_soil_environment_perturbation,
    compare_native_visit_soil_perturbation_to_irf,
    compare_native_visit_soil_to_python,
    linearize_visit_soil_environment_trajectory,
)
from control_carbon.visit_soil import VISITSoilParameters
from control_carbon.visit_source_map import VISIT_SOURCE_COMMIT


VISIT_SOURCE_ROOT = Path(__file__).resolve().parents[2] / "VISIT-matrix" / "visit_local"


@pytest.fixture(scope="module")
def native_bridge(tmp_path_factory):
    if not VISIT_SOURCE_ROOT.is_dir():
        pytest.skip("authoritative VISIT source checkout is unavailable")
    return build_native_visit_soil_bridge(
        VISIT_SOURCE_ROOT, tmp_path_factory.mktemp("native") / "visit_soil_bridge"
    )


@pytest.fixture
def parameters():
    return VISITSoilParameters(
        sr_lf=1.1,
        sr_lc=0.3,
        sr_lr=0.75,
        sr_ha=0.18,
        sr_hi=0.08,
        sr_hp=0.025,
        f_co2_lf=0.5,
        f_co2_lc=0.5,
        f_co2_lr=0.5,
        f_hm_a=0.6,
        f_hm_i=0.3,
        f_hm_p=0.1,
    )


@pytest.fixture
def decomposition_parameters():
    return VISITDecompositionParameters(kml=0.22, kmh=0.1, kmsl=0.21, kmsh=0.07)


@pytest.fixture
def environment():
    return VISITDecompositionEnvironment(
        soil_temperature_10cm=10.0,
        soil_temperature_deep=8.0,
        soil_water_upper=50.0,
        soil_water_lower=100.0,
        soil_aperture_upper=0.8,
        soil_aperture_whole=0.9,
        field_capacity_30cm=120.0,
        field_capacity=300.0,
    )


def test_native_visit_one_day_matches_python(
    native_bridge, parameters, decomposition_parameters, environment
):
    comparison = compare_native_visit_soil_to_python(
        native_bridge,
        VISIT_SOURCE_ROOT,
        initial_state=np.arange(10.0, 100.0, 10.0),
        litter_inputs=np.ones((1, 6)),
        soil_parameters=parameters,
        decomposition_parameters=decomposition_parameters,
        environment=environment,
    )

    assert np.allclose(comparison.state_error, 0.0, atol=1e-13)
    assert np.allclose(
        comparison.heterotrophic_respiration_error, 0.0, atol=1e-13
    )
    assert comparison.native.source_commit == VISIT_SOURCE_COMMIT
    assert all(len(value) == 64 for value in comparison.native.source_sha256.values())


def test_native_visit_litter_pulse_trajectory_matches_python(
    native_bridge, parameters, decomposition_parameters, environment
):
    inputs = np.ones((40, 6))
    inputs[0, 0] += 1.0

    comparison = compare_native_visit_soil_to_python(
        native_bridge,
        VISIT_SOURCE_ROOT,
        initial_state=np.arange(10.0, 100.0, 10.0),
        litter_inputs=inputs,
        soil_parameters=parameters,
        decomposition_parameters=decomposition_parameters,
        environment=environment,
    )

    assert np.max(np.abs(comparison.state_error)) < 1e-12
    assert np.max(np.abs(comparison.heterotrophic_respiration_error)) < 1e-12
    assert np.all(comparison.native.litter_scalars > 0.0)
    assert np.all(comparison.native.humus_scalars > 0.0)


def test_native_visit_litter_pulse_difference_matches_irf(
    native_bridge, parameters, decomposition_parameters, environment
):
    perturbations = np.zeros((100, 6))
    perturbations[0, 0] = 1.0

    comparison = compare_native_visit_soil_perturbation_to_irf(
        native_bridge,
        VISIT_SOURCE_ROOT,
        baseline_litter_inputs=np.ones(6),
        input_perturbations=perturbations,
        soil_parameters=parameters,
        decomposition_parameters=decomposition_parameters,
        environment=environment,
    )

    assert np.max(np.abs(comparison.state_error)) < 2e-10
    assert np.max(np.abs(comparison.heterotrophic_respiration_error)) < 1e-13
    assert comparison.native_heterotrophic_respiration_difference[0] == pytest.approx(0.0)
    assert comparison.native_heterotrophic_respiration_difference[1] > 0.0


def test_native_visit_daily_environment_trajectory_matches_python(
    native_bridge, parameters, decomposition_parameters, environment
):
    n_steps = 120
    days = np.arange(n_steps, dtype=float)
    environments = tuple(
        replace(
            environment,
            soil_temperature_10cm=10.0 + 8.0 * np.sin(2.0 * np.pi * day / n_steps),
            soil_temperature_deep=8.0 + 4.0 * np.sin(2.0 * np.pi * (day - 10.0) / n_steps),
            soil_water_upper=50.0 + 10.0 * np.cos(2.0 * np.pi * day / n_steps),
            soil_water_lower=100.0 + 15.0 * np.cos(2.0 * np.pi * (day - 15.0) / n_steps),
        )
        for day in days
    )

    comparison = compare_native_visit_soil_to_python(
        native_bridge,
        VISIT_SOURCE_ROOT,
        initial_state=np.arange(10.0, 100.0, 10.0),
        litter_inputs=np.ones((n_steps, 6)),
        soil_parameters=parameters,
        decomposition_parameters=decomposition_parameters,
        environment=environments,
    )

    assert np.max(np.abs(comparison.state_error)) < 2e-12
    assert np.max(np.abs(comparison.heterotrophic_respiration_error)) < 1e-15
    assert np.ptp(comparison.native.litter_scalars) > 0.1
    assert np.ptp(comparison.native.humus_scalars) > 0.1
    assert comparison.native.environments.shape == (n_steps, 8)


def test_native_temperature_pulse_tangent_error_is_second_order(
    native_bridge, parameters, decomposition_parameters, environment
):
    n_steps = 80
    inputs = np.ones((n_steps, 6))
    environments = (environment,) * n_steps
    temperature_column = DECOMPOSITION_ENVIRONMENT_NAMES.index(
        "soil_temperature_10cm"
    )
    errors = []
    for amplitude in (2.0, 1.0, 0.5):
        perturbations = np.zeros((n_steps, 8))
        perturbations[5, temperature_column] = amplitude
        comparison = compare_native_visit_soil_environment_perturbation(
            native_bridge,
            VISIT_SOURCE_ROOT,
            initial_state=np.arange(10.0, 100.0, 10.0),
            litter_inputs=inputs,
            soil_parameters=parameters,
            decomposition_parameters=decomposition_parameters,
            baseline_environments=environments,
            environment_perturbations=perturbations,
        )
        assert np.max(np.abs(comparison.native_python_state_error)) < 2e-12
        assert (
            np.max(
                np.abs(comparison.native_python_heterotrophic_respiration_error)
            )
            < 1e-15
        )
        errors.append(np.max(np.abs(comparison.tangent_state_error)))

    assert errors[1] < 0.3 * errors[0]
    assert errors[2] < 0.3 * errors[1]


def test_environment_linearization_rejects_temperature_threshold(
    parameters, decomposition_parameters, environment
):
    threshold_environment = replace(environment, soil_temperature_10cm=-20.0)

    with pytest.raises(ValueError, match="nonsmooth decomposition branch"):
        linearize_visit_soil_environment_trajectory(
            initial_state=np.arange(10.0, 100.0, 10.0),
            litter_inputs=np.ones((1, 6)),
            soil_parameters=parameters,
            decomposition_parameters=decomposition_parameters,
            environments=(threshold_environment,),
        )