import os
from pathlib import Path

import numpy as np
import pytest

from control_carbon.visitc_hydrology import (
    VISITCHydrologyForcing,
    VISITCHydrologyParameters,
    VISITCHydrologyState,
    VISITCCanopyRadiationParameters,
    VISITCCanopyConductanceParameters,
    VISITCPenmanMonteithEnvironment,
)
from control_carbon.visitc_native import (
    build_native_visitc_hydrology_bridge,
    build_native_visitc_pm_bridge,
    compare_native_visitc_hydrology_to_python,
    compare_native_visitc_canopy_conductance_to_python,
    compare_native_visitc_extinction_to_python,
    compare_native_visitc_lai_to_python,
    compare_native_visitc_pm_to_python,
)


def _source_root() -> Path:
    candidates = (
        Path(os.environ.get("VISITC_SOURCE_ROOT", "")),
        Path(__file__).resolve().parents[2] / "VISITc" / "point",
        Path("/tmp/control_carbon_visitc_source/point"),
    )
    for candidate in candidates:
        if candidate.is_dir() and (candidate / "hydro_balance.c").is_file():
            return candidate
    pytest.skip("pinned visit-manager/VISITc source checkout is unavailable")


@pytest.fixture(scope="module")
def source_root():
    return _source_root()


@pytest.fixture(scope="module")
def hydrology_bridge(source_root, tmp_path_factory):
    return build_native_visitc_hydrology_bridge(
        source_root, tmp_path_factory.mktemp("visitc-native") / "hydrology"
    )


@pytest.fixture(scope="module")
def pm_bridge(source_root, tmp_path_factory):
    return build_native_visitc_pm_bridge(
        source_root, tmp_path_factory.mktemp("visitc-native") / "pm"
    )


@pytest.mark.parametrize(
    "state,forcing,site",
    (
        (
            VISITCHydrologyState(2.0, 50.0, 100.0),
            VISITCHydrologyForcing(5.0, 15.0, 8.0, 0.3, 0.2, 0.1, 0.4, 0.8, 0.3, 0.2),
            "",
        ),
        (
            VISITCHydrologyState(0.0, 50.0, 100.0),
            VISITCHydrologyForcing(0.0, -5.0, 0.0),
            "",
        ),
        (
            VISITCHydrologyState(0.0, 1.0, 100.0),
            VISITCHydrologyForcing(0.0, 10.0, 5.0, potential_transpiration_c3=100.0, potential_transpiration_c4=100.0),
            "QHB",
        ),
    ),
)
def test_native_hydrology_all_exposed_outputs_match_python(
    hydrology_bridge, source_root, state, forcing, site
):
    parameters = VISITCHydrologyParameters(
        field_capacity_upper=100.0,
        field_capacity_whole=200.0,
        hydraulic_conductivity=1.0e-8,
        tree_leaf_area_index=2.0,
        c3_leaf_area_index=1.0,
        c4_leaf_area_index=0.5,
        c3_understory_fraction=0.4,
        c4_understory_fraction=0.2,
        site_id=site,
    )
    comparison = compare_native_visitc_hydrology_to_python(
        hydrology_bridge, source_root, state, forcing, parameters
    )
    np.testing.assert_allclose(comparison.error, 0.0, atol=2e-13)
    assert comparison.native.source_commit == "5202debd96df6f88beb7d61f8688fff02ace964a"


@pytest.mark.parametrize("air_temperature,wind,conductance", ((15.0, 2.5, 120.0), (-5.0, 0.0, 0.0)))
def test_native_pm_and_radiation_match_python(
    pm_bridge, source_root, air_temperature, wind, conductance
):
    radiation = VISITCCanopyRadiationParameters(
        leaf_area_index=(2.0, 1.0, 0.5),
        extinction_initial=(0.5, 0.6, 0.7),
        extinction_radiation=(0.45, 0.55, 0.65),
        albedo=(0.12, 0.18, 0.20, 0.15),
        c3_understory_fraction=0.4,
        c4_understory_fraction=0.2,
    )
    environment = VISITCPenmanMonteithEnvironment(
        air_temperature=air_temperature,
        surface_temperature=air_temperature + 2.0,
        air_pressure=1000.0,
        vapor_pressure=5.0 if air_temperature < 0.0 else 12.0,
        vapor_pressure_deficit=8.0,
        wind_speed=wind,
        day_length=12.0,
        incoming_shortwave=220.0,
        cloud_fraction=0.4,
        canopy_conductance_tree=conductance,
        canopy_conductance_c3=conductance / 2.0,
        canopy_conductance_c4=conductance / 4.0,
        radiation=radiation,
    )
    comparison = compare_native_visitc_pm_to_python(
        pm_bridge,
        source_root,
        environment,
        upper_soil_water=50.0,
        field_capacity_upper=100.0,
    )
    np.testing.assert_allclose(comparison.pm_error, 0.0, atol=2e-13)
    np.testing.assert_allclose(comparison.radiation_error, 0.0, atol=3e-13)


def test_native_canopy_conductance_matches_python(pm_bridge):
    parameters = VISITCCanopyConductanceParameters(
        photosynthetic_capacity=15.0,
        radiation_extinction=0.5,
        light_use_efficiency=0.05,
        canopy_top_ppfd=1000.0,
        leaf_area_index=2.0,
        atmospheric_co2=410.0,
        co2_compensation_point=40.0,
        vapor_pressure_deficit=10.0,
        minimum_stomatal_conductance=10.0,
        ball_berry_slope=9.0,
        vpd_scale=15.0,
    )
    assert compare_native_visitc_canopy_conductance_to_python(
        pm_bridge, parameters
    ) == pytest.approx(0.0, abs=1e-14)


def test_native_lai_and_extinction_dependencies_match_python(pm_bridge):
    assert compare_native_visitc_lai_to_python(
        pm_bridge, 2.0, 20.0
    ) == pytest.approx(0.0, abs=1e-14)
    assert compare_native_visitc_extinction_to_python(
        pm_bridge, 0.5, 35.0
    ) == pytest.approx(0.0, abs=1e-14)
