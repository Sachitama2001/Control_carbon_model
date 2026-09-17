"""Build native VISITc bridges and report one-day Python/C agreement."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import tempfile

import numpy as np

from control_carbon import (
    VISITCHydrologyForcing,
    VISITCHydrologyParameters,
    VISITCHydrologyState,
    VISITCCanopyConductanceParameters,
    VISITCCanopyRadiationParameters,
    VISITCPenmanMonteithEnvironment,
    build_native_visitc_hydrology_bridge,
    build_native_visitc_pm_bridge,
    compare_native_visitc_canopy_conductance_to_python,
    compare_native_visitc_hydrology_to_python,
    compare_native_visitc_extinction_to_python,
    compare_native_visitc_lai_to_python,
    compare_native_visitc_pm_to_python,
)


def run(source_root: Path) -> dict[str, object]:
    state = VISITCHydrologyState(2.0, 50.0, 100.0)
    forcing = VISITCHydrologyForcing(
        precipitation=5.0,
        air_temperature=15.0,
        deep_soil_temperature=8.0,
        potential_interception_tree=0.3,
        potential_interception_c3=0.2,
        potential_interception_c4=0.1,
        potential_soil_evaporation=0.4,
        potential_transpiration_tree=0.8,
        potential_transpiration_c3=0.3,
        potential_transpiration_c4=0.2,
    )
    hydrology_parameters = VISITCHydrologyParameters(
        field_capacity_upper=100.0,
        field_capacity_whole=200.0,
        hydraulic_conductivity=1.0e-8,
        tree_leaf_area_index=2.0,
        c3_leaf_area_index=1.0,
        c4_leaf_area_index=0.5,
        c3_understory_fraction=0.4,
        c4_understory_fraction=0.2,
    )
    radiation = VISITCCanopyRadiationParameters(
        leaf_area_index=(2.0, 1.0, 0.5),
        extinction_initial=(0.5, 0.6, 0.7),
        extinction_radiation=(0.45, 0.55, 0.65),
        albedo=(0.12, 0.18, 0.20, 0.15),
        c3_understory_fraction=0.4,
        c4_understory_fraction=0.2,
    )
    pm_environment = VISITCPenmanMonteithEnvironment(
        air_temperature=15.0,
        surface_temperature=17.0,
        air_pressure=1000.0,
        vapor_pressure=12.0,
        vapor_pressure_deficit=8.0,
        wind_speed=2.5,
        day_length=12.0,
        incoming_shortwave=220.0,
        cloud_fraction=0.4,
        canopy_conductance_tree=120.0,
        canopy_conductance_c3=60.0,
        canopy_conductance_c4=30.0,
        radiation=radiation,
    )
    conductance_parameters = VISITCCanopyConductanceParameters(
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
    with tempfile.TemporaryDirectory(prefix="visitc-bridge-") as temporary:
        temporary_path = Path(temporary)
        hydrology_executable = build_native_visitc_hydrology_bridge(
            source_root, temporary_path / "hydrology"
        )
        pm_executable = build_native_visitc_pm_bridge(
            source_root, temporary_path / "pm"
        )
        hydrology = compare_native_visitc_hydrology_to_python(
            hydrology_executable,
            source_root,
            state,
            forcing,
            hydrology_parameters,
        )
        pm = compare_native_visitc_pm_to_python(
            pm_executable,
            source_root,
            pm_environment,
            upper_soil_water=state.upper_soil_water,
            field_capacity_upper=hydrology_parameters.field_capacity_upper,
        )
        conductance_error = compare_native_visitc_canopy_conductance_to_python(
            pm_executable, conductance_parameters
        )
        lai_error = compare_native_visitc_lai_to_python(pm_executable, 2.0, 20.0)
        extinction_error = compare_native_visitc_extinction_to_python(
            pm_executable, 0.5, 35.0
        )
    return {
        "source_commit": hydrology.native.source_commit,
        "hydrology_max_abs_error": float(np.max(np.abs(hydrology.error))),
        "pm_max_abs_error": float(np.max(np.abs(pm.pm_error))),
        "radiation_max_abs_error": float(np.max(np.abs(pm.radiation_error))),
        "canopy_conductance_abs_error": abs(conductance_error),
        "lai_abs_error": abs(lai_error),
        "irradiance_extinction_abs_error": abs(extinction_error),
        "hydrology_source_sha256": hydrology.native.source_sha256,
        "pm_source_sha256": pm.source_sha256,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "source_root",
        type=Path,
        help="path to visit-manager/VISITc/point at the pinned commit",
    )
    print(json.dumps(run(parser.parse_args().source_root), indent=2))


if __name__ == "__main__":
    main()
