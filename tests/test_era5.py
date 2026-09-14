from datetime import UTC, date, datetime, timedelta
import hashlib
from io import BytesIO
import json
from pathlib import Path
from zipfile import ZipFile

import h5py
import numpy as np
import pytest

from control_carbon.era5 import (
    CDS_ERA5_HOURLY_VARIABLES,
    ERA5DailyPointData,
    ERA5HourlyPointData,
    ERA5PointRequest,
    VISIT_WEATHER_SOURCES,
    aggregate_era5_hourly_for_visit,
    build_cds_era5_point_requests,
    decode_cds_era5_archives,
    format_visit_ncep_site_climate,
    nearest_era5_grid_point,
    retrieve_cds_era5_point,
    write_visit_ncep_site_climate,
)


def hourly_data(n_hours=24):
    start = datetime(2020, 1, 1)
    times = tuple(
        start + timedelta(hours=hour)
        for hour in range(n_hours)
    )
    values = {
        "temperature_2m": np.arange(n_hours, dtype=float),
        "dew_point_2m": np.zeros(n_hours),
        "precipitation": np.ones(n_hours),
        "shortwave_radiation": np.full(n_hours, 100.0),
        "surface_pressure": np.full(n_hours, 1000.0),
        "cloud_cover": np.full(n_hours, 50.0),
        "wind_speed_10m": np.full(n_hours, 2.0),
        "wind_direction_10m": np.full(n_hours, 90.0),
        "soil_temperature_0_to_7cm": np.full(n_hours, 10.0),
        "soil_temperature_7_to_28cm": np.full(n_hours, 20.0),
        "soil_temperature_28_to_100cm": np.full(n_hours, 30.0),
        "soil_temperature_100_to_255cm": np.full(n_hours, 40.0),
        "soil_moisture_0_to_7cm": np.full(n_hours, 0.1),
        "soil_moisture_7_to_28cm": np.full(n_hours, 0.2),
        "soil_moisture_28_to_100cm": np.full(n_hours, 0.3),
        "soil_moisture_100_to_255cm": np.full(n_hours, 0.4),
        "snow_depth": np.full(n_hours, np.nan),
        "surface_temperature": np.full(n_hours, 12.0),
    }
    request = ERA5PointRequest(
        latitude=36.146,
        longitude=137.423,
        start_date=date(2020, 1, 1),
        end_date=date(2020, 1, 1),
    )
    return ERA5HourlyPointData(
        request=request,
        times=times,
        variables=values,
        units={name: "test-unit" for name in values},
        provider="synthetic fixture",
        model="ERA5",
        source_url="fixture",
        resolved_latitude=36.25,
        resolved_longitude=137.5,
        grid_elevation_m=None,
        utc_offset_seconds=0,
    )


def test_daily_aggregation_matches_hand_computable_values():
    daily = aggregate_era5_hourly_for_visit(hourly_data())
    vapour_pressure_hpa = 6.1078
    expected_specific_humidity = 0.622 * vapour_pressure_hpa / (
        1000.0 - 0.378 * vapour_pressure_hpa
    )

    assert daily.dates == (date(2020, 1, 1),)
    assert daily.variables["air_temperature_mean_k"] == pytest.approx([284.65])
    assert daily.variables["air_temperature_max_k"] == pytest.approx([296.15])
    assert daily.variables["air_temperature_min_k"] == pytest.approx([273.15])
    assert daily.variables["specific_humidity_kg_kg"] == pytest.approx(
        [expected_specific_humidity]
    )
    assert daily.variables["precipitation_rate_kg_m2_s"] == pytest.approx(
        [24.0 / 86400.0]
    )
    assert daily.variables["surface_temperature_k"] == pytest.approx([285.15])
    assert daily.variables["soil_temperature_0_10cm_k"] == pytest.approx(
        [286.15]
    )
    assert daily.variables["soil_temperature_10_200cm_k"] == pytest.approx(
        [273.15 + (18.0 * 20.0 + 72.0 * 30.0 + 100.0 * 40.0) / 190.0]
    )
    assert daily.variables["soil_water_0_10cm_m3_m3"] == pytest.approx([0.13])
    assert daily.variables["soil_water_10_200cm_m3_m3"] == pytest.approx(
        [(18.0 * 0.2 + 72.0 * 0.3 + 100.0 * 0.4) / 190.0]
    )
    assert daily.variables["u_wind_10m_m_s"] == pytest.approx([-2.0])
    assert daily.variables["v_wind_10m_m_s"] == pytest.approx([0.0], abs=1e-15)
    assert daily.variables["snow_depth_m"] == pytest.approx([0.0])
    assert daily.variables["surface_pressure_pa"] == pytest.approx([100000.0])


def test_daily_aggregation_rejects_missing_utc_hour():
    with pytest.raises(ValueError, match="every UTC hour"):
        aggregate_era5_hourly_for_visit(hourly_data(n_hours=23))


def test_visit_writer_has_source_column_order_and_complete_year_guard(tmp_path):
    daily = aggregate_era5_hourly_for_visit(hourly_data())

    text = format_visit_ncep_site_climate(daily, require_complete_years=False)
    fields = text.split()

    assert len(fields) == 19
    assert fields[:2] == ["2020", "1"]
    assert float(fields[2]) == pytest.approx(284.65)
    assert float(fields[6]) * 86400.0 == pytest.approx(24.0)
    assert float(fields[-1]) == pytest.approx(100000.0)
    with pytest.raises(ValueError, match="complete calendar years"):
        format_visit_ncep_site_climate(daily)

    manifest_path = tmp_path / "climate.manifest.json"
    write_visit_ncep_site_climate(
        daily,
        tmp_path / "climate.txt",
        manifest_path=manifest_path,
        require_complete_years=False,
    )
    visit_manifest = json.loads(manifest_path.read_text())["visit_output"]
    assert "init_site.c::f_init_site" in visit_manifest["source_compatibility_warning"]
    assert "location_proc.c::f_loct_proc" in visit_manifest[
        "source_compatibility_warning"
    ]
    assert {source.path for source in VISIT_WEATHER_SOURCES} == {
        "visit_local/init_site.c",
        "visit_local/location_proc.c",
        "visit_local/experiment.c",
    }
    assert all(source.commit is not None for source in VISIT_WEATHER_SOURCES)


def test_cds_requests_split_months_without_cartesian_date_expansion():
    request = ERA5PointRequest(
        36.146, 137.423, date(2020, 1, 30), date(2020, 2, 2)
    )

    jobs = build_cds_era5_point_requests(request)

    assert len(jobs) == 2
    assert jobs[0].request["month"] == ["01"]
    assert jobs[0].request["day"] == ["30", "31"]
    assert jobs[1].request["month"] == ["02"]
    assert jobs[1].request["day"] == ["01", "02"]
    assert jobs[0].request["variable"] == list(CDS_ERA5_HOURLY_VARIABLES)
    assert len(jobs[0].request["time"]) == 24
    assert jobs[0].request["area"] == [36.25, 137.5, 36.25, 137.5]
    assert jobs[0].request["download_format"] == "zip"
    assert jobs[0].target_name.endswith("20200130_20200131.zip")
    assert jobs[1].target_name.endswith("20200201_20200202.zip")


def test_nearest_era5_grid_point_uses_half_up_ties_and_validates_bounds():
    assert nearest_era5_grid_point(36.125, 137.375) == (36.25, 137.5)
    assert nearest_era5_grid_point(-36.125, -137.375) == (-36.0, -137.25)
    assert nearest_era5_grid_point(90.0, 180.0) == (90.0, 180.0)
    with pytest.raises(ValueError, match="latitude"):
        nearest_era5_grid_point(90.1, 0.0)


def test_cds_archive_decoder_converts_si_units_and_uses_skin_temperature(tmp_path):
    instant_values = {
        "t2m": 280.0,
        "d2m": 275.0,
        "tcc": 0.5,
        "skt": 285.0,
        "stl1": 281.0,
        "stl2": 282.0,
        "stl3": 283.0,
        "stl4": 284.0,
        "swvl1": 0.1,
        "swvl2": 0.2,
        "swvl3": 0.3,
        "swvl4": 0.4,
        "u10": -2.0,
        "v10": 0.0,
        "sd": 0.1,
        "sp": 90000.0,
    }
    archive_path = tmp_path / "era5.zip"
    with ZipFile(archive_path, "w") as archive:
        archive.writestr("instant.nc", cds_netcdf_bytes(instant_values))
        archive.writestr("accum.nc", cds_netcdf_bytes({"tp": 0.001, "ssrd": 360000.0}))
    request = ERA5PointRequest(
        36.146, 137.423, date(2020, 1, 1), date(2020, 1, 1)
    )

    hourly = decode_cds_era5_archives(request, [archive_path])
    daily = aggregate_era5_hourly_for_visit(hourly)

    assert len(hourly.times) == 24
    assert hourly.resolved_latitude == 36.25
    assert hourly.resolved_longitude == 137.5
    assert hourly.variables["temperature_2m"] == pytest.approx([6.85] * 24)
    assert hourly.variables["precipitation"] == pytest.approx([1.0] * 24)
    assert hourly.variables["shortwave_radiation"] == pytest.approx([100.0] * 24)
    assert hourly.variables["surface_pressure"] == pytest.approx([900.0] * 24)
    assert hourly.variables["cloud_cover"] == pytest.approx([50.0] * 24)
    assert hourly.variables["wind_speed_10m"] == pytest.approx([2.0] * 24)
    assert hourly.variables["wind_direction_10m"] == pytest.approx([90.0] * 24)
    assert daily.variables["surface_temperature_k"] == pytest.approx([285.0])
    assert daily.variables["precipitation_rate_kg_m2_s"] == pytest.approx(
        [24.0 / 86400.0]
    )
    artifact = hourly.to_manifest()["source_artifacts"][0]
    assert artifact["size_bytes"] == archive_path.stat().st_size
    assert len(artifact["sha256"]) == 64


def test_cds_retrieval_reuses_valid_exact_range_archive(tmp_path):
    source_archive = tmp_path / "source.zip"
    with ZipFile(source_archive, "w") as archive:
        archive.writestr("data.nc", b"netcdf fixture")

    class RecordingClient:
        def __init__(self):
            self.calls = 0

        def retrieve(self, dataset, request, target):
            self.calls += 1
            Path(target).write_bytes(source_archive.read_bytes())

    client = RecordingClient()
    request = ERA5PointRequest(
        36.146, 137.423, date(2020, 1, 1), date(2020, 1, 1)
    )

    first = retrieve_cds_era5_point(request, tmp_path / "cache", client=client)
    second = retrieve_cds_era5_point(request, tmp_path / "cache", client=client)
    retrieve_cds_era5_point(
        request, tmp_path / "cache", client=client, use_cache=False
    )

    assert first == second
    assert client.calls == 2
    assert first[0].name.endswith("20200101_20200101.zip")


def test_cds_decoder_joins_month_archives_in_time_order(tmp_path):
    values = {
        "t2m": 280.0,
        "d2m": 275.0,
        "tcc": 0.5,
        "skt": 285.0,
        "stl1": 281.0,
        "stl2": 282.0,
        "stl3": 283.0,
        "stl4": 284.0,
        "swvl1": 0.1,
        "swvl2": 0.2,
        "swvl3": 0.3,
        "swvl4": 0.4,
        "u10": -2.0,
        "v10": 0.0,
        "sd": 0.1,
        "sp": 90000.0,
        "tp": 0.001,
        "ssrd": 360000.0,
    }
    first = tmp_path / "january.zip"
    second = tmp_path / "february.zip"
    with ZipFile(first, "w") as archive:
        archive.writestr(
            "all.nc", cds_netcdf_bytes(values, start=datetime(2020, 1, 31))
        )
    with ZipFile(second, "w") as archive:
        archive.writestr(
            "all.nc", cds_netcdf_bytes(values, start=datetime(2020, 2, 1))
        )
    request = ERA5PointRequest(
        36.146, 137.423, date(2020, 1, 31), date(2020, 2, 1)
    )

    hourly = decode_cds_era5_archives(request, [second, first])
    daily = aggregate_era5_hourly_for_visit(hourly)

    assert len(hourly.times) == 48
    assert hourly.times[0] == datetime(2020, 1, 31)
    assert hourly.times[-1] == datetime(2020, 2, 1, 23)
    assert daily.dates == (date(2020, 1, 31), date(2020, 2, 1))
    assert len(hourly.to_manifest()["source_artifacts"]) == 2


def test_complete_leap_year_output_records_coverage_and_checksum(tmp_path):
    one_day = aggregate_era5_hourly_for_visit(hourly_data())
    dates = tuple(date(2020, 1, 1) + timedelta(days=offset) for offset in range(366))
    daily = ERA5DailyPointData(
        request=ERA5PointRequest(36.146, 137.423, dates[0], dates[-1]),
        dates=dates,
        variables={
            name: np.repeat(values, len(dates))
            for name, values in one_day.variables.items()
        },
        units=one_day.units,
        source_manifest=one_day.source_manifest,
        transformations=one_day.transformations,
    )
    output_path = tmp_path / "visit-2020.txt"
    manifest_path = tmp_path / "visit-2020.manifest.json"

    write_visit_ncep_site_climate(
        daily, output_path, manifest_path=manifest_path
    )

    rows = output_path.read_text().splitlines()
    manifest = json.loads(manifest_path.read_text())
    assert len(rows) == 366
    assert rows[59].split()[:2] == ["2020", "60"]
    assert rows[-1].split()[:2] == ["2020", "366"]
    assert manifest["calendar_years"] == [
        {"year": 2020, "day_count": 366, "gregorian_leap_year": True}
    ]
    assert manifest["visit_output"]["artifact"]["sha256"] == hashlib.sha256(
        output_path.read_bytes()
    ).hexdigest()


def cds_netcdf_bytes(values, *, start=datetime(2020, 1, 1), n_hours=24):
    buffer = BytesIO()
    with h5py.File(buffer, "w") as dataset:
        start_timestamp = int(start.replace(tzinfo=UTC).timestamp())
        valid_time = dataset.create_dataset(
            "valid_time",
            data=np.arange(
                start_timestamp, start_timestamp + n_hours * 3600, 3600
            ),
        )
        valid_time.attrs["units"] = "seconds since 1970-01-01"
        dataset.create_dataset("latitude", data=[36.25])
        dataset.create_dataset("longitude", data=[137.5])
        for name, value in values.items():
            dataset.create_dataset(name, data=np.full((n_hours, 1, 1), value))
    return buffer.getvalue()