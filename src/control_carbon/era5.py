"""ERA5 point-data acquisition with explicit source metadata."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import UTC, date, datetime, time, timedelta
from io import BytesIO
import calendar
import hashlib
import json
from pathlib import Path
from typing import Any, Final, Mapping
from zipfile import ZipFile

import numpy as np

from .provenance import SourceRef
from .visit_source_map import VISIT_SOURCE_COMMIT, VISIT_SOURCE_REPOSITORY


ERA5_REQUIRED_HOURLY_VARIABLES: Final[tuple[str, ...]] = (
    "temperature_2m",
    "dew_point_2m",
    "precipitation",
    "shortwave_radiation",
    "surface_pressure",
    "cloud_cover",
    "wind_speed_10m",
    "wind_direction_10m",
    "soil_temperature_0_to_7cm",
    "soil_temperature_7_to_28cm",
    "soil_temperature_28_to_100cm",
    "soil_temperature_100_to_255cm",
    "soil_moisture_0_to_7cm",
    "soil_moisture_7_to_28cm",
    "soil_moisture_28_to_100cm",
    "soil_moisture_100_to_255cm",
    "snow_depth",
    "surface_temperature",
)
CDS_ERA5_SINGLE_LEVELS_DATASET: Final = "reanalysis-era5-single-levels"
CDS_ERA5_HOURLY_VARIABLES: Final[tuple[str, ...]] = (
    "2m_temperature",
    "2m_dewpoint_temperature",
    "total_precipitation",
    "surface_solar_radiation_downwards",
    "total_cloud_cover",
    "skin_temperature",
    "soil_temperature_level_1",
    "soil_temperature_level_2",
    "soil_temperature_level_3",
    "soil_temperature_level_4",
    "volumetric_soil_water_layer_1",
    "volumetric_soil_water_layer_2",
    "volumetric_soil_water_layer_3",
    "volumetric_soil_water_layer_4",
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "snow_depth",
    "surface_pressure",
)
VISIT_NCEP_DAILY_VARIABLES: Final[tuple[str, ...]] = (
    "air_temperature_mean_k",
    "air_temperature_max_k",
    "air_temperature_min_k",
    "specific_humidity_kg_kg",
    "precipitation_rate_kg_m2_s",
    "downward_shortwave_w_m2",
    "cloud_cover_percent",
    "surface_temperature_k",
    "soil_temperature_0_10cm_k",
    "soil_temperature_10_200cm_k",
    "soil_temperature_300cm_k",
    "soil_water_0_10cm_m3_m3",
    "soil_water_10_200cm_m3_m3",
    "u_wind_10m_m_s",
    "v_wind_10m_m_s",
    "snow_depth_m",
    "surface_pressure_pa",
)
VISIT_NCEP_COLUMNS: Final[tuple[str, ...]] = (
    "year",
    "day_of_year_1_based",
    *VISIT_NCEP_DAILY_VARIABLES,
)
VISIT_WEATHER_SOURCES: Final[tuple[SourceRef, ...]] = (
    SourceRef(
        path="visit_local/init_site.c",
        symbol="f_init_site",
        role="site climate format selection, column parsing, and unit conversion",
        repository=VISIT_SOURCE_REPOSITORY,
        commit=VISIT_SOURCE_COMMIT,
        function="f_init_site",
        approximation_level="native-source",
    ),
    SourceRef(
        path="visit_local/location_proc.c",
        symbol="f_loct_proc",
        role="daily climate selection, overrides, and atmospheric derivations",
        repository=VISIT_SOURCE_REPOSITORY,
        commit=VISIT_SOURCE_COMMIT,
        function="f_loct_proc",
        approximation_level="native-source",
    ),
    SourceRef(
        path="visit_local/experiment.c",
        symbol="f_experiment",
        role="daily experiment calendar and climate-array indexing",
        repository=VISIT_SOURCE_REPOSITORY,
        commit=VISIT_SOURCE_COMMIT,
        function="f_experiment",
        approximation_level="native-source",
    ),
)


@dataclass(frozen=True)
class ERA5PointRequest:
    """A UTC, inclusive-date request for one ERA5 grid point."""

    latitude: float
    longitude: float
    start_date: date
    end_date: date

    def __post_init__(self) -> None:
        if not np.isfinite(self.latitude) or not -90.0 <= self.latitude <= 90.0:
            raise ValueError("latitude must lie in [-90, 90]")
        if not np.isfinite(self.longitude) or not -180.0 <= self.longitude <= 180.0:
            raise ValueError("longitude must lie in [-180, 180]")
        if self.end_date < self.start_date:
            raise ValueError("end_date must not precede start_date")


@dataclass(frozen=True)
class CDSERA5MonthlyRequest:
    """One calendar month's direct Copernicus CDS retrieval job."""

    dataset: str
    request: Mapping[str, Any]
    target_name: str


def build_cds_era5_point_requests(
    point_request: ERA5PointRequest,
) -> tuple[CDSERA5MonthlyRequest, ...]:
    """Build monthly direct-CDS jobs without selecting dates outside the range."""
    requests = []
    grid_latitude, grid_longitude = nearest_era5_grid_point(
        point_request.latitude, point_request.longitude
    )
    month_start = date(point_request.start_date.year, point_request.start_date.month, 1)
    while month_start <= point_request.end_date:
        month_end = date(
            month_start.year,
            month_start.month,
            calendar.monthrange(month_start.year, month_start.month)[1],
        )
        selected_start = max(point_request.start_date, month_start)
        selected_end = min(point_request.end_date, month_end)
        days = [
            f"{day:02d}"
            for day in range(selected_start.day, selected_end.day + 1)
        ]
        request = {
            "product_type": ["reanalysis"],
            "variable": list(CDS_ERA5_HOURLY_VARIABLES),
            "year": [f"{month_start.year:04d}"],
            "month": [f"{month_start.month:02d}"],
            "day": days,
            "time": [f"{hour:02d}:00" for hour in range(24)],
            "data_format": "netcdf",
            "download_format": "zip",
            "area": [
                grid_latitude,
                grid_longitude,
                grid_latitude,
                grid_longitude,
            ],
        }
        requests.append(
            CDSERA5MonthlyRequest(
                dataset=CDS_ERA5_SINGLE_LEVELS_DATASET,
                request=request,
                target_name=(
                    f"era5_point_{point_request.latitude:+.3f}_"
                    f"{point_request.longitude:+.3f}_"
                    f"{selected_start:%Y%m%d}_{selected_end:%Y%m%d}.zip"
                ),
            )
        )
        if month_start.month == 12:
            month_start = date(month_start.year + 1, 1, 1)
        else:
            month_start = date(month_start.year, month_start.month + 1, 1)
    return tuple(requests)


def nearest_era5_grid_point(latitude: float, longitude: float) -> tuple[float, float]:
    """Return the nearest point on the regular 0.25-degree ERA5 CDS grid."""
    if not np.isfinite(latitude) or not -90.0 <= latitude <= 90.0:
        raise ValueError("latitude must lie in [-90, 90]")
    if not np.isfinite(longitude) or not -180.0 <= longitude <= 180.0:
        raise ValueError("longitude must lie in [-180, 180]")

    spacing = 0.25
    latitude_index = int(np.floor((latitude + 90.0) / spacing + 0.5))
    longitude_index = int(np.floor((longitude + 180.0) / spacing + 0.5))
    grid_latitude = min(90.0, max(-90.0, -90.0 + spacing * latitude_index))
    grid_longitude = min(
        180.0, max(-180.0, -180.0 + spacing * longitude_index)
    )
    return grid_latitude, grid_longitude


def retrieve_cds_era5_point(
    point_request: ERA5PointRequest,
    output_directory: str | Path,
    *,
    client: Any | None = None,
    use_cache: bool = True,
) -> tuple[Path, ...]:
    """Run official CDS jobs using credentials managed by ``cdsapi``."""
    retrieval_client = client
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    targets = []
    for job in build_cds_era5_point_requests(point_request):
        target = output_directory / job.target_name
        if not use_cache or not _is_valid_cds_archive(target):
            if retrieval_client is None:
                try:
                    import cdsapi
                except ImportError as error:
                    raise RuntimeError(
                        "direct CDS retrieval requires the optional 'cdsapi' package"
                    ) from error
                retrieval_client = cdsapi.Client()
            temporary_target = target.with_suffix(target.suffix + ".part")
            temporary_target.unlink(missing_ok=True)
            try:
                retrieval_client.retrieve(
                    job.dataset, dict(job.request), str(temporary_target)
                )
                if not _is_valid_cds_archive(temporary_target):
                    raise ValueError(
                        f"CDS retrieval did not produce a valid NetCDF ZIP: {target}"
                    )
                temporary_target.replace(target)
            finally:
                temporary_target.unlink(missing_ok=True)
        targets.append(target)
    return tuple(targets)


def _is_valid_cds_archive(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        with ZipFile(path) as archive:
            members = archive.namelist()
            return (
                archive.testzip() is None
                and any(name.lower().endswith(".nc") for name in members)
            )
    except (OSError, ValueError):
        return False


def decode_cds_era5_archives(
    point_request: ERA5PointRequest,
    archive_paths: tuple[str | Path, ...] | list[str | Path],
) -> ERA5HourlyPointData:
    """Decode monthly CDS ZIP archives into the canonical hourly point schema."""
    try:
        import h5py
    except ImportError as error:
        raise RuntimeError(
            "decoding CDS NetCDF archives requires the optional 'h5py' package"
        ) from error
    if not archive_paths:
        raise ValueError("at least one CDS archive path is required")

    decoded_parts: list[tuple[tuple[datetime, ...], dict[str, np.ndarray]]] = []
    resolved_point: tuple[float, float] | None = None
    source_paths: list[str] = []
    source_artifacts: list[dict[str, Any]] = []
    for archive_path_value in archive_paths:
        archive_path = Path(archive_path_value)
        source_paths.append(str(archive_path.resolve()))
        source_artifacts.append(_file_artifact(archive_path))
        archive_variables: dict[str, np.ndarray] = {}
        archive_times: tuple[datetime, ...] | None = None
        with ZipFile(archive_path) as archive:
            netcdf_members = [
                name for name in archive.namelist() if name.lower().endswith(".nc")
            ]
            if not netcdf_members:
                raise ValueError(f"CDS archive contains no NetCDF members: {archive_path}")
            for member_name in netcdf_members:
                with h5py.File(BytesIO(archive.read(member_name)), "r") as dataset:
                    member_times = _decode_cds_valid_times(dataset)
                    if archive_times is None:
                        archive_times = member_times
                    elif archive_times != member_times:
                        raise ValueError(
                            f"CDS NetCDF members have different time axes: {archive_path}"
                        )
                    member_point = _decode_cds_grid_point(dataset)
                    if resolved_point is None:
                        resolved_point = member_point
                    elif resolved_point != member_point:
                        raise ValueError("CDS archives contain different grid points")
                    for name in _CDS_SHORT_NAME_TO_CANONICAL:
                        if name in dataset:
                            if name in archive_variables:
                                raise ValueError(f"duplicate CDS variable {name}")
                            archive_variables[name] = _read_cds_point_variable(
                                dataset[name], len(member_times), name
                            )
        if archive_times is None:
            raise ValueError(f"CDS archive contains no time axis: {archive_path}")
        missing = set(_CDS_SHORT_NAME_TO_CANONICAL) - set(archive_variables)
        if missing:
            raise ValueError(
                "CDS archive is missing variables: " + ", ".join(sorted(missing))
            )
        decoded_parts.append((archive_times, archive_variables))

    decoded_parts.sort(key=lambda part: part[0][0])
    times = tuple(timestamp for part_times, _ in decoded_parts for timestamp in part_times)
    if len(set(times)) != len(times):
        raise ValueError("CDS archives contain duplicate timestamps")
    raw = {
        name: np.concatenate([part_variables[name] for _, part_variables in decoded_parts])
        for name in _CDS_SHORT_NAME_TO_CANONICAL
    }
    u_wind = raw["u10"]
    v_wind = raw["v10"]
    variables = {
        "temperature_2m": raw["t2m"] - 273.15,
        "dew_point_2m": raw["d2m"] - 273.15,
        "precipitation": raw["tp"] * 1000.0,
        "shortwave_radiation": raw["ssrd"] / 3600.0,
        "surface_pressure": raw["sp"] / 100.0,
        "cloud_cover": raw["tcc"] * 100.0,
        "wind_speed_10m": np.hypot(u_wind, v_wind),
        "wind_direction_10m": (
            np.rad2deg(np.arctan2(-u_wind, -v_wind)) % 360.0
        ),
        "soil_temperature_0_to_7cm": raw["stl1"] - 273.15,
        "soil_temperature_7_to_28cm": raw["stl2"] - 273.15,
        "soil_temperature_28_to_100cm": raw["stl3"] - 273.15,
        "soil_temperature_100_to_255cm": raw["stl4"] - 273.15,
        "soil_moisture_0_to_7cm": raw["swvl1"],
        "soil_moisture_7_to_28cm": raw["swvl2"],
        "soil_moisture_28_to_100cm": raw["swvl3"],
        "soil_moisture_100_to_255cm": raw["swvl4"],
        "snow_depth": raw["sd"],
        "surface_temperature": raw["skt"] - 273.15,
    }
    units = {
        "temperature_2m": "°C",
        "dew_point_2m": "°C",
        "precipitation": "mm",
        "shortwave_radiation": "W/m²",
        "surface_pressure": "hPa",
        "cloud_cover": "%",
        "wind_speed_10m": "m/s",
        "wind_direction_10m": "°",
        "soil_temperature_0_to_7cm": "°C",
        "soil_temperature_7_to_28cm": "°C",
        "soil_temperature_28_to_100cm": "°C",
        "soil_temperature_100_to_255cm": "°C",
        "soil_moisture_0_to_7cm": "m³/m³",
        "soil_moisture_7_to_28cm": "m³/m³",
        "soil_moisture_28_to_100cm": "m³/m³",
        "soil_moisture_100_to_255cm": "m³/m³",
        "snow_depth": "m",
        "surface_temperature": "°C",
    }
    if resolved_point is None:
        raise ValueError("CDS archives contain no grid coordinates")
    return ERA5HourlyPointData(
        request=point_request,
        times=times,
        variables=variables,
        units=units,
        provider="Copernicus Climate Data Store",
        model="ERA5",
        source_url=";".join(source_paths),
        resolved_latitude=resolved_point[0],
        resolved_longitude=resolved_point[1],
        grid_elevation_m=None,
        utc_offset_seconds=0,
        source_artifacts=tuple(source_artifacts),
    )


def fetch_cds_era5_point(
    point_request: ERA5PointRequest,
    archive_directory: str | Path,
    *,
    client: Any | None = None,
    use_cache: bool = True,
) -> ERA5HourlyPointData:
    """Retrieve official CDS archives and decode them to canonical hourly data."""
    archive_paths = retrieve_cds_era5_point(
        point_request, archive_directory, client=client, use_cache=use_cache
    )
    return decode_cds_era5_archives(point_request, archive_paths)


_CDS_SHORT_NAME_TO_CANONICAL: Final[Mapping[str, str]] = {
    "t2m": "temperature_2m",
    "d2m": "dew_point_2m",
    "tp": "precipitation",
    "ssrd": "shortwave_radiation",
    "tcc": "cloud_cover",
    "skt": "surface_temperature",
    "stl1": "soil_temperature_0_to_7cm",
    "stl2": "soil_temperature_7_to_28cm",
    "stl3": "soil_temperature_28_to_100cm",
    "stl4": "soil_temperature_100_to_255cm",
    "swvl1": "soil_moisture_0_to_7cm",
    "swvl2": "soil_moisture_7_to_28cm",
    "swvl3": "soil_moisture_28_to_100cm",
    "swvl4": "soil_moisture_100_to_255cm",
    "u10": "wind_speed_10m",
    "v10": "wind_direction_10m",
    "sd": "snow_depth",
    "sp": "surface_pressure",
}


def _decode_cds_valid_times(dataset: Any) -> tuple[datetime, ...]:
    if "valid_time" not in dataset:
        raise ValueError("CDS NetCDF member is missing valid_time")
    units = _decode_hdf5_attribute(dataset["valid_time"].attrs.get("units"))
    if units != "seconds since 1970-01-01":
        raise ValueError(f"unsupported CDS time units: {units}")
    return tuple(
        datetime.fromtimestamp(int(value), tz=UTC).replace(tzinfo=None)
        for value in dataset["valid_time"][...]
    )


def _decode_cds_grid_point(dataset: Any) -> tuple[float, float]:
    if "latitude" not in dataset or "longitude" not in dataset:
        raise ValueError("CDS NetCDF member is missing point coordinates")
    latitudes = np.asarray(dataset["latitude"][...], dtype=float).reshape(-1)
    longitudes = np.asarray(dataset["longitude"][...], dtype=float).reshape(-1)
    if latitudes.size != 1 or longitudes.size != 1:
        raise ValueError("CDS decoder requires exactly one grid point")
    return float(latitudes[0]), float(longitudes[0])


def _read_cds_point_variable(variable: Any, n_times: int, name: str) -> np.ndarray:
    values = np.asarray(variable[...], dtype=float)
    if values.shape != (n_times, 1, 1):
        raise ValueError(
            f"CDS variable {name} must have shape ({n_times}, 1, 1), got {values.shape}"
        )
    return values[:, 0, 0]


def _decode_hdf5_attribute(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("ascii")
    if isinstance(value, np.ndarray) and value.size == 1:
        return _decode_hdf5_attribute(value.reshape(-1)[0])
    return str(value)


@dataclass(frozen=True)
class ERA5HourlyPointData:
    """Hourly ERA5 values and the metadata needed for reproducibility."""

    request: ERA5PointRequest
    times: tuple[datetime, ...]
    variables: Mapping[str, np.ndarray]
    units: Mapping[str, str]
    provider: str
    model: str
    source_url: str
    resolved_latitude: float
    resolved_longitude: float
    grid_elevation_m: float | None
    utc_offset_seconds: int
    source_artifacts: tuple[Mapping[str, Any], ...] = ()

    def __post_init__(self) -> None:
        if not self.times:
            raise ValueError("ERA5 response contains no hourly timestamps")
        expected_size = len(self.times)
        for name in ERA5_REQUIRED_HOURLY_VARIABLES:
            if name not in self.variables:
                raise ValueError(f"ERA5 response is missing {name}")
            if self.variables[name].shape != (expected_size,):
                raise ValueError(f"ERA5 variable {name} has the wrong length")
        if self.utc_offset_seconds != 0:
            raise ValueError("ERA5 hourly data must use UTC")

    def to_manifest(self) -> dict[str, Any]:
        """Return JSON-compatible acquisition metadata without data arrays."""
        return {
            "provider": self.provider,
            "model": self.model,
            "source_url": self.source_url,
            "requested_point": {
                "latitude": self.request.latitude,
                "longitude": self.request.longitude,
            },
            "resolved_grid_point": {
                "latitude": self.resolved_latitude,
                "longitude": self.resolved_longitude,
                "elevation_m": self.grid_elevation_m,
            },
            "start_date": self.request.start_date.isoformat(),
            "end_date": self.request.end_date.isoformat(),
            "timezone": "UTC",
            "hourly_sample_count": len(self.times),
            "variables": list(self.variables),
            "units": dict(self.units),
            "source_artifacts": [dict(artifact) for artifact in self.source_artifacts],
        }


@dataclass(frozen=True)
class ERA5DailyPointData:
    """Daily meteorology normalized to the VISIT NCEP site's physical fields."""

    request: ERA5PointRequest
    dates: tuple[date, ...]
    variables: Mapping[str, np.ndarray]
    units: Mapping[str, str]
    source_manifest: Mapping[str, Any]
    transformations: Mapping[str, str]

    def __post_init__(self) -> None:
        if not self.dates:
            raise ValueError("daily ERA5 data contains no dates")
        expected_size = len(self.dates)
        for name in VISIT_NCEP_DAILY_VARIABLES:
            if name not in self.variables:
                raise ValueError(f"daily ERA5 data is missing {name}")
            values = self.variables[name]
            if values.shape != (expected_size,):
                raise ValueError(f"daily variable {name} has the wrong length")
            if not np.all(np.isfinite(values)):
                raise ValueError(f"daily variable {name} contains non-finite values")

    def to_manifest(self) -> dict[str, Any]:
        """Return JSON-compatible source and transformation metadata."""
        return {
            "source": dict(self.source_manifest),
            "adapter_approximation_level": "source-grounded input adapter",
            "visit_sources": [source.to_dict() for source in VISIT_WEATHER_SOURCES],
            "daily_sample_count": len(self.dates),
            "calendar_years": [
                {
                    "year": year,
                    "day_count": sum(current.year == year for current in self.dates),
                    "gregorian_leap_year": calendar.isleap(year),
                }
                for year in sorted({current.year for current in self.dates})
            ],
            "variables": list(VISIT_NCEP_DAILY_VARIABLES),
            "units": dict(self.units),
            "transformations": dict(self.transformations),
        }


def aggregate_era5_hourly_for_visit(
    hourly: ERA5HourlyPointData,
) -> ERA5DailyPointData:
    """Aggregate complete UTC days and derive fields used by VISIT."""
    expected_times = tuple(
        datetime.combine(hourly.request.start_date, time()) + timedelta(hours=hour)
        for hour in range(
            ((hourly.request.end_date - hourly.request.start_date).days + 1) * 24
        )
    )
    if hourly.times != expected_times:
        raise ValueError(
            "ERA5 timestamps must contain every UTC hour in the requested date range"
        )

    non_finite = {
        name
        for name, values in hourly.variables.items()
        if name != "snow_depth" and not np.all(np.isfinite(values))
    }
    if non_finite:
        raise ValueError(
            "ERA5 hourly variables contain non-finite values: "
            + ", ".join(sorted(non_finite))
        )

    n_days = len(expected_times) // 24

    def daily(name: str) -> np.ndarray:
        return hourly.variables[name].reshape(n_days, 24)

    air_temperature_c = daily("temperature_2m")
    dew_point_c = daily("dew_point_2m")
    surface_pressure_hpa = daily("surface_pressure")
    vapour_pressure_hpa = _saturation_vapour_pressure_hpa(dew_point_c)
    specific_humidity = 0.622 * vapour_pressure_hpa / (
        surface_pressure_hpa - 0.378 * vapour_pressure_hpa
    )

    wind_speed = daily("wind_speed_10m")
    wind_direction_radians = np.deg2rad(daily("wind_direction_10m"))
    u_wind = -wind_speed * np.sin(wind_direction_radians)
    v_wind = -wind_speed * np.cos(wind_direction_radians)

    shallow_soil_temperature_c = daily("soil_temperature_0_to_7cm")
    surface_temperature_c = daily("surface_temperature")
    soil_temperature_7_28cm_c = daily("soil_temperature_7_to_28cm")
    soil_temperature_28_100cm_c = daily("soil_temperature_28_to_100cm")
    deep_soil_temperature_c = daily("soil_temperature_100_to_255cm")
    shallow_soil_water = daily("soil_moisture_0_to_7cm")
    soil_water_7_28cm = daily("soil_moisture_7_to_28cm")
    soil_water_28_100cm = daily("soil_moisture_28_to_100cm")
    deep_soil_water = daily("soil_moisture_100_to_255cm")

    soil_temperature_0_10cm_c = (
        7.0 * shallow_soil_temperature_c + 3.0 * soil_temperature_7_28cm_c
    ) / 10.0
    soil_temperature_10_200cm_c = (
        18.0 * soil_temperature_7_28cm_c
        + 72.0 * soil_temperature_28_100cm_c
        + 100.0 * deep_soil_temperature_c
    ) / 190.0
    soil_water_0_10cm = (7.0 * shallow_soil_water + 3.0 * soil_water_7_28cm) / 10.0
    soil_water_10_200cm = (
        18.0 * soil_water_7_28cm
        + 72.0 * soil_water_28_100cm
        + 100.0 * deep_soil_water
    ) / 190.0

    snow_depth = np.nan_to_num(daily("snow_depth"), nan=0.0)
    kelvin_offset = 273.15
    variables = {
        "air_temperature_mean_k": air_temperature_c.mean(axis=1) + kelvin_offset,
        "air_temperature_max_k": air_temperature_c.max(axis=1) + kelvin_offset,
        "air_temperature_min_k": air_temperature_c.min(axis=1) + kelvin_offset,
        "specific_humidity_kg_kg": specific_humidity.mean(axis=1),
        "precipitation_rate_kg_m2_s": daily("precipitation").sum(axis=1)
        / 86400.0,
        "downward_shortwave_w_m2": daily("shortwave_radiation").mean(axis=1),
        "cloud_cover_percent": daily("cloud_cover").mean(axis=1),
        "surface_temperature_k": surface_temperature_c.mean(axis=1)
        + kelvin_offset,
        "soil_temperature_0_10cm_k": soil_temperature_0_10cm_c.mean(axis=1)
        + kelvin_offset,
        "soil_temperature_10_200cm_k": soil_temperature_10_200cm_c.mean(axis=1)
        + kelvin_offset,
        "soil_temperature_300cm_k": deep_soil_temperature_c.mean(axis=1)
        + kelvin_offset,
        "soil_water_0_10cm_m3_m3": soil_water_0_10cm.mean(axis=1),
        "soil_water_10_200cm_m3_m3": soil_water_10_200cm.mean(axis=1),
        "u_wind_10m_m_s": u_wind.mean(axis=1),
        "v_wind_10m_m_s": v_wind.mean(axis=1),
        "snow_depth_m": snow_depth.mean(axis=1),
        "surface_pressure_pa": surface_pressure_hpa.mean(axis=1) * 100.0,
    }
    units = {
        "air_temperature_mean_k": "K",
        "air_temperature_max_k": "K",
        "air_temperature_min_k": "K",
        "specific_humidity_kg_kg": "kg kg-1",
        "precipitation_rate_kg_m2_s": "kg m-2 s-1",
        "downward_shortwave_w_m2": "W m-2",
        "cloud_cover_percent": "%",
        "surface_temperature_k": "K",
        "soil_temperature_0_10cm_k": "K",
        "soil_temperature_10_200cm_k": "K",
        "soil_temperature_300cm_k": "K",
        "soil_water_0_10cm_m3_m3": "m3 m-3",
        "soil_water_10_200cm_m3_m3": "m3 m-3",
        "u_wind_10m_m_s": "m s-1",
        "v_wind_10m_m_s": "m s-1",
        "snow_depth_m": "m",
        "surface_pressure_pa": "Pa",
    }
    transformations = {
        "calendar": "strict 24-hour UTC days; no missing-hour interpolation",
        "specific_humidity_kg_kg": (
            "derived hourly from dew-point vapour pressure and surface pressure, "
            "then averaged"
        ),
        "precipitation_rate_kg_m2_s": (
            "sum of 24 preceding-hour precipitation amounts divided by 86400 s"
        ),
        "surface_temperature_k": (
            "ERA5 skin temperature from the official CDS"
        ),
        "soil_layers": (
            "ERA5 layers thickness-weighted to VISIT 0-10 and "
            "10-200 cm intervals"
        ),
        "soil_temperature_300cm_k": (
            "deepest available soil temperature; VISIT reads but does not use this field"
        ),
        "wind_components": (
            "meteorological direction and speed converted hourly to u/v, then averaged"
        ),
        "snow_depth_m": (
            "missing hourly values filled with zero; VISIT reads but does not use this field"
        ),
    }
    dates = tuple(
        hourly.request.start_date + timedelta(days=offset)
        for offset in range(n_days)
    )
    return ERA5DailyPointData(
        request=hourly.request,
        dates=dates,
        variables=variables,
        units=units,
        source_manifest=hourly.to_manifest(),
        transformations=transformations,
    )


def _saturation_vapour_pressure_hpa(temperature_c: np.ndarray) -> np.ndarray:
    """VISIT's water/ice saturation formula evaluated at dew point."""
    water = 6.1078 * np.power(10.0, 7.5 * temperature_c / (237.3 + temperature_c))
    ice = 6.1078 * np.power(10.0, 9.5 * temperature_c / (265.3 + temperature_c))
    return np.where(temperature_c >= 0.0, water, ice)


def format_visit_ncep_site_climate(
    daily: ERA5DailyPointData, *, require_complete_years: bool = True
) -> str:
    """Format the 19-column climate file read by VISIT ``WMODE=1``."""
    _validate_visit_dates(daily.dates, require_complete_years=require_complete_years)
    lines = []
    for index, current_date in enumerate(daily.dates):
        fields = [
            str(current_date.year),
            str(current_date.timetuple().tm_yday),
            *(f"{daily.variables[name][index]:.12g}" for name in VISIT_NCEP_DAILY_VARIABLES),
        ]
        lines.append(" ".join(fields))
    return "\n".join(lines) + "\n"


def write_visit_ncep_site_climate(
    daily: ERA5DailyPointData,
    output_path: str | Path,
    *,
    manifest_path: str | Path | None = None,
    require_complete_years: bool = True,
) -> None:
    """Write a VISIT site climate file and optional provenance manifest."""
    climate_text = format_visit_ncep_site_climate(
        daily, require_complete_years=require_complete_years
    )
    output_path = Path(output_path)
    output_path.write_text(climate_text, encoding="ascii")
    if manifest_path is not None:
        manifest = daily.to_manifest()
        manifest["visit_output"] = {
            "working_mode": 1,
            "site_format": "intended default NCEP/NCAR branch",
            "columns": list(VISIT_NCEP_COLUMNS),
            "day_of_year_convention": "1-based (VISIT reads but ignores this field)",
            "complete_years_required": require_complete_years,
            "artifact": _file_artifact(output_path),
            "source_compatibility_warning": (
                "visit_local/init_site.c::f_init_site and "
                "visit_local/location_proc.c::f_loct_proc in source commit "
                "3285bd8e131a932e338b59892751648fd9edcc7b have a missing '== 0' "
                "in the CEAMIP_TMK strcmp term. Correct both selectors before "
                "using this file with a non-CEAMIP site."
            ),
        }
        Path(manifest_path).write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )


def _file_artifact(path: Path) -> dict[str, Any]:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def _validate_visit_dates(
    dates: tuple[date, ...], *, require_complete_years: bool
) -> None:
    expected = tuple(
        dates[0] + timedelta(days=offset) for offset in range(len(dates))
    )
    if dates != expected:
        raise ValueError("VISIT climate dates must be consecutive")
    incompatible_years = {
        current_date.year
        for current_date in dates
        if (current_date.year % 4 == 0) != calendar.isleap(current_date.year)
    }
    if incompatible_years:
        years = ", ".join(str(year) for year in sorted(incompatible_years))
        raise ValueError(
            f"VISIT's year%4 leap-year rule differs from the Gregorian calendar: {years}"
        )
    if require_complete_years and (
        dates[0] != date(dates[0].year, 1, 1)
        or dates[-1] != date(dates[-1].year, 12, 31)
    ):
        raise ValueError("VISIT climate files must span complete calendar years")

