# ERA5 to VISIT weather framework

## Scope and authority

This adapter acquires hourly ERA5 point data, normalizes it to complete UTC
days, and writes the climate-file layout intended by the default NCEP/NCAR
branch of VISIT `WMODE=1`.

The VISIT source authority is `Sachitama2001/VISIT-matrix` commit
`3285bd8e131a932e338b59892751648fd9edcc7b`. The controlling source paths are:

- `visit_local/setting.c::f_setting`: `WMODE` and site file paths;
- `visit_local/initialize.c::f_initialize`: initialization dispatch;
- `visit_local/init_site.c::f_init_site`: site metadata and climate formats;
- `visit_local/init_region.c::f_init_region`: regional NCEP grid files;
- `visit_local/location_init.c::f_loct_init`: climatologies derived from input;
- `visit_local/location_proc.c::f_loct_proc`: daily selection and overrides;
- `visit_local/experiment.c::f_experiment` and `spinup.c::f_spinup`: calendars.

The implementation is a source-grounded input adapter. It does not change
VISIT's equations and is not itself a VISIT model implementation.

## Acquisition paths

| Path | Role | Authentication | Result |
| --- | --- | --- | --- |
| Copernicus CDS | Authoritative bulk/publication path | CDS account and `cdsapi` configuration | Monthly ZIP archives containing NetCDF point data |

Direct CDS requests are split by month. This avoids the Cartesian expansion
that can occur if years, months, and days from an arbitrary interval are put
into one request. Point coordinates are snapped to the nearest regular
0.25-degree ERA5 grid point before constructing the zero-area MARS selection;
the original requested coordinate remains in the caller's `ERA5PointRequest`.
The requested CDS variables include ERA5 skin temperature and native 10 m u/v
wind components. CDS separates instantaneous and accumulated variables into
NetCDF members within the ZIP archive. The decoder validates the common time
axis and single grid point, joins those members, converts ERA5 SI units to the
canonical hourly schema, and retains skin temperature for VISIT's surface
temperature field.

Each archive name includes the exact selected start and end dates. A non-empty
ZIP with readable members and at least one NetCDF member is reused by default;
new downloads are validated and atomically moved into place. Use
`--refresh-cds` when an ERA5T update or another provenance requirement calls
for replacing a valid cache. Cache-only decoding does not initialize a CDS
client and therefore does not require network access. Archive and final
VISIT-file SHA-256 values are recorded in the manifest.

## VISIT configuration paths

`f_setting` defines only `WMODE=1` (site) and `WMODE=2` (region). There is no
`WMODE=3` in this source snapshot.

### Site mode formats

`f_init_site` reads `grid->site_id` from the SITEF file before selecting the
climate format. That value overwrites the `SITE` value previously read from
`setting.txt`.

| Selector in `f_init_site` | Records | Values after `year day` | Source behavior |
| --- | --- | --- | --- |
| CEAMIP site-ID list | 365 daily rows/year | mean/max/min 2 m temperature in C, precipitation, relative humidity, VPD in Pa, shortwave, longwave, wind | Converts mean temperature to K and VPD to hPa; sets cloud to 50 percent and all surface/soil temperatures to air temperature; max/min and longwave are discarded |
| `LUCMIP0` to `LUCMIP3` | 4 subdaily rows/day, 365 days/year | longwave, shortwave, 2 m temperature in K, precipitation, specific humidity, pressure in Pa, u wind, v wind | Averages state variables, sums precipitation, derives wind speed and VPD, fixes surface/soil temperature to air temperature |
| `EX_ASIAMIP == 1` | 365 daily rows/year | Same columns as CEAMIP | Same conversions and fixed substitutes as CEAMIP |
| Intended default NCEP/NCAR branch | 365 or 366 daily rows/year | The 17 fields listed below | Uses most fields directly and derives wind speed and VPD |

The intended default NCEP/NCAR row has 19 whitespace-separated fields:

1. year;
2. day of year (read but not validated or used);
3. mean 2 m air temperature, K;
4. maximum 2 m air temperature, K (discarded);
5. minimum 2 m air temperature, K (discarded);
6. 2 m specific humidity, kg/kg;
7. precipitation rate, kg m-2 s-1, multiplied by 86400 on input;
8. downward shortwave radiation, W m-2;
9. total cloud cover, percent;
10. surface temperature, K;
11. 0-10 cm soil temperature, K;
12. 10-200 cm soil temperature, K;
13. 300 cm soil temperature, K (discarded);
14. 0-10 cm soil water content (discarded);
15. 10-200 cm soil water content (discarded);
16. 10 m u wind, m/s;
17. 10 m v wind, m/s;
18. snow depth (discarded);
19. surface pressure, Pa.

### Regional mode format

`f_init_region` expects ten files named for downward shortwave radiation, cloud
cover, 2 m air temperature, skin temperature, two soil temperatures, specific
humidity, precipitation rate, and 10 m u/v wind. Every date record contains a
fixed 16 by 19 NCEP grid. The reader is hard-coded to 1948-2005 and to fixed
Japanese-domain latitude/longitude arrays. It computes VPD using a fixed 500 m
altitude rather than each grid cell's topography.

This regional reader is not on the active initialization path in the pinned
snapshot: `f_initialize` calls `f_init_site` unconditionally, and there is no
call to `f_init_region`. Supporting `WMODE=2` therefore requires a source
dispatch correction before an ERA5 regional writer would be useful.

## Source compatibility findings

The pinned source must not consume the generated default-NCEP file unchanged:

1. The CEAMIP selector in `init_site.c::f_init_site` contains
   `strcmp(grid->site_id, "CEAMIP_TMK")` without `== 0`. It is true for nearly
   every non-`CEAMIP_TMK` ID, so the intended NCEP, LUCMIP, and ASIAMIP branches
   are normally pre-empted. The analogous CEAMIP humidity selector in
   `location_proc.c::f_loct_proc` has the same omission. Correct both selectors
   before using a generic site ID with this writer.
2. VISIT uses `year % 4 == 0`, not the Gregorian leap-year rule. The writer
   rejects century years where those definitions differ.
3. Normal experiments can use 366 days, but `f_spinup` forces 365 days. CEAMIP,
   LUCMIP, and ASIAMIP readers also force 365 days.
4. `f_loct_init` constructs 1980-2009 climatologies using direct
   `f-BYR` indexing. The source assumes the loaded climate period covers those
   years; a later isolated year can lead to invalid indexing.
5. The LUCMIP VPD expression divides the already daily-derived deficit by four
   once more. This appears inconsistent with the surrounding averaging, but
   the adapter does not silently correct that source behavior.

The generated JSON manifest repeats the first compatibility warning next to
the output column contract.

## Hourly to daily transformations

All non-snow variables must be finite for every requested UTC hour. Timestamps
must exactly cover 00:00 through 23:00 for each date; no missing-hour
interpolation is performed.

| VISIT field | ERA5/CDS transformation |
| --- | --- |
| Mean/max/min air temperature | Mean/max/min of 24 hourly 2 m temperatures, C to K |
| Specific humidity | Derive hourly from dew-point vapour pressure and surface pressure, then average |
| Precipitation rate | Sum 24 preceding-hour amounts in mm and divide by 86400 s |
| Downward shortwave | Mean of 24 preceding-hour mean fluxes |
| Cloud cover | Mean of hourly percent cover |
| Surface temperature | Mean ERA5 skin temperature |
| 0-10 cm soil values | Thickness weighting of 0-7 and 7-28 cm layers: `(7*x1 + 3*x2)/10` |
| 10-200 cm soil values | Thickness weighting: `(18*x2 + 72*x3 + 100*x4)/190` |
| 300 cm soil temperature | Deepest available layer; retained only to fill a discarded VISIT column |
| u/v wind | Convert every hourly meteorological speed/direction pair to u/v, then average |
| Snow depth | Missing values become zero only because VISIT discards this input field |
| Surface pressure | Mean hourly hPa converted to Pa |

The layer weighting assumes each ERA5 layer value represents its
documented layer mean. It is an approximation at VISIT's requested boundaries,
not an observation at an exact depth.

For the 2020-01-01 Takayama CDS check, weighted soil temperatures differed from
using one nearest native layer by 0.12 K for 0-10 cm and -1.97 K for 10-200 cm.
These are one-day diagnostics, not general uncertainty estimates.

## Production period design

The unmodified `f_loct_init` computes its climatology with direct accesses for
1980-2009. A native VISIT input intended to use this path must therefore include
at least complete calendar years 1980 through 2009, even if the scientific
experiment begins later. Generate complete years, retain all monthly CDS ZIPs
and their checksums, and set `BYR/EYR` consistently. A shorter ERA5 file is
valid for adapter diagnostics but not for that unmodified initialization path.

## Soil-temperature and soil-water coupling policy

ERA5 soil temperature is mapped to the temperature fields already read by the
VISIT NCEP site format: the ERA5 layers are thickness-weighted to 0-10 cm and
10-200 cm and become `tmp10_soil` and `tmp200_soil` after VISIT converts K to C.

ERA5 volumetric soil water is not assigned directly to decomposition water
variables. In native VISIT:

- `soilwtr_l = mass->sw30`, the internally evolved 0-30 cm bucket water in mm;
- `soilwtr_h = mass->sww`, the internally evolved deeper/whole bucket water in mm;
- field capacities are water depths in mm, derived from soil texture/site data;
- `soilappr_l/w` are diagnosed as clipped `1 - soil_water / field_capacity`.

ERA5 reports volumetric water in m3/m3 over different layer boundaries. A
conversion would require compatible layer thickness, rooting depth, porosity,
field capacity, and a decision about whether ERA5 replaces or only validates
VISIT's prognostic hydrology. Until that contract is implemented and tested,
ERA5 soil water is retained as diagnostic input metadata and is not injected
into `soilwtr_l/h`.

## Runtime overrides after file reading

The climate file is only the first stage. `f_loct_proc` can subsequently replace
or alter values in this order:

1. sequential data for `BYR <= climy <= EYR`;
2. average-climate or LARS spin-up (`EX_SPINUP=2` or `4`);
3. GCM anomaly prediction (`PREDICT >= 1` beyond `EYR`);
4. sensitivity controls for temperature, precipitation, and radiation;
5. decomposition warming experiment (`EX_DECTMP`);
6. hard-coded site corrections for TKY, TMK, LSH, TSE, SKT, MKL, and QHB;
7. snow-pack damping of soil temperatures.

`FIX_CLIM` also changes the year/day selected by the experiment and spin-up
loops. A run configuration must therefore record both the climate manifest and
the compile-time values in `setting.h`.

## Usage

Production VISIT files omit `--allow-partial-years`; the CLI then requires
complete calendar years. It always writes a provenance manifest beside the
climate file.

Official CDS request construction and retrieval:

```python
from datetime import date
from control_carbon import ERA5PointRequest, retrieve_cds_era5_point

request = ERA5PointRequest(
    36.146, 137.423, date(2020, 1, 1), date(2020, 1, 3)
)
paths = retrieve_cds_era5_point(request, "era5-cds")
```

Official CDS retrieval, decoding, daily aggregation, and VISIT output can also
be run in one command:

```bash
control-carbon-era5 \
   --latitude 36.146 --longitude 137.423 \
   --start-date 2020-01-01 --end-date 2020-01-01 \
   --output /tmp/era5-cds-takayama.visit.txt \
   --allow-partial-years
```

Install the optional dependency with `pip install -e '.[era5]'`, configure the
CDS API by its official instructions, and accept the dataset terms. Credentials
remain under `cdsapi` management and are never accepted by this package.

## Numerical guarantees

`tests/test_era5.py` verifies:

- ERA5, UTC, and nearest-grid-point request parameters;
- rejection of a missing hourly timestamp;
- hand-computable temperature, humidity, precipitation, soil-layer, wind, and
  pressure transformations;
- precipitation mass preservation through the VISIT rate conversion;
- the exact 19-column VISIT output order and complete-year guard;
- monthly CDS splitting with no out-of-range date combinations.

The live CDS check retrieved 24 hourly records for 2020-01-01 near Takayama and
resolved `(36.146, 137.423)` to ERA5 grid point `(36.25, 137.5)`. The complete
CDS ZIP, decoded daily values, VISIT output, and SHA-256 manifest were validated.