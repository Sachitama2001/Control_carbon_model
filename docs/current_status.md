# Current implementation status

This file is the short checkpoint after the source-grounded soil, ERA5, and
minimal native-soil validation cycles. The Japanese progress narrative and
P0-P5 record is in `progress_and_next_plan.md`. The active research order is
now `deterministic_tipping_plan.md`; A-K phases remain a longer-term backlog.

## Completed vertical slices

### Generic mathematics

- deterministic continuous nonlinear integration with one varying driver;
- frozen root residuals, analytic Jacobians, stability and conditioning;
- signed compartment algebra with nonnegative transfer/input structure;
- hypothetical saturating-feedback branches and analytically known fold;
- finite-ramp B/R benchmarks, final holds, basin checks and solver refinement;
- temperature-only VISIT soil control and reproducible PNG/NPZ/JSON results;
- continuous and discrete LTI validation;
- discrete poles, stability, transfer/frequency response, impulse, step, and
  forced response;
- constant-input discrete fixed points with residual and conditioning data.

### VISIT soil carbon

- authoritative source pinned to `Sachitama2001/VISIT-matrix` commit
  `3285bd8e131a932e338b59892751648fd9edcc7b`;
- source-order 9-pool daily carbon update;
- equivalent fixed-coefficient discrete matrices;
- carbon balance, clipping, fixed point, Rh IRF, and direct-trajectory
  comparison;
- `frl` and `frh` environmental scalars with undefined native `stype=1/2`
  behavior rejected and the inferred correction separately labelled.
- native C bridge compiling pinned `soil_proc.c` and `decomposition.c`;
- one-day, 40-day direct trajectory, and 100-day native perturbation-versus-IRF
  validation for the 9-pool carbon subsystem.
- daily varying decomposition environments and a tangent environmental model;
- grouped discrete modes, Rh residues, and fixed-point environmental sensitivity;
- reproducible native/nonlinear/tangent temperature-pulse comparison figure.

### VISIT plant process slices

- the first P5 slice reproduces normal foliage/stem/root turnover;
- TKY deciduous shedding, final shedding, crop-stage 70% turnover, and drought
  leaf shedding are explicit source-labelled branches;
- immediate structural-pool subtraction and litter flux are covered by
  hand-computable tests;
- photosynthesis, respiration, allocation, NSC, and survival reallocation are
  not yet integrated into one daily plant map.
- Q10 acclimation, size-dependent maintenance respiration, and growth
  respiration are implemented and matched to native C functions;
- positive/nonpositive EPP, LAI, dormant-season, and crop-grain allocation
  branches are implemented and matched to native `f_allocation`;
- NSC storage, leaf emergence, survival reallocation, and GPP remain to be
  integrated.

### ERA5 to VISIT weather

- official Copernicus CDS point acquisition;
- nearest 0.25-degree CDS grid selection;
- monthly CDS ZIP retrieval with exact-range names, cache validation, and
  atomic replacement;
- instantaneous/accumulated NetCDF decoding;
- strict complete-UTC-day aggregation;
- intended VISIT `WMODE=1` NCEP/NCAR 19-column output;
- source archive and output SHA-256 records in the JSON manifest;
- synthetic multi-month and complete leap-year regression coverage;
- successful live CDS retrieval and end-to-end one-day VISIT output.

## Supplementary parameter workbook

The user-supplied source
`/mnt/d/CT/VISIT-matrix/visit_local/INPUT/parameter_VISITc_16.xlsx` is pinned by
SHA-256
`f8748b9dca3a7e7e38a2aa7ce93fd54c22653948f2467fd2ed7c1f6181edbad3`.

It contains 16 MOD12 land-cover columns and sheets for Tree, Herb3, Herb4,
Soil, Canopy, and the runtime text layout. It is supplementary because it is
not part of the pinned Git commit. `parameter.c::set_parameter` remains the
authority for runtime ordering. The referenced `parameter_S1b.txt` is absent
from the inspected snapshot, so workbook values must not yet be described as
values proven to have been used by a native run.

The Soil sheet provides all parameters used by the current Python soil and
decomposition slices: `sr_lf/lc/lr`, `sr_ha/hi/hp`, `kml/kmh/kmsl/kmsh`,
`f_co2_lf/lc/lr`, and `f_hm_a/i/p`. Workbook spellings `f_hm_I` and
`root_spt_a/b` correspond to C symbols `f_hm_i` and `root_dpt_a/b`.

## Reproducible checks

```bash
python -m pip install -e '.[dev,era5]'
pytest -q
git diff --check
```

Official CDS smoke test requires `~/.cdsapirc` and accepted dataset terms:

```bash
control-carbon-era5 \
  --latitude 36.146 --longitude 137.423 \
  --start-date 2020-01-01 --end-date 2020-01-01 \
  --output /tmp/era5-cds.visit.txt \
  --allow-partial-years
```

Live network checks are intentionally not part of the default pytest suite.
The suite tests provider parsing and conversion with deterministic fixtures.

## Known blockers before native VISIT validation

- CEAMIP `strcmp` conditions in `init_site.c` and `location_proc.c` omit
  `== 0` for one term and pre-empt the intended generic NCEP branch;
- `initialize.c` calls the site initializer unconditionally, leaving the
  regional initializer unreachable;
- `f_loct_init` directly constructs a 1980-2009 climatology, so a production
  native run needs that coverage or a documented source correction;
- spin-up forces 365 days while normal experiments can use 366;
- the runtime global parameter text and complete native input bundle are not
  present.

## Next scientific milestone

Connect one source-informed plant feedback to the deterministic branch/basin
analysis, initially varying temperature only and fixing other meteorology.
The new scalar benchmarks verify mathematical tools, not VISIT/TKY tipping.
Generic continuation, critical-rate estimation, return experiments, integrated
plant dynamics and higher-dimensional basins remain pending. Stochastic and
N-tipping tools are explicitly excluded from the current stage. See
`deterministic_tipping_plan.md` for acceptance criteria and reproducible checks.