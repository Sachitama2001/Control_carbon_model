# Current implementation status

## Temperature–NSC autonomous checkpoint (2026-09-22)

See `temperature_nsc_first_results.md` and the Japanese TeX report. An independent
eight-carbon-state theoretical model and source/assumption ledger are implemented.
Forty-two standard spinups converge; a supported-canopy closure creates bistability
and a fold, but 8 rate-family runs and 24 controls do not establish R-tipping.
Near-fold transients are checked to 2000 years. NSC mortality is not necessary
for the observed bistability. Numeric parameter constraints, disconnected
branches/cycles and seasonal extension remain open; no ecological calibration
or universal no-tipping result is claimed. Existing water work is unchanged.

## T/P v2 explicit spinup and configurable processes (2026-09-22)

- Plan revision 1.1 adds only initialization/nonnegativity §4.4; exact previous
  text is `temperature_precipitation_tipping_research_plan_v1.md`. v1 results
  and source files are unchanged, with hashes verified using that archive.
- New `tp_experiments_v2.py` integrates spinup from a non-equilibrium seed.
  Ten models converge in 375–400 years; 40 factorial experiments follow.
  The original four mortality variants match v1 trajectories within 2.13e-7.
- JSON configuration: `configs/tp_v2_spinup.json`. `tp_processes.py` supplies
  selectable old-VISIT temperature/canopy/soil-scalar/ET-supply/baseflow parts,
  with explicit ODE mappings and native C component comparisons.
- Every accepted step plus dense interior samples is guarded against negative
  pools; retry smaller steps or fail, never silently clip. Audited 1,095,865
  points, no negative pools/capacity excess, zero clipping mass.
- The negative v1 plot quantity was I_C, not a stock. Equation choice changes
  its sign in some v2 comparisons; no rate/basin or ecological claim follows.
- Results/configuration/provenance: `docs/tp_v2_results.md` and
  `docs/tp_v2_spinup_and_functions.md`; artifacts in
  `artifacts/temperature_precipitation/v2/spinup_comparison/`.
- Validation: 190 passed / 7 skipped, including 14 new dedicated tests.

## Temperature–precipitation first experiments (2026-09-22)

- New independent `temperature_precipitation.py`: eight physical states, T/P
  forcing, mass budgets, moving water-capacity bounds, carbon capacity,
  equilibrium and Jacobian. Defaults are uncalibrated assumptions.
- WP1 passed; first WP2 controls run for none/heat/water/additive mortality.
  Each hierarchy has 0/T/P/TP, two rainfall endpoints, and 208 total
  single-driver continuation points including reverse traversal.
- P=4 mm/day remains above field capacity and yields a carbon-response
  plateau; P=2 produces quantitative interactions without an explicit product
  mortality term. All sampled positive branches/endpoints are stable.
- Water mortality: I_C is -9.901 Mg C/ha after 20 years, crosses zero around
  year 33.1, and reaches +37.380 after 200 years, approaching the frozen
  equilibrium contrast +37.416. This is a transient interaction, not R-tipping.
- Fixed absolute-water interventions leave the moving-capacity domain after
  670–1050 days; rejected for ecological attribution after violation. Other
  normal eight-state trajectories pass bounds and budgets.
- Tests: 176 passed / 7 skipped; 13 dedicated tests. Cross-solver normalized
  differences <=2.41e-8. WP2 attribution still needs a valid clamp; WP3–5,
  boundary equilibria, seasonal reference, and rate/basin tests remain pending.
- Results and reproducibility: `docs/temperature_precipitation_first_results.md`
  and companion TeX; artifacts under `artifacts/temperature_precipitation/`.

## Temperature–precipitation theory plan (2026-09-22)

- `temperature_precipitation_tipping_research_plan.md` defines a separate
  low-dimensional track with matching leaf, stem, root, and soil carbon and
  water stores (eight states);
- temperature and precipitation are the only varying external inputs, with
  continuous ODE dynamics and monthly forcing resolution;
- the experiment design contains baseline, temperature-only,
  precipitation-only, and combined runs plus pathway-clamped attribution;
- the plan distinguishes quantitative non-additivity, frozen-system
  bifurcation, and rate-induced basin transitions;
- an eastern/southern Amazon seasonal-forest-type idealized site is the first
  ecological context, not a calibrated site model;
- no eight-state implementation, parameter calibration, or R-tipping result is
  claimed yet.

## Native summer-state checkpoint (2026-09-21)

- Latest workbook: `artifacts/visit_matrix_workbook/VISIT_TKY_matrix_equations_native_20000715.xlsx`.
- Pinned 2013 source run from spinup with supplied TKY data; all 37 states and
  15 history variables populated for the start of 2000-07-15.
- Missing inputs: 70 -> 2; both remaining values are undefined source diagnostics.
- Native versus workbook one-day update: 37/37 states, max scaled error 5.44e-16.
- Corrected workbook runtime capacities: upper 104.810272811 mm and deep
  704.059727189 mm after Saxton overwrite, rather than raw Config 64.21/808.87.
- Native vapor pressure follows a source `strcmp` comparison omission; retained
  for exactness and documented rather than silently repaired.
- Full repository verification: 163 passed, 7 skipped.  The skipped cases are
  optional native/dependency paths, not failures of the generated workbook.

## Summer forcing import checkpoint (2026-09-17)

- New dated workbook: `artifacts/visit_matrix_workbook/VISIT_TKY_matrix_equations_20000715.xlsx`.
- Reference 2000-07-15, actual supplied visitb meteorology/deposition/GHG;
  24 inputs resolved, missing scalar inputs 70 → 46. Forty-one sheets.
- Ten history values reconstructed from prior temperature with unknown initial
  memory; only values fixed by resets are admitted. No arbitrary zero fill.
- Summer stocks unavailable: daily output contains fluxes, not pool states.
  Twenty-one spinup-end restart stocks are shown separately, never adopted as summer X.
- Pinned Config/model retained. Source differences and humidity inconsistency
  are documented in the workbook and `docs/visit_matrix_workbook.md`.
- New import tests and the existing workbook tests: 16 passed. Excel application
  recalculation remains untested; generation checks formula caches and structure.
- Full test suite after import changes: 162 passed, 7 skipped.

## TKY C/N/water Excel checkpoint (2026-09-17)

- Delivered workbook target: `artifacts/visit_matrix_workbook/VISIT_TKY_matrix_equations.xlsx`.
- Japanese guide and exact scope: `docs/visit_matrix_workbook.md`.
- Selected old `visit_local` daily pathway (`FLUX_SCHEME=0`), TKY DBF,
  37 states / 132 fluxes / B(37×78), plus phenology and previous-day gc memory.
- Includes real, acyclic Excel formulas, cached results, Config_TKY defaults,
  source function/line provenance, C/N/water budgets and explicit missing inputs.
- Seven artificial cases / 1,379 native-C comparisons agree at rounding scale;
  13 dedicated tests pass. Full driver, all branches, unresolved diagnostics,
  and actual Excel-engine recalculation are outside this verification.
- This transcribes the daily plant--NSC loop for inspection. It does not itself
  establish an attracting alternative trajectory or validate an ODE embedding.

## Capacity/QSE exploration checkpoint (2026-09-17)

- `capacity_qse_attractor_exploration.tex` derives the nonlinear capacity map
  \(\widehat X_c(X)=-M(X)^{-1}u(X)\), the self-consistency condition
  \(X^*=\widehat X_c(X^*)\), and the equilibrium Jacobian
  \(J=M(I-D_X\widehat X_c)\);
- `capacity_feedback.py` gives a deliberately model-agnostic positive-carbon
  counterexample with two nonzero stable QSE branches and a moving basin
  boundary; fast and slow ramps select different branches while the tracked
  branch remains stable;
- this establishes mathematical possibility, not multiple attractors in VISIT;
- `visit_r_tipping_process_screening.tex` records that the explanatory slides'
  individual photosynthesis, stomatal, respiration, and decomposition formulas
  do not alone establish R-tipping;
- the pinned executable source instead prioritizes the coupled
  LAI--GPP--piecewise-allocation--NSC loop, whose full daily update is not yet
  implemented;
- a second equilibrium is still a QSE. A literal non-QSE attracting trajectory
  under constant forcing would require a periodic or more complex attractor and
  remains unproven.
- verification at this checkpoint: 146 tests passed and 7 source-dependent
  tests skipped; both Japanese reports compile without layout warnings.

This file is the short checkpoint after the first explicit carbon-water
implementation. The historical soil, ERA5, plant-slice, and tipping work
remains available. The active research order is now
`carbon_water_research_plan.md`; `deterministic_tipping_plan.md` is a deferred
analysis track rather than the project driver.

## Active carbon-water checkpoint

- current `visit-manager/VISITc` source is pinned independently at commit
  `5202debd96df6f88beb7d61f8688fff02ace964a`;
- native `snwa`, `sw30`, and `sww` stores plus the source-order hydrology fluxes
  are inventoried in `visitc_carbon_water_source_ledger.md`;
- `visitc_hydrology.py` transcribes `f_hydrology` and the seven Penman-Monteith
  potentials, including atmosphere, aerodynamic resistance, LAI/radiation,
  soil resistance, and canopy-conductance dependencies;
- minimal native bridges compare every exposed one-day hydrology output and
  directly execute `hydro_flows.c`, `radiation.c`, and `f_canopy_cond`;
- the apparent baseflow double subtraction is reproduced and exposed through
  a budget residual, not silently corrected;
- `coupled_carbon_water.py` implements a four-carbon/four-water continuous ODE
  with exact internal mass cancellation and water-potential-gradient transport;
- normalized transport storage was replaced by the native VISITc soil
  retention relation and a reduced SurEau-Ecos pressure-volume/capacitance relation;
- instantaneous frozen-water carbon capacity is kept separate from the full
  coupled equilibrium;
- productivity, transit-time, and interaction effects on capacity are
  decomposed exactly;
- full Jacobian blocks and the water Schur complement are implemented;
- the idealized dry-down/re-wetting example records budget residuals and solver
  refinement.

The storage-to-potential functions are now source-labelled. Conductance
magnitudes, plant saturated water/allometry, carbon process parameters, and
external-flux closures remain illustrative. It is not a calibrated VISITc or
species model.

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
python examples/run_coupled_carbon_water.py
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

Verified after the native hydrology/PM and hydraulic-constitutive milestone:
`118 passed` when both pinned source checkouts are available. Source-dependent
native tests skip when their separate checkout is absent. The updated
carbon-water example reported a stable baseline, maximum step-refinement state
difference `4.35e-7`, carbon budget residual `2.17e-18`, and water budget
residual `4.44e-16` for the illustrative run.

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

For the paper-oriented capacity question, the next source-grounded milestone is
to complete the pinned VISIT plant daily map through GPP, respiration, turnover,
piecewise allocation, NSC, and survival reallocation. Test constant-forcing
fixed points first, then the annual Poincare map, before applying rate ramps.
In the broader carbon-water track, dynamic water storage versus its quasi-steady
Schur reduction remains the central comparison. Neither track should tune
parameters merely to force tipping.
