# Current implementation status

## Scalar equilibria and hydraulic-failure capacity (2026-10-14)

The latest [user plan](research_plan_hydraulic_failure_capacity_20261014.md)
supersedes earlier next-stage ramp/cycle/basin priorities.
[Derivation, proofs and results](spac_carbon_scalar_equilibrium_results.md)
reduce all positive E0/E1/E4 frozen roots to D_curve(L), with full state recovery.
Analytic E0/E1 single maxima and outward-checked E4 beta1 monotonicity
classify root counts globally on the positive conducting domain.
E0/E1 high-L roots are stable; E1 uses a transformed irreducible Metzler core.
E4 stability remains sampled/numerical: known .002664752905 MPa crossing
and new formal .004791823126 MPa restabilization, outside routine/25 C forcing.
No new l1/cycle campaign or continuous E4 stability certificate.
Positive equilibria satisfy K=C*, regardless of attraction.
gs=0 implies immediate K=0; both closed imply frozen water and C=exp(Mt)C0 ->0.
gx-only closure is invariant and gives C/K ->0 under baseline unsafe R0 bound;
E4 may retain positive gs. Formal zero-carbon continua and kinks kept separate.
262 full-root records;8 fixed-.00244 MPa examples;24 new tests;
236 related regression cases passed. Old A/B/C sources/configs unchanged.
Adjacent-float zero-leaf endpoint tests preserve strictly positive roots.
Artifacts: [scalar_equilibria_final_v5_20261014](../artifacts/spac_carbon/scalar_equilibria_final_v5_20261014/).
Current stop condition completed: scalar reduction, root classification,
capacity-loss propositions. Further ramps/global basins/new mortality deferred.

## WP-C2 and bounded WP-C3 (2026-10-14)

The [tracking / local-section report](spac_carbon_stage_c23_results.md) follows
the capacity-first plan. New independent modules provide analytical F_v and
full-root tangents, exact/linearized tracking, frozen non-normal gain and moving
Lyapunov metrics, plus phase-fixed E4 Poincare geometry.
Six predefined routine-MPa ramps: .00266/10day exits L_min at~15.78year;
the five remaining cases stay positive through300year hold but are not classified
as recovered. Large-path linearization is inaccurate; small paths verify quadratic
accuracy. All sampled big-path metric drift tests fail the sufficient guarantee,
not the necessary stability condition. No continuous/uniform certificate.
E4 r=.002 cycle has~474year radial e-folding and a codim1 local stable tangent
candidate; this is a DIFFERENT VPD slice from the rate endpoint, not a proved
nonlinear/global recovery-exit boundary. Scalar20/120year benchmark reproduces
low/high transition; linear periodic benchmark is nonautonomous.
13 new tests;212 related regression passed. A/B/C1 models unchanged.
Artifacts: `artifacts/spac_carbon/stage_c23_tracking_final_v2_20261014/`.
Next: same-slice cycle continuation and nonlinear local manifold/outcome sections.
C4 new feedback laws, C5 manuscript, autonomous stable-cycle attractors and true
E4 R-tipping remain unverified/deferred.

## Stage C capacity-first: WP-C0/C1 (2026-10-13)

The [new plan](SPAC_nonlinear_carbon_storage_capacity_research_plan.md) takes
precedence over standalone global tipping exploration. The
[definition ledger](spac_carbon_capacity_definition_ledger.md) and
[Stage C results](spac_carbon_stage_c_results.md) fix K/P signs and distinguish
the map, conditional roots, full roots, attracting sets and actual trajectories.
New pure diagnostics preserve Stage A/B; 27 new tests and 199 directly related
regression cases pass. 232 bounded frozen-root records, 15 counterexamples,
E4 transient and periodic diagnostics show P=0 does not imply H=0, pool
cancellation does not imply P=0, and ||P|| can increase on the repelling cycle.
Numerical Floquet modes imply a local codimension-one stable-manifold candidate;
the ~1657-year radial e-folding time precludes casual short-horizon classification.
The global manifold/recovery-exit boundary, speed experiment, stable non-equilibrium
attractor and true R-tipping remain unverified. A conditional analytical tracking
bound is derived, but its uniform assumptions are not certified for E4.
Artifacts: `artifacts/spac_carbon/stage_c_capacity_final_v2_20261013/`.
Next WP-C2/C3: capacity-target tracking and bounded Floquet-directed sections.

## E4 local nonlinear milestone (2026-10-13)

The [B3 E4 results](spac_carbon_stage_b3_e4_results.md) follow the latest MPa
unit/scale review rather than treating high-VPD E3 cycles as routine climate.
Exact four-state tensors, transversality and positive l1 numerically resolve
E4's candidate as subcritical Hopf at .002664752905 MPa. Independent periodic
collocation, full R/S periodic lift, IVP comparison, Floquet multipliers and
tolerance refinement verify local repelling cycles, NOT attractors.
Seven bounded synthetic sensitivities show that wood exponent/loss assumptions
can alter climate admissibility or produce a fold and an already-unstable
complex-pair crossing. Original physiology remains unchanged.
Targeted/neighboring tests:163 passed. Artifacts:
`artifacts/spac_carbon/stage_b3_e4_final_v2_20261013/`.
Next: E4 global boundary geometry and coefficient provenance; post-exit
physiology requires a separate justified choice before mortality/R-tipping.

## SPAC units and exploration scales (2026-10-13)

The [unit/scale review](spac_carbon_units_and_scales.md) introduces an explicit
MPa-facing boundary without changing original coefficients or equations.
VPD remains vapor-pressure deficit, NOT atmospheric water potential.
User-approved synthetic assumptions: p_atm=.101325 MPa, diagnostic 25 C,
routine VPD 0--.003 MPa. Existing canonical APIs/configs/runners remain
historical mol/mol reproduction interfaces. New tests: 27 passed;
related baseline six-file regression: 127 passed.
Pressure curves, derivative/shift conversion, shooting and outward certificate
inclusion are verified; E2/E3 11/26-cell certificates regenerated.
E2/E3 high-VPD Hopf exceeds saturation at 25 C; mark formal exploration.
E4 .002664753 MPa remains in range. Initial next-stage work checks exact
four-state/passive-RS factorization and local quartic Hurwitz crossing;
E4 l1/nondegenerate Hopf and cycles remain open.
Artifacts: `artifacts/spac_carbon/pressure_review_final_v2_20261013/`,
`artifacts/spac_carbon/phase3d_mpa_verified_20261013/`.

## SPAC–carbon Stage B2 (2026-10-12)

The [short B2 plan](SPAC_carbon_stage_B2_plan.md) is implemented through the
[first B2 milestone](spac_carbon_stage_b2_results.md). E2/E3's exact three-state
dynamic core gives a cubic Hurwitz test; outward-rounded interval certificates
cover the entire positive frozen branch for D=.01-.03. This does not imply
nonautonomous tracking/no-R. Analytic higher derivatives, eigenvalue
transversality and positive first Lyapunov coefficients numerically identify
subcritical Hopf. Amplitude-fixed shooting, positive full-six-state periodic
lifts and Floquet multipliers verify small repelling cycles on the stable
equilibrium side, not additional attracting states. Formal zero-leaf limits
and conservation-aware post-exit requirements are documented; no mortality
regime or new physiological feedback was invented.
Targeted/neighboring tests: 129 passed. Stage A/B1 sources are unchanged.
Final artifacts: `artifacts/spac_carbon/stage_b2_verified_20261012/`.
Next: global recovery/exit boundary geometry, then user-guided physical regime
selection; E4 nonlinear analysis and genuine R-tipping remain open.

## SPAC–carbon Stage B1 (2026-10-09)

The [new user plan](SPAC_carbon_stage_B_plan.md) prioritizes structural
flux feedback and defers Leuning. The independent
[B1 results](spac_carbon_stage_b_results.md) implement E0–E4 with fixed capacities
and fixed Stage A references, correct series conductance, unchanged GPP/budgets,
full Jacobians and positive-equilibrium pseudo-arclength continuation.
E1 shifts the synthetic fold from .02439716 to .01526558 mol/mol; E2/E3/E4 have
oscillatory stability boundary candidates rather than detected folds in the
stated scan. Full-six-state multistart roots, basin slices, sensitivities and
same-path 10-day/10-year comparisons are saved. E2/E3 show rate-dependent
vegetation-domain exit versus verified recovery, NOT proven R-tipping:
no physical alternate attractor is defined beyond L_min. An analytic E2
invariant-box non-exit bound is rate-independent but is not a convergence theorem.
Targeted + neighboring tests: 115 passed. Stage A is unchanged and revalidated.
Outputs: `artifacts/spac_carbon/stage_b1_final_20261009/`.
Next: continuous-path stability proof, nonlinear oscillatory boundary analysis,
and physical zero-leaf/domain definition; E4 transients and physiological closure
remain deferred.

## SPAC–carbon Stage A (2026-10-09)

The [user-plan Stage A analysis](spac_carbon_stage_a.md) adds an independent
six-state one-way Manzoni/four-carbon model. It proves the concave-canopy R0
threshold and conditional unique positive equilibrium, input-off zero capacity
and stock decay, block-triangular equilibrium lifting, and no irreversible
finite-temperature-pulse extinction when final R0>1 and water remains safe.
VPD fold crossing is not R-tipping; persistent-temperature extinction crosses
the carbon threshold. Synthetic GPP/loss/temperature coefficients and zero-leaf
hydraulic limitations are explicit. No VISIT or existing model changed.
Dedicated + neighboring regressions: 88 passed.
Outputs: `artifacts/spac_carbon/stage_a_final_20261009/`.
Next-step priority is superseded by the user Stage B plan above: first structural
flux feedback in a restricted positive-leaf domain; Leuning is postponed.

## Manzoni hydraulic reproduction (2026-10-07)

The user-directed [first reproduction milestone](manzoni2014_reproduction.md)
is complete: local paper and Appendix A read; equations (1)--(4), admissible
active equilibria, full Jacobian, and Fig. 2 supply/demand, phase and branch plots
implemented independently. Wet-root zero-limit and -1 MPa root boundaries give
critical D = 0.02439716312 and 0.01307440031 mol/mol. Capacitances are explicitly
assumed numerical fixtures. The fully nonconducting equilibrium is a continuum,
not an isolated full-system sink. No carbon coupling, mortality or R-tipping claim.
Dedicated and neighboring tests: 38 passed; equilibrium residual < 3e-15 MPa/day;
sampled trajectory tolerance-refinement error < 1.1e-8 MPa. Artifacts persist under
`artifacts/manzoni2014/figure2_20261007/` (git-ignored).
Follow-up: [remaining figure results](manzoni2014_remaining_results.md).
Fig. 3 thresholds and Appendix A's 0--3 finite-equilibrium counts are implemented,
with full Jacobians and unbounded log-domain root isolation. Soil-root continuity
and the fixed damage factor 0.6 are implemented; user-approved assumed soil
coefficients produce Fig. 4 sensitivity, not original numerical replication.
Final artifacts: `artifacts/manzoni2014/remaining_final_20261007/`.
All 3364 common coarse/fine A2 counts agree. Regressions now total 66 passed.
Original observational percentiles/points, soil coefficients and capacitances
remain unresolved; no carbon/cover coupling or R-tipping experiment was added.
Existing cover/carbon and native VISIT results are unchanged.

## Cover, temperature and humidity plan (2026-10-02)

The [current plan](forest_grass_cover_climate_research_plan.md) adds a grass-cover
fraction to ten carbon stores as a new analytical baseline. Its next tasks are
the carbon-equilibrium extension proof, sourced climate response laws, and
no-R/R conditions under warming and decreasing relative humidity. This is a
documentation-led direction. The 11-state baseline is implemented below;
climate experiments remain unimplemented.

### Stage A baseline (2026-10-02)

The 11-state area-exclusive cover/carbon baseline is now implemented in
`src/control_carbon/cover_carbon_baseline.py`, with its derivation and assumptions
in `docs/cover_carbon_baseline.md`. It reuses the existing published cover vector
field, adds ten carbon pools with area-conversion budgets, and provides carbon
capacity, potential, an analytic Jacobian, and budget diagnostics. The block
triangular Jacobian proves that each stable interior cover equilibrium lifts to a
stable full equilibrium; distinct interior cover states have distinct grass
living-carbon pools, while total-carbon ordering is not generally determined.
The explicit carbon coefficients used in dedicated tests are uncalibrated fixtures.
No temperature/relative-humidity mapping or climate-driven tipping result is
implemented or claimed. Next: source-constrain one climate-to-process response
before considering forcing-rate experiments.

The baseline is visualized in the
[cover/carbon phase atlas](figures/cover_carbon_baseline/cover_carbon_baseline_atlas.pdf):
it includes the parameter phase diagram, scalar cover phase line, adopted
function shapes, and fixed-cover carbon capacities. Reproduce it with
`python examples/plot_cover_carbon_baseline.py`.

### Climate-response literature ledger (2026-10-02)

The [response ledger](cover_climate_response_ledger.md) maps \(T,R\) to
atmospheric VPD and separates candidate effects on NPP, recruitment, fire/fuel
moisture, and soil decomposition. Only the physical \(T,R\to\mathrm{VPD}\)
transform is ready for a pure-function implementation. Climate coefficients for
the 11-state model remain uncalibrated; \(\alpha\), the cover sigmoid, and soil
decomposition stay fixed until appropriate site measurements are available.

## Tree–grass carbon–climate plan (2026-09-30)

The [research plan](tree_grass_carbon_climate_research_plan.md) now has a first
[process ledger](tree_grass_process_ledger.md), [provisional ten-pool equations](tree_grass_carbon_equations.md),
and [conditional proofs](tree_grass_bistability_analysis.md). These are independent,
uncalibrated theoretical equations, not a native VISIT or SPITFIRE implementation.
Shade-only positive-equilibrium uniqueness, full-system stability and a conditional
no-R result are derived. With endogenous fire, equilibrium elimination yields a
scalar self-consistency equation, with full-system stability still to be checked.
A reproducible checker validates carbon budgets, boundary signs, frozen capacity,
equilibrium elimination and the rank-one Jacobian identity, including step refinement.
See [first exploration record](tree_grass_first_exploration.md) and
[check results](tree_grass_carbon_checks.json). Climate forcing maps, fuel moisture,
resprouting and fire bistability/R-tipping remain unresolved.

## Forest–grass analytic checkpoint (2026-09-30)

Visualization: [Japanese five-page phase atlas](forest_grass_phase_atlas.md)
and PNG/vector-PDF figures in docs/figures/forest_grass_phase_atlas/.
The atlas distinguishes frozen basins, stochastic mean escape times, ramp
endpoints, and asymptotic rate-dependent destinations; it does not change the model.

Read the two user-provided PDFs; see forest_grass_nbr_tipping.md and its
Japanese TeX/PDF. Kumar K & Dutta (2026), (2.3)/(4.1), supplies a fire-mediated
cover model, NOT native VISIT. Kéfi et al. (2010) is vegetation–water, not
forest–grass; its steady water budget explains why retaining feedback matters.
Analytic results: two interior stable communities; nondegenerate B-folds;
reflected-noise first-passage formula for N; no low-to-high R-tipping for
alpha-only forcing confined to a bistable interval; a constructed coupled
alpha/phi1 path has low-speed tracking, high-speed tipping, and a unique
critical duration. The constructed path is an existence witness, NOT a
calibrated temperature response or a reproduction of the paper's rate value.
New independent module/tests/runner: forest_grass_tipping.py,
test_forest_grass_tipping.py, analyze_forest_grass_tipping.py.
Artifacts: artifacts/forest_grass_tipping/analytic_20260930/ (git-ignored).
VISIT source audit confirms top-down shade but does not establish native
bistability or a grass–fire closure. No native code or older experiment changed.
Validation: 7 dedicated tests; full suite 257 passed / 7 skipped.
Next: constrain an ecologically justified reciprocal feedback and forcing path,
then attach mass-conserving carbon pools before making flux/capacity claims.

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
