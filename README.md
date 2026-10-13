# Control Carbon Model

2026-10-13: [B3 E4 nonlinear results](docs/spac_carbon_stage_b3_e4_results.md)
use the latest MPa-facing scale review to prioritize E4 within routine VPD.
Exact four-state derivatives and independent periodic collocation numerically
verify local subcritical Hopf and repelling full-six-state cycles.
Seven predefined sensitivities expose dependence on wood allometry/loss
assumptions; no empirical mortality/global basin/R-tipping claim.
163 targeted/neighboring tests pass.
Run `python examples/analyze_spac_carbon_e4.py --output NEW_PATH`.

2026-10-13: Current SPAC analyses use **VPD in MPa**, still a positive gas
pressure deficit, not atmospheric water potential. Start with the Japanese
[units and scales review](docs/spac_carbon_units_and_scales.md),
[MPa configuration](configs/spac_carbon_pressure_review.json) and
[PressureSPAC interface](src/control_carbon/spac_carbon_mpa.py).
Explicit assumptions: atmospheric pressure 0.101325 MPa, diagnostic 25 C,
routine VPD 0--0.003 MPa. High-VPD E2/E3 Hopf results are formal exploration,
not valid forcing at 25 C. Archived paper/Stage A/B1/B2 APIs, configs and
runners below retain mol/mol for historical reproduction.
New unit tests: 27 passed; six related baseline files: 127 passed.
The review additionally checks E4's exact four-state local Hurwitz crossing;
E4 nonlinear Hopf classification remains open.
Run `python examples/review_spac_carbon_units.py --output NEW_PATH`.

2026-10-13: Start with the Japanese
[research logic and methods report](docs/spac_carbon_research_logic_report.md)
for an accessible account of why each Stage A/B1/B2 method was chosen.
It separates analytic deductions, continuous interval certificates,
numerical verification and unresolved global/physiological hypotheses.
This is a documentation-only synthesis, not a new experiment or tipping claim.

2026-10-12: [Stage B2 plan](docs/SPAC_carbon_stage_B2_plan.md) and
[results](docs/spac_carbon_stage_b2_results.md) add continuous-path
Routh–Hurwitz interval certificates for E2/E3, analytic higher derivatives
and numerical subcritical-Hopf verification. Nonzero-amplitude shooting
and full-six-state Floquet analysis confirm small **repelling** periodic
orbits, not alternate attractors. Zero-leaf/domain-exit physical requirements
are documented; Stage A/B1 are unchanged. 129 targeted/neighboring tests pass.
Run `python examples/analyze_spac_carbon_bifurcation.py --output NEW_PATH`.

2026-10-09: User-directed [Stage B1 results](docs/spac_carbon_stage_b_results.md)
add independent E0–E4 structural water-flux feedback, fixed reference biomass
and capacities, exact series conductance, full 6x6 Jacobians, pseudo-arclength,
multistart roots, basin slices and bounded sensitivities. Leaf supply alone
lowers the synthetic fold; demand feedback can instead produce oscillatory
stability candidates. Same-path slow/fast experiments demonstrate
vegetation-domain exit versus positive-equilibrium recovery, **not verified
R-tipping or mortality**. Stage A and native VISIT are unchanged; 115 targeted/
neighboring tests pass. Leuning is deferred under the
[new user plan](docs/SPAC_carbon_stage_B_plan.md).
Run `python examples/analyze_spac_carbon_two_way.py --output NEW_PATH`.

2026-10-09: User-directed [SPAC–carbon Stage A](docs/spac_carbon_stage_a.md)
independently couples checked Manzoni hydraulics to four carbon pools.
Analytic R0 thresholds, input-off decay, equilibrium lifting and a conditional
no-temperature-pulse-extinction theorem are checked alongside VPD fold and
temperature controls (88 dedicated/neighboring tests). This is synthetic theory,
not VISIT, Leuning, calibrated mortality, or a genuine R-tipping demonstration.
Run `python examples/analyze_spac_carbon.py --output NEW_PATH`.

2026-10-07: User-directed [Manzoni (2014) hydraulic reproduction](docs/manzoni2014_reproduction.md).
An independent two-storage model reproduces Fig. 2 mechanisms and active-branch
fold thresholds, with full-Jacobian and integration checks (38 dedicated/neighboring
tests passed). VPD is mol/mol; numerical capacitances are not paper-calibrated.
This is not a native VISIT, carbon-coupling, or R-tipping result.
Run `python examples/reproduce_manzoni2014.py --output NEW_PATH`.

The [remaining-figure follow-up](docs/manzoni2014_remaining_results.md) adds
Fig. 3 thresholds, sigmoid finite-equilibrium counts (0--3), and root continuity.
User-approved assumed soil coefficients support a Fig. 4 sensitivity study, not
original numerical replication. Coarse/fine common counts agree; 66 dedicated/
neighboring tests pass. Observational points/percentiles and original soil values
remain unresolved. Run `python examples/reproduce_manzoni2014_remaining.py --shape-step .25 --output NEW_PATH`.

2026-10-02: Start the current tree–grass track with the
[forest/grass cover and climate research plan](docs/forest_grass_cover_climate_research_plan.md).
It specifies one cover fraction plus ten carbon pools, analytical transfer of
bistability to carbon storage, temperature/relative-humidity forcing, no-R and
R-tipping conditions, and carbon capacity diagnostics. Documentation only;
climate-driven tipping remains a research question.

2026-09-30: The current tree–grass research direction is specified in the
[Japanese carbon–climate research plan](docs/tree_grass_carbon_climate_research_plan.md).
It starts from ten carbon pools (tree/grass × leaf/stem/root/litter/humus),
prioritizes analytical bistability and climate-driven rate-induced tipping,
and distinguishes frozen-coefficient carbon capacity from nonlinear equilibria.
First exploration: [process/source ledger](docs/tree_grass_process_ledger.md),
[provisional ten-pool equations](docs/tree_grass_carbon_equations.md), and
[conditional analytic results](docs/tree_grass_bistability_analysis.md).
The shade-only control has at most one positive equilibrium and a conditional
no-R result in the full ten-pool system. Fire equilibria reduce to a scalar
self-consistency equation, but fire bistability and climate-rate tipping remain
unproved. Climate responses and post-topkill recovery are not calibrated.
Existing cover-model results remain separate.

2026-09-22: Added a separate
[temperature–NSC R-tipping research plan](docs/temperature_nsc_r_tipping_research_plan.md).
It specifies a homogeneous, area-based four-structural-carbon/four-mobile-carbon
model informed by the organ-level NSC structure of Ninomiya et al. (2023).
Only temperature is varied; leaf-dependent photosynthesis, temperature-dependent
respiration, NSC-funded leaf flush, and NSC-concentration-dependent mortality are
tested for basin geometry and genuine rate-induced transitions. No R-tipping result
is claimed by the plan itself.

2026-09-22: Added a
[temperature–precipitation tipping research plan](docs/temperature_precipitation_tipping_research_plan.md)
for a deliberately low-dimensional theoretical track. It starts from matching
leaf/stem/root/soil carbon and water stores (eight states), treats temperature and
precipitation as the only varying external drivers, and separates quantitative
compound effects, frozen-system bifurcations, and genuine rate-induced basin
transitions. It complements rather than replaces the native VISIT audits and
the eight-state carbon–water ODE.

2026-09-17: Added a [VISIT process screening report](docs/visit_r_tipping_process_screening.md)
and a [capacity/QSE/attractor exploration](docs/capacity_qse_attractor_exploration.md).
The latter proves that a state-dependent Luo-style capacity is a map whose
self-consistent fixed points need not be unique, and gives a reproducible
nonzero-to-nonzero R-tipping counterexample. It does **not** claim that VISIT
has multiple attractors. The next source-grounded test is the complete plant
daily map, especially the LAI--GPP--allocation--NSC loop.

2026-09-15: Added a separate [minimal carbon–nitrogen R-tipping study](docs/nitrogen_tipping_minimal.md)
with a Japanese TeX/PDF report, elemental budgets, a monotone-limitation no-tipping
result, and a conditional excess-nitrogen benchmark. This is not a calibrated forest model.

A research project for analyzing terrestrial ecosystem carbon-cycle models
through matrix equations, explicit water dynamics, control theory, impulse
responses, quasi-static equilibria (QSEs), and source-sink diagnostics.

The active research question is now:

> How does adding explicit water storage and SPAC transport to the land-carbon
> matrix framework change carbon storage capacity, disequilibrium, and response
> time?

The first implementation has two deliberately separate tracks:

- an auditable, native-validated transcription of current VISITc hydrology and
  its Penman-Monteith dependency chain;
- a reduced continuous model with leaf/stem/root/soil carbon and water stores.

R-tipping is a later, conditional analysis rather than a required result. See
the [carbon-water research plan](docs/carbon_water_research_plan.md) and
[equations](docs/coupled_carbon_water_equations.md). The older
[deterministic tipping plan](docs/deterministic_tipping_plan.md) is retained as
a deferred analysis track and benchmark record.

The long-term goal is **not to reproduce one specific land model**, but to build a model-agnostic framework that can compare many terrestrial ecosystem models through common state-space and input-output dynamical quantities.

## Current reference models

The established nine-pool soil-carbon adapter is sourced from
`Sachitama2001/VISIT-matrix/visit_local` at commit
`3285bd8e131a932e338b59892751648fd9edcc7b`.

The new water audit targets [`visit-manager/VISITc`](https://github.com/visit-manager/VISITc)
at commit `5202debd96df6f88beb7d61f8688fff02ace964a`. These snapshots are
not interchangeable; their provenance registries and documentation are kept
separate.

The existing Python matrix implementation in `VISIT-matrix/visit_matrix` is useful prior work, but the VISIT C source is treated as the authority for state definitions, update order, process equations, and provenance.

For the current paper-oriented capacity question, start with
[the Japanese TeX report](docs/capacity_qse_attractor_exploration.tex) or its
[compiled PDF](docs/capacity_qse_attractor_exploration.pdf). The companion
[VISIT slide/source screening PDF](docs/visit_r_tipping_process_screening.pdf)
records why allocation and plant carbon memory are prioritized over individual
photosynthesis or respiration response curves.

## Start here

Japanese mathematical report (2026-09-15):
[PDF](docs/carbon_water_matrix_report.pdf) /
[TeX source](docs/carbon_water_matrix_report.tex) /
[build instructions](docs/carbon_water_matrix_report_build.md).
It develops the model-agnostic carbon–water matrix framework and records
corrections and outstanding unit/boundary issues in the earlier equation notes.

Minimal R-tipping follow-up: [Japanese TeX](docs/water_tipping_minimal.tex) /
[PDF](docs/water_tipping_minimal.pdf) /
[reproduction and scope](docs/water_tipping_minimal.md).
This is an independent, uncalibrated theoretical benchmark, not evidence of
R-tipping in VISIT. It includes a rainfall-only no-tipping result, an analytical
finite-rate threshold for a paired environmental path, and a two-state water-memory example.

For coding agents and new contributors:

1. `AGENTS.md` - coding-agent rules and the active task order.
2. `docs/carbon_water_research_plan.md` - scientific story, hypotheses, experiment ladder, and acceptance criteria.
3. `docs/coupled_carbon_water_equations.md` - complete equations and diagnostics for the eight-state model.
4. `docs/visitc_carbon_water_source_ledger.md` - current VISITc stores, fluxes, call order, source links, and audit findings.
5. `docs/literature_index.md` - local PDFs, stable links, and evidence-to-equation routing.
6. `HANDOFF.md` and `docs/current_status.md` - implementation checkpoint and immediate next work.
7. `docs/architecture_and_conventions.md` - package, time, sign, unit, and testing conventions.
8. `docs/implementation_roadmap.md` - full engineering and research backlog.

## Active goals

- Reproduce VISITc water bookkeeping without changing its source order.
- Couple water stores and transport to a small carbon compartment system.
- Distinguish instantaneous carbon capacity, coupled equilibrium, and periodic tracking trajectories.
- Decompose capacity change into productivity, transit-time, and interaction effects.
- Test when dynamic water storage differs from a quasi-steady hydraulic reduction.
- Preserve model provenance and mass balance at every abstraction level.
- Retain the existing IRF, modal, native-soil, ERA5, and nonlinear benchmark infrastructure.

## Reproduce the first carbon-water experiment

```bash
python -m pip install -e '.[dev,era5]'
python examples/run_coupled_carbon_water.py
pytest -q
```

The example starts at a stable coupled equilibrium, applies an idealized
precipitation dry-down and re-wetting, and writes NPZ data plus a JSON manifest
under `artifacts/coupled_carbon_water/`. It includes carbon/water budget
residuals and a maximum-step refinement comparison. The parameter set is
illustrative, not calibrated VISITc.

## Carbon-water code

`src/control_carbon/visitc_source_map.py` pins current VISITc source provenance.

`src/control_carbon/visitc_hydrology.py` follows
`point/hydro_balance.c::f_hydrology` in source order while taking seven
Penman-Monteith potentials as effective inputs. It also transcribes the native
atmospheric, LAI/radiation, canopy-conductance, and PM dependency chain. It deliberately exposes the
apparent double subtraction of baseflow instead of silently correcting it.

`src/control_carbon/visitc_native.py` builds minimal bridges against the pinned
current VISITc source. One-day hydrology stores/fluxes, PM potentials,
radiation partitioning, and canopy conductance are compared directly with C.

`src/control_carbon/hydraulic_relations.py` implements the VISITc soil
storage-to-potential relation and an explicitly reduced SurEau-Ecos
pressure-volume/capacitance relation for plant stores. The full mathematical
specification is in
[`docs/carbon_water_differential_equations.md`](docs/carbon_water_differential_equations.md).

`src/control_carbon/coupled_carbon_water.py` implements the eight-state reduced
ODE with potential-gradient SPAC transport, incidence-based water balance, frozen carbon capacity, Wei-style capacity
decomposition, coupled equilibria, finite-difference Jacobian blocks, the water
Schur complement, and time integration.

## Existing code

The new [nonlinear core](src/control_carbon/nonlinear.py) provides single-driver
continuous integration, frozen-equilibrium diagnostics, and signed compartment
algebra. [Deterministic benchmarks](src/control_carbon/tipping.py) add analytic
saturating-feedback branches, a fold experiment, and finite-ramp basin tests.
These hypothetical models are not native VISIT equations or TKY predictions.

Reproduce the temperature-only VISIT soil control and synthetic B/R figures:

```bash
python -m pip install -e '.[dev,analysis]'
PYTHONPATH=src python examples/plot_deterministic_tipping.py
```

The script writes a PNG, numerical NPZ arrays, and a JSON experiment manifest
under `artifacts/`. No CDS access or native compiler is needed for this example.

`src/control_carbon/state_space.py` contains generic continuous- and discrete-time linear state-space utilities, including impulse responses, forced simulation, and constant-input fixed points.

`src/control_carbon/provenance.py` and `src/control_carbon/visit_source_map.py` provide machine-readable source provenance pinned to the authoritative VISIT source revision.

`src/control_carbon/visit_decomposition.py` reconstructs the VISIT litter and humus temperature/moisture multipliers. The undefined `frh()` branches in the source snapshot are rejected by default; the apparent `stype=0/1/2` correction is available only as an explicitly labelled inferred behavior.

`src/control_carbon/visit_soil.py` implements the source-derived daily 9-pool soil carbon balance, constant-environment fixed points, heterotrophic-respiration impulse responses, and common-condition trajectory data for comparing IRF predictions with direct daily simulation.

`src/control_carbon/visit_native.py` builds a minimal harness against the pinned VISIT `soil_proc.c` and `decomposition.c` sources and compares native C baseline/perturbation trajectories with the Python 9-pool implementation and its IRF. See `docs/native_visit_soil_bridge.md` for scope and exclusions.

Grouped discrete soil modes, basis-independent participation, Rh residues, and implicit fixed-point environmental sensitivities are documented in `docs/visit_soil_modal_qse.md`.

`src/control_carbon/visit_plant.py` begins the P5 plant reconstruction with source-grounded foliage/stem/root turnover, Q10 and size-dependent respiration, and piecewise allocation. Each slice has hand calculations and direct native C checks, but they are not yet integrated into a complete plant daily update.

`src/control_carbon/era5.py` provides point-scale ERA5 acquisition, strict UTC daily aggregation, an intended VISIT `WMODE=1` NCEP/NCAR climate writer, and monthly request construction for direct Copernicus CDS retrieval. See `docs/era5_visit_weather.md` for source-format differences, known compatibility defects in the pinned VISIT snapshot, transformation assumptions, and the tested command-line workflow.

The native soil bridge supports day-varying environments, local environmental
Jacobians, modal analysis, QSE sensitivity, and a native/nonlinear/tangent
comparison figure. Plant turnover, respiration, and allocation slices have
also been individually matched to native C. Their full daily integration
remains pending. Immediate work now follows the
[carbon-water plan](docs/carbon_water_research_plan.md); deterministic tipping
benchmarks remain available but are not the main implementation driver.
