# Control Carbon Model

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
