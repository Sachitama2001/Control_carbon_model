# Control Carbon Model — IDE Coding Agent Handoff

> 2026-09-21 native-state follow-up: pinned 2013 VISIT was compiled and run
> with supplied TKY forcing/deposition/GHG and reconstructed text parameters.
> `VISIT_TKY_matrix_equations_native_20000715.xlsx` now has all 37 start states,
> 15 histories and post-location forcing for 2000-07-15. Missing inputs fell
> 70 -> 2 (only two source-undefined diagnostics). All 37 one-day endpoints
> match native C, maximum scaled error 5.44e-16. Runtime Saxton capacities
> replace raw Config capacities (104.8103 / 704.0597 mm). Native execution also
> exposed a missing `==0` in a `strcmp` vapor-pressure branch; actual TKY vp is
> retained and flagged. See `docs/visit_native_snapshot_20000715.md`.

> 2026-09-17 actual-data follow-up: summer snapshot workbook
> `artifacts/visit_matrix_workbook/VISIT_TKY_matrix_equations_20000715.xlsx`.
> `visit_workbook_data.py` imports supplied visitb data read-only, retaining the
> pinned equations/Config. 24 of 70 missing inputs resolved (14 forcing/aggregate
> inputs, 10 histories reconstructed with initially unknown memory); 46 remain.
> All 37 summer states remain unknown. The 21 spinup-end restart stocks are a
> separate reference sheet, not substituted into summer X. New sheets 06–08
> record conversion, provenance, restart and SHA256 comparisons. Older template
> remains unchanged. Source-version differences in mean-temperature period,
> annual-rain accumulation and VPD branch are documented; pinned rules win.

> 2026-09-17 TKY matrix workbook: see `docs/visit_matrix_workbook.md`.
> User-requested Excel audit covers the selected old `visit_local` daily
> C/N/water path with 37 states, 132 flux channels and a 37×78 B matrix.
> Config_TKY defaults are imported; actual states/weather/history remain blank.
> `visit_workbook.py`, `_export.py`, `_validation.py` generate real Excel formulas,
> cached values and source/unknown inventories. Seven artificial cases match
> 1,379 native-C values to rounding error; 13 dedicated tests pass. It is a
> source-order daily effective-rate representation, not a continuous VISIT ODE
> or the 30-minute pathway. Native baseflow/N-unit/undefined-gas issues are
> exposed, not silently fixed. Excel-application recalculation remains untested.

> 2026-09-17 capacity/QSE exploration: `docs/capacity_qse_attractor_exploration.*`
> derives the state-dependent capacity map, its self-consistency equation, and
> the full-equilibrium Jacobian `J=M(I-D capacity)`. A model-agnostic positive
> scalar counterexample has two nonzero stable capacity branches and gives
> different final branches for 20-year and 120-year ramps; it is not a VISIT
> result or calibration. `docs/visit_r_tipping_process_screening.*` audits the
> 2020 explanatory slides and the pinned old VISIT source. The source-grounded
> next target is the complete daily plant map, especially LAI--GPP--allocation--NSC.
> New module/tests/example: `capacity_feedback`; nine dedicated tests pass.
> Full suite at this checkpoint: 146 passed, 7 source-dependent skips.

> 2026-09-15 carbon–nitrogen follow-up: `docs/nitrogen_tipping_minimal.md` and TeX/PDF.
> Independent 2-state fixed-stoichiometry benchmark: monotone N limitation has no
> permanent R-tipping when the final positive equilibrium exists. An explicitly
> hypothetical excess-N inhibition law admits rate-dependent basin changes;
> the fast-step basin entry has an analytic sufficient proof. Equal-dose centered
> ramps and 3-state delayed recycling are tested. No VISIT/water coefficients changed.
> New module/tests/example: `minimal_nitrogen_tipping`; full suite: 144 passed.

> 2026-09-15 theory follow-up: see `docs/water_tipping_minimal.md` and its TeX/PDF.
> User-requested minimal R-tipping analysis is a separate uncalibrated benchmark
> (`minimal_water_tipping.py`), not a completion of the native calibration path.
> It proves a scalar rainfall-only no-tipping result and a finite-rate threshold
> for a paired environmental path, plus a conditional dynamic-water example.
> VISITa was inspected at `5c513196f21c1b1efb9ead540e3d40865b15e07e` but not run.
> Existing hydraulic unit/boundary caveats remain unresolved. Full suite: 129 passed.

> Active priority (2026-09-14): explicit water storage and SPAC transport in
> the land-carbon matrix framework. Read `docs/carbon_water_research_plan.md`,
> `docs/coupled_carbon_water_equations.md`, and
> `docs/visitc_carbon_water_source_ledger.md` before using the historical
> skeleton below. Deterministic tipping code remains a tested, deferred track.

## 0. Purpose of this document

This repository is intended to become a **model-agnostic control-theoretic framework for terrestrial ecosystem carbon-cycle models**. VISIT is the first source-grounded implementation and validation target, not the final scope.

The immediate coding-agent task is **not** to produce a polished VISIT emulator as quickly as possible. The task is to build a rigorous chain from process-model source code to state-space/IRF representations while preserving provenance and clearly separating exact equations from reductions and local approximations.

The active chain is now:

```text
current VISITc water source
  -> source-order water ledger and native comparison
  -> four-component carbon-water balances
  -> coupled equilibrium/capacity/Jacobian diagnostics
  -> dynamic-vs-quasi-steady water experiments
  -> periodic tracking and conditional tipping analysis
```

Implemented in the current working tree:

- `visitc_source_map.py`: current `visit-manager/VISITc` provenance at
  `5202debd96df6f88beb7d61f8688fff02ace964a`;
- `visitc_hydrology.py`: `f_hydrology`, atmospheric/resistance, LAI/radiation,
  canopy-conductance, and PM transcriptions with an explicit baseflow audit;
- `visitc_native.py`: direct native-C checks for all exposed one-day hydrology,
  radiation, PM, and canopy-conductance outputs;
- `hydraulic_relations.py`: native VISITc soil retention and reduced
  SurEau-Ecos plant pressure-volume/capacitance relations;
- `coupled_carbon_water.py`: eight-state reduced ODE, budgets, carbon capacity,
  capacity-change decomposition, coupled equilibrium, Jacobian blocks, Schur
  reduction, and simulation;
- `examples/run_coupled_carbon_water.py`: idealized dry-down/re-wetting
  experiment with refinement and provenance outputs;
- dedicated tests and carbon-water documentation.

Do not call the eight-state continuous model “VISITc.” It is a synthesis. The
implemented hydrology/PM slices are native-validated, but they are not a
complete source-faithful VISITc model: initialization, all site branches, and
the complete carbon system remain outside the bridge.

The research objective is to make it possible to compare many terrestrial ecosystem models in a common language:

- state variables and flux topology,
- quasi-static equilibrium (QSE),
- disequilibrium from QSE,
- impulse response functions (IRFs),
- poles / modal time scales / residues,
- transfer functions and frequency responses,
- NEP/NEE/NECB input-output behavior,
- local linearization and, later, periodic/LTV response,
- forcing-rate response and analytical source/sink criteria.

## 1. Authoritative repositories

### Primary implementation repository

- `Sachitama2001/Control_carbon_model`

### Source/reference repository for the first model adapter

- `Sachitama2001/VISIT-matrix`
- pinned commit `3285bd8e131a932e338b59892751648fd9edcc7b`
- authoritative VISIT site-model source is under `visit_local/`

Important: `VISIT-matrix/visit_matrix/` already contains a reduced 18-pool matrix implementation. Use it as prior work and a cross-check, **not as the authority**. The C source in `visit_local/` is the authority when deciding whether an equation, state, input, or parameter is actually present in VISIT.

### Source/reference repository for the active water audit

- `visit-manager/VISITc`
- pinned commit `5202debd96df6f88beb7d61f8688fff02ace964a`
- inspected source is under `point/`

Keep this provenance separate from the older VISIT-matrix adapter. A formula
found in one snapshot is not automatically authoritative for the other.

## 2. Current repository status

> Status note: this section describes the original skeleton. The implemented
> checkpoint has advanced substantially; read `docs/current_status.md` before
> using the older task ordering below.

Current package skeleton:

```text
Control_carbon_model/
├── README.md
├── HANDOFF.md
├── pyproject.toml
├── docs/
│   └── visit_state_space_source_map.md
├── src/control_carbon/
│   ├── __init__.py
│   ├── state_space.py
│   └── visit_source_map.py
└── tests/
    └── test_state_space.py
```

Existing functionality:

1. `state_space.py`
   - generic continuous-time LTI representation,
   - poles,
   - transfer-function evaluation,
   - frequency response,
   - output IRF `C exp(A t) B`,
   - state IRF `exp(A t) B`.

2. `visit_source_map.py`
   - source-provenance registry for VISIT variables/processes.

3. `docs/visit_state_space_source_map.md`
   - first source-grounded decomposition of VISIT into state/input/output/process layers.

4. `tests/test_state_space.py`
   - basic analytical consistency checks for generic LTI utilities.

## 3. Scientific design rule

Every VISIT-derived model object must belong to one explicit approximation level.

### Level 0 — native nonlinear discrete map

```text
x[k+1] = F_VISIT(x[k], z[k], q[k]; theta)
y[k]   = H_VISIT(x[k], z[k], q[k]; theta)
```

This is the conceptual reference. VISIT is run daily and contains sequential updates and hybrid/piecewise rules.

### Level 1 — source-grounded reduced compartment model

```text
dx/dt = M(t) x + B(t) mu(t)
```

or, where the native daily discretization matters,

```text
x[k+1] = A_d[k] x[k] + B_d[k] mu[k].
```

The soil carbon subsystem is especially close to this form because, conditional on environmental scalars, decomposition is first-order in carbon stocks.

### Level 2 — local linearization of the native model

```text
delta x[k+1] = A_d delta x[k] + B_d delta z[k]
delta y[k]   = C_d delta x[k] + D_d delta z[k]
```

with Jacobians evaluated at an explicitly defined operating point.

### Level 3 — periodic linear time-varying representation

Because phenology and meteorology are strongly seasonal, an annual periodic trajectory is likely a more physically meaningful reference than a single annual-mean equilibrium.

### Level 4 — model-agnostic comparison layer

Different terrestrial ecosystem models should expose the same abstract interfaces even when their native state dimensions and process structures differ.

Never collapse these levels into a single class called simply `VISITModel` without metadata indicating which approximation is being used.

## 4. Central mathematical targets

For a nonlinear model

```text
x_dot = f(x, z, theta)
y     = h(x, z, theta)
```

or its discrete equivalent, implement tools to obtain a local state-space system

```text
delta x_dot = A delta x + B delta z
delta y     = C delta x + D delta z
```

and then

```text
G(s) = C (sI - A)^(-1) B + D.
```

For carbon stock total

```text
Ctot = w^T x
```

and NEP-like outputs, the long-term target is to connect:

```text
QSE geometry -> disequilibrium -> forcing rate -> NEP sign
```

with IRF / transfer-function quantities.

A recurring object will be the moving QSE

```text
x_star(z) : f(x_star, z) = 0
```

and its sensitivity

```text
dx_star/dz.
```

In an LTI/local-linear setting, the low-frequency response to forcing rate should recover the QSE sensitivity. This identity should eventually become a tested theorem/example in the repository rather than remain only in notes.

## 5. First model: VISIT

### 5.1 Minimal structural-carbon state

Use the following 18-pool ordering as the first reduced carbon state:

```text
0  tree_fol
1  tree_stm
2  tree_rot
3  c3_fol
4  c3_stm
5  c3_rot
6  c4_fol
7  c4_stm
8  c4_rot
9  ltr_tf
10 ltr_tc
11 ltr_tr
12 ltr_gf
13 ltr_gc
14 ltr_gr
15 msl_a
16 msl_i
17 msl_p
```

Do not describe this as the complete native VISIT state. `nsch_storage`, crop `grain`, sub-compartments, phenology memory, water stores, and potentially other variables can carry dynamical memory.

### 5.2 External forcing candidates

At minimum track:

```text
Tair, Tsurface, Tsoil10, Tsoil200,
SW radiation, precipitation, humidity/VPD, wind, atmospheric CO2
```

For early IRF work, effective inputs are allowed:

```text
GPP or EPP by PFT,
f_tm_l,
f_tm_h or pool-specific humus scalars,
allocation fractions,
phenology/regime variables.
```

But label them explicitly as **effective/internal drivers**, not meteorological forcing.

### 5.3 Outputs

At minimum:

```text
GPP, NPP, heterotrophic respiration, NEP, total carbon
```

Keep `NEP` distinct from `NECB`. Harvest, BVOC, DOC export, disturbance, and other lateral/non-respiratory losses can cause `dCtot/dt`, NEP, and NECB to differ.

## 6. Source-provenance requirement

Every VISIT adapter equation or matrix element must have provenance metadata with:

- repository,
- source path,
- function name,
- source variable names,
- equation/rule summary,
- approximation notes,
- regime assumptions,
- unit conversion.

A matrix coefficient with no provenance should not be merged.

Preferred pattern:

```python
SourceRef(
    repository="Sachitama2001/VISIT-matrix",
    path="visit_local/soil_proc.c",
    function="f_cycle_soil",
    symbol="degrade_tf",
    note="ltr_tf * sr_lf / 1000 * f_tm_l",
)
```

## 7. Implementation order

Follow the roadmap in `docs/implementation_roadmap.md`. The short version is:

1. Harden generic state-space API.
2. Complete provenance representation.
3. Implement VISIT soil subsystem exactly from C equations.
4. Verify one-step mass balance algebraically.
5. Implement reduced plant structural-carbon subsystem.
6. Make every approximation explicit.
7. Build discrete-time local-linear tools.
8. Build QSE solver and QSE sensitivity.
9. Implement IRF/step/ramp/frequency diagnostics.
10. Derive and test NEP transfer-function identities.
11. Compare reduced model IRFs against native VISIT perturbation experiments.
12. Add periodic/LTV treatment.
13. Generalize adapters for additional ecosystem models.

## 8. Definition of done for the first major milestone

Milestone M1 is complete when all of the following are true:

- the 18-pool reduced VISIT carbon model can be built from documented parameters and drivers;
- each nonzero matrix/flux relation is traceable to VISIT C source;
- one-day soil and plant-litter mass-balance unit tests pass;
- NEP from the reduced model is explicitly computed and tested;
- constant-forcing QSE can be solved;
- the local LTI/discrete state-space system around QSE is constructible;
- impulse and step responses are available for at least EPP, temperature/decomposition-scalar perturbations;
- direct time integration agrees with convolution/IRF in the linearized regime;
- README contains a reproducible example;
- no result is labelled "VISIT" without its approximation level.

## 9. Things the coding agent must not do

- Do not silently copy numerical defaults from `VISIT-matrix/visit_matrix` and call them VISIT parameters unless verified against C parameter initialization/input files.
- Do not assume `LAI` is always an independent state; in normal code it is recalculated from foliage mass, but fixed-phenology modes can overwrite it.
- Do not assume the native code is a simultaneous continuous ODE; it performs sequential daily updates.
- Do not treat all forcing as temperature-only. Temperature is the first scientific focus, but soil water, radiation, CO2, and phenology can dominate particular paths.
- Do not infer exact native Jacobians from the reduced 18-pool model.
- Do not optimize performance before establishing source fidelity and tests.
- Do not merge large undocumented refactors.

## 10. Recommended working style

For each new VISIT process:

1. locate the exact C routine;
2. write the equation in a source note;
3. identify state/input/parameter/output roles;
4. implement the smallest pure Python function reproducing the algebra;
5. add a unit test with hand-computable values;
6. only then insert it into a matrix or nonlinear adapter;
7. add provenance metadata;
8. document any approximation or branch excluded.

## 11. Expected eventual architecture

```text
src/control_carbon/
├── core/
│   ├── continuous.py
│   ├── discrete.py
│   ├── ltv.py
│   ├── irf.py
│   ├── qse.py
│   ├── modal.py
│   └── reduction.py
├── provenance/
│   └── source.py
├── adapters/
│   ├── base.py
│   └── visit/
│       ├── schema.py
│       ├── source_map.py
│       ├── soil.py
│       ├── plant.py
│       ├── outputs.py
│       ├── reduced18.py
│       ├── linearize.py
│       └── native_bridge.py
└── experiments/
    ├── irf.py
    ├── qse_tracking.py
    └── forcing_rate.py
```

Do not move files merely to match this tree immediately; migrate incrementally with tests.

## 12. Long-term research direction

The final scientific contribution should not be "a matrix version of VISIT". The intended contribution is a framework able to answer, across many terrestrial ecosystem models:

- what are the dominant carbon-cycle modes and time scales?
- which forcing frequencies project onto which modes?
- which internal modes are visible in NEP and which are hidden by cancellation?
- when does NEP change sign because of modal competition?
- how does QSE motion map to forcing-rate thresholds?
- how do model parameters shift poles, residues, QSE sensitivities, and source/sink boundaries?
- which features are robust across model structures?

VISIT is the first testbed for proving that this program can be carried out rigorously from real ecosystem-model source code.
